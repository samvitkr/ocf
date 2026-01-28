import torch
import torch.nn as nn
import numpy as np
import h5py
import scipy.io as sio
import os

# ================= CONFIGURATION =================
# Paths (Adjust if needed)
DATA_DIR = '/users/1/kuma0458/open_channel_ret180/data'
DATA_FILENAME = 'filter_calc_input.mat'
CODE_DIR = '/users/1/kuma0458/open_channel_ret180/codes'
MODEL_FILENAME = 'siren_filter_model.pth'
OUTPUT_FILENAME = 'siren_filters.mat'

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# =================================================

# --- 1. Re-define the Model Architecture ---
# (Must match training script exactly)
class SineLayer(nn.Module):
    def __init__(self, in_f, out_f, bias=True, is_first=False, omega_0=30):
        super().__init__()
        self.omega_0 = omega_0
        self.linear = nn.Linear(in_f, out_f, bias=bias)
        
    def forward(self, x):
        return torch.sin(self.omega_0 * self.linear(x))

class DualFilterSiren(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            SineLayer(3, 256, is_first=True, omega_0=1.0),
            SineLayer(256, 256, is_first=False, omega_0=1.0),
            SineLayer(256, 256, is_first=False, omega_0=1.0),
            nn.Linear(256, 2)
        )
    def forward(self, x): return self.net(x)

def generate_filters():
    print(f"Running Inference on: {DEVICE}")

    # --- 2. Load Coordinates from Input File ---
    input_path = os.path.join(DATA_DIR, DATA_FILENAME)
    print(f"Loading grid from: {input_path}")
    
    with h5py.File(input_path, 'r') as f:
        # Load grid and check shape
        # Recall: h5py reads as (Nx, Ny, Nz), so we transpose to (Nz, Ny, Nx)
        # to match Python's Z-first convention.
        Kxn = np.array(f['Kxn']).transpose()
        Kyn = np.array(f['Kyn']).transpose()
        Zn  = np.array(f['Zn']).transpose()
        
        original_shape = Kxn.shape
        print(f"Grid Shape (Z, Y, X): {original_shape}")

    # Flatten and Normalize (EXACTLY as done in training)
    coords_flat = np.stack([Kxn.flatten(), Kyn.flatten(), Zn.flatten()], axis=1)
    
    # --- CRITICAL: Normalize Coords to [-1, 1] ---
    # We must use the exact same logic as training
    min_vals = coords_flat.min(axis=0)
    max_vals = coords_flat.max(axis=0)
    denom = max_vals - min_vals
    denom[denom == 0] = 1.0
    coords_flat = 2 * (coords_flat - min_vals) / denom - 1
    
    # Convert to Tensor
    coords_tensor = torch.tensor(coords_flat, dtype=torch.float32)

    # --- 3. Load Trained Model ---
    model = DualFilterSiren().to(DEVICE)
    model_path = os.path.join(CODE_DIR, MODEL_FILENAME)
    
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")
        
    model.load_state_dict(torch.load(model_path, map_location=DEVICE))
    model.eval()
    print("Model loaded successfully.")

    # --- 4. Batch Inference (Prevent Out of Memory) ---
    batch_size = 65536
    num_points = coords_tensor.shape[0]
    probabilities = []
    
    print("Starting batch inference...")
    with torch.no_grad():
        for i in range(0, num_points, batch_size):
            batch = coords_tensor[i : i+batch_size].to(DEVICE)
            
            # Forward pass
            logits = model(batch)
            
            # Softmax to get probabilities (0 to 1)
            probs = torch.nn.functional.softmax(logits, dim=1)
            
            probabilities.append(probs.cpu().numpy())
            
            if i % (batch_size * 10) == 0:
                print(f"Processed {i}/{num_points} points...", flush=True)

    # Concatenate all batches
    full_probs = np.concatenate(probabilities, axis=0)

    # --- 5. Reshape and Save ---
    # full_probs has shape (N_points, 2). 
    # Col 0 = Probability of Class 0 (Positive Region)
    # Col 1 = Probability of Class 1 (Negative Region)
    
    mask_pos_flat = full_probs[:, 0]
    mask_neg_flat = full_probs[:, 1]
    
    # Reshape back to 3D (Nz, Ny, Nx)
    mask_pos_3d = mask_pos_flat.reshape(original_shape)
    mask_neg_3d = mask_neg_flat.reshape(original_shape)

    # Transpose BACK to (Nx, Ny, Nz) for MATLAB compatibility
    # MATLAB expects x-first
    mask_pos_matlab = mask_pos_3d.transpose()
    mask_neg_matlab = mask_neg_3d.transpose()

    # Save to .mat
    save_path = os.path.join(DATA_DIR, OUTPUT_FILENAME)
    print(f"Saving to {save_path}...")
    
    sio.savemat(save_path, {
        'mask_pos': mask_pos_matlab,
        'mask_neg': mask_neg_matlab
    })
    print("✅ Done! File ready for MATLAB.")

if __name__ == "__main__":
    generate_filters()
