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

# MATCH THE FILENAME FROM YOUR TRAINING LOG
MODEL_FILENAME = 'siren_plateau_filter_model.pth' 
OUTPUT_FILENAME = 'siren_filters_plateau.mat'

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
        # omega_0=1.0 matches your "smooth blob" training settings
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
        Kxn = np.array(f['Kxn']).transpose()
        Kyn = np.array(f['Kyn']).transpose()
        Zn  = np.array(f['Zn']).transpose()
        
        original_shape = Kxn.shape
        print(f"Grid Shape (Z, Y, X): {original_shape}")

    # --- 3. Flatten Coordinates (NO ABSOLUTE VALUE) ---
    # We use the raw coordinates to respect specific quadrants
    coords_flat = np.stack([Kxn.flatten(), Kyn.flatten(), Zn.flatten()], axis=1)
    
    # --- Normalize Coords to [-1, 1] ---
    min_vals = coords_flat.min(axis=0)
    max_vals = coords_flat.max(axis=0)
    denom = max_vals - min_vals
    denom[denom == 0] = 1.0
    coords_flat = 2 * (coords_flat - min_vals) / denom - 1
    
    # Convert to Tensor
    coords_tensor = torch.tensor(coords_flat, dtype=torch.float32)

    # --- 4. Load Trained Model ---
    model = DualFilterSiren().to(DEVICE)
    model_path = os.path.join(CODE_DIR, MODEL_FILENAME)
    
    if not os.path.exists(model_path):
        # Fallback check
        alt_path = os.path.join(CODE_DIR, 'siren_plateau_filter_model.pth')
        if os.path.exists(alt_path):
             model_path = alt_path
        else:
             raise FileNotFoundError(f"Model file not found: {model_path}")
        
    print(f"Loading model from: {model_path}")
    model.load_state_dict(torch.load(model_path, map_location=DEVICE))
    model.eval()

    # --- 5. Batch Inference ---
    batch_size = 65536
    num_points = coords_tensor.shape[0]
    probabilities = []
    
    print("Starting batch inference...")
    with torch.no_grad():
        for i in range(0, num_points, batch_size):
            batch = coords_tensor[i : i+batch_size].to(DEVICE)
            
            # Forward pass
            logits = model(batch)
            
            # --- CRITICAL: SIGMOID ---
            # Allows independent plateaus (max integral) and independent zeros (clean background)
            probs = torch.sigmoid(logits)
            
            probabilities.append(probs.cpu().numpy())
            
            if i % (batch_size * 10) == 0:
                print(f"Processed {i}/{num_points} points...", flush=True)

    # Concatenate all batches
    full_probs = np.concatenate(probabilities, axis=0)

    # --- 6. Reshape and Save ---
    # Col 0 = Positive Filter, Col 1 = Negative Filter
    mask_pos_flat = full_probs[:, 0]
    mask_neg_flat = full_probs[:, 1]
    
    # Reshape back to 3D (Nz, Ny, Nx)
    mask_pos_3d = mask_pos_flat.reshape(original_shape)
    mask_neg_3d = mask_neg_flat.reshape(original_shape)

    # Transpose BACK to (Nx, Ny, Nz) for MATLAB compatibility
    mask_pos_matlab = mask_pos_3d.transpose()
    mask_neg_matlab = mask_neg_3d.transpose()

    # Optional Cleanup: Snap very small values to pure 0
    mask_pos_matlab[mask_pos_matlab < 0.01] = 0
    mask_neg_matlab[mask_neg_matlab < 0.01] = 0

    # Save to .mat
    save_path = os.path.join(DATA_DIR, OUTPUT_FILENAME)
    print(f"Saving to {save_path}...")
    
    sio.savemat(save_path, {
        'mask_pos': mask_pos_matlab,
        'mask_neg': mask_neg_matlab
    })
    print("Done! Standard Plateau filters ready for MATLAB.")

if __name__ == "__main__":
    generate_filters()
