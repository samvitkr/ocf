import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import h5py
import os
import sys

# ================= CONFIGURATION =================
DATA_DIR = '/users/1/kuma0458/open_channel_ret180/data'
CODE_DIR = '/users/1/kuma0458/open_channel_ret180/codes'
DATA_FILENAME = 'filter_calc_input.mat'
MODEL_FILENAME = 'siren_plateau_filter_model.pth'

# ==========================================
# 1. DATA LOADING (PLATEAU FILTER STRATEGY)
# ==========================================
def load_matlab_data():
    input_path = os.path.join(DATA_DIR, DATA_FILENAME)
    print(f"Loading data from {input_path}...")
    
    with h5py.File(input_path, 'r') as f:
        # Load arrays and transpose to (Nz, Ny, Nx)
        Kxn = np.array(f['Kxn']).transpose()
        Kyn = np.array(f['Kyn']).transpose()
        Zn  = np.array(f['Zn']).transpose()
        phin = np.array(f['phin']).transpose() # Flux Data
        
        # Load Weights (W) if available, otherwise derive from phin
        if 'W' in f:
            W = np.array(f['W']).transpose()
        else:
            W = np.abs(phin) 

    # --- Flatten & Normalize Coordinates ---
    coords_flat = np.stack([Kxn.flatten(), Kyn.flatten(), Zn.flatten()], axis=1)
    
    min_vals = coords_flat.min(axis=0)
    max_vals = coords_flat.max(axis=0)
    denom = max_vals - min_vals
    denom[denom == 0] = 1.0
    coords_flat = 2 * (coords_flat - min_vals) / denom - 1

    # --- THE FIX: BINARY TARGETS + AMPLITUDE WEIGHTS ---
    print("Generating Binary Targets for Integral Preservation...")
    
    # 1. BINARY TARGETS (0 or 1)
    # This forces the model to predict 1.0 (max filter) whenever the sign is correct.
    # It prevents the filter from dimming the signal in weak regions.
    
    # Positive Filter Target: 1 if flux > 0, else 0
    target_pos = (phin.flatten() > 0).astype(np.float32)
    
    # Negative Filter Target: 1 if flux < 0, else 0
    target_neg = (phin.flatten() < 0).astype(np.float32)
    
    targets_flat = np.stack([target_pos, target_neg], axis=1)

    # 2. CONTINUOUS WEIGHTS (Softness Control)
    # We use the normalized flux magnitude as the "importance".
    # - High Flux: Weight ~ 1. Model MUST predict 1.0.
    # - Low Flux: Weight ~ 0. Model is allowed to transition smoothly.
    mag = np.abs(phin.flatten())
    if mag.max() > 0:
        weights_flat = mag / mag.max()
    else:
        weights_flat = mag

    # Add a small floor (0.05) to ensure the background (Target 0) is actively suppressed
    weights_flat = weights_flat + 0.05

    print(f"Targets prepared. Positive Ratio: {target_pos.mean():.1%}, Negative Ratio: {target_neg.mean():.1%}")

    return (torch.tensor(coords_flat, dtype=torch.float32),
            torch.tensor(targets_flat, dtype=torch.float32),
            torch.tensor(weights_flat, dtype=torch.float32))

# ==========================================
# 2. MODEL DEFINITION
# ==========================================
class SineLayer(nn.Module):
    def __init__(self, in_f, out_f, bias=True, is_first=False, omega_0=30):
        super().__init__()
        self.omega_0 = omega_0
        self.linear = nn.Linear(in_f, out_f, bias=bias)
        
        with torch.no_grad():
            if is_first:
                self.linear.weight.uniform_(-1/in_f, 1/in_f)
            else:
                self.linear.weight.uniform_(-np.sqrt(6/in_f)/omega_0, np.sqrt(6/in_f)/omega_0)
        
    def forward(self, x):
        return torch.sin(self.omega_0 * self.linear(x))

class DualFilterSiren(nn.Module):
    def __init__(self):
        super().__init__()
        # omega_0=1.0 ensures smooth, large-scale blobs (prevents high-freq noise)
        self.net = nn.Sequential(
            SineLayer(3, 256, is_first=True, omega_0=1.0),
            SineLayer(256, 256, is_first=False, omega_0=1.0),
            SineLayer(256, 256, is_first=False, omega_0=1.0),
            nn.Linear(256, 2)
        )
    def forward(self, x): return self.net(x)

# ==========================================
# 3. TRAINING LOOP
# ==========================================
def train_model():
    # Setup Device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Running on: {device}")

    # 1. Load Data
    coords, targets, weights = load_matlab_data()
    
    # 2. Move to GPU
    print("Moving dataset to GPU...")
    coords = coords.to(device)
    targets = targets.to(device)
    weights = weights.to(device)

    # Initialize Model
    model = DualFilterSiren().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
    
    # Loss: Binary Cross Entropy with Logits
    # Suitable for binary targets (0/1) + sigmoid activation
    # reduction='none' allows manual weighting per pixel
    criterion = nn.BCEWithLogitsLoss(reduction='none')

    # Batching Setup
    N = coords.shape[0]
    batch_size = 65536
    num_batches = int(np.ceil(N / batch_size))
    
    epochs = 200 # Sufficient for convergence
    print(f"Starting Training (Integral Maximization Mode)...")

    model.train()

    for epoch in range(epochs):
        epoch_loss = 0
        
        # Shuffle indices
        indices = torch.randperm(N, device=device)
        
        for i in range(num_batches):
            start = i * batch_size
            end = min(start + batch_size, N)
            idx = indices[start:end]

            # Inputs (Already on GPU)
            batch_coords = coords[idx]
            batch_targets = targets[idx]          # Binary (0 or 1)
            batch_weights = weights[idx].unsqueeze(1) # Amplitude Weights

            optimizer.zero_grad()
            
            # Forward Pass
            logits = model(batch_coords)
            
            # Weighted Loss Calculation
            # Forces model to 1.0 where Target=1 & Weight=High
            # Forces model to 0.0 where Target=0 & Weight>0
            raw_loss = criterion(logits, batch_targets)
            loss = (raw_loss * batch_weights).mean()
            
            loss.backward()
            optimizer.step()
            
            epoch_loss += loss.item()

        if epoch % 10 == 0:
            avg_loss = epoch_loss / num_batches
            print(f"Epoch {epoch} | Loss: {avg_loss:.6f}", flush=True)

    # Save Model
    if not os.path.exists(CODE_DIR):
        os.makedirs(CODE_DIR)
        
    save_path = os.path.join(CODE_DIR, MODEL_FILENAME)
    torch.save(model.state_dict(), save_path)
    print(f" Model saved to {save_path}")

if __name__ == "__main__":
    train_model()
