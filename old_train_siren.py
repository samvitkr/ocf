import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
import numpy as np
import h5py
import os

# ==========================================
# 0. CONFIGURATION
# ==========================================
# Paths based on your MSI structure
DATA_DIR = '/users/1/kuma0458/open_channel_ret180/data'
DATA_FILENAME = 'filter_calc_input.mat'
CODE_DIR = '/users/1/kuma0458/open_channel_ret180/codes'
MODEL_FILENAME = 'siren_filter_model.pth'

# ==========================================
# 1. DATA LOADER
# ==========================================
def load_matlab_data():
    filepath = os.path.join(DATA_DIR, DATA_FILENAME)
    print(f"Loading data from: {filepath}")
    
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    with h5py.File(filepath, 'r') as f:
        # Load and transpose variables
        if 'phin' not in f.keys():
             raise KeyError(f"'phin' not found. Keys: {list(f.keys())}")

        # Transpose to fix MATLAB -> Python dimension flip
        phin = np.array(f['phin']).transpose()
        W    = np.array(f['W']).transpose()
        Kxn  = np.array(f['Kxn']).transpose()
        Kyn  = np.array(f['Kyn']).transpose()
        Zn   = np.array(f['Zn']).transpose()
        
    print(f"Data Shapes - Phi: {phin.shape}, W: {W.shape}")

    # Flatten for MLP input
    coords_flat = np.stack([Kxn.flatten(), Kyn.flatten(), Zn.flatten()], axis=1)
    
    # Target: 0 for Pos, 1 for Neg
    targets_flat = (phin.flatten() < 0).astype(np.longlong) 
    
    # Weights: Absolute normalized flux
    weights_flat = W.flatten()
    
    return (torch.tensor(coords_flat, dtype=torch.float32),
            torch.tensor(targets_flat, dtype=torch.long),
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

    # 1. Load Data (CPU)
    coords, targets, weights = load_matlab_data()

    # 2. APPLY SOFT AMPLITUDE WEIGHTING WITH FLOOR
    print("applying Soft Magnitude Weighting (with background suppression)...")
    
    # A. Calculate Magnitude (Absolute Flux)
    # This captures the natural fade: Bright Red/Blue -> Dim -> White
    mag = torch.abs(weights)
    
    # B. Normalize to [0, 1]
    # So the brightest red/blue spot has weight 1.0
    if mag.max() > 0:
        weights = mag / mag.max()
    
    # C. ADD THE FLOOR (Crucial for Clean White Background)
    # We add 0.05 to EVERYTHING.
    # - In the white region (flux~0), weight becomes 0.05. Target is 0. 
    #   The model is forced to learn "Silence" (White) here.
    # - In the red/blue lobes (flux~1), weight becomes 1.05.
    #   The structure remains the priority, but edges fade smoothly.
    background_floor = 0.05
    weights = weights + background_floor
    
    print(f"Weights scaled. Floor: {background_floor}. Range: [{weights.min():.4f}, {weights.max():.4f}]")
    # 3. MOVE EVERYTHING TO GPU NOW (The Speed Fix)
    print("Moving entire dataset to GPU VRAM... (This takes a few seconds)")
    coords = coords.to(device)
    targets = targets.to(device)
    weights = weights.to(device)
    print("Data successfully loaded on GPU.")

    # Initialize Model
    model = DualFilterSiren().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)

    # Manual Batching Setup
    N = coords.shape[0]
    batch_size = 65536
    num_batches = int(np.ceil(N / batch_size))

    epochs = 100
    print(f"Starting High-Speed Training for {epochs} epochs...")

    model.train() # Set mode

    for epoch in range(epochs):
        epoch_loss = 0

        # Shuffle indices directly on GPU
        indices = torch.randperm(N, device=device)

        # --- Manual Loop (No DataLoader Overhead) ---
        for i in range(num_batches):
            start = i * batch_size
            end = min(start + batch_size, N)
            idx = indices[start:end]

            # Data is ALREADY on GPU - Zero transfer time
            batch_coords = coords[idx]
            batch_targets = targets[idx]
            batch_weights = weights[idx]

            optimizer.zero_grad()
            logits = model(batch_coords)

            raw_loss = F.cross_entropy(logits, batch_targets, reduction='none')
            loss = (raw_loss * batch_weights).mean()

            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()

        # Print progress
        avg_loss = epoch_loss / num_batches
        print(f"Epoch {epoch} | Loss: {avg_loss:.6f}", flush=True)

    # Save
    save_path = os.path.join(CODE_DIR, MODEL_FILENAME)
    torch.save(model.state_dict(), save_path)
    print(f"Model saved to {save_path}")

if __name__ == "__main__":
    train_model()
