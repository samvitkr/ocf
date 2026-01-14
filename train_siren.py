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
            SineLayer(3, 256, is_first=True, omega_0=30),
            SineLayer(256, 256, is_first=False, omega_0=30),
            SineLayer(256, 256, is_first=False, omega_0=30),
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

    # Load Data
    coords, targets, weights = load_matlab_data()
    
    # --- FIX 2: Scale weights so gradients aren't tiny ---
    # If the mean weight is 0.0001, the network learns very slowly.
    # We boost them so the average weight is closer to 1.0.
    mean_weight = weights.mean()
    if mean_weight > 0:
        print(f"Original Mean Weight: {mean_weight:.6f}. Scaling weights by {1.0/mean_weight:.2f}...")
        weights = weights / mean_weight
    
    dataset = TensorDataset(coords, targets, weights)
    # num_workers=0 to prevent deadlocks
    dataloader = DataLoader(dataset, batch_size=65536, shuffle=True, num_workers=0)

    # Initialize Model
    model = DualFilterSiren().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4) # Standard SIREN LR
    
    epochs = 2000
    print(f"Starting training for {epochs} epochs...")

    for epoch in range(epochs):
        epoch_loss = 0
        count = 0
        
        for batch_coords, batch_targets, batch_weights in dataloader:
            batch_coords = batch_coords.to(device)
            batch_targets = batch_targets.to(device)
            batch_weights = batch_weights.to(device)

            optimizer.zero_grad()
            logits = model(batch_coords)
            
            # Weighted Loss
            raw_loss = F.cross_entropy(logits, batch_targets, reduction='none')
            loss = (raw_loss * batch_weights).mean()
            
            loss.backward()
            optimizer.step()
            
            epoch_loss += loss.item()
            count += 1
            
        # --- FIX 1: DEDENT THIS BLOCK (Move it LEFT) ---
        # It must align with 'for batch...', NOT be inside it.
        #if epoch % 1 == 0:
            avg_loss = epoch_loss / count
            print(f"Epoch {epoch} | Loss: {avg_loss:.6f}", flush=True)

    # Save Model
    save_path = os.path.join(CODE_DIR, MODEL_FILENAME)
    torch.save(model.state_dict(), save_path)
    print(f"Model saved to {save_path}")
if __name__ == "__main__":
    train_model()
