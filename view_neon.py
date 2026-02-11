import os
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from rasterio.plot import show

def view_data():
    base_dir = "neon_data"
    
    # 1. Setup Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # --- Part A: RGB ---
    rgb_file = glob.glob(f"{base_dir}/RGB/*.tif")[0]
    with rasterio.open(rgb_file) as src:
        axes[0].set_title("RGB (0.1m Resolution)")
        show(src, ax=axes[0])

    # --- Part B: LiDAR (Canopy Height Model) ---
    chm_file = glob.glob(f"{base_dir}/CHM/*.tif")[0]
    with rasterio.open(chm_file) as src:
        chm_data = src.read(1)
        axes[1].set_title("LiDAR Canopy Height (1m)")
        # Use terrain-style colormap
        im2 = axes[1].imshow(chm_data, cmap='viridis')
        plt.colorbar(im2, ax=axes[1], label='Height (m)')

    # --- Part C: Hyperspectral (Clipped to False Color) ---
    hs_file = glob.glob(f"{base_dir}/HS/*.h5")[0]
    with h5py.File(hs_file, 'r') as f:
        # Navigate the NEON H5 structure
        site = list(f.keys())[0]
        reflectance = f[site]['Reflectance']
        data = reflectance['Reflectance_Data']
        
        # NEON Hyperspectral has 426 bands. 
        # Let's grab bands for a False Color Infrared (CIR) view:
        # NIR (~860nm), Red (~650nm), Green (~550nm)
        # Indices are roughly: Band 90, Band 58, Band 34
        r_band = data[:, :, 90]
        g_band = data[:, :, 58]
        b_band = data[:, :, 34]
        
        rgb_stack = np.stack([r_band, g_band, b_band], axis=-1)
        # Clean up data (scale 0-10000 to 0-1)
        rgb_stack = rgb_stack.astype(float) / 10000.0
        rgb_stack = np.clip(rgb_stack, 0, 1)
        
        axes[2].set_title("Hyperspectral False Color (1m)")
        axes[2].imshow(rgb_stack)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    view_data()