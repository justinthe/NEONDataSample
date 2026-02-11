import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

class SpectralExplorer:
    def __init__(self, h5_path):
        self.h5_path = h5_path
        if not os.path.exists(h5_path):
            raise FileNotFoundError(f"Could not find H5 file at {h5_path}")
            
    def get_pixel_signature(self, row, col):
        """Extracts and cleans a single pixel's spectral curve."""
        with h5py.File(self.h5_path, 'r') as f:
            site = list(f.keys())[0]
            refl = f[site]['Reflectance']['Reflectance_Data']
            wavelengths = f[site]['Reflectance']['Metadata']['Spectral_Data']['Wavelength'][()]
            
            # Extract curve and scale (NEON data is stored as integers * 10000)
            signature = refl[row, col, :].astype(float) / 10000.0
            
            # Mask out NoData and atmospheric absorption windows (water vapor)
            signature[signature < 0] = np.nan
            # Common water vapor absorption bands to hide for cleaner plots
            signature[190:210] = np.nan # ~1350-1450nm
            signature[280:315] = np.nan # ~1800-2000nm
            
            return wavelengths, signature

    def plot_comparison(self, pixels_dict):
        """
        Plots multiple signatures for comparison.
        pixels_dict format: {"Label": (row, col), ...}
        """
        plt.figure(figsize=(12, 7))
        
        for label, coords in pixels_dict.items():
            row, col = coords
            waves, sig = self.get_pixel_signature(row, col)
            plt.plot(waves, sig, label=label, linewidth=2)

        # Scientific Plot Styling
        plt.title(f"BART Site: Hyperspectral Species Comparison", fontsize=14)
        plt.xlabel("Wavelength (nm)", fontsize=12)
        plt.ylabel("Reflectance (%)", fontsize=12)
        
        # Highlight Red Edge and NIR Plateau
        plt.axvspan(680, 750, color='gray', alpha=0.1, label='Red Edge Transition')
        
        plt.legend()
        plt.grid(True, linestyle=':', alpha=0.7)
        plt.xlim(400, 2500)
        plt.ylim(0, 0.6) # Standard reflectance range 0-60%
        
        plt.tight_layout()
        plt.show()

# --- Execution ---
if __name__ == "__main__":
    # Update this path to your BART H5 file
    target_h5 = "neon_data_BART/HS/NEON_D01_BART_DP3_315000_4878000_reflectance.h5"
    
    explorer = SpectralExplorer(target_h5)
    
    # Let's define some interesting points in the BART forest
    # (These coordinates are samples; try changing them to explore your specific tile!)
    comparison_points = {
        "Sugar Maple Stand": (520, 480),
        "Conifer Pocket": (210, 815),
        "Bare Soil/Path": (100, 100)
    }
    
    explorer.plot_comparison(comparison_points)