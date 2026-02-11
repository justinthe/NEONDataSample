import os
import h5py
import numpy as np
import rasterio
from rasterio.transform import from_origin
from glob import glob

class NEONTileProcessor:
    def __init__(self, data_folder):
        """
        Initializes the processor for a specific NEON data directory.
        """
        self.base_dir = data_folder
        self.output_dir = os.path.join(data_folder, "qgis_ready")
        os.makedirs(self.output_dir, exist_ok=True)

    def _get_h5_metadata(self, h5_file_path):
        """Extracts spatial metadata (UTM coordinates and CRS) from NEON H5."""
        with h5py.File(h5_file_path, 'r') as f:
            site = list(f.keys())[0]
            coord_sys = f[site]['Reflectance']['Metadata']['Coordinate_System']
            map_info = coord_sys['Map_Info'][()].decode('utf-8').split(',')
            
            # Index 3 is Easting (X), Index 4 is Northing (Y)
            x_min = float(map_info[3])
            y_max = float(map_info[4])
            
            crs_str = "EPSG:32611" # Default fallback for SJER
            if 'Coordinate_System_String' in coord_sys:
                crs_str = coord_sys['Coordinate_System_String'][()].decode('utf-8')
            
            return x_min, y_max, crs_str

    def _save_as_tif(self, data, output_name, x_min, y_max, crs_str):
        """Saves 2D arrays with type-specific NoData handling."""
        transform = from_origin(x_min, y_max, 1, 1)
        output_path = os.path.join(self.output_dir, output_name)
        
        # Ensure 2D
        if data.ndim > 2:
            data = data[:, :, 0] if data.shape[2] == 1 else data[0, :, :]

        dtype = data.dtype
        data_to_save = np.copy(data)

        # Robust NoData Assignment based on Data Type
        if np.issubdtype(dtype, np.floating):
            # For NDVI, Biomass, etc.
            nodata_val = -9999.0
            data_to_save[np.isnan(data_to_save)] = nodata_val
            data_to_save[np.isinf(data_to_save)] = nodata_val
        elif dtype == np.uint8:
            # For Masks (Shadow, Haze) - 0 and 1 are valid data, use 255 for empty
            nodata_val = 255
        else:
            # For Reflectance bands (int16)
            nodata_val = -9999

        with rasterio.open(
            output_path, 'w',
            driver='GTiff',
            height=data_to_save.shape[0],
            width=data_to_save.shape[1],
            count=1,
            dtype=dtype,
            crs=crs_str,
            transform=transform,
            nodata=nodata_val
        ) as dst:
            dst.write(data_to_save, 1)
        print(f"  [Created] {output_name} (Type: {dtype})")

    def process_h5(self):
        """Recursively crawls H5 for ancillary masks and reference RGB bands."""
        h5_files = glob(os.path.join(self.base_dir, "HS", "*.h5"))
        for h5_path in h5_files:
            print(f"Extracting H5 Layers: {os.path.basename(h5_path)}")
            x_min, y_max, crs_str = self._get_h5_metadata(h5_path)
            
            with h5py.File(h5_path, 'r') as f:
                site = list(f.keys())[0]
                layers_to_find = ['Cast_Shadow', 'Data_Selection_Index', 'Haze_Cloud_Water_Map']
                
                def scanner(name, obj):
                    if isinstance(obj, h5py.Dataset):
                        basename = name.split('/')[-1]
                        if basename in layers_to_find:
                            self._save_as_tif(obj[()], f"Layer_{basename}.tif", x_min, y_max, crs_str)
                
                f[site].visititems(scanner)
                
                # Extract RGB Bands for reference
                refl = f[site]['Reflectance']['Reflectance_Data']
                for name, idx in [("Red_B58", 58), ("Green_B34", 34), ("Blue_B19", 19)]:
                    self._save_as_tif(refl[:, :, idx], f"Band_{name}.tif", x_min, y_max, crs_str)

    def calculate_ndvi(self):
        """Calculates quantitative NDVI from the spectral cube."""
        h5_files = glob(os.path.join(self.base_dir, "HS", "*.h5"))
        for h5_path in h5_files:
            print(f"Calculating NDVI: {os.path.basename(h5_path)}")
            x_min, y_max, crs_str = self._get_h5_metadata(h5_path)
            
            with h5py.File(h5_path, 'r') as f:
                site = list(f.keys())[0]
                refl = f[site]['Reflectance']['Reflectance_Data']
                nir = refl[:, :, 90].astype(np.float32)
                red = refl[:, :, 58].astype(np.float32)
                
                # Clean NoData
                nir[nir == -9999] = np.nan
                red[red == -9999] = np.nan
                
                with np.errstate(divide='ignore', invalid='ignore'):
                    ndvi = (nir - red) / (nir + red)
                
                self._save_as_tif(ndvi.astype(np.float32), "Product_NDVI.tif", x_min, y_max, crs_str)

    def calculate_cleaned_tree_health(self, height_threshold=5.0): # Set to 5m for BART
        """Fuses LiDAR and Shadow data to isolate sunlit tree health."""
        print(f"Cleaning Tree Health (Height > {height_threshold}m)...")
        ndvi_path = os.path.join(self.output_dir, "Product_NDVI.tif")
        shadow_path = os.path.join(self.output_dir, "Layer_Cast_Shadow.tif")
        chm_files = glob(os.path.join(self.base_dir, "LiDAR", "*.tif"))
        
        with rasterio.open(ndvi_path) as ndvi_src, \
             rasterio.open(shadow_path) as shadow_src, \
             rasterio.open(chm_files[0]) as chm_src:
            
            ndvi = ndvi_src.read(1).astype(np.float32)
            shadow = shadow_src.read(1)
            chm = chm_src.read(1)
            
            # 1. Perform Masking
            # In BART, shadow==0 is sun, chm > 5.0m is forest canopy
            is_valid_tree = (chm >= height_threshold) & (shadow == 0)
            
            # 2. FORCE float32 and explicit NaN placement
            cleaned_ndvi = np.full(ndvi.shape, np.nan, dtype=np.float32)
            cleaned_ndvi[is_valid_tree] = ndvi[is_valid_tree]
            
            # 3. Final safety check: ensure no Inf values escaped
            cleaned_ndvi[np.isinf(cleaned_ndvi)] = np.nan
            
            # 4. Save using our standardized helper
            self._save_as_tif(
                cleaned_ndvi, 
                "Product_Tree_Health_Sunlit.tif", 
                ndvi_src.transform.c, 
                ndvi_src.transform.f, 
                ndvi_src.crs
            )

    def calculate_biomass_proxy(self):
        """Calculates Biomass Proxy = Height * Cleaned NDVI."""
        print("Calculating Tree Biomass Proxy...")
        # health_path = os.path.join(self.output_dir, "Product_Tree_Health_Sunlit.tif")
        health_path = os.path.join(self.output_dir, "Product_NDVI.tif")
        chm_files = glob(os.path.join(self.base_dir, "LiDAR", "*.tif"))
        
        print(f"Checking NDVI path: {health_path}")
        print(f"Searching for LiDAR in: {chm_files}")
        
        
        with rasterio.open(health_path) as ndvi_src, rasterio.open(chm_files[0]) as chm_src:
            ndvi = ndvi_src.read(1)
            chm = chm_src.read(1)
            biomass = np.where(np.isnan(ndvi), np.nan, chm * ndvi)
            biomass = np.where(biomass < 0, 0, biomass)
            
            self._save_as_tif(biomass.astype(np.float32), "Product_Biomass_Proxy.tif", 
                              ndvi_src.transform.c, ndvi_src.transform.f, ndvi_src.crs)

    def run_all(self):
        """Executes the full scientific pipeline."""
        print("\n--- NEON EXPERT PIPELINE STARTING ---")
        self.process_h5()      
        self.calculate_ndvi()  
        self.calculate_cleaned_tree_health()
        self.calculate_biomass_proxy()
        print("\n--- PIPELINE COMPLETE: Load 'qgis_ready' folder in QGIS ---")

if __name__ == "__main__":
    processor = NEONTileProcessor(data_folder="neon_data_BART")
    processor.run_all()