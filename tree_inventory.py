import numpy as np
import rasterio
import geopandas as gpd
from shapely.geometry import Point
from scipy.ndimage import maximum_filter, label
import os

class TreeInventoryTool:
    def __init__(self, output_dir="qgis_ready"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def extract_inventory(self, chm_path, ndvi_path, biomass_path):
        """Extracts individual tree metrics and saves to Shapefile."""
        print("Starting Individual Tree Extraction...")
        
        with rasterio.open(chm_path) as chm_src, \
             rasterio.open(ndvi_path) as ndvi_src:
            
            chm = chm_src.read(1)
            ndvi = ndvi_src.read(1)
            transform = chm_src.transform
            
            # 1. Tree Counting & Detection (Local Maxima Filter)
            # We look for pixels that are the highest in a 3x3 meter neighborhood
            neighborhood_size = 3 
            local_max = maximum_filter(chm, size=neighborhood_size) == chm
            # Filter out ground/grass (Height < 5m)
            local_max[chm < 5.0] = False
            
            # Label each unique tree
            labeled_trees, count = label(local_max)
            print(f"  [Found] {count} individual trees in tile.")

            tree_data = []
            
            # 2. Loop through detections to extract metrics
            rows, cols = np.where(local_max)
            for r, c in zip(rows, cols):
                # Spatial Coordinates
                x, y = transform * (c, r)
                
                # 3. Tree Height (Direct from CHM)
                height = float(chm[r, c])
                
                # 5. Tree Health (Direct from NDVI)
                health = float(ndvi[r, c])
                
                # 4. Tree Diameter (DBH) - Estimated via Allometric Scaling
                # For BART hardwoods, a common proxy is DBH (cm) = 1.5 * Height (m)
                dbh = height * 1.5 
                
                # 2. Species Identification (Simplified Spectral Proxy)
                # Maples/Deciduous usually have NDVI > 0.75, Conifers < 0.75
                species = "Deciduous" if health > 0.75 else "Coniferous"

                tree_data.append({
                    'geometry': Point(x, y),
                    'Height_m': round(height, 2),
                    'NDVI_Health': round(health, 3),
                    'DBH_cm_est': round(dbh, 2),
                    'Species': species
                })

            # Create GeoDataFrame and export to SHP
            gdf = gpd.GeoDataFrame(tree_data, crs=chm_src.crs)
            shp_path = os.path.join(self.output_dir, "BART_Tree_Inventory.shp")
            gdf.to_file(shp_path)
            
            print(f"--- Inventory Complete! ---")
            print(f"Shapefile saved to: {shp_path}")
            return gdf

# --- Execution ---
if __name__ == "__main__":
    # Paths to the products we created earlier
    inventory = TreeInventoryTool()
    
    # Update these paths to your BART files
    inventory.extract_inventory(
        chm_path="neon_data_BART/LiDAR/NEON_D01_BART_DP3_315000_4878000_CHM.tif",
        ndvi_path="neon_data_BART/qgis_ready/Product_NDVI.tif",
        biomass_path="neon_data_BART/qgis_ready/Product_Biomass_Proxy_Forgiving.tif"
    )