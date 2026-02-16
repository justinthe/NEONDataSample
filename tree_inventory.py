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
            tree_id_counter = 1
            
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
                
                # Approximate bounding box size based on tree diameter (simulating CV algorithm)
                # Convert DBH to approximate crown diameter (rough estimation: crown_diameter ~= DBH * 0.5)
                # Then convert to meters and calculate bbox in coordinate space
                crown_radius_m = (dbh / 100) * 0.5 / 2  # Convert cm to m, apply scaling factor, get radius
                pixel_size = abs(transform[0])  # Get pixel resolution from transform
                bbox_offset = max(crown_radius_m, pixel_size * 2)  # Minimum 2 pixel offset
                
                bbox_xmin = x - bbox_offset
                bbox_ymin = y - bbox_offset
                bbox_xmax = x + bbox_offset
                bbox_ymax = y + bbox_offset

                tree_data.append({
                    'tree_id': tree_id_counter,
                    'site_id': 'BART',
                    'geometry': Point(x, y),
                    'height_m': round(height, 2),
                    'health': round(health, 3),
                    'diameter_cm': round(dbh, 2),
                    'species': species,
                    'bbox_xmin': round(bbox_xmin, 2),
                    'bbox_ymin': round(bbox_ymin, 2),
                    'bbox_xmax': round(bbox_xmax, 2),
                    'bbox_ymax': round(bbox_ymax, 2)
                })
                tree_id_counter += 1

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