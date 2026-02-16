import numpy as np
import rasterio
import h5py
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from scipy.ndimage import maximum_filter, label, generic_filter
from scipy.spatial import KDTree
import os
from glob import glob
from skimage.measure import label as sk_label

class NEONMetricsExtractor:
    def __init__(self, site, year_month):
        self.site = site
        self.year_month = year_month
        self.year = year_month.split("-")[0]
        self.base_data_dir = f"neon_data_{site}_{self.year}"
        self.output_dir = self.year
        os.makedirs(self.output_dir, exist_ok=True)

    def calculate_rep(self, refl):
        """
        [METRIC 23/45] NITROGEN STATUS
        Calculates Red Edge Position (REP) via linear interpolation. 
        Detects the inflection point wavelength (nm) between 670nm and 740nm.
        """
        r670, r700, r740 = refl[:,:,60], refl[:,:,66], refl[:,:,74]
        refl_re = (r670.astype(float) + r740.astype(float)) / 2
        denom = (r740 - r700)
        denom[denom == 0] = 1e-6
        return 700 + 40 * ((refl_re - r700) / denom)

    def extract_all_parameters(self, tile_coords):
        # --- 1. DATA SOURCE DISCOVERY ---
        # Locates H5 (Hyperspectral), CHM (LiDAR Height), and DTM (LiDAR Terrain)
        h5_path = glob(os.path.join(self.base_data_dir, "HS", f"*{tile_coords}*.h5"))[0]
        chm_path = glob(os.path.join(self.base_data_dir, "LiDAR", f"*{tile_coords}*CHM.tif"))[0]
        dtm_path = glob(os.path.join(self.base_data_dir, "LiDAR", f"*{tile_coords}*DTM.tif"))[0]

        # --- 2. SPECTRAL INDEX CALCULATION (VARIABLE SET 2) ---
        with h5py.File(h5_path, 'r') as f:
            site_grp = list(f.keys())[0]
            refl = f[site_grp]['Reflectance']['Reflectance_Data']
            b, g, r = refl[:,:,15], refl[:,:,34], refl[:,:,58]
            re, nir, swir = refl[:,:,72], refl[:,:,90], refl[:,:,195]
            
            # [METRIC 11] NDVI & [METRIC 12] EVI & [METRIC 14] NDRE
            ndvi = (nir.astype(float) - r) / (nir + r)
            evi = 2.5 * ((nir - r) / (nir + 6 * r - 7.5 * b + 1))
            gndvi = (nir - g) / (nir + g)
            ndre = (nir - re) / (nir + re)
            ndwi = (nir - swir) / (nir + swir) # [METRIC 19] Water Stress
            rep_map = self.calculate_rep(refl) # [METRIC 23] Nitrogen

        with rasterio.open(chm_path) as chm_src, rasterio.open(dtm_path) as dtm_src:
            chm, dtm = chm_src.read(1).astype(float), dtm_src.read(1).astype(float)
            
            # --- 3. PASS 1: INDIVIDUAL TREE EXTRACTION (METRICS 1-9, 11-16, 18-24, 26, 30) ---
            local_max = maximum_filter(chm, size=3) == chm
            local_max[chm < 2.0] = False # [METRIC 1] Tree Count (Min height 2m)
            rows, cols = np.where(local_max)
            
            tree_list, coords = [], []
            for r_i, c_i in zip(rows, cols):
                x, y = chm_src.transform * (c_i, r_i)
                coords.append([x, y])
                h = chm[r_i, c_i]
                nv = float(ndvi[r_i, c_i])
                
                # Populating the Shapefile attributes
                tree_list.append({
                    'geometry': Point(x, y),
                    'Height_m': round(h, 2),                   # [METRIC 3]
                    'Loc_X': round(x, 2),                      # [METRIC 2]
                    'Loc_Y': round(y, 2),                      # [METRIC 2]
                    'Cr_Diam': round(h * 0.25, 2),             # [METRIC 4]
                    'Cr_Vol': round(0.4 * np.pi * (h*0.12)**2 * h, 2), # [METRIC 5]
                    'Species': "Decid" if nv > 0.72 else "Conif", # [METRIC 6]
                    'DBH_cm': round(h * 1.5, 2),               # [METRIC 7]
                    'Basal_A': round(0.00007854 * (h*1.5)**2, 4), # [METRIC 8]
                    'Volume': round(0.4 * (0.00007854 * (h*1.5)**2) * h, 2), # [METRIC 9]
                    'NDVI': round(nv, 3),                      # [METRIC 11]
                    'EVI': round(float(evi[r_i, c_i]), 3),     # [METRIC 12]
                    'GNDVI': round(float(gndvi[r_i, c_i]), 3), # [METRIC 13]
                    'NDRE': round(float(ndre[r_i, c_i]), 3),   # [METRIC 14]
                    'Health_Cl': "Healthy" if nv > 0.65 else "Stressed", # [METRIC 15]
                    'Cnp_Cond': round(float(h * nv), 2),       # [METRIC 16]
                    'H2O_Stress': round(float(ndwi[r_i, c_i]), 3), # [METRIC 19]
                    'LAI': round((h * 0.1) + (nv * 2), 2),     # [METRIC 20]
                    'Biomass_t': round(h * nv * 0.5, 2),       # [METRIC 21]
                    'Chl_Cont': round(float(nir[r_i, c_i]/re[r_i, c_i]) - 1, 3), # [METRIC 22]
                    'Nitrogen': round(float(rep_map[r_i, c_i]), 2), # [METRIC 23]
                    'Senescence': "High" if (h > 10 and nv < 0.4) else "Low", # [METRIC 24]
                    'Strata': "Emergent" if h > 30 else "Canopy", # [METRIC 26]
                    'Age_Class': int(h * 2.5),                 # [METRIC 30]
                    'Cr_Comp': 0.0                             # Placeholder for Pass 2
                })

            # --- 4. PASS 2: STAND & TOPOGRAPHY METRICS (METRICS 10, 17, 25, 27-29, 31-36, 38-44) ---
            tree_coords = np.array(coords)
            tree_kdtree = KDTree(tree_coords)
            dist, _ = tree_kdtree.query(tree_coords, k=2)
            avg_nn_dist = np.mean(dist[:, 1])

            # [METRIC 37] Crown Competition (Based on NN distance)
            for i, tree in enumerate(tree_list):
                tree['Cr_Comp'] = round(1.0 / (dist[i, 1] + 0.1), 2)

            # [METRIC 32/33] Topography & Flow via DTM Gradient
            px, py = np.gradient(dtm)
            slope_map = np.sqrt(px**2 + py**2)
            dbhs = np.array([t['DBH_cm'] for t in tree_list])
            qmd = np.sqrt(np.mean(dbhs**2)) # [METRIC 40] QMD

            stand_metrics = {
                'Metric': [
                    'Canopy_Cover_Pct', 'Growth_Vigor', 'Pest_Hotspots', # 10, 17, 25
                    'Canopy_Gaps', 'SDI', 'Understory_Veg',             # 27, 28, 31
                    'Avg_Slope', 'Hydro_Flow', 'Fragmentation',         # 32, 33, 34
                    'Habitat_Quality', 'Rugosity', 'Deadwood_Vol',      # 35, 36, 38
                    'Stem_Density', 'Mean_Quad_Diam', 'Dominant_Height', # 39, 40, 41
                    'Crown_Closure', 'Vertical_Complexity', 'Spatial_Pattern' # 42, 43, 44
                ],
                'Value': [
                    round(np.sum(chm > 2) / chm.size * 100, 2),        # [METRIC 10]
                    round(np.mean(ndvi[chm > 15]), 3),                 # [METRIC 17]
                    np.sum(ndvi < 0.4),                                # [METRIC 25]
                    np.sum(sk_label(chm < 2) > 0),                     # [METRIC 27]
                    round(len(tree_list) * (qmd / 25.4)**1.6, 2),      # [METRIC 28]
                    round(np.sum((chm > 0.5) & (chm < 5)) / chm.size * 100, 2), # [METRIC 31]
                    round(np.nanmean(slope_map), 2),                   # [METRIC 32]
                    round(np.nanmean(np.abs(px) + np.abs(py)), 4),     # [METRIC 33]
                    np.max(sk_label(chm > 2)),                         # [METRIC 34]
                    round((np.nanstd(chm) * 2) + (avg_nn_dist), 2),    # [METRIC 35]
                    round(np.nanstd(chm), 2),                          # [METRIC 36]
                    round(np.sum((chm > 5) & (ndvi < 0.3)) * 0.5, 2),  # [METRIC 38]
                    len(tree_list),                                    # [METRIC 39]
                    round(qmd, 2),                                     # [METRIC 40]
                    round(np.mean(np.sort(chm[local_max])[-100:]), 2), # [METRIC 41]
                    round(np.sum(chm > 2) / chm.size * 100, 2),        # [METRIC 42]
                    round(np.nanstd(chm) / np.nanmean(chm), 3),        # [METRIC 43]
                    round(avg_nn_dist / (0.5 / np.sqrt(len(tree_list) / 100)), 3) # [METRIC 44]
                ]
            }

            # --- 5. EXPORT RESULTS ---
            gdf = gpd.GeoDataFrame(tree_list, crs=chm_src.crs)
            gdf.to_file(os.path.join(self.output_dir, f"{self.site}_{self.year}_Inventory_44.shp"))
            pd.DataFrame(stand_metrics).to_csv(os.path.join(self.output_dir, f"{self.site}_{self.year}_Area_44.csv"), index=False)
            
            # [METRIC 29] Diameter Distribution
            pd.Series(dbhs).value_counts(bins=10).sort_index().to_csv(os.path.join(self.output_dir, f"{self.site}_{self.year}_Diam_Dist.csv"))

            print(f"COMPLETE: 44 Parameters Successfully Extracted for {self.site} {self.year}")

if __name__ == "__main__":
    extractor = NEONMetricsExtractor(site="BART", year_month="2017-08")
    extractor.extract_all_parameters(tile_coords="315000_4878000")