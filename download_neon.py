import requests
import os
from tqdm import tqdm

class NEONDownloader:
    def __init__(self, site, year_month, tile_coords):
        self.site = site
        self.year_month = year_month
        self.tile_coords = tile_coords
        self.year = year_month.split("-")[0]
        self.base_dir = f"neon_data_{site}_{self.year}"
        
        # Product IDs for all 44 metrics 
        self.products = {
            "RGB": "DP3.30010.001",    # Visual/Species Identification [cite: 13]
            "LiDAR_CHM": "DP3.30015.001", # Height/Volume [cite: 7, 18]
            "LiDAR_DTM": "DP3.30024.001", # Topography/Hydrology [cite: 63, 64]
            "HS": "DP3.30006.001"      # Health/NDVI/LAI/Chemistry [cite: 23, 40, 44]
        }

    def download_all(self):
        print(f"--- Fetching Data for {self.site} ({self.year_month}) ---")
        for name, p_id in self.products.items():
            api_url = f"https://data.neonscience.org/api/v0/data/{p_id}/{self.site}/{self.year_month}"
            res = requests.get(api_url)
            
            if res.status_code == 200:
                files = res.json()['data']['files']
                target_files = [f for f in files if self.tile_coords in f['name'] and f['name'].lower().endswith(('.h5', '.tif'))]
                
                if not target_files:
                    print(f"  [!] No matching tile {self.tile_coords} for {name}")
                    continue

                for f_meta in target_files:
                    subfolder = "LiDAR" if "LiDAR" in name else name
                    dest_folder = os.path.join(self.base_dir, subfolder)
                    self._request_download(f_meta['url'], dest_folder, f_meta['name'])
            else:
                print(f"  [!] Site/Date unavailable for {name}")

    def _request_download(self, url, folder, filename):
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, filename)
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        
        with open(path, 'wb') as f, tqdm(total=total_size, unit='iB', unit_scale=True) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                bar.update(f.write(chunk))

if __name__ == "__main__":
    # CHANGE THESE PARAMETERS HERE
    downloader = NEONDownloader(site="BART", year_month="2017-08", tile_coords="315000_4878000")
    downloader.download_all()