import requests
import os
from tqdm import tqdm

def download_file(url, folder, filename):
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, filename)
    head = requests.head(url)
    total_size = int(head.headers.get('content-length', 0))
    
    if os.path.exists(path) and os.path.getsize(path) == total_size:
        print(f"  [Skipping] {filename}")
        return

    print(f"  [Target] {filename} ({total_size / (1024*1024):.2f} MB)")
    response = requests.get(url, stream=True)
    with open(path, 'wb') as f, tqdm(total=total_size, unit='iB', unit_scale=True) as bar:
        for chunk in response.iter_content(chunk_size=8192):
            bar.update(f.write(chunk))

def get_bart_sample():
    # Site: Bartlett Experimental Forest (BART)
    # Tile: 315000 E, 4878000 N (UTM Zone 19N)
    # Date: August 2019
    site = "BART"
    year_month = "2019-08" 
    tile_coords = "315000_4878000"
    
    products = {
        "RGB": "DP3.30010.001",
        "LiDAR": "DP3.30015.001",
        "HS": "DP3.30006.001"
    }

    base_dir = "neon_data_BART"
    
    print(f"--- Searching NEON Data: {site} {year_month} ---")
    for name, p_id in products.items():
        api_url = f"https://data.neonscience.org/api/v0/data/{p_id}/{site}/{year_month}"
        res = requests.get(api_url)
        if res.status_code == 200:
            files = res.json()['data']['files']
            target_files = [f for f in files if tile_coords in f['name'] and f['name'].endswith(('.h5', '.tif'))]
            for f_meta in target_files:
                download_file(f_meta['url'], os.path.join(base_dir, name), f_meta['name'])
        else:
            print(f"  [!] No data for {name} in {year_month}")

if __name__ == "__main__":
    get_bart_sample()