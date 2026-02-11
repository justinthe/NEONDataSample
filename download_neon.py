import requests
import os
from tqdm import tqdm

def download_file(url, folder, filename):
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, filename)
    
    # Get the file size from headers without downloading yet
    head = requests.head(url)
    total_size = int(head.headers.get('content-length', 0))
    
    if os.path.exists(path):
        if os.path.getsize(path) == total_size:
            print(f"  [Skipping] {filename} (Already exists)")
            return
        else:
            print(f"  [Resume] {filename} (Size mismatch, re-downloading)")

    # Download with progress bar
    print(f"  [Target] {filename} ({total_size / (1024*1024):.2f} MB)")
    
    response = requests.get(url, stream=True)
    
    with open(path, 'wb') as f, tqdm(
        desc=filename,
        total=total_size,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as bar:
        for chunk in response.iter_content(chunk_size=8192):
            size = f.write(chunk)
            bar.update(size)

def get_neon_sample():
    # Site: San Joaquin Experimental Range (SJER)
    # Tile: 257000 Easting, 4112000 Northing (UTM Zone 11N)
    site = "SJER"
    year_month = "2021-03" 
    tile_coords = "257000_4112000"
    
    products = {
        "RGB": "DP3.30010.001",    # Camera Mosaic
        "LiDAR": "DP3.30015.001",  # Canopy Height Model
        "HS": "DP3.30006.001"      # Hyperspectral Reflectance
    }

    base_dir = "neon_data"
    
    print("--- NEON AOP Data Finder ---")
    
    for name, p_id in products.items():
        api_url = f"https://data.neonscience.org/api/v0/data/{p_id}/{site}/{year_month}"
        
        try:
            res = requests.get(api_url)
            res.raise_for_status()
            files = res.json()['data']['files']
            
            # Filter for the specific 1km tile and correct extension
            target_files = [f for f in files if tile_coords in f['name'] and f['name'].endswith(('.h5', '.tif'))]
            
            if not target_files:
                print(f"\n[!] No matching tile found for {name}")
                continue

            for f_meta in target_files:
                download_file(f_meta['url'], os.path.join(base_dir, name), f_meta['name'])
                
        except Exception as e:
            print(f"\n[!] Error fetching {name}: {e}")

if __name__ == "__main__":
    get_neon_sample()
    print("\nAll downloads complete!")