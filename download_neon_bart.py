import requests
import os
from tqdm import tqdm

def download_file(url, folder, filename):
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, filename)
    
    # Check if file exists and has correct size
    try:
        head = requests.head(url)
        total_size = int(head.headers.get('content-length', 0))
    except:
        total_size = 0
    
    if os.path.exists(path) and os.path.getsize(path) == total_size:
        print(f"  [Skipping] {filename}")
        return

    print(f"  [Target] {filename} ({total_size / (1024*1024):.2f} MB)")
    response = requests.get(url, stream=True)
    with open(path, 'wb') as f, tqdm(total=total_size, unit='iB', unit_scale=True) as bar:
        for chunk in response.iter_content(chunk_size=8192):
            bar.update(f.write(chunk))

def get_bart_comprehensive_data_2017():
    """
    Downloads all necessary NEON AOP products for 44 Forest Metrics at BART (August 2017):
    - RGB: DP3.30010.001 (Basic Inventory Metric 6) [cite: 13]
    - CHM: DP3.30015.001 (Basic Inventory Metric 3) 
    - DTM: DP3.30024.001 (Structure & Ecology Metric 32) 
    - HS:  DP3.30006.001 (Advanced Health Metric 11) [cite: 23]
    """
    site = "BART"
    year_month = "2017-08" # Switched to August to ensure tile availability
    tile_coords = "315000_4878000"
    
    products = {
        "RGB": "DP3.30010.001",    
        "LiDAR_CHM": "DP3.30015.001", 
        "LiDAR_DTM": "DP3.30024.001", 
        "HS": "DP3.30006.001"      
    }

    base_dir = "neon_data_BART_2017"
    
    print(f"--- Searching NEON Data for BART {year_month} ---")
    for name, p_id in products.items():
        api_url = f"https://data.neonscience.org/api/v0/data/{p_id}/{site}/{year_month}"
        res = requests.get(api_url)
        
        if res.status_code == 200:
            files = res.json()['data']['files']
            # Search for the 1km tile coordinates
            target_files = [f for f in files if tile_coords in f['name'] and f['name'].lower().endswith(('.h5', '.tif'))]
            
            if not target_files:
                print(f"  [!] No matching tile {tile_coords} found for {name} in {year_month}.")
                continue

            for f_meta in target_files:
                subfolder = "LiDAR" if "LiDAR" in name else name
                download_file(f_meta['url'], os.path.join(base_dir, subfolder), f_meta['name'])
        else:
            print(f"  [!] No data found for {name} ({p_id}) at {site} in {year_month}.")

if __name__ == "__main__":
    get_bart_comprehensive_data_2017()