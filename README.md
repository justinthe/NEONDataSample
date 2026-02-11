### **README.md (Updated)**

# NEON AOP Expert Multi-Sensor Pipeline

This repository provides an end-to-end workflow for processing **NEON (National Ecological Observatory Network)** Airborne Observation Platform data. It specializes in **Sensor Fusion**, combining LiDAR, Hyperspectral, and RGB data to move from raw sensor counts to individual tree-level metrics.

## 🚀 Key Features

### 1. Individual Tree Inventory (ITC)

The pipeline now performs automated **Individual Tree Crown (ITC)** detection. It identifies treetops and extracts a digital inventory including:

* **Tree Count:** Automated counting via local maxima filtering.
* **Species Identification:** Spectral classification (Deciduous vs. Coniferous).
* **Tree Height:** LiDAR-derived vertical structure.
* **Estimated Diameter (DBH):** Allometric scaling based on regional forest metrics.
* **Health (NDVI):** Tree-top physiological vigor.

### 2. Multi-Sensor Biomass Fusion

Calculates a **Biomass Proxy** by integrating structural volume (LiDAR) with spectral activity (Hyperspectral). The pipeline includes a "Forgiving" mode to handle the dense, complex shadows found in closed-canopy forests like New Hampshire's BART site.

### 3. Spectral Explorer

A class-based tool to extract and plot **Spectral Signatures** ( to ). This allows for the visual comparison of "Red Edge" and "NIR Plateau" differences between hardwood maples and softwood hemlocks.

---

## 📂 Output Data Products

All outputs are saved to `qgis_ready/` with full CRS anchoring (**EPSG:32619** for BART):

| File | Type | Description |
| --- | --- | --- |
| **BART_Tree_Inventory.shp** | **Vector** | Points for every tree with height, DBH, and species attributes. |
| **Product_Biomass_Proxy.tif** | **Raster** | Continuous map of forest carbon potential. |
| **Product_Tree_Health_Sunlit.tif** | **Raster** | Pure NDVI signal masked for canopy and sunlit leaves. |
| **Layer_Cast_Shadow.tif** | **Raster** | Binary mask of solar shadows for quality control. |

---

## 🛠️ Installation & Setup

### 1. Environment Requirements

Ensure your `neon_env` is configured with the necessary scientific stack:

```bash
conda activate neon_env
pip install requests numpy rasterio h5py tqdm scipy geopandas shapely matplotlib

```

### 2. Execution Flow

1. **Download:** `python download_neon.py` (configured for BART 2019-08).
2. **Process:** `python neon_expert_tool.py` (generates the raster stack).
3. **Inventory:** `python tree_inventory.py` (generates the Shapefile).
4. **Explore:** `python spectral_explorer.py` (visualize species curves).

---

## 📊 QGIS Visualization Guide

1. **Map Layout:** Use the RGB tile as the base layer.
2. **Tree Inventory:** Drag in the `.shp` file. Use **Categorized Symbology** on the `Species` column.
3. **Labeling:** Enable labels for the `Height_m` attribute to see tree heights in the map view.
4. **Biomass Heatmap:** Apply a `Viridis` or `YlGn` ramp to the Biomass Proxy to identify high-productivity forest stands.

---

### **Next Step**

You now have a complete toolkit! Would you like me to add a **Statistical Summary** script that creates a bar chart showing the distribution of tree heights (e.g., how many trees are , , etc.) for your BART site?