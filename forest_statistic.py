import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import os

class ForestStatsReporter:
    def __init__(self, shp_path):
        self.shp_path = shp_path
        self.output_dir = os.path.dirname(shp_path)
        self.df = gpd.read_file(shp_path)

    def generate_summary(self):
        """Prints a text summary of the forest inventory."""
        total_trees = len(self.df)
        avg_height = self.df['Height_m'].mean()
        species_counts = self.df['Species'].value_counts()
        
        print("\n--- BART Forest Inventory Summary ---")
        print(f"Total Trees Detected: {total_trees}")
        print(f"Average Tree Height:  {avg_height:.2f} meters")
        print(f"Species Breakdown:\n{species_counts}")
        print("-------------------------------------\n")

    def plot_height_distribution(self):
        """Creates a histogram showing the distribution of tree heights."""
        plt.figure(figsize=(10, 6))
        
        # Plotting the histogram
        n, bins, patches = plt.hist(self.df['Height_m'], bins=15, color='skyblue', edgecolor='black', alpha=0.7)
        
        plt.title('BART Site: Tree Height Distribution (Structure)', fontsize=14)
        plt.xlabel('Height (meters)', fontsize=12)
        plt.ylabel('Number of Trees', fontsize=12)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Save the plot
        plot_path = os.path.join(self.output_dir, "BART_Height_Distribution.png")
        plt.savefig(plot_path)
        print(f"  [Saved] Height Distribution Chart: {plot_path}")
        plt.show()

    def plot_species_biomass(self):
        """Compare average estimated DBH (size) by species."""
        plt.figure(figsize=(8, 6))
        
        # Group by species and calculate mean DBH
        species_stats = self.df.groupby('Species')['DBH_cm_est'].mean()
        species_stats.plot(kind='bar', color=['darkgreen', 'orange'], alpha=0.8)
        
        plt.title('Average Tree Diameter (DBH) by Species Type', fontsize=14)
        plt.ylabel('Mean Estimated DBH (cm)', fontsize=12)
        plt.xlabel('Species Group', fontsize=12)
        plt.xticks(rotation=0)
        
        plot_path = os.path.join(self.output_dir, "BART_Species_Comparison.png")
        plt.savefig(plot_path)
        print(f"  [Saved] Species Comparison Chart: {plot_path}")
        plt.show()

if __name__ == "__main__":
    shp_file = "qgis_ready/BART_Tree_Inventory.shp"
    reporter = ForestStatsReporter(shp_file)
    reporter.generate_summary()
    reporter.plot_height_distribution()
    reporter.plot_species_biomass()