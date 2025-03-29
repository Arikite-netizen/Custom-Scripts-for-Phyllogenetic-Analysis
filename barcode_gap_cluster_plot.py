import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.widgets import Cursor

# Define input files
distance_file = "Full_Distance_Matrix.csv"
pca_file = "pca_scores_cox1_with_clusters.csv"
output_csv_file = "Barcode_Gap_Plot_Data.csv"

# Initialize dictionaries
intrak2p_dict = {}  # Intraspecific distances
interk2p_dict = {}  # Interspecific distances
barcode_gap_dict = {}  # Barcode gap dictionary

# Read PCA file with clusters
pca_data = pd.read_csv(pca_file)

# Standardize species names (replace underscores with spaces)
pca_data["Species"] = pca_data["Species"].str.replace("_", " ")

# Extract cluster information
species_clusters = dict(zip(pca_data["Species"], pca_data["Cluster"]))
cluster_names = dict(zip(pca_data["Cluster"], pca_data["Cluster_Name"]))

# Define cluster colors
cluster_colors = {
    "Cluster 1": "red",
    "Cluster 2": "blue",
    "Cluster 3": "green"
}

# Define cluster representatives
representatives = {
    "Cluster 1": "Hydrophis schistosus",
    "Cluster 2": "Harpadon nehereus",
    "Cluster 3": "Ablennes hians"
}

# Read the genetic distance matrix
def read_distance_matrix(file):
    """
    Reads a tab-separated distance matrix and returns a cleaned dictionary.
    """
    global species_list

    # Read the file
    df = pd.read_csv(file, sep="\t", dtype=str)

    # Strip spaces from headers
    df.columns = df.columns.str.strip()
    df.iloc[:, 0] = df.iloc[:, 0].str.strip()

    # Extract species names
    species_list = df.columns[1:].tolist()

    # Convert matrix into a dictionary
    distance_dict = {}
    for i in range(len(df)):
        species_name = df.iloc[i, 0].strip()
        distances = df.iloc[i, 1:].fillna(0).astype(float).tolist()

        distance_dict[species_name] = distances

    return distance_dict

# Compute intraspecific and interspecific distances
def compute_barcode_gap(distance_dict):
    global intrak2p_dict, interk2p_dict, barcode_gap_dict

    for i, species1 in enumerate(species_list):
        for j, species2 in enumerate(species_list):
            dist = distance_dict[species1][j]

            if species1 == species2:  # Intraspecific distance
                if species1 not in intrak2p_dict:
                    intrak2p_dict[species1] = []
                intrak2p_dict[species1].append(dist)

            else:  # Interspecific distance
                if species1 not in interk2p_dict:
                    interk2p_dict[species1] = []
                interk2p_dict[species1].append(dist)

    # Compute average intraspecific and interspecific distances
    for species in species_list:
        intra_avg = sum(intrak2p_dict.get(species, [0])) / max(len(intrak2p_dict.get(species, [])), 1)
        inter_avg = sum(interk2p_dict.get(species, [0])) / max(len(interk2p_dict.get(species, [])), 1)

        barcode_gap_dict[species] = (intra_avg, inter_avg)

# Generate barcode gap plot with clusters
def barcode_gap_plot_fun(barcode_gap_dict):
    """
    Generates a barcode gap plot with clusters, highlights representatives, and adds clear legends.
    """
    global species_clusters, cluster_colors, representatives

    # Filter out zero values
    barcode_gap_dict = {k: v for k, v in barcode_gap_dict.items() if v[0] > 0 and v[1] > 0}

    if not barcode_gap_dict:
        print("Error: No valid barcode gap data found after filtering out zero values!")
        return

    # Extract values
    x_values = []
    y_values = []
    species_labels = []
    colors = []
    cluster_legend_added = set()

    fig, ax = plt.subplots(figsize=(8, 6))

    # Scatter plot with cluster colors
    for species, (intra, inter) in barcode_gap_dict.items():
        cluster_name = cluster_names.get(species_clusters.get(species, -1), "Unknown")
        color = cluster_colors.get(cluster_name, "gray")

        # Add to scatter plot
        sc = ax.scatter(intra, inter, color=color, alpha=0.6, edgecolors="black", label=cluster_name if cluster_name not in cluster_legend_added else "")
        cluster_legend_added.add(cluster_name)

        x_values.append(intra)
        y_values.append(inter)
        species_labels.append(species)

    # Restore the Equal Line (Diagonal)
    max_val = max(max(x_values), max(y_values))
    ax.plot([0, max_val], [0, max_val], linestyle="--", color="black", label="Equal Line")

    # Highlight cluster representatives separately
    for cluster, species in representatives.items():
        if species in barcode_gap_dict:
            intra, inter = barcode_gap_dict[species]
            ax.scatter(intra, inter, color=cluster_colors[cluster], edgecolors="black", s=150, marker="*", label=f"{species} (Representative)")

    # Hover labels
    annot = ax.annotate("", xy=(0,0), xytext=(10,10), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        index = ind["ind"][0]
        annot.xy = sc.get_offsets()[index]
        annot.set_text(species_labels[index])
        annot.set_visible(True)
        fig.canvas.draw_idle()

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)

    # Labels and title
    ax.set_xlabel("Intraspecific Distance")
    ax.set_ylabel("Interspecific Distance")
    ax.set_title("Barcode Gap Analysis with Clusters")
    ax.legend()
    ax.grid()

    plt.savefig("Barcode_Gap_Plot_with_Clusters.png")
    plt.show()

    # Save CSV output
    plot_data = pd.DataFrame({
        "Species": species_labels,
        "Intraspecific_Distance": x_values,
        "Interspecific_Distance": y_values,
        "Cluster": [species_clusters.get(species, -1) for species in species_labels],
        "Cluster_Name": [cluster_names.get(species_clusters.get(species, -1), "Unknown") for species in species_labels]
    })
    plot_data.to_csv(output_csv_file, index=False)
    print(f"✅ Plot data saved to: {output_csv_file}")

# Main function
def main():
    print("Loading genetic distance matrix...")
    distance_dict = read_distance_matrix(distance_file)

    print("Computing barcode gap values...")
    compute_barcode_gap(distance_dict)

    print("Generating barcode gap plot with clusters...")
    barcode_gap_plot_fun(barcode_gap_dict)

if __name__ == "__main__":
    main()
