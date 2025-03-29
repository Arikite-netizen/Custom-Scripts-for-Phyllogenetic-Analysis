import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Load network data
network_file = "cox1_combined.nexus_network.txt"

# Read network file
edges = []
mutation_steps = {}  # Store mutation steps per species

with open(network_file, "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 3:
            species1, mutations, species2 = parts
            mutations = int(mutations)
            edges.append((species1, species2, mutations))

            # Track mutation steps for coloring nodes
            mutation_steps[species1] = mutation_steps.get(species1, 0) + mutations
            mutation_steps[species2] = mutation_steps.get(species2, 0) + mutations

# Create a network graph
G = nx.Graph()
for edge in edges:
    G.add_edge(edge[0], edge[1], weight=edge[2])

# Identify high mutation rate species (Top 15%)
sorted_mutations = sorted(mutation_steps.values(), reverse=True)
high_mutation_threshold = sorted_mutations[int(len(sorted_mutations) * 0.15)]  # Top 15%

# Identify major connectors (nodes with high degree)
degree_centrality = nx.degree_centrality(G)
high_degree_threshold = np.percentile(list(degree_centrality.values()), 85)  # Top 15% connectors

# Assign 3D positions for layout
pos = nx.spring_layout(G, dim=3, seed=42)

# Extract node positions & properties
node_x, node_y, node_z = [], [], []
node_colors = []
node_sizes = []
node_labels = {}

mutation_values = list(mutation_steps.values())
min_mutation = min(mutation_values)
max_mutation = max(mutation_values)

# Fix: Use new Matplotlib colormap function
cmap = plt.get_cmap("coolwarm")  # Updated fix
plotly_colorscale = [[i / 100, mcolors.rgb2hex(cmap(i / 100))] for i in range(101)]

# Define key species to highlight
highlight_species = {"Ablennes hians", "Harpadon nehereus", "Hydrophis schistosus"}

for node in G.nodes():
    x, y, z = pos[node]
    node_x.append(x)
    node_y.append(y)
    node_z.append(z)

    # Normalize mutation steps to color scale
    norm_mutation = (mutation_steps[node] - min_mutation) / (max_mutation - min_mutation)
    node_colors.append(norm_mutation)  # This will be mapped to the colorscale

    # Determine node size:
    if mutation_steps[node] >= high_mutation_threshold or degree_centrality[node] >= high_degree_threshold:
        node_sizes.append(12)  # Larger size for high mutation & major connectors
        if node in highlight_species or mutation_steps[node] >= high_mutation_threshold:
            node_labels[node] = node  # Label key nodes
    else:
        node_sizes.append(5)  # Default small size

# Create edge traces (gene flow visualization)
edge_x, edge_y, edge_z = [], [], []
edge_width = []  # Fix: Use a single float value for width

for edge in G.edges():
    x0, y0, z0 = pos[edge[0]]
    x1, y1, z1 = pos[edge[1]]
    edge_x.extend([x0, x1, None])  # None for line breaks
    edge_y.extend([y0, y1, None])
    edge_z.extend([z0, z1, None])

    mutation_distance = G[edge[0]][edge[1]]['weight']
    edge_width.append(max(1, 5 - (mutation_distance / 10)))  # Fix: Use a single width value

# Fix: Assign a single width value for edges
edge_trace = go.Scatter3d(
    x=edge_x, y=edge_y, z=edge_z,
    line=dict(width=2, color='gray'),  # Fixed width instead of list
    hoverinfo='none',
    mode='lines'
)

# Create 3D node trace
node_trace = go.Scatter3d(
    x=node_x, y=node_y, z=node_z,
    mode='markers+text',
    marker=dict(
        size=node_sizes,
        cmin=0,
        cmax=1,
        color=node_colors,
        opacity=0.9,
        line=dict(width=1, color="black"),
        colorscale=plotly_colorscale,  # Fixed colorscale issue
        showscale=True,
        colorbar=dict(title="Mutation Steps (Low → High)")
    ),
    text=[node_labels.get(node, "") for node in G.nodes()],
    hoverinfo="text"
)

# Create figure
fig = go.Figure(data=[edge_trace, node_trace])

# Update layout for better 3D view
fig.update_layout(
    title="3D Gene Flow Network with High Mutation Nodes",
    margin=dict(l=0, r=0, b=0, t=40),
    showlegend=False,
    scene=dict(
        xaxis=dict(title="X-axis"),
        yaxis=dict(title="Y-axis"),
        zaxis=dict(title="Z-axis"),
    )
)

# Show interactive plot
fig.show()
