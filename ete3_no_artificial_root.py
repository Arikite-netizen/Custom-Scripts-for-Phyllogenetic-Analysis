from ete3 import Tree, TreeStyle, NodeStyle
import matplotlib.pyplot as plt

def load_trees(tree_files):
    """Load trees from files and return a list of Tree objects."""
    trees = []
    for file in tree_files:
        with open(file, "r") as f:
            trees.append(Tree(f.readline().strip()))
    return trees

def color_branches(tree, color):
    """Apply a specific color to the branches of a tree."""
    for node in tree.traverse():
        style = NodeStyle()
        style["fgcolor"] = color  # Set branch color
        style["hz_line_color"] = color
        style["vt_line_color"] = color
        style["size"] = 0  # Hide node circles
        node.set_style(style)

def plot_legend(labels, colors, legend_file="legend.png"):
    """Generate a separate legend image using Matplotlib."""
    fig, ax = plt.subplots(figsize=(4, 3))
    for label, color in zip(labels, colors):
        ax.plot([], [], color=color, label=label, linewidth=4)

    ax.legend(loc="center", fontsize=10, frameon=False)
    ax.axis("off")  # Hide axes
    plt.savefig(legend_file, bbox_inches="tight", dpi=300)  # Save as a separate image
    plt.close()
    print(f"Legend saved as {legend_file}")

def plot_separate_trees(trees, colors, labels, output_prefix="phylogenetic_tree"):
    """Generate separate tree images and save them without merging."""
    for i, tree in enumerate(trees):
        color_branches(tree, colors[i])

        ts = TreeStyle()
        ts.mode = "c"  # Circular layout
        ts.show_leaf_name = True

        output_file = f"{output_prefix}_{labels[i].replace(' ', '_')}.png"
        tree.render(output_file, tree_style=ts)
        print(f"Tree saved as {output_file}")

# Define file names, colors, and labels
tree_files = ["ChainA_aligned.fas.treefile", "ChainB_aligned.fas.treefile",
              "ChainC_aligned.fas.treefile", "ChainD_aligned.fas.treefile", "ChainE_aligned.fas.treefile"]
colors = ["skyblue", "orange", "green", "darkred", "purple"]
labels = ["Chain A", "Chain B", "Chain C", "Chain D", "Chain E"]

trees = load_trees(tree_files)

# Generate individual trees and separate legend
plot_separate_trees(trees, colors, labels)
plot_legend(labels, colors)  # Creates a separate legend image
