import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Load the combined entropy landscape data
entropy_data = pd.read_csv('combined_entropy_landscape.csv')

# Extract positions and chain names
positions = entropy_data['Position']
chain_names = entropy_data.columns[1:]

# Dynamic figure size based on number of chains
num_chains = len(chain_names)
fig, axes = plt.subplots(num_chains, 1, figsize=(18, 4 * num_chains), sharex=True)

# Ensure axes is iterable
if num_chains == 1:
    axes = [axes]

# Color palette
palette = sns.color_palette('tab10', num_chains)

for idx, chain in enumerate(chain_names):
    ax = axes[idx]
    entropy_values = entropy_data[chain]

    # Thresholds for highlights
    low_entropy_threshold = np.percentile(entropy_values.dropna(), 25)
    high_entropy_threshold = np.percentile(entropy_values.dropna(), 75)

    # Plot entropy
    ax.plot(positions, entropy_values, color=palette[idx], linewidth=1.5)

    # Highlight regions
    ax.fill_between(positions, entropy_values, low_entropy_threshold,
                    where=(entropy_values <= low_entropy_threshold),
                    color=palette[idx], alpha=0.15, label='Low Entropy')

    ax.fill_between(positions, entropy_values, high_entropy_threshold,
                    where=(entropy_values >= high_entropy_threshold),
                    color=palette[idx], alpha=0.25, label='High Entropy')

    ax.set_ylabel('Entropy', fontsize=12)
    ax.set_title(f"{chain} Entropy Landscape", fontsize=14, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend()

# Shared X-axis label
plt.xlabel('Position in Alignment', fontsize=14)

# Adjust spacing to prevent overlap
plt.subplots_adjust(hspace=0.4)

plt.show()
