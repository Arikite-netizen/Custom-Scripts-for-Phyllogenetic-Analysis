import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO
import math
import os
from collections import Counter

# ===== Load phylogenetic organism data =====
phylo_data = pd.read_csv("./phyllogenetic_analysis.csv")

# Flatten the organism names and remove any NaNs
organisms = phylo_data.values.flatten()
organisms = [org for org in organisms if pd.notnull(org)]

# Count organism frequencies
counts = Counter(organisms)
total = sum(counts.values())
probabilities = [count / total for count in counts.values()]

# Shannon entropy from organism frequencies
phylo_entropy = -np.sum([p * np.log2(p) for p in probabilities])
print(f"Shannon entropy from organism frequencies: {phylo_entropy:.4f}")

# ===== Function to calculate positional entropy =====
def calculate_positional_entropy(alignment):
    alignment_length = alignment.get_alignment_length()
    entropy_scores = np.zeros(alignment_length)

    for column in range(alignment_length):
        column_bases = alignment[:, column]
        base_frequencies = {base: column_bases.count(base) / len(column_bases) for base in set(column_bases)}
        entropy = -sum([freq * math.log(freq, 2) for freq in base_frequencies.values()])
        entropy_scores[column] = entropy

    return entropy_scores

# ===== Process aligned FASTA files =====
fasta_files = ["ChainA_aligned.fas", "ChainB_aligned.fas", "ChainC_aligned.fas", "ChainD_aligned.fas", "ChainE_aligned.fas"]

entropy_dict = {}
max_length = 0

for fasta_file in fasta_files:
    chain_name = os.path.splitext(fasta_file)[0]
    alignment = AlignIO.read(fasta_file, "fasta")
    entropy_scores = calculate_positional_entropy(alignment)

    entropy_dict[chain_name] = entropy_scores
    if len(entropy_scores) > max_length:
        max_length = len(entropy_scores)

# ===== Create DataFrame for entropy across chains =====
entropy_df = pd.DataFrame()

for chain, entropy_scores in entropy_dict.items():
    padded_entropy = np.pad(entropy_scores, (0, max_length - len(entropy_scores)), constant_values=np.nan)
    entropy_df[chain] = padded_entropy

entropy_df.index = range(1, max_length + 1)

# Save combined entropy (optional)
entropy_df.to_csv("combined_entropy_landscape.csv", index_label="Position")

# ===== Plot the heatmap with phylogenetic entropy bar =====
fig = plt.figure(figsize=(18, 8))
grid = plt.GridSpec(1, 10, wspace=0.3, hspace=0.1)

# Heatmap for chain entropy (no black lines)
ax_main = fig.add_subplot(grid[0, :-1])
sns.heatmap(entropy_df.T, cmap="coolwarm", cbar_kws={'label': 'Entropy'}, ax=ax_main,
            xticklabels=50, yticklabels=True)

ax_main.set_xlabel('Position in Alignment')
ax_main.set_ylabel('Protein Chain')
ax_main.set_title('Entropy Landscape Heatmap Across Chains (A–E)')

# Vertical bar for phylogenetic entropy
ax_entropy = fig.add_subplot(grid[0, -1])
entropy_bar = np.full((len(entropy_df.columns), 1), phylo_entropy)
sns.heatmap(entropy_bar, cmap="viridis", cbar_kws={'label': 'Phylogenetic Entropy'},
            ax=ax_entropy, xticklabels=False, yticklabels=False)

ax_entropy.set_xlabel('')
ax_entropy.set_ylabel('')

plt.show()
