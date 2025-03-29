import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the positional entropy data
entropy_data = pd.read_csv('positional_entropy.csv')

# Extract position and entropy values
positions = entropy_data['Position']
entropy_values = entropy_data['Entropy']

# Define thresholds for low and high entropy
low_entropy_threshold = np.percentile(entropy_values, 25)  # Lower 25% as low entropy
high_entropy_threshold = np.percentile(entropy_values, 75)  # Upper 25% as high entropy

# Identify low and high entropy positions
low_entropy_positions = positions[entropy_values <= low_entropy_threshold]
high_entropy_positions = positions[entropy_values >= high_entropy_threshold]

# Create figure
plt.figure(figsize=(12, 6))

# Plot the entropy curve
plt.plot(positions, entropy_values, linestyle='-', color='black', linewidth=1.5, label="Entropy", alpha=0.9)

# Highlight low entropy regions in cyan
plt.fill_between(positions, entropy_values, low_entropy_threshold,
                 where=(entropy_values <= low_entropy_threshold),
                 color='cyan', alpha=0.5, label="Low Entropy Regions (Conserved)")

# Highlight high entropy regions in red
plt.fill_between(positions, entropy_values, high_entropy_threshold,
                 where=(entropy_values >= high_entropy_threshold),
                 color='red', alpha=0.5, label="High Entropy Regions (Variable)")

# Add vertical grid lines at low entropy positions
for pos in low_entropy_positions[::5]:  # Plot every 5th position to reduce clutter
    plt.axvline(x=pos, color='blue', linestyle='--', alpha=0.2)

# Add labels at key low entropy positions, placed near the x-axis
for i, pos in enumerate(low_entropy_positions[::10]):  # Adjust step size to reduce clutter
    plt.text(pos, min(entropy_values) - 0.1, f"{int(pos)}",
             fontsize=6, ha='right', color='blue')

# Titles and Labels
plt.title('Entropy Analysis of Conserved & Variable Regions in COX1 Gene', fontsize=14, fontweight='bold')
plt.xlabel('Position in COX1 Gene Sequence', fontsize=12)
plt.ylabel('Entropy Score', fontsize=12)

# Enable grid and legend
plt.grid(True, linestyle='--', alpha=0.3)
plt.legend()

# Show plot
plt.show()
