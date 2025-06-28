import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load aspect ratio data initial
#df1 = pd.read_csv("../results/aspect_ratios_quad_blossom_initial.csv")
df1 = pd.read_csv("../results/aspect_ratios_quad_medians_initial.csv")
max_ratio1 = df1['ratio'].max()
print(f"Max aspect ratio initial: {max_ratio1:.4f}")
bins1 = np.linspace(1, max_ratio1, 50)
hist1, bin_edges1 = np.histogram(df1['ratio'], bins=bins1)
bin_centers1 = (bin_edges1[:-1] + bin_edges1[1:]) / 2

# Load aspect ratio data optimized
#df2 = pd.read_csv("../results/aspect_ratios_quad_blossom_optimized_angles.csv")
df2 = pd.read_csv("../results/aspect_ratios_quad_medians_optimized_angles.csv")
max_ratio2 = df2['ratio'].max()
print(f"Max aspect ratio optimization: {max_ratio2:.4f}")
bins2 = np.linspace(1, max_ratio2, 50)
hist2, bin_edges2 = np.histogram(df2['ratio'], bins=bins2)
bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

# Load aspect ratio data optimized
#df3 = pd.read_csv("../results/aspect_ratios_quad_blossom_optimized.csv")
df3 = pd.read_csv("../results/aspect_ratios_quad_medians_optimized.csv")
max_ratio3 = df3['ratio'].max()
print(f"Max aspect ratio optimization: {max_ratio3:.4f}")
bins3 = np.linspace(1, max_ratio3, 50)
hist3, bin_edges3 = np.histogram(df3['ratio'], bins=bins3)
bin_centers3 = (bin_edges3[:-1] + bin_edges3[1:]) / 2

# Plot as line chart
plt.figure(figsize=(10, 6))
plt.plot(bin_centers1, hist1, linestyle='-', linewidth=1.5, color='blue', label=r'$\mathbf{Initial\ Medians}$' + f'\nmax = {max_ratio1:.2f}')
plt.plot(bin_centers2, hist2, linestyle='-', linewidth=1.5, color='#fc03b1', label=r'$\mathbf{g_3}$' + f'\nmax = {max_ratio2:.2f}')
plt.plot(bin_centers3, hist3, linestyle='-', linewidth=1.5, color='green', label=r'$\mathbf{g_4}$' + f'\nmax = {max_ratio3:.2f}')

# Labels and ticks
plt.xlabel("Aspect Ratio")
plt.ylabel("# Quads")

# Dynamically choose fewer xticks
xmax = max(max_ratio1, max_ratio2)
xticks = np.round(np.linspace(1, xmax + 0.05, num=8), 2)
plt.xticks(xticks)

plt.grid(alpha=0.3, linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()

# Save plot
plt.savefig("aspect_ratio_medians.png", dpi=600, bbox_inches='tight')
