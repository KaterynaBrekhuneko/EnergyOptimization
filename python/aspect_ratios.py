import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load aspect ratio data initial
df1 = pd.read_csv("../results/aspect_ratios_initial_simple-polygon_20_4bd3c2e5.csv")

max_ratio1 = df1['ratio'].max()
print(f"Max aspect ratio initial: {max_ratio1:.4f}")

bins1 = np.linspace(0.57, max_ratio1, 50)
hist1, bin_edges1 = np.histogram(df1['ratio'], bins=bins1)

bin_centers1 = (bin_edges1[:-1] + bin_edges1[1:]) / 2

# Load aspect ratio data sigmoid or lloyd
df2 = pd.read_csv("../results/aspect_ratios_lloyd_simple-polygon_20_4bd3c2e5.csv")

max_ratio2 = df2['ratio'].max()
print(f"Max aspect ratio optimization: {max_ratio2:.4f}")

bins2 = np.linspace(0.57, max_ratio2, 50)
hist2, bin_edges2 = np.histogram(df2['ratio'], bins=bins2)

bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

# Load aspect ratio data quadratic
df3 = pd.read_csv("../results/aspect_ratios_sq_penalty_simple-polygon_20_4bd3c2e5.csv")

max_ratio3 = df3['ratio'].max()
print(f"Max aspect ratio equilateral: {max_ratio3:.4f}")

bins3 = np.linspace(0.57, max_ratio3, 50)
hist3, bin_edges3 = np.histogram(df3['ratio'], bins=bins3)

bin_centers3 = (bin_edges3[:-1] + bin_edges3[1:]) / 2

# Load aspect ratio data quadratic with penalty
"""df4 = pd.read_csv("../results/aspect_ratios_equilateral_mod.csv")

max_ratio4 = df4['ratio'].max()
print(f"Max aspect ratio equilateral with penalty: {max_ratio4:.4f}")

bins4 = np.linspace(0.57, max_ratio4, 50)
hist4, bin_edges4 = np.histogram(df4['ratio'], bins=bins4)

count = (df4['ratio'] <= 0.7).sum()
total = len(df4['ratio'])
percentage = (count / total) * 100
print(f"Percentage of ratios ≤ 0.7: {percentage:.2f}%")

bin_centers4 = (bin_edges4[:-1] + bin_edges4[1:]) / 2"""

# Plot as line chart
plt.figure(figsize=(10, 6))
plt.plot(bin_centers1, hist1, linestyle='-', linewidth=1.5, color='blue', label=r'$\mathbf{Terminator\ Algorithm}$' + f'\nmax = {max_ratio1:.2f}')
plt.plot(bin_centers2, hist2, linestyle='-', linewidth=1.5, color='#fc03b1', label=r'$\mathbf{Batch\ Delaunay\ +\ Lloyd}$' + f'\nmax = {max_ratio2:.2f}')
plt.plot(bin_centers3, hist3, linestyle='-', linewidth=1.5, color='green', label=r'$\mathbf{Batch\ Delaunay\ +\ g_2}$' + f'\nmax = {max_ratio3:.2f}')
#plt.plot(bin_centers4, hist4, linestyle='-', linewidth=1.5, color='orange', label=r'$\mathbf{Quadratic\ with\ inverse\ penalty}$' + f'\nmax = {max_ratio4:.2f}')

# Labels and ticks
plt.xlabel("Aspect Ratio")
plt.ylabel("# Triangles")
max = max(max_ratio1, max_ratio2)
plt.xticks(np.arange(0.57, max, 0.1))
plt.legend()
# Save plot
plt.savefig("aspect_ratio_simple-polygon_20_4bd3c2e5.png", dpi=600, bbox_inches='tight')