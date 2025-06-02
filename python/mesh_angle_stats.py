import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
df1 = pd.read_csv("../results/angle_stats_initial_simple-polygon_20_4bd3c2e5.csv")

min_angle1 = df1['angle'].min()
max_angle1 = df1['angle'].max()
print(f"Min angle before opt: {min_angle1:.2f}°")
print(f"Max angle before opt: {max_angle1:.2f}°")

bins1 = np.arange(0, 181, 1)
hist1, bin_edges1 = np.histogram(df1['angle'], bins=bins1)

bin_centers1 = (bin_edges1[:-1] + bin_edges1[1:]) / 2

# after optimization sigmoid or lloyd
df2 = pd.read_csv("../results/angle_stats_lloyd_simple-polygon_20_4bd3c2e5.csv")

min_angle2 = df2['angle'].min()
max_angle2 = df2['angle'].max()
print(f"Min angle after sigmoid: {min_angle2:.2f}°")
print(f"Max angle after sigmoid: {max_angle2:.2f}°")

bins2 = np.arange(0, 181, 1)
hist2, bin_edges2 = np.histogram(df2['angle'], bins=bins2)

bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

# after optimization quadratic
df3 = pd.read_csv("../results/angle_stats_sq_penalty_simple-polygon_20_4bd3c2e5.csv")

min_angle3 = df3['angle'].min()
max_angle3 = df3['angle'].max()
print(f"Min angle quadratic: {min_angle3:.2f}°")
print(f"Max angle quadratic: {max_angle3:.2f}°")

bins3 = np.arange(0, 181, 1)
hist3, bin_edges3 = np.histogram(df3['angle'], bins=bins3)

bin_centers3 = (bin_edges3[:-1] + bin_edges3[1:]) / 2

# after optimization quadratic with penalty
"""df4 = pd.read_csv("../results/angle_data_equilateral_mod.csv")

min_angle4 = df4['angle'].min()
max_angle4 = df4['angle'].max()
print(f"Min angle quadratic with penalty: {min_angle4:.2f}°")
print(f"Max angle quadratic with penalty: {max_angle4:.2f}°")

bins4 = np.arange(0, 181, 1)
hist4, bin_edges4 = np.histogram(df4['angle'], bins=bins4)

bin_centers4 = (bin_edges4[:-1] + bin_edges4[1:]) / 2"""

# Plot as line chart
plt.figure(figsize=(10, 6))
plt.plot(bin_centers1, hist1, linestyle='-', linewidth=1.5, color='blue', label=r'$\mathbf{Delaunay\ Refinement}$' + f'\nangles in [{min_angle1:.2f}°, {max_angle1:.2f}°]')
plt.plot(bin_centers2, hist2, linestyle='-', linewidth=1.5, color='#fc03b1', label=r'$\mathbf{Batch\ Delaunay\ +\ Lloyd}$' + f'\n angles in [{min_angle2:.2f}°, {max_angle2:.2f}°]')
plt.plot(bin_centers3, hist3, linestyle='-', linewidth=1.5, color='green', label=r'$\mathbf{Batch\ Delaunay\ +\ g_2}$' + f'\n angles in [{min_angle3:.2f}°, {max_angle3:.2f}°]')
#plt.plot(bin_centers4, hist4, linestyle='-', linewidth=1.5, color='orange', label=r'$\mathbf{Quadratic\ with\ inverse\ penalty}$' + f'\n angles in [{min_angle4:.2f}°, {max_angle4:.2f}°]')
#plt.plot(bin_centers3, hist3, linestyle='-', linewidth=1.5, color='green')

#plt.fill_between(bin_centers1, hist1, color='blue', alpha=0.1)
#plt.fill_between(bin_centers2, hist2, color='red', alpha=0.1)
#plt.fill_between(bin_centers3, hist3, color='green', alpha=0.1)
#plt.fill_between(bin_centers4, hist4, color='orange', alpha=0.1)

# Labels and grid
plt.xlabel("Angle (degrees)")
plt.ylabel("# Angles")
# Set x-axis ticks every 30 degrees
plt.xticks(np.arange(0, 181, 30))
plt.legend()
#plt.tight_layout()
plt.savefig("angle_distribution_simple-polygon_20_4bd3c2e5.png", dpi=600, bbox_inches='tight')
