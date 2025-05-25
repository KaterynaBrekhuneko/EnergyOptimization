import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
df1 = pd.read_csv("../results/angle_data_before.csv")

# Print min and max angle
min_angle1 = df1['angle'].min()
max_angle1 = df1['angle'].max()
print(f"Min angle before opt: {min_angle1:.2f}°")
print(f"Max angle before opt: {max_angle1:.2f}°")

# Histogram bin edges (1° bins from 0° to 180°)
bins1 = np.arange(0, 181, 1)
hist1, bin_edges1 = np.histogram(df1['angle'], bins=bins1)

# X = bin centers
bin_centers1 = (bin_edges1[:-1] + bin_edges1[1:]) / 2

# after optimization lloyd
# Load data
df2 = pd.read_csv("../results/angle_data_lloyd.csv")

# Print min and max angle
min_angle2 = df2['angle'].min()
max_angle2 = df2['angle'].max()
print(f"Min angle after lloyd: {min_angle2:.2f}°")
print(f"Max angle fter lloyd: {max_angle2:.2f}°")

# Histogram bin edges (1° bins from 0° to 180°)
bins2 = np.arange(0, 181, 1)
hist2, bin_edges2 = np.histogram(df2['angle'], bins=bins2)

# X = bin centers
bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

# after optimization energy
# Load data
df3 = pd.read_csv("../results/angle_data_equilateral.csv")

# Print min and max angle
min_angle3 = df3['angle'].min()
max_angle3 = df3['angle'].max()
print(f"Min angle energy: {min_angle3:.2f}°")
print(f"Max angle energy: {max_angle3:.2f}°")

# Histogram bin edges (1° bins from 0° to 180°)
bins3 = np.arange(0, 181, 1)
hist3, bin_edges3 = np.histogram(df3['angle'], bins=bins3)

# X = bin centers
bin_centers3 = (bin_edges3[:-1] + bin_edges3[1:]) / 2

# Plot as line chart
plt.figure(figsize=(10, 6))
plt.plot(bin_centers1, hist1, linestyle='-', linewidth=1.5, color='blue')
plt.plot(bin_centers2, hist2, linestyle='-', linewidth=1.5, color='red')
plt.plot(bin_centers3, hist3, linestyle='-', linewidth=1.5, color='green')

# Labels and grid
plt.title("Triangle Mesh Angle Distribution")
plt.xlabel("Angle (degrees)")
plt.ylabel("# Angles")
#plt.tight_layout()
plt.savefig("angle_distribution_line_chart.png", dpi=300)