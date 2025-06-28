import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df1 = pd.read_csv("../results/fully_nonobtuse_ortho.csv", index_col=False)
#df1 = pd.read_csv("../results/fully_nonobtuse_simple_copy.csv", index_col=False)
#df1 = pd.read_csv("../results/fully_nonobtuse_pointset.csv", index_col=False)
names = df1['name'].str.removeprefix("ortho_")
#names = df1['name'].str.removeprefix("simple-polygon_")
#names = df1['name'].str.removeprefix("point-set_")
#names = df1['name'].str.removeprefix("simple-polygon-exterior_")
#names = df1['name'].str.removeprefix("simple-polygon-exterior-")
steiner1 = df1['#steiner']

df2 = pd.read_csv("../results/fully_nonobtuse_ortho_offcenter.csv", index_col=False)
#df2 = pd.read_csv("../results/fully_nonobtuse_simple_offcenter.csv", index_col=False)
#df2 = pd.read_csv("../results/fully_nonobtuse_pointset_offcenter.csv", index_col=False)
names = df2['name'].str.removeprefix("ortho_")
#names = df2['name'].str.removeprefix("simple-polygon_")
#names = df2['name'].str.removeprefix("point-set_")
steiner2 = df2['#steiner']

# Optimal
steiner3 = df1['#steiner_opt']

# Plot
plt.figure(figsize=(12, 6))
plt.plot(names, steiner1, marker='o', linestyle='-', color='red', label=f'Traditional Delaunay + ' + r'$f_3$')
plt.plot(names, steiner2, marker='o', linestyle='-', color='green', label=f'Off-Center Delaunay + ' + r'$f_3$')
plt.plot(names, steiner3, marker='o', linestyle='-', color='blue', label=f'CG:SHOP Best Score')

# Formatting
plt.title("Number of Steiner Points (Orthogonal Polygons)")
plt.xlabel("Instance")
plt.ylabel("# Steiner")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.legend()

# Save or show
plt.savefig("fully_nonobtuse_ortho.png", dpi=600, bbox_inches='tight')