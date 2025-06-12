import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# My
#df = pd.read_csv("../results/fully_nonobtuse_ortho.csv", index_col=False)
#df = pd.read_csv("../results/fully_nonobtuse_simple.csv", index_col=False)
df = pd.read_csv("../results/fully_nonobtuse_pointset.csv", index_col=False)
#names = df['name'].str.removeprefix("ortho_")
#names = df['name'].str.removeprefix("simple-polygon_")
names = df['name'].str.removeprefix("point-set_")
#names = df['name'].str.removeprefix("simple-polygon-exterior_")
#names = df['name'].str.removeprefix("simple-polygon-exterior-")
steiner1 = df['#steiner']

# Optimal
steiner2 = df['#steiner_opt']

# Plot
plt.figure(figsize=(12, 6))
plt.plot(names, steiner1, marker='o', linestyle='-', color='red', label=f'Traditional Delaunay +' + r'$f_3$')
plt.plot(names, steiner2, marker='o', linestyle='-', color='blue', label=f'CG:SHOP Best Score')

# Formatting
plt.title("Number of steiner (Point-sets)")
plt.xlabel("Instance")
plt.ylabel("# Steiner")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.legend()

# Save or show
plt.savefig("fully_nonobtuse_pointset.png", dpi=600, bbox_inches='tight')