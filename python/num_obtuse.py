import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Initial
#df1 = pd.read_csv("../results/modified3_initial_ln_ortho.csv", index_col=False)
#df1 = pd.read_csv("../results/modified3_initial_ln_simple.csv", index_col=False)
#df1 = pd.read_csv("../results/modified3_initial_ln_point_set.csv", index_col=False)
#df1 = pd.read_csv("../results/modified3_initial_ln_simple_exterior_pt1.csv", index_col=False)
df1 = pd.read_csv("../results/modified3_initial_ln_simple_exterior_pt2.csv", index_col=False)
#names1 = df1['instance'].str.removeprefix("ortho_")
#names1 = df1['instance'].str.removeprefix("simple-polygon_")
#names1 = df1['instance'].str.removeprefix("point-set_")
#names1 = df1['instance'].str.removeprefix("simple-polygon-exterior_")
names1 = df1['instance'].str.removeprefix("simple-polygon-exterior-")
obtuse1 = df1['obtuse_meshing']

# f1
#df2 = pd.read_csv("../results/modified3_initial_ln_ortho.csv", index_col=False)
#df2 = pd.read_csv("../results/modified3_initial_ln_simple.csv", index_col=False)
#df2 = pd.read_csv("../results/modified3_initial_ln_point_set.csv", index_col=False)
#df2 = pd.read_csv("../results/modified3_initial_ln_simple_exterior_pt1.csv", index_col=False)
df2 = pd.read_csv("../results/modified3_initial_ln_simple_exterior_pt2.csv", index_col=False)
#names2 = df2['instance'].str.removeprefix("ortho_")
#names2 = df2['instance'].str.removeprefix("simple-polygon_")
#names2 = df2['instance'].str.removeprefix("point-set_")
#names2 = df2['instance'].str.removeprefix("simple-polygon-exterior_")
names2 = df2['instance'].str.removeprefix("simple-polygon-exterior-")
obtuse2 = df2['obtuse_optimized']

# f2
#df3 = pd.read_csv("../results/modified3_initial_sigmoid_ortho.csv", index_col=False)
#df3 = pd.read_csv("../results/modified3_initial_sigmoid_simple.csv", index_col=False)
#df3 = pd.read_csv("../results/modified3_initial_sigmoid_point_set.csv", index_col=False)
#df3 = pd.read_csv("../results/modified3_initial_sigmoid_simple_exterior_pt1.csv", index_col=False)
df3 = pd.read_csv("../results/modified3_initial_sigmoid_simple_exterior_pt2.csv", index_col=False)
#names3 = df3['instance'].str.removeprefix("ortho_")
#names3 = df3['instance'].str.removeprefix("simple-polygon_")
#names3 = df3['instance'].str.removeprefix("point-set_")
#names3 = df3['instance'].str.removeprefix("simple-polygon-exterior_")
names3 = df3['instance'].str.removeprefix("simple-polygon-exterior-")
obtuse3 = df3['obtuse_optimized']

# f3
#df4 = pd.read_csv("../results/modified3_initial_refined_sigmoid_ortho.csv", index_col=False)
#df4 = pd.read_csv("../results/modified3_initial_refined_sigmoid_simple.csv", index_col=False)
#df4 = pd.read_csv("../results/modified3_initial_refined_sigmoid_point_set.csv", index_col=False)
#df4 = pd.read_csv("../results/modified3_initial_refined_sigmoid_simple_exterior_pt1.csv", index_col=False)
df4 = pd.read_csv("../results/modified3_initial_refined_sigmoid_simple_exterior_pt2.csv", index_col=False)
#names4 = df4['instance'].str.removeprefix("ortho_")
#names4 = df4['instance'].str.removeprefix("simple-polygon_")
#names4 = df4['instance'].str.removeprefix("point-set_")
#names4 = df4['instance'].str.removeprefix("simple-polygon-exterior_")
names4 = df4['instance'].str.removeprefix("simple-polygon-exterior-")
obtuse4 = df4['obtuse_optimized']

# Plot
plt.figure(figsize=(12, 6))
#plt.plot(names1, obtuse1, marker='o', linestyle='-', color='orange', label=r'$\mathbf{Before\ Optimization}$')
plt.plot(names2, obtuse2, marker='o', linestyle='-', color='green', label=r'$\mathbf{f_1}$')
plt.plot(names3, obtuse3, marker='o', linestyle='-', color='blue', label=r'$\mathbf{f_2}$')
plt.plot(names4, obtuse4, marker='o', linestyle='-', color='red', label=r'$\mathbf{f_3}$')

# Formatting
plt.title("Number of obtuse triangles before and after different optimizations (Simple Polygons with Constraints (Part 2))")
plt.xlabel("Instance")
plt.ylabel("# Obtuse")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.legend()

# Save or show
plt.savefig("2_steiner_initial_simple_exterior_pt2.png", dpi=600, bbox_inches='tight')

# --- Prepare for Bar Plot ---
"""x = np.arange(len(names1))  # Instance indices
width = 0.2  # Width of each bar

# --- Create Bar Plot ---
fig, ax = plt.subplots(figsize=(14, 6))

ax.bar(x - 1.5 * width, obtuse1, width, label=r'$\mathbf{Before}$', color='orange')
ax.bar(x - 0.5 * width, obtuse2, width, label=r'$\mathbf{f_1}$', color='green')
ax.bar(x + 0.5 * width, obtuse3, width, label=r'$\mathbf{f_2}$', color='blue')
ax.bar(x + 1.5 * width, obtuse4, width, label=r'$\mathbf{f_3}$', color='red')

# --- Formatting ---
ax.set_xlabel("Instance")
ax.set_ylabel("# Steiner")
ax.set_title("Number of Steiner Points Before and After Optimization (Orthogonal Polygons)")
ax.set_xticks(x)
ax.set_xticklabels(names1, rotation=45, ha='right')
ax.legend()
ax.grid(axis='y', linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig("steiner_initial_ortho_bar.png", dpi=600, bbox_inches='tight')"""