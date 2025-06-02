import pandas as pd
import matplotlib.pyplot as plt

# Initial
#df1 = pd.read_csv("../results/mean_deviations_initial_ortho.csv", index_col=False)
df1 = pd.read_csv("../results/mean_deviations_initial_simple.csv", index_col=False)
#df1 = pd.read_csv("../results/mean_deviations_initial_point_set.csv", index_col=False)
#df1 = pd.read_csv("../results/mean_deviations_initial_simple_exterior_pt1.csv", index_col=False)
#df1 = pd.read_csv("../results/mean_deviations_initial_simple_exterior_pt2.csv", index_col=False)
#names1 = df1['name'].str.removeprefix("ortho_")
names1 = df1['name'].str.removeprefix("simple-polygon_")
#names1 = df1['name'].str.removeprefix("point-set_")
#names1 = df1['name'].str.removeprefix("simple-polygon-exterior_")
#names1 = df1['name'].str.removeprefix("simple-polygon-exterior-")
steiner1 = df1['num_steiner']

# Lloyd
#df2 = pd.read_csv("../results/mean_deviations_better_mesh_lloyd_ortho.csv", index_col=False)
df2 = pd.read_csv("../results/mean_deviations_better_mesh_lloyd_simple.csv", index_col=False)
#df2 = pd.read_csv("../results/mean_deviations_better_mesh_lloyd_point_set.csv", index_col=False)
#df2 = pd.read_csv("../results/mean_deviations_better_mesh_lloyd_simple_exterior_pt1.csv", index_col=False)
#df2 = pd.read_csv("../results/mean_deviations_better_mesh_lloyd_simple_exterior_pt2.csv", index_col=False)
#df2 = pd.read_csv("../results/mean_deviations_lloyd_ortho.csv", index_col=False)
#df2 = pd.read_csv("../results/mean_deviations_lloyd_simple.csv", index_col=False)
#df2 = pd.read_csv("../results/mean_deviations_lloyd_point_set.csv", index_col=False)
#df2 = pd.read_csv("../results/mean_deviations_lloyd_simple_exterior_pt1.csv", index_col=False)
#df2 = pd.read_csv("../results/mean_deviations_lloyd_simple_exterior_pt2.csv", index_col=False)
#names2 = df2['name'].str.removeprefix("ortho_")
names2 = df2['name'].str.removeprefix("simple-polygon_")
#names2 = df2['name'].str.removeprefix("point-set_")
#names2 = df2['name'].str.removeprefix("simple-polygon-exterior_")
#names2 = df2['name'].str.removeprefix("simple-polygon-exterior-")
steiner2 = df2['num_steiner']

# Quadratic
#df3 = pd.read_csv("../results/mean_deviations_better_mesh_sq_penalty_ortho.csv", index_col=False)
df3 = pd.read_csv("../results/mean_deviations_better_mesh_sq_penalty_simple.csv", index_col=False)
#df3 = pd.read_csv("../results/mean_deviations_better_mesh_sq_penalty_point_set.csv", index_col=False)
#df3 = pd.read_csv("../results/mean_deviations_better_mesh_sq_penalty_simple_exterior_pt1.csv", index_col=False)
#df3 = pd.read_csv("../results/mean_deviations_better_mesh_sq_penalty_simple_exterior_pt2.csv", index_col=False)
#df3 = pd.read_csv("../results/mean_deviations_energy_sq_ortho.csv", index_col=False)
#df3 = pd.read_csv("../results/mean_deviations_energy_sq_simple.csv", index_col=False)
#df3 = pd.read_csv("../results/mean_deviations_energy_sq_point_set.csv", index_col=False)
#df3 = pd.read_csv("../results/mean_deviations_sq_simple_exterior_pt1.csv", index_col=False)
#df3 = pd.read_csv("../results/mean_deviations_sq_simple_exterior_pt2.csv", index_col=False)
#names3 = df3['name'].str.removeprefix("ortho_")
names3 = df3['name'].str.removeprefix("simple-polygon_")
#names3 = df3['name'].str.removeprefix("point-set_")
#names3 = df3['name'].str.removeprefix("simple-polygon-exterior_")
#names3 = df3['name'].str.removeprefix("simple-polygon-exterior-")
steiner3 = df3['num_steiner']

# Plot
plt.figure(figsize=(12, 6))
#plt.bar(names, deviations, color='skyblue')
plt.plot(names1, steiner1, marker='o', linestyle='-', color='blue', label=r'$\mathbf{Delaunay\ Refinement}$')
plt.plot(names2, steiner2, marker='o', linestyle='-', color='#fc03b1', label=r'$\mathbf{Lloyd}$')
plt.plot(names3, steiner3, marker='o', linestyle='-', color='green', label=r'$\mathbf{Quadratic\ with\ inverse\ penalty}$')

# Formatting
plt.title("Steiner Points per Instance (Simple Polygons)")
plt.xlabel("Instance")
plt.ylabel("# Steiner Points")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.legend()

# Save or show
plt.savefig("steiner_batch_mesh_simple.png", dpi=600, bbox_inches='tight')
# plt.show()

# calculate average deviation decrease

# Merge all on 'name'
merged = df1[['name', 'num_steiner']].rename(columns={'num_steiner': 'steiner1'}) \
    .merge(df2[['name', 'num_steiner']].rename(columns={'num_steiner': 'steiner2'}), on='name') \
    .merge(df3[['name', 'num_steiner']].rename(columns={'num_steiner': 'steiner3'}), on='name')
    #.merge(df4[['name', 'deviation']].rename(columns={'deviation': 'deviation4'}), on='name')

# Compute differences
merged['diff2'] = merged['steiner2'] - merged['steiner1']
merged['diff3'] = merged['steiner3'] - merged['steiner1']
#merged['diff4'] = merged['deviation4'] - merged['deviation1']

# Compute average differences (signed)
mean_diff2 = merged['diff2'].mean()
mean_diff3 = merged['diff3'].mean()
#mean_diff4 = merged['diff4'].mean()

print(f"Average num steiner (Lloyd): {mean_diff2:.4f}")
print(f"Average num steiner (Quadratic): {mean_diff3:.4f}")
#print(f"Average deviation change (Quadratic + penalty): {mean_diff4:.4f}")