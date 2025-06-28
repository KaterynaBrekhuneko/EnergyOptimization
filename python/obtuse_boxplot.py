import pandas as pd
import matplotlib.pyplot as plt

# Load or define your three DataFrames
df1 = pd.concat([
    pd.read_csv("../results/fully_nonobtuse_ortho.csv", index_col=False),
    pd.read_csv("../results/fully_nonobtuse_simple_copy.csv", index_col=False),
    pd.read_csv("../results/fully_nonobtuse_pointset.csv", index_col=False)
], ignore_index=True)
df2 = pd.concat([
    pd.read_csv("../results/fully_nonobtuse_ortho_offcenter.csv", index_col=False),
    pd.read_csv("../results/fully_nonobtuse_simple_offcenter.csv", index_col=False),
    pd.read_csv("../results/fully_nonobtuse_pointset_offcenter.csv", index_col=False)
], ignore_index=True)


# Rename 'steiner' column to identify the algorithm
df1 = df1.rename(columns={"#steiner": "Traditional Delaunay"})
df1 = df1.rename(columns={"#steiner_opt": "CG:SHOP Best Solution"})
df2 = df2.rename(columns={"#steiner": "Off-Center Delaunay"})

# Merge all dataframes on 'name' (i.e., instance name)
df_merged = df1.merge(df2, on="name")

# Collect data
data = [
    df_merged["CG:SHOP Best Solution"],
    df_merged["Traditional Delaunay"],
    df_merged["Off-Center Delaunay"]
]
labels = ["CG:SHOP Best Solution", "Traditional Delaunay", "Off-Center Delaunay"]

# Plot
fig, ax = plt.subplots()
box = ax.boxplot(data, patch_artist=True, tick_labels=labels)

# Colors
for patch in box['boxes']:
    patch.set_facecolor('lightgreen')

for median in box['medians']:
    median.set_color('mediumseagreen')
    median.set_linewidth(2)

# Horizontal grid only
ax.yaxis.grid(True)
ax.xaxis.grid(False)

# Labels and layout
ax.set_ylabel("Number of Steiner Points")
ax.set_title("Comparison of Steiner Point Counts")
plt.tight_layout()
plt.savefig("steiner_box_plot_fully_non_obtuse.png", dpi=600, bbox_inches='tight')
