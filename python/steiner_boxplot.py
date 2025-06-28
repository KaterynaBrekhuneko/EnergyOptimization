import pandas as pd
import matplotlib.pyplot as plt

# Load or define your three DataFrames
df1 = pd.read_csv("../results/steiner_merged_initial.csv", index_col=False)
df2 = pd.read_csv("../results/steiner_merged_better_mesh_lloyd.csv", index_col=False)
df3 = pd.read_csv("../results/steiner_merged_better_mesh_sq_penalty.csv", index_col=False)

# Rename 'steiner' column to identify the algorithm
df1 = df1.rename(columns={"num_steiner": "Delaunay Refinement"})
df2 = df2.rename(columns={"num_steiner": "Batch Refinement + Lloyd"})
df3 = df3.rename(columns={"num_steiner": "Batch Refinement + "+r"$g_2$"})

# Merge all dataframes on 'name' (i.e., instance name)
df_merged = df1.merge(df2, on="name").merge(df3, on="name")

# Collect data
data = [
    df_merged["Delaunay Refinement"],
    df_merged["Batch Refinement + Lloyd"],
    df_merged["Batch Refinement + "+r"$g_2$"]
]
labels = ["Terminator Algorithm", "Batch Refinement + Lloyd", "Batch Refinement + $g_2$"]

# Plot
fig, ax = plt.subplots()
box = ax.boxplot(data, patch_artist=True, tick_labels=labels)

# Colors
for patch in box['boxes']:
    patch.set_facecolor('lightblue')

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
plt.savefig("steiner_box_plot.png", dpi=600, bbox_inches='tight')
