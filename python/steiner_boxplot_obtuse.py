import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- Load All Datasets ---

# Function to load and extract steiner counts
def load_steiner_data(initial_path, f1_path, f2_path, f3_path):
    df_init = pd.read_csv(initial_path)
    df_f1   = pd.read_csv(f1_path)
    df_f2   = pd.read_csv(f2_path)
    df_f3   = pd.read_csv(f3_path)
    return (df_init['steiner_meshing'],
            df_f1['obtuse_optimized'],
            df_f2['obtuse_optimized'],
            df_f3['obtuse_optimized'])

# --- File paths for each dataset ---
paths = {
    "ortho": (
        "../results/modified3_initial_ln_ortho.csv",
        "../results/modified3_initial_ln_ortho.csv",
        "../results/modified3_initial_sigmoid_ortho.csv",
        "../results/modified3_initial_refined_sigmoid_ortho.csv",
    ),
    "simple": (
        "../results/modified3_initial_ln_simple.csv",
        "../results/modified3_initial_ln_simple.csv",
        "../results/modified3_initial_sigmoid_simple.csv",
        "../results/modified3_initial_refined_sigmoid_simple.csv",
    ),
    "point_set": (
        "../results/modified3_initial_ln_point_set.csv",
        "../results/modified3_initial_ln_point_set.csv",
        "../results/modified3_initial_sigmoid_point_set.csv",
        "../results/modified3_initial_refined_sigmoid_point_set.csv",
    ),
    "simple_exterior": (
        "../results/modified3_initial_ln_simple_exterior_pt1.csv",
        "../results/modified3_initial_ln_simple_exterior_pt1.csv",
        "../results/modified3_initial_sigmoid_simple_exterior_pt1.csv",
        "../results/modified3_initial_refined_sigmoid_simple_exterior_pt1.csv",
    )
}

# --- Aggregate Data ---
all_initial, all_f1, all_f2, all_f3 = [], [], [], []

for dataset in paths.values():
    initial, f1, f2, f3 = load_steiner_data(*dataset)
    all_initial.extend(initial)
    all_f1.extend(f1)
    all_f2.extend(f2)
    all_f3.extend(f3)

# --- Plot Box Plot ---
data = [all_initial, all_f1, all_f2, all_f3]
labels = [r'$\mathbf{Before}$', r'$\mathbf{f_1}$', r'$\mathbf{f_2}$', r'$\mathbf{f_3}$']
colors = ['orange', 'green', 'blue', 'red']

fig, ax = plt.subplots(figsize=(10, 6))
box = ax.boxplot(data, patch_artist=True, labels=labels)

for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)

ax.set_ylabel("Number of Obtuse Triangles")
ax.set_title("Distribution of Obtuse Triangles Before and After Optimization Across All Instances")
ax.grid(axis='y', linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig("steiner_boxplot_all.png", dpi=600)

med_before = np.median(all_initial)
med_f1 = np.median(all_f1)
med_f2 = np.median(all_f2)
med_f3 = np.median(all_f3)

print(f"Median Steiner Points (Before Optimization): {med_before:.2f}")
print(f"Median Steiner Points (f1): {med_f1:.2f}")
print(f"Median Steiner Points (f2): {med_f2:.2f}")
print(f"Median Steiner Points (f3): {med_f3:.2f}")
