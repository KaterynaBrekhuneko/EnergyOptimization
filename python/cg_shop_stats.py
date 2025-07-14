import pandas as pd

# Load all CSV files
cg_shop = pd.read_csv("cg_shop_simple_exterior.csv", sep=";")
ln = pd.read_csv("../results/modified3_initial_ln_simple_exterior.csv")
sigmoid = pd.read_csv("../results/modified3_initial_sigmoid_simple_exterior.csv")
refined_sigmoid = pd.read_csv("../results/modified3_initial_refined_sigmoid_simple_exterior.csv")

# Convert num obtuse to dictionaries for fast lookup
obtuse_ln = ln.set_index("instance")["obtuse_optimized"].to_dict()
obtuse_sigmoid = sigmoid.set_index("instance")["obtuse_optimized"].to_dict()
obtuse_refined_sigmoid = refined_sigmoid.set_index("instance")["obtuse_optimized"].to_dict()

steiner_ln = ln.set_index("instance")["steiner_optimized"].to_dict()
steiner_sigmoid = sigmoid.set_index("instance")["steiner_optimized"].to_dict()
steiner_refined_sigmoid = refined_sigmoid.set_index("instance")["steiner_optimized"].to_dict()

# Prepare result list
result_rows = []

# Process each name from fully_nonobtuse_ortho.csv
for name in cg_shop["instance"]:

    row = cg_shop[cg_shop["instance"] == name].iloc[0]

    if name not in ln["instance"].values:
        steiner_dortmund = int(row["steiner_dortmund"])
        obtuse_dortmund = int(row["obtuse_dortmund"])
        # Assemble new row
        result_rows.append({
            "instance": name,
            "obtuse_best": row["obtuse_best"],
            "steiner_best": row["steiner_best"],
            "obtuse_dortmund": obtuse_dortmund,
            "steiner_dortmund": steiner_dortmund
        })
        continue

    # Compute min steiner_opt from both sources
    o_opt1 = obtuse_ln.get(name, float('inf'))
    o_opt2 = obtuse_sigmoid.get(name, float('inf'))
    o_opt3 = obtuse_refined_sigmoid.get(name, float('inf'))

    s_opt1 = steiner_ln.get(name, float('inf'))
    s_opt2 = steiner_sigmoid.get(name, float('inf'))
    s_opt3 = steiner_refined_sigmoid.get(name, float('inf'))
    options = [
        (o_opt1, s_opt1),
        (o_opt2, s_opt2),
        (o_opt3, s_opt3)
    ]
    min_obtuse = min(o for o, _ in options)
    min_steiner_opt = min(s for o, s in options if o == min_obtuse)
    min_obtuse_opt = min_obtuse

    #if name == "simple-polygon-exterior_10_40642b31":
        #print(name)
        #print(options)

    #print(name)
    #print(min_obtuse_opt)
    #print(min_steiner_opt)

    # Decide steiner_dortmund
    if row["obtuse_dortmund"] == 0:
        if min_obtuse_opt == 0:
            steiner_dortmund = min(int(row["steiner_dortmund"]), min_steiner_opt)
        else:
            steiner_dortmund = int(row["steiner_dortmund"])    
        obtuse_dortmund = 0
    else:
        if row["obtuse_dortmund"] < min_obtuse_opt:
            steiner_dortmund = int(row["steiner_dortmund"]) 
            obtuse_dortmund = row["obtuse_dortmund"]
        elif row["obtuse_dortmund"] == min_obtuse_opt:
            steiner_dortmund = min(int(row["steiner_dortmund"]), min_steiner_opt)
            obtuse_dortmund = min_obtuse_opt
        else:
            steiner_dortmund = min_steiner_opt
            obtuse_dortmund = min_obtuse_opt

    # Assemble new row
    result_rows.append({
        "instance": name,
        "obtuse_best": row["obtuse_best"],
        "steiner_best": row["steiner_best"],
        "obtuse_dortmund": obtuse_dortmund,
        "steiner_dortmund": steiner_dortmund
    })

# Write to output CSV
pd.DataFrame(result_rows).to_csv("cg_shop_result.csv", sep=";", index=False)
