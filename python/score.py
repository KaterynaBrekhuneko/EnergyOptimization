import pandas as pd
import math

# Load all CSV files
cg_shop = pd.read_csv("cg_shop.csv", sep=";")
cg_shop_result = pd.read_csv("cg_shop_result.csv", sep=";")

score = 0.0
num_non_obtuse = 0
num_better = 0

for name in cg_shop["instance"]:

    row = cg_shop[cg_shop["instance"] == name].iloc[0]
    row_result = cg_shop_result[cg_shop_result["instance"] == name].iloc[0]

    if int(row["obtuse_dortmund"]) != 0 and int(row_result["obtuse_dortmund"]) == 0:
        num_non_obtuse = num_non_obtuse + 1
        num_better = num_better + 1
    elif (int(row["obtuse_dortmund"]) == 0 and int(row_result["obtuse_dortmund"]) == 0) and (int(row["steiner_dortmund"]) > int(row_result["steiner_dortmund"])):
        num_better = num_better + 1
    elif (int(row["obtuse_dortmund"]) != 0 and int(row_result["obtuse_dortmund"]) != 0) and (int(row["obtuse_dortmund"]) > int(row_result["obtuse_dortmund"])):
        num_better = num_better + 1

    if int(row["obtuse_dortmund"]) == 0:
        steiner_dortmund = int(row["steiner_dortmund"])
        steiner_best = int(row["steiner_best"])
        score = score + 0.5 + steiner_best/(2*steiner_dortmund)
    else:
        obtuse_dortmund = int(row["obtuse_dortmund"])
        score = score + 0.5*math.pow(0.97, obtuse_dortmund)

print(score)

print("Num additional non-obtuse triangulations: %n", num_non_obtuse)
print("Num better triangulations: %n", num_better)