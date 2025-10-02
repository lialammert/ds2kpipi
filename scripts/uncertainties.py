import os
import json
import ast
import numpy as np

yields = {}
araw = {}
araw_sig = {}

number_of_w = 9
strat = "Vega_sWeights"

labels = ["KPIPI"] + [f"KSK_{i}" for i in range(number_of_w)] + [f"Ds_{i}" for i in range(number_of_w)]

with open(f"./{strat}/0_iter_kpipi/weighted_final_results.txt", "r") as f:
    lines = f.readlines()
yields["KPIPI"] = ast.literal_eval(lines[1].strip())
araw["KPIPI"] = ast.literal_eval(lines[2].strip())
araw_sig["KPIPI"] = ast.literal_eval(lines[3].strip())

for w in range(number_of_w):      
    with open(f"./{strat}/{w}_iter/weighted_final_results.txt", "r") as f:
        lines = f.readlines()
    label = f"KSK_{w}"
    yields[label] = ast.literal_eval(lines[1].strip())
    araw[label] = ast.literal_eval(lines[2].strip())   
    araw_sig[label] = ast.literal_eval(lines[3].strip())
    
    with open(f"./Ds_prod_asym/{strat}/{w}_iter/Ds_Prod_Asym.txt", "r") as f:
        lines = f.readlines()
    label = f"Ds_{w}"
    # if w == 0 :
    #     # araw[label] = [ast.literal_eval(lines[-2].strip())[0]] * 20
    #     # araw_sig[label] = [ast.literal_eval(lines[-1].strip())[0]] * 20
    #     araw[label] = ast.literal_eval(lines[-2].strip())[1:]  
    #     araw_sig[label] = ast.literal_eval(lines[-1].strip())[1:]
    # else:
    araw[label] = ast.literal_eval(lines[-2].strip())   
    araw_sig[label] = ast.literal_eval(lines[-1].strip())

full_uncertainty = {}
internal_uncertainty = {}
difference = [{} for _ in range(number_of_w)]

hi = np.sqrt((araw_sig["KPIPI"][5])**2 + (araw_sig["KSK_1"][5])**2)
print(hi)
hii = np.sqrt((araw_sig["KPIPI"][5])**2 + (araw_sig["KSK_1"][5])**2 +  araw_sig[f"Ds_1"][5]**2) #araw_sig[f"KPIPI_Ds"][5]**2 +
print(hii)

for w in range(number_of_w):
    for i in range(20):
        internal_uncertainty[f"{w} - Bin {i}"] = np.sqrt((araw_sig["KPIPI"][i])**2+(araw_sig[f"KSK_{w}"][i])**2)
        full_uncertainty[f"{w} - Bin {i}"] = np.sqrt((araw_sig["KPIPI"][i])**2 + (araw_sig[f"KSK_{w}"][i])**2 + araw_sig[f"Ds_{w}"][i]**2) # + araw_sig[f"KPIPI_Ds"][i]**2
        difference[w][i] = full_uncertainty[f"{w} - Bin {i}"] - internal_uncertainty[f"{w} - Bin {i}"]
    
print(internal_uncertainty)
print(full_uncertainty)
print(difference)

with open("uncertainty_Vega_sWeights.json", "w") as f:
    json.dump({
        # "internal_uncertainty": internal_uncertainty,
        # "full_uncertainty": full_uncertainty,
        "difference": difference
    }, f, indent=4)
    # json.dump(difference, f,  indent=4)
    

