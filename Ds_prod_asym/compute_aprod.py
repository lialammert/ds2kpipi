import numpy as np
import os
import json
import ast  # for safely parsing the list strings
from per_bin_prod_asym import compute_per_bin_prod_asym

"""
Using measured values and uncertainties of the Ds production asymmetry in PT_ETA (values to change in the code if needed),
compute the Ds production asymmetry using the entries of the ksk and kpipi PT_ETA histogram
save values as vector 
Give as arguments the list of kpipi files, ksk files, and the weights files
"""

def compute_prod_asym(kpipi_files, ksk_files, kpipi_weights_file, ksk_weights_file):
    prod_asym = compute_per_bin_prod_asym(kpipi_files, ksk_files, kpipi_weights_file, ksk_weights_file)
    results = {}
    
    folder = f"./Ds_prod_asym/{strat}/"
    os.system(f"mkdir -p {folder}")

    with open(f"{folder}/Ds_Prod_Asym.json", "w") as f:
        pass
        
    for w in range(8):
        w = int(w)
        results[f"{w}w"] = {}

        a_prod = [1.19, 0.65, -1.34, 0.86, -0.54, 1.64, 1.00, 0.70, 0.20, -0.56, -0.37, 0, 0.72, -1.81, 0.71, 0]
        stat = [1.17, 0.84, 0.83, 1.41, 0.83, 0.76, 0.88, 1.91, 0.85, 0.87, 1.27, 0, 0.72, 0.97, 2.59, 0]
        syst = [0.35, 0.45, 0.18, 0.27, 0.23, 0.21, 0.13, 0.43, 0.38, 0.15, 0.21, 0, 0.11, 0.18, 0.53, 0]

        diff = {}
        central_values_sum = []
        prod_asym_ksk_unweighted = prod_asym["KsK_0w"]

        for i in range(len(kpipi_files)):
            prod_asym_kpipi = prod_asym[f"Kpipi_Bin{i}"]
            if w == 0 : 
                vector_kpipi = prod_asym_kpipi
                vector_ksk = prod_asym_ksk_unweighted
            else :
                prod_asym_ksk = prod_asym[f"{w}w-KsK_Bin{i}"]
                vector_kpipi = prod_asym_kpipi
                vector_ksk = prod_asym_ksk
            diff[i] = [y - x for x, y in zip(vector_kpipi, vector_ksk)]

        stat_error = []
        syst_error = []

        for i in range(len(diff)):
            label = f"Bin {i}"
            val = sum(x * y for x, y in zip(diff[i], a_prod))
            central_values_sum.append(round(val, 5))
            stat_val = np.sqrt(np.sum((np.array(diff[i])**2) * (np.array(stat)**2)))
            stat_error.append(round(stat_val, 5))
            syst_val = np.sqrt(np.sum((np.array(diff[i])**2) * (np.array(syst)**2)))
            syst_error.append(round(syst_val,5))

        error = np.zeros(len(diff))

        for i in range(len(stat_error)):
            error[i] = round(np.sqrt(stat_error[i]**2 + syst_error[i]**2), 4)

        for i in range(len(central_values_sum)):
            results[f"{w}w"][f"bin_{i}"] = {
                "a_prod": central_values_sum[i],
                "a_prod_error": error[i]
            }

    with open(f"{folder}/Ds_Prod_Asym.json", "w") as f:
        json.dump(results, f)

strat = "Vega_sWeights"
kpipi_files = [f"dalitz_bin_{i}_hlt1_kp.root" for i in range(20)]  # List of Kpipi files
ksk_files = ["fcut_hlt1_kp_sweights.root" for i in range(20)]  # Assuming the same file for all bins, adjust as necessary
kpipi_weights_file = "weights_kpipi.root"  # Path to Kpipi weights file
ksk_weights_file = "weights_ksk.root"  # Path to KSK weights
aprod = compute_prod_asym(kpipi_files, ksk_files, kpipi_weights_file, ksk_weights_file)