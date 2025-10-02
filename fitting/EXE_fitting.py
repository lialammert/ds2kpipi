import ROOT
from ROOT import RooFit as RF
from ROOT import TMath
import numpy as np
from fitting import Exponential, ExecuteFit, CrystalBall_Gaussian
from simultaneous_fitting import SimultaneousFitter
from weighting.weighting_strategies import STRATEGIES
from array import array
import json
import time
import sys
import os
import argparse
from tqdm import tqdm
import re
ROOT.gROOT.ProcessLine(".L extras/lhcbStyle.C");
RDF = ROOT.ROOT.RDataFrame

"""
script that calls the SimultaneousFitter class for fitting
parameters to give : `root_files_directory`, range for the mass `mass_range`, number of bins `bins`
raw asymmetries are saved in `a_raw.json` and `weighted_ksk_fit_results.txt` files.
"""

start_time = time.time()

def get_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--strategy', required="0")
    parser.add_argument('-na', '--name', required="0")
    parser.add_argument('-mi', '--max_iter', required="0", type=int)
    parser.add_argument('-d', '--decay', required="0")
    parser.add_argument('-b', '--bin', required="0")
    return parser

options = get_argparser().parse_args()

root_files_directory = "Dalitz_Bins_Kp_TOS"

weighting_strategy = STRATEGIES[options.strategy](name=options.name, max_iter=options.max_iter)

bin = int(options.bin)

if bin != -1 : folder = f"{weighting_strategy.name}/BIN{bin}"
else : folder = f"{weighting_strategy.name}/KSK"
os.system(f'mkdir -p ./{folder}')

fit_results_file = f"./{folder}/weighted_ksk_fit_results.txt"
final_results_file = f"./{folder}/a_raw.json"
output_files = [fit_results_file, final_results_file]
for file in output_files:
    with open(file, "w") as f:
        pass

if bin == -1 :
    range_of_w = [0]
    plot_file = "ksk"
else :
    range_of_w = list(range(1, weighting_strategy.max_iter+1))
    plot_file = f"Bin{bin}"
    
results_dict = {}  
    
canvas = ROOT.TCanvas()
canvas.Print(f'./{folder}/{plot_file}.pdf[')

for w in range_of_w:
    decay = options.decay
    param_dict = {}

    kpipi_input_files = [f"dalitz_bin_{i}_hlt1_kp.root" for i in range(20)]
    ksk_input_file = "fcut_hlt1_kp_sweights.root"
    ksk_weights_files = f"./{weighting_strategy.name}/{w}_iter/weights_ksk_{w}.root"
    kpipi_weights_files = f"./{weighting_strategy.name}/{w}_iter/weights_kpipi_{w}.root"

    mass_range = [1920,2025]
    bins = 200
        
    final = True # True to save initial parameters for fits in the .json file

    if w == 0 and decay == "ksk" :
        print(f"#### FITTING ENTIRE KSK SAMPLE WITH NO WEIGHTING #### \n")
        simultaneous_fitter = SimultaneousFitter(decay, 0, bins, folder, 
                                                 root_files_directory, canvas, output_files, plot_file, results_dict,
                                                 w, mass_range, 
                                                 kpipi_input_files[0], ksk_input_file, 
                                                 ksk_weights_files, kpipi_weights_files)
        simultaneous_fitter.fit()
    else : 
        print(f"#### FITTING BIN {bin}, WEIGHTED {w} TIMES #### \n")
        simultaneous_fitter = SimultaneousFitter(decay,  
                                                 bin, bins, folder, 
                                                 root_files_directory, canvas, output_files, plot_file, results_dict,
                                                 w, mass_range, 
                                                 kpipi_input_files[bin], ksk_input_file, 
                                                 ksk_weights_files, kpipi_weights_files)
        simultaneous_fitter.fit()

    # clean_yields = [int(x) for x in total_yields]

    end_time = time.time()  
    elapsed_time = end_time - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"All fits performed in {minutes} min {seconds} sec")
    print("\n")
    print("\n")
    print("========================================================== \n")
    print("\n")
    print("\n")

canvas.Print(f'./{folder}/{plot_file}.pdf]')

with open(final_results_file, "a") as file:
    json_str = json.dumps(results_dict, indent=2)
    json_str = re.sub(
        r'\[\s+([^\[\]]+?)\s+\]',
        lambda m: "[" + " ".join(m.group(1).strip().split()) + "]",
        json_str
    )
    file.write(json_str)