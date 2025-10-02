import ROOT
import numpy as np
import math
from tqdm import tqdm
import time
from array import array
import os
import copy
from .weighting_strategies import STRATEGIES
from .step1_create_datasets import Create_Datasets
from .step2_weights_per_bin import Get_Weights_Per_Bin
from .step3_weights_per_val import Get_Weights_Per_Val
from .step4_save_weights import Save_Weights
from .step5_plot_histograms import Plot_Histograms

RDF = ROOT.ROOT.RDataFrame

"""
ROOT data files are taken from folder `Dalitz_Bins_Kp_TOS/` -> to change in the code
weights are saved in a folder with the name of the strategy, for each iteration
"""

total_start_time = time.time()
start_time = time.time()

def get_argparser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required="0")
    parser.add_argument('-n', '--bin_number', required="0")
    parser.add_argument('-s', '--strategy', required="0")
    parser.add_argument('-na', '--name', required="0")
    parser.add_argument('-mi', '--max_iter', required="0", type=int)
    parser.add_argument('-fc', '--first_coarse', default=False)
    return parser

options = get_argparser().parse_args()
kpipi_file = ROOT.TFile("Dalitz_Bins_Kp_TOS/"+options.input_file)
bin_number = int(options.bin_number)
ksk_file = ROOT.TFile("Dalitz_Bins_Kp_TOS/fcut_hlt1_kp_sweights.root")

tos = "kp"

''' can give binning as an argument, as [bin_edges_weigh, bin_edges_plot] '''

weighting_strategy = STRATEGIES[options.strategy](name=options.name, max_iter=options.max_iter, first_coarse=options.first_coarse)  # Use the strategy from the command line argument
first_coarse = weighting_strategy.first_coarse  # If True, use coarse binning for the first iteration
vars_to_reweigh = weighting_strategy.get_variables() * (weighting_strategy.max_iter)

if weighting_strategy.max_iter == 0:
    max_iteration = 0
    ordering = weighting_strategy.order
else :
    max_iteration = len(weighting_strategy.get_variables())*(weighting_strategy.max_iter)
    if first_coarse : 
        max_iteration += 1
        vars_to_reweigh.insert(0, vars_to_reweigh[0])
    if isinstance(weighting_strategy.order, int) : ordering = np.repeat(weighting_strategy.order,max_iteration)  # Create an array with the same order for all iterations
    else : ordering = weighting_strategy.order * (weighting_strategy.max_iter)
    if first_coarse :
        ordering.insert(0, ordering[0])
        # if all(x == weighting_strategy.order[0] for x in ordering): different_orders = False
        # else : different_orders = True  # If the order is a list, it means different orders for each iteration
weight_vars_flat = [x for pair in weighting_strategy.get_variables() for x in pair]
if first_coarse: 
    weight_vars_flat.insert(0, weight_vars_flat[1])
    weight_vars_flat.insert(0, weight_vars_flat[1])  # Add the first two variables for the coarse binning

print(f"total number of iterations : {max_iteration} (+1 if coarse 1st step)")
print(f"ordering : {ordering}")
print(f"reweighting scheme : {vars_to_reweigh}")

kinematic_vars, formula_vars_tree_kpipi, formula_vars_tree_ksk, variable_ranges, bin_edges_weigh, bin_edges_plot = weighting_strategy.binning  # Get variable ranges and bin edges from the weighting strategy
if first_coarse : variable_ranges_coarse, bin_edges_weigh_coarse, bin_edges_plot_coarse = weighting_strategy.coarse_binning  # Get coarse binning for weights
bins_weigh, bins_plot = len(bin_edges_weigh[0])-1, len(bin_edges_plot[0])-1

all_variables = copy.deepcopy(kinematic_vars)
all_variables += [x for x in weight_vars_flat if x not in kinematic_vars]  # Ensure no duplicates

start_time = time.time()
create_datasets = Create_Datasets(kpipi_file, ksk_file, kinematic_vars, weight_vars_flat, all_variables, formula_vars_tree_kpipi, formula_vars_tree_ksk, bin_edges_weigh, variable_ranges, bins_plot, bin_edges_plot, bins_weigh, weighting_strategy)
ds_kpipi, ds_ksk, kpipi_indices, ksk_indices, weights_kpipi, weights_ksk = create_datasets.create_datasets()
bin_edges_weigh, bin_edges_plot = create_datasets.get_bin_edges()

end_time = time.time()  
elapsed_time = end_time - start_time
minutes = int(elapsed_time // 60)
seconds = int(elapsed_time % 60)
print(f"Datasets created successfully in {minutes} min {seconds} sec")

if max_iteration == 0:
    folder = weighting_strategy.name+"/0_iter"
    os.system(f'mkdir -p ./{folder}/weighted_plots')
    save_weights = Save_Weights(weights_kpipi, weights_ksk, bin_number, folder, 0)
    save_weights.save_weights()
    histograms = Plot_Histograms(ds_kpipi, ds_ksk, kpipi_indices, ksk_indices, weights_kpipi, weights_ksk, all_variables, bin_edges_plot, bins_plot, ordering[0], folder, bin_number, max_iteration, first_coarse,0)
    histograms.plot_variable()

for i in range(max_iteration+1):
    start_time = time.time()
    print("--------------------------")
    n = len(weighting_strategy.get_variables())
    folder = weighting_strategy.name+f"/{i}_iter"
    
    if i == 0 :
        os.system(f'mkdir -p ./{folder}/weighted_plots')
        print("Saving the unweighted histograms.")
        save_weights = Save_Weights(weights_kpipi, weights_ksk, bin_number, folder, i)
        save_weights.save_weights()
        histograms = Plot_Histograms(ds_kpipi, ds_ksk, kpipi_indices, ksk_indices, weights_kpipi, weights_ksk, all_variables, bin_edges_plot, bins_plot, ordering[0], folder, bin_number, max_iteration, first_coarse,i)
        histograms.plot_variable()
        if first_coarse :
            print(f"First step is coarse reweighting.")
            weights_per_bin = Get_Weights_Per_Bin(ds_kpipi, ds_ksk, kpipi_indices, ksk_indices, weights_kpipi, weights_ksk, all_variables, vars_to_reweigh[0], bin_edges_weigh_coarse, 1, i) 
            weights_per_bin = weights_per_bin.get_weights()
            weights_per_val = Get_Weights_Per_Val(ds_kpipi, ds_ksk, kpipi_indices, ksk_indices, weights_kpipi, weights_ksk,all_variables, vars_to_reweigh[0], bin_edges_weigh_coarse, 1, max_iteration,1, weights_per_bin)
            weights_kpipi, weights_ksk = weights_per_val.get_weights()
        continue
    print("i only print this if i > 0")
        # if not first_coarse : continue

    round_of_w = int((i)/n)
    if i%n == 0: print(f"Round {round_of_w} of reweighting.")
    variables = vars_to_reweigh[i-1]
    print(f"Iteration number : {i}, variables being reweighted : {variables}")
    weights_per_bin = Get_Weights_Per_Bin(ds_kpipi, ds_ksk, kpipi_indices, ksk_indices, weights_kpipi, weights_ksk, all_variables, variables, bin_edges_weigh, ordering[i-1], i)
    weights_per_bin = weights_per_bin.get_weights()
    end_time = time.time()  
    elapsed_time = end_time - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"Weights per bin calculated successfully in {minutes} min {seconds} sec")

    start_time = time.time()
    weights_per_val = Get_Weights_Per_Val(ds_kpipi, ds_ksk, kpipi_indices, ksk_indices, weights_kpipi, weights_ksk,all_variables, variables, bin_edges_weigh, i+1, max_iteration, ordering[i-1], weights_per_bin)
    weights_kpipi, weights_ksk = weights_per_val.get_weights()
    end_time = time.time()  
    elapsed_time = end_time - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"Weights per value calculated successfully for iteration {i+1} in {minutes} min {seconds} sec.")
    
    # if i+1 == max_iteration :
    if i%n == 0:
        folder = weighting_strategy.name+f"/{round_of_w}_iter"
        os.system(f'mkdir -p ./{folder}/weighted_plots')
        print(f"Round {round_of_w} of reweighting completed.")
        start_time = time.time()
        save_weights = Save_Weights(weights_kpipi, weights_ksk, bin_number, folder, round_of_w)
        save_weights.save_weights()
        end_time = time.time()  
        elapsed_time = end_time - start_time
        minutes = int(elapsed_time // 60)
        seconds = int(elapsed_time % 60)
        print(f"Weights saved successfully in {minutes} min and {seconds} sec.")
        
        start_time = time.time()
        histograms = Plot_Histograms(ds_kpipi, ds_ksk, kpipi_indices, ksk_indices, weights_kpipi, weights_ksk, all_variables, bin_edges_plot, bins_plot, ordering[i-1], folder, bin_number, max_iteration, first_coarse, round_of_w)
        histograms.plot_variable()
        end_time = time.time()
        elapsed_time = end_time - start_time
        minutes = int(elapsed_time // 60)
        seconds = int(elapsed_time % 60)
        print(f"Histograms plotted successfully in {minutes} min and {seconds} sec.")
        
total_end_time = time.time()
total_elapsed_time = total_end_time - total_start_time
total_minutes = int(total_elapsed_time // 60)
total_seconds = int(total_elapsed_time % 60)
print(f"Total runtime: {total_minutes} min {total_seconds} sec")
print("All steps completed successfully.")