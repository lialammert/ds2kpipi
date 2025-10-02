import ROOT
import numpy as np
from tqdm import tqdm
import time
from array import array
import copy

RDF = ROOT.ROOT.RDataFrame


class Create_Datasets :
    
    """
    Class to create datasets from ROOT files.
    Defines kinematic variables in datasets.
    Creates weight tuples, initialised to 1 everywhere.
    """
    
    def __init__(self, kpipi_file, ksk_file, kinematic_vars, weight_vars, all_variables, formula_kpipi, formula_ksk, bin_edges_weigh, variable_ranges, bins_plot, bin_edges_plot, bins_weigh, weighting_strategy):
        self.kpipi_file = kpipi_file
        self.ksk_file = ksk_file
        self.kinematic_vars = kinematic_vars
        self.weight_vars = weight_vars
        self.all_variables = all_variables
        self.formula_kpipi = formula_kpipi
        self.formula_ksk = formula_ksk
        self.bin_edges_weigh = bin_edges_weigh
        self.variable_ranges = variable_ranges
        self.bins_plot = bins_plot
        self.bin_edges_plot = bin_edges_plot
        self.bins_weigh = bins_weigh
        self.weighting_strategy = weighting_strategy
        
    def create_datasets(self):
        self.tree_kpipi = self.kpipi_file.Get("DecayTree")
        self.tree_ksk = self.ksk_file.Get("DecayTree")
        
        #  creating dataframes then datasets with all kinematic variables of interest + sWeights
        self.df_kpipi = RDF(self.tree_kpipi)
        self.df_ksk = RDF(self.tree_ksk)
        for k,kin in enumerate(self.kinematic_vars) :
            if kin in self.df_kpipi.GetColumnNames() : self.df_kpipi = self.df_kpipi.Redefine(kin, self.formula_kpipi[k])
            else : self.df_kpipi = self.df_kpipi.Define(kin, self.formula_kpipi[k])
            if kin in self.df_ksk.GetColumnNames() : self.df_ksk = self.df_ksk.Redefine(kin, self.formula_ksk[k])
            else : self.df_ksk = self.df_ksk.Define(kin, self.formula_ksk[k])
        
        self.ds_kpipi = self.df_kpipi.AsNumpy(columns=self.all_variables+["sWeight"])
        self.ds_ksk = self.df_ksk.AsNumpy(columns=self.all_variables+["sWeight"])
        
        # initialise weights to 1 everywhere
        self.weights_kpipi_final = np.ones(len(self.ds_kpipi[self.all_variables[0]])) 
        self.weights_ksk_final = np.ones(len(self.ds_ksk[self.all_variables[0]])) 
        
        print("Before cut")
        print("Number of entries in ds_kpipi:", len(self.ds_kpipi[self.all_variables[0]]))
        print("Number of entries in ds_ksk:", len(self.ds_ksk[self.all_variables[0]]))
        print("sum of weights ksk : ", np.sum(self.weights_ksk_final))
        print(self.weights_ksk_final[:10])
        
        # we want the weights on all values so we cut after creating the weight tuples
        self.perform_cuts() # apply (mass cut) and Hlt1 cuts
        
        print("After cut")
        print("Number of entries in ds_kpipi:", len(self.ds_kpipi[self.all_variables[0]]))
        print("Number of entries in ds_ksk:", len(self.ds_ksk[self.all_variables[0]]))
        print("sum of weights ksk : ", np.sum(self.weights_ksk_final))
        print(self.weights_ksk_final[:10])
        
        return self.ds_kpipi, self.ds_ksk, self.indices_kpipi, self.indices_ksk, self.weights_kpipi_final, self.weights_ksk_final
    
    def get_bin_edges(self):
        for ikk in range (len(self.kinematic_vars), len(self.all_variables)) :
            self.variable_ranges.append([min(self.df_kpipi.Min(self.all_variables[ikk]).GetValue(), self.df_ksk.Min(self.all_variables[ikk]).GetValue()), max(self.df_kpipi.Max(self.all_variables[ikk]).GetValue(), self.df_ksk.Max(self.all_variables[ikk]).GetValue())])
            self.bin_edges_plot[ikk]=np.array(np.linspace(self.variable_ranges[ikk][0], self.variable_ranges[ikk][1], self.bins_plot+1))
            self.bin_edges_weigh[ikk]=np.array(np.linspace(self.variable_ranges[ikk][0], self.variable_ranges[ikk][1], self.bins_weigh+1))
        
        return self.bin_edges_weigh, self.bin_edges_plot 
    
    def perform_cuts(self):
        # self.mass_cut = "(Dp_M > 1920) && (Dp_M < 2025)" # don't apply since using sWeights
        # Hlt1 cuts
        kpipi_cut = " (Kp_Hlt1TrackMVADecision_TOS)"
        ksk_cut = "(hp_Hlt1TrackMVADecision_TOS)"
        self.ds_kpipi = self.df_kpipi.Filter(kpipi_cut).AsNumpy(columns=self.all_variables+["sWeight"])
        self.ds_ksk = self.df_ksk.Filter(ksk_cut).AsNumpy(columns=self.all_variables+["sWeight"])
    
        # Get indices
        self.indices_kpipi = self.df_kpipi.Filter(kpipi_cut).Define("entry_index_kpipi", "rdfentry_").AsNumpy()["entry_index_kpipi"]
        self.indices_ksk = self.df_ksk.Filter(ksk_cut).Define("entry_index_ksk", "rdfentry_").AsNumpy()["entry_index_ksk"]