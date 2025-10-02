import ROOT
import numpy as np
from tqdm import tqdm
import time
from array import array
import copy
import logging

RDF = ROOT.ROOT.RDataFrame

class Get_Weights_Per_Bin :
    
    """
    Class to compute weights per bin of the variable(s) being reweighted.
    1D or 2D reweighting.
    Order = 0, hist_kpipi / hist_ksk, weights for ksk
    Order = 1, hist_ksk / hist_kpipi, weights for kpipi
    Returns : 1 weight value per bin
    """
    
    def __init__(self, ds_kpipi,ds_ksk, kpipi_indices, ksk_indices, weights_kpipi, weights_ksk, all_variables, variables, bin_edges, order, iter):
        self.ds_kpipi = ds_kpipi
        self.ds_ksk = ds_ksk
        self.kpipi_indices = kpipi_indices
        self.ksk_indices = ksk_indices
        self.weights_kpipi_init = weights_kpipi
        self.weights_ksk_init = weights_ksk
        self.all_variables = all_variables
        self.variables = variables
        self.bin_edges = bin_edges
        self.order = order
        self.iter = iter
        if len(variables) == 1:
            print(f"Order : {self.order}")
            print(f"Performing 1D reweighting for variable: {variables[0]} and bins = {bin_edges[self.all_variables.index(variables[0])]}")
            self.weights = self.compute_1d_weights(variables[0], bin_edges[self.all_variables.index(variables[0])])
        elif len(variables) == 2:
            print(f"Order : {self.order}")
            print(f"Performing 2D reweighting for variables: {variables} and bins = {bin_edges[self.all_variables.index(variables[0])]}, {bin_edges[self.all_variables.index(variables[1])]}")
            self.weights = self.compute_2d_weights(variables[0], variables[1], bin_edges[self.all_variables.index(variables[0])], bin_edges[self.all_variables.index(variables[1])])

    def compute_1d_weights(self, var_name, bin_edge):
        self.var_name = var_name
        self.bin_edge = bin_edge
        """Compute the weights based on histograms of the variable."""
        hist_kpipi = ROOT.TH1F(f"kpp_{self.var_name}", "hist1", len(self.bin_edge)-1, np.array(self.bin_edge))
        hist_ksk = ROOT.TH1F(f"ksk_{self.var_name}", "hist2", len(self.bin_edge)-1, np.array(self.bin_edge))
        
        self.weights_kpipi = copy.deepcopy(self.weights_kpipi_init)
        self.weights_ksk = copy.deepcopy(self.weights_ksk_init)
        
        # Fill the histograms with entries and weights
        # weights = previously computed weights (or 1 if it's the first iteration) * sweights
        
        for i in range(len(self.ds_kpipi[var_name])):
            entry = self.ds_kpipi[var_name][i]
            sweight = self.ds_kpipi["sWeight"][i]
            hist_kpipi.Fill(entry, self.weights_kpipi[self.kpipi_indices[i]] * sweight)

        for i in range(len(self.ds_ksk[var_name])):
            entry = self.ds_ksk[var_name][i]
            sweight = self.ds_ksk["sWeight"][i]
            hist_ksk.Fill(entry, self.weights_ksk[self.ksk_indices[i]] * sweight)

        hist_kpipi.Scale(1.0 / hist_kpipi.Integral())
        hist_ksk.Scale(1.0 / hist_ksk.Integral())  
        
        if self.order == 0 :
            weights = np.ones(hist_kpipi.GetNbinsX())
            for bin in range(1,hist_kpipi.GetNbinsX()+1):
                if hist_ksk.GetBinContent(bin) != 0:
                    ratio = hist_kpipi.GetBinContent(bin) / hist_ksk.GetBinContent(bin)
                else:
                    ratio = 1.0
                weights[bin-1] = ratio
        if self.order == 1 :
            weights = np.ones(hist_ksk.GetNbinsX())
            for bin in range(1,hist_ksk.GetNbinsX()+1):
                if hist_kpipi.GetBinContent(bin) != 0:
                    ratio = hist_ksk.GetBinContent(bin) / hist_kpipi.GetBinContent(bin)
                else:
                    ratio = 1.0
                weights[bin-1] = ratio
                
        return weights

    def compute_2d_weights(self, var1_name, var2_name, bin1_edge, bin2_edge):
        self.var1_name = var1_name
        self.var2_name = var2_name
        self.bin1_edge = bin1_edge
        self.bin2_edge = bin2_edge
        """Compute the weights based on histograms of the variable."""
        hist_kpipi = ROOT.TH2F(f"kpp_{self.var1_name}_{self.var2_name}", "", len(self.bin1_edge)-1, np.array(self.bin1_edge),len(self.bin2_edge)-1, np.array(self.bin2_edge))
        hist_ksk = ROOT.TH2F(f"ksk_{self.var1_name}_{self.var2_name}", "", len(self.bin1_edge)-1, np.array(self.bin1_edge),len(self.bin2_edge)-1, np.array(self.bin2_edge))
        
        self.weights_kpipi = copy.deepcopy(self.weights_kpipi_init)
        self.weights_ksk = copy.deepcopy(self.weights_ksk_init)
        
        # Fill the histograms with entries and weights
        # weights = previously computed weights (or 1 if it's the first iteration) * sweights

        for i in tqdm(range(len(self.ds_kpipi[self.var1_name]))):
            entry1 = self.ds_kpipi[self.var1_name][i]
            entry2 = self.ds_kpipi[self.var2_name][i]
            if i == 0 : print("sweights applied to kpipi")
            sweight = self.ds_kpipi["sWeight"][i]
            hist_kpipi.Fill(entry1, entry2, self.weights_kpipi[self.kpipi_indices[i]] * sweight)

        for i in tqdm(range(len(self.ds_ksk[self.var1_name]))):
            entry1 = self.ds_ksk[self.var1_name][i]
            entry2 = self.ds_ksk[self.var2_name][i] 
            if i == 0 : print("sweights applied to ksk")
            sweight = self.ds_ksk["sWeight"][i]
            hist_ksk.Fill(entry1, entry2, self.weights_ksk[self.ksk_indices[i]] * sweight)
        
        hist_kpipi.Scale(1.0 / hist_kpipi.Integral())
        hist_ksk.Scale(1.0 / hist_ksk.Integral())
        
        if self.order == 0 :
            weights = [[1.0 for _ in range(hist_kpipi.GetNbinsY())] for _ in range(hist_kpipi.GetNbinsX())]
            for xbin in range(1,hist_kpipi.GetNbinsX()+1) :
                for ybin in range(1,hist_kpipi.GetNbinsY()+1) :
                    if hist_ksk.GetBinContent(xbin,ybin) != 0 :
                        ratio = hist_kpipi.GetBinContent(xbin,ybin) / hist_ksk.GetBinContent(xbin,ybin)
                    else : 
                        ratio = 1.0
                    weights[xbin-1][ybin-1] = ratio
        if self.order == 1 : 
            weights = [[1.0 for _ in range(hist_ksk.GetNbinsY())] for _ in range(hist_ksk.GetNbinsX())]
            for xbin in range(1,hist_ksk.GetNbinsX()+1) :
                for ybin in range(1,hist_ksk.GetNbinsY()+1) :
                    if hist_kpipi.GetBinContent(xbin,ybin) != 0 :
                        ratio = hist_ksk.GetBinContent(xbin,ybin) / hist_kpipi.GetBinContent(xbin,ybin)
                    else : 
                        ratio = 1.0
                    weights[xbin-1][ybin-1] = ratio
                    
        return weights

    def get_weights(self):
        if hasattr(self, 'weights'):
            return self.weights
        else:
            raise ValueError("Weights have not been computed yet.")
        