import ROOT
import numpy as np
from tqdm import tqdm
import time
from array import array
import copy
import logging

RDF = ROOT.ROOT.RDataFrame

class Get_Weights_Per_Val :
    
    """
    Class to compute weights per value in the ntuple of the variable(s) being reweighted.
    Loops over dataset and finds ntuple index for each weight (from weight_per_bin class).
    Returns : 1 weight value per value in the ntuple.
    """

    def __init__(self, dataset_kpipi, dataset_ksk, indices_kpipi, indices_ksk, weights_kpipi, weights_ksk, all_variables,variables, bin_edges, iteration, max_iteration, order, weights_per_bin):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.ds_kpipi = dataset_kpipi
        self.ds_ksk = dataset_ksk
        self.indices_kpipi = indices_kpipi
        self.indices_ksk = indices_ksk
        self.weights_kpipi_final = weights_kpipi
        self.weights_ksk_final = weights_ksk
        self.all_variables = all_variables
        self.variables = variables
        self.bin_edges = bin_edges
        self.num_iter = iteration
        self.max_iter = max_iteration
        self.order = order
        self.weights_per_bin = weights_per_bin
        if len(variables) == 1:
            self.apply_1d_weights(variables[0], bin_edges[0])
        elif len(variables) == 2:
            self.apply_2d_weights(variables[0], variables[1], bin_edges[self.all_variables.index(variables[0])], bin_edges[self.all_variables.index(variables[1])])
        
    def get_bin_id( self, datapoint, bin1_edge, bin2_edge ):
        if datapoint[0] > bin1_edge[-1]:
            self.logger.warning("Datapoint is overflowing the x-axis! ({},{})".format(datapoint[0], datapoint[1]))
        
        if datapoint[0] < bin1_edge[0]:
            self.logger.warning("Datapoint is underflowing the x-axis! ({},{})".format(datapoint[0], datapoint[1]))

        x_bin_id = np.searchsorted( bin1_edge, datapoint[0] ) 

        if datapoint[1] > bin2_edge[-1]:
            self.logger.warning("Datapoint is overflowing the y-axis! ({},{})".format(datapoint[0], datapoint[1]))
        
        if datapoint[1] < bin2_edge[0]:
            self.logger.warning("Datapoint is underflowing the y-axis! ({},{})".format(datapoint[0], datapoint[1]))

        y_bin_id = np.searchsorted( bin2_edge, datapoint[1])
        return [ x_bin_id - 1, y_bin_id - 1]
    
        
    def apply_1d_weights(self, var_name, bin_edge):
        """Apply the weights to dataframe based on the kinematic variable."""
        self.bin_edge = bin_edge
        self.weighting_var1 = var_name
        self.weights_to_apply = copy.deepcopy(self.weights_per_bin)
        
        if self.order == 0 :
            for i in tqdm(range(len(self.ds_ksk[self.weighting_var1]))):
                xval = self.ds_ksk[self.weighting_var1][i]
                xbin = np.searchsorted(bin_edge, xval)-1
                if xbin < 0 or xbin >= len(self.weights_to_apply):
                    self.logger.warning(f"out-of-range bin index: {xbin} for value {xval}")
                    xbin = max(0, min(xbin, len(self.weights_to_apply) - 1))
                self.weights_ksk_final[self.indices_ksk[i]] *= self.weights_to_apply[xbin]
                    
        if self.order == 1 : 
            for i in tqdm(range(len(self.ds_kpipi[self.weighting_var1]))):
                xval = self.ds_kpipi[self.weighting_var1][i]
                xbin = np.searchsorted(bin_edge, xval)-1
                if xbin < 0 or xbin >= len(self.weights_to_apply):
                    self.logger.warning(f"out-of-range bin index: {xbin} for value {xval}")
                    xbin = max(0, min(xbin, len(self.weights_to_apply) - 1))
                self.weights_kpipi_final[self.indices_kpipi[i]] *= self.weights_to_apply[xbin]
        
        # Normalize the weights
        one_over_norm_factor_ksk = np.sum([w for w in self.weights_ksk_final]) / np.sum([w**2 for w in self.weights_ksk_final])
        one_over_norm_factor_kpipi = np.sum([w for w in self.weights_kpipi_final]) / np.sum([w**2 for w in self.weights_kpipi_final])
        self.weights_ksk_final *= one_over_norm_factor_ksk
        self.weights_kpipi_final *= one_over_norm_factor_kpipi
        
        print("APPLIED THE WEIGHTS")
        print(sum(self.weights_kpipi_final))
        print(sum(self.weights_ksk_final))


    def apply_2d_weights(self, var1_name, var2_name, bin1_edge, bin2_edge):
        """Apply the weights to dataframe based on the kinematic variable."""
        self.bin1_edge = bin1_edge
        self.bin2_edge = bin2_edge
        self.weighting_var1 = var1_name
        self.weighting_var2 = var2_name
        self.weights_to_apply = copy.deepcopy(self.weights_per_bin)
    
        if self.order == 0 :
            for i in tqdm(range(len(self.ds_ksk[self.weighting_var1]))):
            # for index in tqdm((self.indices_kpipi)):
                xval = self.ds_ksk[self.weighting_var1][i]
                yval = self.ds_ksk[self.weighting_var2][i]
                [xbin,ybin] = self.get_bin_id([xval,yval], bin1_edge, bin2_edge)
                if xbin < 0 or xbin >= len(self.weights_to_apply):
                    self.logger.warning(f"out-of-range bin index: {xbin} for value {xval}")
                    xbin = max(0, min(xbin, len(self.weights_to_apply) - 1))
                if ybin < 0 or ybin >= len(self.weights_to_apply[0]):
                    self.logger.warning(f"out-of-range bin index: {ybin} for value {yval}")
                    ybin = max(0, min(ybin, len(self.weights_to_apply[0]) - 1))
                self.weights_ksk_final[self.indices_ksk[i]] *= self.weights_to_apply[xbin][ybin]

        if self.order == 1 : 
            for i in tqdm(range(len(self.ds_kpipi[self.weighting_var1]))):
            # for index in tqdm((self.indices_ksk)):
                xval = self.ds_kpipi[self.weighting_var1][i]
                yval = self.ds_kpipi[self.weighting_var2][i]
                [xbin,ybin] = self.get_bin_id([xval,yval], bin1_edge, bin2_edge)
                if xbin < 0 or xbin >= len(self.weights_to_apply):
                    self.logger.warning(f"out-of-range bin index: {xbin} for value {xval}")
                    xbin = max(0, min(xbin, len(self.weights_to_apply) - 1))
                if ybin < 0 or ybin >= len(self.weights_to_apply[0]):
                    self.logger.warning(f"out-of-range bin index: {ybin} for value {yval}")
                    ybin = max(0, min(ybin, len(self.weights_to_apply[0]) - 1))
                self.weights_kpipi_final[self.indices_kpipi[i]] *= self.weights_to_apply[xbin][ybin]
        
        if self.num_iter == self.max_iter : 
            one_over_norm_factor_ksk = np.sum([w for w in self.weights_ksk_final]) / np.sum([w**2 for w in self.weights_ksk_final])
            one_over_norm_factor_kpipi = np.sum([w for w in self.weights_kpipi_final]) / np.sum([w**2 for w in self.weights_kpipi_final])
            self.weights_ksk_final *= one_over_norm_factor_ksk
            self.weights_kpipi_final *= one_over_norm_factor_kpipi
        
        print("APPLIED THE WEIGHTS")
        print(sum(self.weights_kpipi_final))
        print(sum(self.weights_ksk_final))
    
    def get_weights(self):
        """Return the final weights."""
        return self.weights_kpipi_final, self.weights_ksk_final