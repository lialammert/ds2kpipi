import ROOT
import numpy as np
from tqdm import tqdm
import time
from array import array
import json

RDF = ROOT.ROOT.RDataFrame

class Save_Weights:
    
    """
    Class to save the final weights for kpp and ksk in ROOT files.
    Individual files for each bin.
    """
    
    def __init__(self, weights_kpipi_final, weights_ksk_final, bin_number, folder, iter):
        self.weights_kpipi_final = weights_kpipi_final
        self.weights_ksk_final = weights_ksk_final
        self.bin_number = bin_number
        self.folder = folder
        self.iter = iter
    
    def save_weights (self) :
        print(f"SAVING WEIGHTS to folder {self.folder}")
        file1 = ROOT.TFile(f"./weights_kpipi_{self.iter}/weights_kpipi_{self.bin_number}.root", "UPDATE")
        weight_kpipi = array("f", [0])
        new_tree_kpipi = ROOT.TTree(f"weights_kpipi_tree_{self.bin_number}", "Tree with kpp weights")
        new_tree_kpipi.Branch(f"weights_kpipi_{self.bin_number}", weight_kpipi, f"weights_kpipi_{self.bin_number}/F")
        for w in self.weights_kpipi_final:
            weight_kpipi[0] = w
            new_tree_kpipi.Fill()
        new_tree_kpipi.Write()
        file1.Close()
        
        file2 = ROOT.TFile(f"./weights_ksk_{self.iter}/weights_ksk_{self.bin_number}.root", "UPDATE")
        weight_ksk = array("f", [0])
        new_tree_ksk = ROOT.TTree(f"weights_ksk_tree_{self.bin_number}", "Tree with ksk weights")
        new_tree_ksk.Branch(f"weights_ksk_{self.bin_number}", weight_ksk, f"weights_ksk_{self.bin_number}/F")
        for w in self.weights_ksk_final:
            weight_ksk[0] = w
            new_tree_ksk.Fill()
        new_tree_ksk.Write()
        file2.Close()
        