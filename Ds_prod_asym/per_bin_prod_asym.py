import ROOT
from ROOT import TMath
import numpy as np
from array import array
from tqdm import tqdm
import os
import sys
RDF = ROOT.ROOT.RDataFrame

"""
Function to compute the production asymmetry per bin.
This function processes the KSK and Kpipi files, computes the 2D histograms for Dp_PT and Dp_ETA,
and saves the results in the specified folder structure.
"""

def compute_per_bin_prod_asym(kpipi_files, ksk_files, kpipi_weights_file, ksk_weights_file):
    ksk_file = ksk_files[0]  # Assuming a single KSK file for all bins
    results={}
    strat = "Vega_sWeights"  # Change this to the desired strategy
    for w in range(2):
        if w == 0:
            files = [ksk_file] + kpipi_files
            weights_file = kpipi_weights_file
        else : 
            print(len(ksk_files))
            files = ksk_files
            weights_file = ksk_weights_file

        bin_edges_pt = [2.5e3,4.7e3,6.5e3,8.5e3,25e3]
        bin_edges_eta = [2.0,3.0,3.5,4.0,4.5]
        
        hist_2d_pt_eta= {}
        hist_2d_weighted = {}

        for ff,fil in enumerate(files) : 
            sys.stdout.flush()
            print(f"Processing file: {fil}")
            file = ROOT.TFile.Open(fil)
            tree = file.Get('DecayTree')
            data_frame = RDF(tree)
            vars = ["Dp_PT", "Dp_P", "Dp_ETA"]
            formulas = ["TMath::Sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY)",
                    "TMath::Sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY + Dp_PZ*Dp_PZ)",
                    "0.5*TMath::Log((Dp_P + Dp_PZ)/(Dp_P - Dp_PZ))"]
            for e, v in enumerate(vars):
                if v not in data_frame.GetColumnNames():
                    data_frame = data_frame.Define(v, formulas[e])
                else : 
                    data_frame = data_frame.Redefine(v, formulas[e])

            values = []
            weight_file = ROOT.TFile.Open(weights_file)

            if w == 0:
                hist_2d_pt_eta[ff] = ROOT.TH2D("pt_eta","pt_eta", len(bin_edges_pt)-1,  np.array(bin_edges_pt), len(bin_edges_eta)-1,  np.array(bin_edges_eta))
                data_set = data_frame.AsNumpy(columns=["Dp_M", "Dp_PT", "Dp_ETA","sWeight"])
                for i in range(len(data_set["Dp_PT"])):
                    hist_2d_pt_eta[ff].Fill(data_set['Dp_PT'][i], data_set['Dp_ETA'][i], data_set['sWeight'][i])
                
                integral = hist_2d_pt_eta[ff].Integral()
                if integral > 0:
                    hist_2d_pt_eta[ff].Scale(1. / integral)

                for ix in range(1, hist_2d_pt_eta[ff].GetNbinsX() + 1):
                    for iy in range(1, hist_2d_pt_eta[ff].GetNbinsY() + 1):
                        value = hist_2d_pt_eta[ff].GetBinContent(ix, iy)
                        values.append(round(value,5))
                        
                print("writing values to file")
                if ff == 0:
                    results["KsK_0w"] = values
                else : 
                    results[f"Kpipi_Bin{ff-1}"] = values
                    
            else : 
                weights_tree = weight_file.Get(f"weights_ksk_tree_{ff}")
                data_frame_w = RDF(weights_tree)
                print(f"filling histogram with weights from bin {ff}")
                data_set = data_frame.AsNumpy(columns=["Dp_M", "Dp_PT", "Dp_ETA","sWeight"])
                data_set_w = data_frame_w.AsNumpy(columns=[f"weights_ksk_{ff}"])

                hist_2d_weighted[ff] = ROOT.TH2D(f"KSK_{ff}",f"KSK_{ff}", len(bin_edges_pt)-1,  np.array(bin_edges_pt), len(bin_edges_eta)-1,  np.array(bin_edges_eta))
                eta = data_set['Dp_ETA']
                pt = data_set['Dp_PT']
                weight = data_set_w[f"weights_ksk_{ff}"]
                for i in range(len(pt)):
                    hist_2d_weighted[ff].Fill(pt[i], eta[i], weight[i]* data_set['sWeight'][i])

                integral = hist_2d_weighted[ff].Integral()
                if integral > 0:
                    hist_2d_weighted[ff].Scale(1. / integral)
                    
                for ix in range(1, hist_2d_weighted[ff].GetNbinsX() + 1):
                    for iy in range(1, hist_2d_weighted[ff].GetNbinsY() + 1):
                        value = hist_2d_weighted[ff].GetBinContent(ix, iy)
                        values.append(round(value,4))
                        
                results[f"{w}w-KsK_Bin{ff}"] = values

        file.Close()
    return results