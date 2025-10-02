import ROOT
from ROOT import TMath
import sys
import numpy as np
# from fitting import DoubleGaussian, CBGaussian, Cruijff,Exponential, Bernstein, ExecuteFit
from array import array
import json
ROOT.gROOT.ProcessLine(".L extras/lhcbStyle.C");
RDF = ROOT.ROOT.RDataFrame

def get_argparser():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', required="0")
    # parser.add_argument('-i2', '--input_file2', required="0")
    # parser.add_argument('-o', '--output_file', required="0")
    # parser.add_argument('-p', '--plot_file', required="0")

    return parser

options = get_argparser().parse_args()
input_file = options.input_file #'kpipi_small.root' # options.input_file
# input_file2 = options.input_file2
# output_file = options.output_file #'./mass_fit/kpipi_mfit_test.root' # options.output_file
# plot_file = options.plot_file #'./mass_fit/kpipi_mfit_test.pdf' # options.plot_file

# for Dp_M
# mass_range = [1925,2025]
bins = 4
bin_edges_pt = [2.5e3,4.7e3,6.5e3,8.5e3,25e3]
bin_edges_eta = [2.0,3.0,3.5,4.0,4.5]

mass_cut = "(Dp_M > 1920) && (Dp_M < 2025)"
kpipi_cut = "&& (Kp_Hlt1TrackMVADecision_TOS)" #"&& (pippim_Hlt1TwoTrackMVADecision_TOS)" # || pip_Hlt1TrackMVADecision_TOS || pim_Hlt1TrackMVADecision_TOS)"#"&& (pip_Hlt1TrackMVADecision_TOS) && (pim_Hlt1TrackMVADecision_TOS)"
ksk_cut = "&& (hp_Hlt1TrackMVADecision_TOS)"

file = ROOT.TFile.Open(input_file)
tree = file.Get('DecayTree')
data_frame = RDF(tree)
w=1
weights_file = f"./Draco_prod/{w}_iter/weights_ksk_{w}.root"
weight_file = ROOT.TFile.Open(weights_file)
bin_number = np.linspace(0,19,20) #[0,1,10,11,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9]
hist_2d_ksk = {}

for n, nbins in enumerate(bin_number):
    nbins = int(nbins)
    weights_tree = weight_file.Get(f"weights_ksk_tree_{nbins}")
    tree.AddFriend(weights_tree)
    hist_2d_ksk[n] = ROOT.TH2D(f"KSK_{nbins}",f"KSK_{nbins}", bins,  np.array(bin_edges_pt), bins,  np.array(bin_edges_eta))
    
tree.SetAlias("Dp_PT", "TMath::Sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY)")
tree.SetAlias("Dp_ETA", "0.5*TMath::Log((Dp_P + Dp_PZ)/(Dp_P - Dp_PZ))")
data_set = data_frame.Filter(mass_cut+ksk_cut).AsNumpy(columns=["Dp_M", "Dp_PT", "Dp_ETA"])

ROOT.gStyle.SetPalette(ROOT.kBird) #kViridis, kCool, kBird
canvas = ROOT.TCanvas()
canvas.Divide(2,2)
canvas.Print("pt_eta_weighted.pdf[")

for n, nbins in enumerate(bin_number):
    print(f"Processing bin {n}")
    nbins = int(nbins)
    data_set_w = data_frame.Filter(mass_cut+ksk_cut).AsNumpy(columns=[f"weights_ksk_{nbins}"])
    eta = data_set['Dp_ETA']
    pt = data_set['Dp_PT']
    weights_column = f"weights_ksk_{nbins}"
    weight = data_set_w[weights_column]
    hist_2d_ksk[n].FillN(len(pt), array('d', pt), array('d', eta), array('d', weight))
    
    canvas.SetRightMargin(0.15)
    hist_2d_ksk[n].Scale(1./hist_2d_ksk[n].Integral())
        
    file2 = ROOT.TFile.Open(f"FCUT_Dalitz_Bins/fcut_dalitz_bin_{n}_hlt1_kp.root")  # Adjust the file name as needed
    tree2 = file2.Get('DecayTree')
    tree2.SetAlias("Dp_PT", "TMath::Sqrt(Dp_PX*dp_PX + Dp_PY*Dp_PY)")
    tree2.SetAlias("Dp_ETA", "0.5*TMath::Log((Dp_P + Dp_PZ)/(Dp_P - Dp_PZ))")
    hist_2d_kpipi = ROOT.TH2D(f"KPIPI_{n}",f"KPIPI_{n}", bins,  np.array(bin_edges_pt), bins,  np.array(bin_edges_eta))
    tree2.Draw(f"Dp_ETA:Dp_PT>>KPIPI_{n}", mass_cut+kpipi_cut,"goff") #y:x>>histo
    hist_2d_kpipi.Scale(1./hist_2d_kpipi.Integral())

    canvas.cd(1)
    hist_2d_kpipi.GetXaxis().SetTitle("Ds transverse momentum")
    hist_2d_kpipi.GetYaxis().SetTitle("Ds pseudorapidity")
    hist_2d_kpipi.SetStats(0)
    hist_2d_kpipi.Draw("COLZ")  # or "SCAT" for scatter plot if you prefer
    text_kpp = ROOT.TLatex(20e3, 4.2, 'K#pi#pi')
    text_kpp.SetTextSize(0.1)
    text_kpp.Draw("same")

    canvas.cd(2)
    hist_2d_ksk[n].GetXaxis().SetTitle("Ds transverse momentum")
    hist_2d_ksk[n].GetYaxis().SetTitle("Ds pseudorapidity")
    hist_2d_ksk[n].SetStats(0)
    hist_2d_ksk[n].Draw("COLZ")  # or "SCAT" for scatter plot if you prefer
    text_ksk = ROOT.TLatex(20e3, 4.2, 'K_{S}K')
    text_ksk.SetTextSize(0.1)
    text_ksk.Draw("same")

    canvas.cd(3)
    diff = hist_2d_kpipi.Clone(f"diff_{n}")
    diff.Add(hist_2d_ksk[n],-1)
    diff.GetXaxis().SetTitle("Ds transverse momentum")
    diff.GetYaxis().SetTitle("Ds pseudorapidity")
    diff.Draw("COLZ TEXT")
    # print(np.sum(diff.GetEntries()))
    text_d = ROOT.TLatex(20e3, 4.2, 'diff')
    text_d.SetTextSize(0.1)
    text_d.Draw("same")
    
    canvas.cd(0)
    text_pad = ROOT.TLatex(0.7, 0.3, f"Bin {n}")
    text_pad.SetNDC()
    text_pad.SetTextSize(0.05)
    text_pad.Draw("same")
    
    values = []
    for ix in range(1, diff.GetNbinsX() + 1):
        for iy in range(1, diff.GetNbinsY() + 1):
            value = diff.GetBinContent(ix, iy)
            values.append(round(value,7))
    print(values)

    canvas.Update()
    canvas.Print("pt_eta_weighted.pdf")
    
canvas.Print("pt_eta_weighted.pdf]")