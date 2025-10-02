import ROOT
from ROOT import TMath
import sys
import numpy as np
ROOT.gROOT.ProcessLine(".L extras/lhcbStyle.C");
RDF = ROOT.ROOT.RDataFrame

def get_argparser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--plot_file', required=True)
    return parser

options = get_argparser().parse_args()
plot_file = options.plot_file

def makeHistos(plotFileName, variables, bins):
    
    tree = ROOT.TChain("DecayTree")
    files = [i for i in range(20)]
    for i in files:
        tree.Add(f"./Dalitz_Bins_Kp_TOS/dalitz_bin_{i}_hlt1_kp.root")

    mx = variables[0]#"Msq_pp_dtf"
    my = variables[1] #"Msq_kp_dtf"
    df = RDF(tree)
    ds = df.Filter("Kp_Hlt1TrackMVADecision_TOS").AsNumpy(columns=[mx, my, "sWeight"])

    x_range = [0, 3.7]
    y_range = [0, 2.5]

    kp = ROOT.TH2D("dalitz", "dalitz",bins,x_range[0],x_range[1], bins, y_range[0], y_range[1]) #, tree.GetMinimum(mx)) #,tree.GetMaximum(mx),bins,tree.GetMinimum(my),tree.GetMaximum(my))
    # tree.Draw("("+my+")/1e6:("+mx+")/1e6>>dalitz", "TMath::Abs(Dp_M-1969)<10","goff") #TMath::Abs(Dp_M-1969)<10
    for i in range(len(ds[mx])):
        kp.Fill(ds[mx][i]**2/1e6, ds[my][i]**2/1e6, ds["sWeight"][i])

    m1 = ROOT.TH1D("m1", "m1",bins,x_range[0],x_range[1])
    for i in range(len(ds[mx])):
        m1.Fill(ds[mx][i]**2/1e6, ds["sWeight"][i])
    
    m2 = ROOT.TH1D("m2", "m2",bins,y_range[0],y_range[1])
    for i in range(len(ds[my])):
        m2.Fill(ds[my][i]**2/1e6, ds["sWeight"][i])
     
    dalitz_bin_edges=[]   
    with open('./extras/dalitz_bin_edges.txt', 'r') as file:
        for line in file:
            numbers = line.split()
            dalitz_bin_edges.append([float(i) for i in numbers])
    print(dalitz_bin_edges)

    canvas = ROOT.TCanvas("canvas")
    
    canvas.cd(0)
    canvas.SetRightMargin(0.2) 
    canvas.Print(plotFileName+"[")
    kp.SetStats(0) 
    kp.SetLineColor(ROOT.kBlack)
    kp.SetLineWidth(2)  
    kp.GetYaxis().SetTitle("#it{m}^{2}(#it{#pi}^{\pm}#it{#pi}^{\mp}) [GeV^{2}/#it{c}^{4}]")
    kp.GetXaxis().SetTitle("#it{m}^{2}(#it{K}^{\pm}#it{#pi}^{\mp}) [GeV^{2}/#it{c}^{4}]")
    kp.GetZaxis().SetTitle("Candidates")
    kp.GetZaxis().SetTitleOffset(0.9) 
    ROOT.gStyle.SetPalette(ROOT.kViridis)  # Choose colormap 
    kp.Draw("colz")

    canvas.SetLogz(True)
    text_2 = ROOT.TLatex(2, 2.2, "LHCb 2024 Block 1")
    text_2.SetTextSize(0.045)
    text_2.Draw("same")
    text_ = ROOT.TLatex(2, 2.05, "HLT1TrackMVA_TOS(K^{+})")
    text_.SetTextSize(0.03)
    text_.Draw("same")
    text = ROOT.TLatex(2.0, 1.85, "#it{D}_{#it{s}}^{\pm} #rightarrow #it{K}^{\pm}#it{#pi}^{\pm}#it{#pi}^{\mp}")
    text.SetTextSize(0.06)
    text.Draw("same")
    
    boxes, texts = [], []
    
    """To draw Dalitz bins on the plot."""
    for idx, bin in enumerate(dalitz_bin_edges):
        x_min, x_max, y_min, y_max = bin
        print(bin)
    
        # Draw rectangle
        box = ROOT.TBox(x_min, y_min, x_max, y_max)
        box.SetLineColor(ROOT.kPink)
        box.SetLineWidth(2)
        box.SetFillStyle(0)
        box.SetLineStyle(1)
        boxes.append(box)
        box.Draw("same")
    
        # Add bin number text
        x_center = (x_min + x_max) / 2.0
        y_center = (y_min + y_max) / 2.0 
        text = ROOT.TLatex(x_center, y_center, str(idx))
        text.SetTextColor(ROOT.kPink)
        text.SetTextSize(0.035)
        text.SetTextAlign(22)
        texts.append(text)
        text.Draw("same")

        # Update canvas
        canvas.Update()
    
    canvas.cd(1)
    canvas.SetLeftMargin(0.15) 
    canvas.Print(plotFileName+"[")
    m1.SetStats(0)
    m1.SetLineColor(ROOT.kBlue)
    m1.SetLineWidth(2)
    m1.SetMarkerSize(0)
    m1.SetTitle("#it{m}^{2}(#it{K}^{+} #pi^{-})")
    m1.GetYaxis().SetTitle("Number of Events")
    m1.GetYaxis().SetTitleOffset(1.2)  
    m1.GetXaxis().SetTitle("#it{m}^{2}(#it{#pi}^{\pm} #it{#pi}^{\mp}) [GeV^{2}/#it{c}^{4}]")
    m1.Draw("h")
    
    canvas.cd(2)
    canvas.SetLeftMargin(0.15) 
    canvas.Print(plotFileName+"[")
    m2.SetStats(0) 
    m2.SetLineColor(ROOT.kBlue)
    m2.SetLineWidth(2)
    m2.SetMarkerSize(0)
    m2.SetTitle("#it{m}^{2}(#pi^{+} #pi^{-})")
    m2.GetYaxis().SetTitle("Number of Events")
    m2.GetYaxis().SetTitleOffset(1.2) 
    m2.GetXaxis().SetTitle("#it{m}^{2}(#it{K}^{\pm} #it{#pi}^{\mp}) [GeV^{2}/#it{c}^{4}]")
    m2.Draw("h")
    
    canvas.Print(plotFileName)  
    canvas.Print(plotFileName+"]")
    canvas.SaveAs(plotFileName.replace(".pdf",".png"))
    
    
bins = 100 # change number of bins

variables = ["KpPim_M_Constr", "PipPim_M_Constr"]   
makeHistos(plot_file,variables, bins)
    



