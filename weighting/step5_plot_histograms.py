import ROOT
import numpy as np
import math
from tqdm import tqdm
import time
from array import array
import json
import logging
ROOT.gROOT.ProcessLine(".L extras/lhcbStyle.C");

RDF = ROOT.ROOT.RDataFrame

class Plot_Histograms:
    
    """
    Class to plot histograms of the kinematic variables after each reweighting iteration.
    Plots : (P,ETA,PHI,PT) of the Ds and the Kaon, for original (sweighted) and kinematically weighted Kpipi and KSK data
    Plots the ratio or the pull of the weighted to the original histograms.
    """
    
    def __init__(self, ds_kpipi, ds_ksk, indices_kpipi, indices_ksk, weights_kpipi, weights_ksk, variable_names, bin_edges_plot, bins_plot, order, folder, bin_number, max_iter, first_coarse, iter):
        self.ds_kpipi = ds_kpipi
        self.ds_ksk = ds_ksk
        self.indices_df1 = indices_kpipi
        self.indices_df2 = indices_ksk
        self.weights_kpipi_final = weights_kpipi
        self.weights_ksk_final = weights_ksk
        self.var_names = variable_names
        self.bin_edges_plot = bin_edges_plot
        self.bins_plot = bins_plot
        self.order = order  # 0 for KSK, 1 for KPIPI
        self.folder = folder
        self.bin_number = bin_number
        self.max_iter = max_iter
        self.sum_ratios = None
        self.sum_of_sum_ratios = np.zeros(20)
        self.first_coarse = first_coarse
        self.iter = iter
        if max_iter == 0 : 
            self.first_coarse = False  # If no iterations, we don't have different orders
        
    def plot_variable(self) :
        title = f"{int(self.max_iter/3)}_it"
        canvas = ROOT.TCanvas() 
        canvas.Print(f"./{self.folder}/weighted_plots/Bin{self.bin_number}.pdf[")
        
        legends = []
        
        pages = math.ceil(len(self.var_names)/4)
        
        for ii in range(pages) :
            canvas.Clear()
            ii = ii*4   
            self.var_names_temp = self.var_names[ii:ii+4]
            self.bins_temp = self.bin_edges_plot[ii:ii+4]
            print(self.var_names_temp)

            pink, blue, ksk_w, kpipi_w = list(np.ones(4)), list(np.ones(4)), list(np.ones(4)), list(np.ones(4))

            for i,p in enumerate(pink) : 
                self.bin_edge_temp = array('d',self.bins_temp[i])
                pink[i] = ROOT.TH1F(f'kpp_{self.var_names_temp[i]}', title, len(self.bin_edge_temp)-1, np.array(self.bin_edge_temp))
                blue[i] = ROOT.TH1F(f'ksk_{self.var_names_temp[i]}', title, len(self.bin_edge_temp)-1, np.array(self.bin_edge_temp))
                ksk_w[i] = ROOT.TH1F(f'ksk_w_{self.var_names_temp[i]}', title, len(self.bin_edge_temp)-1, np.array(self.bin_edge_temp))
                kpipi_w[i] = ROOT.TH1F(f'kpp_w_{self.var_names_temp[i]}', title, len(self.bin_edge_temp)-1, np.array(self.bin_edge_temp))

            print("STARTING TO PLOT")
            
            print("sum of weights ksk : ", np.sum(self.weights_ksk_final))
            print(self.weights_ksk_final[:10])

            for i,p in enumerate(pink) :
                for x in tqdm(range(len(self.ds_kpipi[self.var_names_temp[i]]))) :
                    entry = self.ds_kpipi[self.var_names_temp[i]][x]
                    sweight = self.ds_kpipi["sWeight"][x]
                    p.Fill(entry, sweight)#, self.weights_kpipi[x])
                    if self.order == 1 or self.first_coarse :
                        # if self.iter == 0: kpipi_w[i].Fill(entry, sweight * self.weights_kpipi_final[self.indices_df1[x]])
                        kpipi_w[i].Fill(entry, self.weights_kpipi_final[self.indices_df1[x]]*sweight)
            
            for i,b in enumerate(blue) :
                for x in tqdm(range(len(self.ds_ksk[self.var_names_temp[i]]))) :
                # for x in tqdm((self.indices_df1)):
                    entry = self.ds_ksk[self.var_names_temp[i]][x]
                    sweight = self.ds_ksk["sWeight"][x]
                    b.Fill(entry, sweight)#, self.weights_kpipi[x])
                    if self.order == 0 or self.first_coarse:
                        # if self.iter == 0 : ksk_w[i].Fill(entry, sweight * self.weights_ksk_final[self.indices_df2[x]])
                        ksk_w[i].Fill(entry, self.weights_ksk_final[self.indices_df2[x]]*sweight)

            for p in pink : 
                p.Scale(1.0/p.Integral())
                p.SetLineColor(ROOT.kMagenta)
                p.SetFillColor(ROOT.kMagenta)
                p.SetFillStyle(3345)    
                p.SetLineWidth(1)
            for b in blue :
                b.Scale(1.0/b.Integral())
                b.SetLineColor(ROOT.kCyan+2)
                b.SetFillColor(ROOT.kCyan+2)
                b.SetFillStyle(3354)
                b.SetLineWidth(1)
            if self.order == 0 or self.first_coarse: 
                for k in ksk_w : 
                    k.Scale(1.0/k.Integral())
            if self.order == 1 or self.first_coarse: 
                for k in kpipi_w : 
                    k.Scale(1.0/k.Integral())
            
            leg1 = ROOT.TLegend(0.3, 0, 0.7, 1)
            leg1.SetHeader(f"Bin {self.bin_number} - 2D Kinematic Reweighting")# : ({self.weighting_var1}, {self.weighting_var2})")
            leg1.SetNColumns(3) 
            leg1.SetTextSize(0.3)
            leg1.AddEntry(pink[0], "K#pi#pi TOS","F")
            leg1.AddEntry(blue[0], "K_{S}K TOS","F")
            if self.first_coarse :
                leg1.AddEntry(ksk_w[0], "K_{S}K weighted", "PL")
                leg1.AddEntry(kpipi_w[0], "K#pi#pi weighted", "PL")
            else : 
                if self.order == 0 : leg1.AddEntry(ksk_w[0], "K_{S}K weighted", "PL")
                else : leg1.AddEntry(kpipi_w[0], "K#pi#pi weighted", "PL")
                     
            canvas.Divide(2,2)
            
            diff_ = []
            text_boxes = []   
            
            for i, p in enumerate(pink) :
                canvas.cd(i+1)   
                if self.first_coarse :
                    # plot the ratio
                    diff = ROOT.TRatioPlot(kpipi_w[i], ksk_w[i]) 
                    # plot the pull
                    # diff = ROOT.TRatioPlot(p, ksk_w[i],"diffsig") 
                else :
                    if self.order == 0 :
                        # plot the ratio
                        diff = ROOT.TRatioPlot(p, ksk_w[i]) 
                        # plot the pull
                        # diff = ROOT.TRatioPlot(p, ksk_w[i],"diffsig")
                    if self.order == 1 :
                        # plot the ratio
                        diff = ROOT.TRatioPlot(blue[i], kpipi_w[i]) 
                        # plot the ratio
                        # diff = ROOT.TRatioPlot(blue[i], kpipi_w[i],"diffsig")
            
                diff.Draw() 
                diff_.append(diff)
                diff.GetUpperRefYaxis().SetTitle('Normalized Entries')
                diff.GetUpperRefYaxis().SetLabelSize(0.05)
                diff.GetUpperRefYaxis().SetTitleOffset(1)
                diff.GetUpperRefYaxis().SetTitleSize(0.07)
                
                diff.GetXaxis().SetTitle(self.var_names_temp[i])
                diff.GetLowerRefXaxis().SetLabelSize(0.05)
                diff.GetLowerRefXaxis().SetTitleSize(10)
                diff.SetUpTopMargin(0.1)
                diff.SetLowBottomMargin(0.2)
                if i == 0 or i == 1 : diff.SetLowBottomMargin(0.25)
                if i == 2 or i == 3 : diff.SetUpTopMargin(0.25)
                
                diff.GetLowYaxis().SetNdivisions(505)
                diff.GetLowerRefYaxis().SetTitle("Ratio")
                diff.GetLowerRefYaxis().SetTitleOffset(1)
                diff.GetLowerRefYaxis().SetLabelSize(0.05)
                diff.GetLowerRefYaxis().SetTitleSize(0.06)
                
                diff.GetLowerRefGraph().SetMarkerSize(0.3)
                
                ratio_hist = diff.GetLowerRefGraph()
                
                if self.order == 0 : max_val = max(h.GetMaximum() for h in [p,blue[i],ksk_w[i]])
                else : max_val = max(h.GetMaximum() for h in [p, blue[i], kpipi_w[i]])
                
                diff.GetUpperPad().cd()
                diff.GetUpperRefYaxis().SetRangeUser(0, max_val*1.1)
                if self.first_coarse :
                    blue[i].GetYaxis().SetTitle("")
                    blue[i].GetYaxis().SetLabelSize(0)
                    blue[i].GetXaxis().SetTitle("")
                    blue[i].GetXaxis().SetLabelSize(0)
                    blue[i].Draw("hist same")
                    p.GetYaxis().SetTitle("")
                    p.GetYaxis().SetLabelSize(0)
                    p.GetXaxis().SetTitle("")
                    p.GetXaxis().SetLabelSize(0)
                    p.Draw("hist same")
                    kpipi_w[i].SetMarkerStyle(5)
                    kpipi_w[i].SetMarkerSize(0.7)
                    kpipi_w[i].SetMarkerColor(ROOT.kRed)
                    kpipi_w[i].Draw("e1 same")
                    ksk_w[i].SetMarkerStyle(5) 
                    ksk_w[i].SetMarkerSize(0.7)
                    ksk_w[i].SetMarkerColor(ROOT.kBlue)
                    ksk_w[i].Draw("e1 same")
                else :
                    if self.order == 0 :
                        blue[i].GetYaxis().SetTitle("")
                        blue[i].GetYaxis().SetLabelSize(0)
                        blue[i].GetXaxis().SetTitle("")
                        blue[i].GetXaxis().SetLabelSize(0)
                        blue[i].Draw("hist same")
                        ksk_w[i].SetMarkerStyle(5) 
                        ksk_w[i].SetMarkerSize(0.7)
                        ksk_w[i].Draw("e1 same")
                    if self.order == 1 : 
                        p.GetYaxis().SetTitle("")
                        p.GetYaxis().SetLabelSize(0)
                        p.GetXaxis().SetTitle("")
                        p.GetXaxis().SetLabelSize(0)
                        p.Draw("hist same")
                        kpipi_w[i].SetMarkerStyle(5)
                        kpipi_w[i].SetMarkerSize(0.7)
                        kpipi_w[i].Draw("e1 same")
                    
                text = ROOT.TPaveText(0.7, 0.6, 0.9, 0.7, "NDC")
                text.SetTextSize(0.08)
                text.SetLineColor(ROOT.kBlack)
                text.SetFillColor(ROOT.kWhite)
                text.SetBorderSize(0)
                text.AddText(self.var_names_temp[i])   
                text.Draw()
                text_boxes.append(text)
                # canvas.Update()
            legends.append(leg1)
                
            canvas.cd()
            legend_pad = ROOT.TPad("legend_pad", "", 0, 0.42, 1, 0.52)  # bottom strip
            legend_pad.SetFillStyle(4000)
            legend_pad.SetBorderSize(0)
            legend_pad.Draw()
            legend_pad.cd()
            if leg1 : leg1.Draw()
                
            canvas.Print(f"./{self.folder}/weighted_plots/Bin{self.bin_number}.pdf")
        
        canvas.Print(f"./{self.folder}/weighted_plots/Bin{self.bin_number}.pdf]")
