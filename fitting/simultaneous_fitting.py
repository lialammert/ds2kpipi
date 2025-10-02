import ROOT
from ROOT import RooFit as RF
from ROOT import TMath
import numpy as np
from fitting import Exponential, ExecuteFit, CrystalBall_Gaussian
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
Simultaneous Fitting for both Ds_plus and Ds_minus invariant mass distributions of the KsK sample before and after weighting.
Fit functions are defined in the `fitting` module and initial parameters are taken from "fitting_parameters.json".
sWeighted distributions are plotted along with the kinematically weighted data.
Yields, asymmetry and chi2 are saved in "a_raw.json".
Fit results, fit parameters are saved in "weighted_ksk_fit_results.txt".
To change : root_files_directory, kinematic variables, binning, and strategy.
"""


class SimultaneousFitter:
    def __init__(self,decay, 
                 n, bins, folder, 
                 root_files_directory, canvas, output_files, plot_file, results_dict,
                 w, mass_range,
                 kpipi_input_files, ksk_input_file, 
                 ksk_weights_files, kpipi_weights_files):
        self.decay = decay
        self.n = n
        self.bins = bins
        self.folder = folder
        self.root_files_directory = root_files_directory
        self.canvas = canvas
        self.output_files = output_files
        self.plot_file = plot_file
        self.results_dict = results_dict
        self.w = w
        self.mass_range = mass_range
        self.kpipi_input_files = kpipi_input_files
        self.ksk_input_file = ksk_input_file
        self.ksk_weights_files = ksk_weights_files
        self.kpipi_weights_files = kpipi_weights_files
        
        
        self.bin_edges_dp = np.linspace(mass_range[0], mass_range[1], bins + 1)
        self.var = "Dp_DTF_PV_MKS_M"
        self.mass = ROOT.RooRealVar("mass", 
                            "Mass of D_{s}^{+}(MeV)", 
                            mass_range[0], mass_range[1])
        self.kinematic_vars = ["Dp_DTF_PV_MKS_M", "Dp_M"] 
        self.formula_vars_tree_kpipi = ["Dp_DTF_PV_M", "Dp_M"] 
        self.formula_vars_tree_ksk = ["Dp_DTF_PV_MKS_M", "Dp_M"] 
        self.particle_id = ["(hp_PARTICLE_ID == 321)", "(hp_PARTICLE_ID == -321)"]
        self.particle_id_kpipi = ["(Kp_PARTICLE_ID == 321)", "(Kp_PARTICLE_ID == -321)"]
        self.kpipi_cut = "(Kp_Hlt1TrackMVADecision_TOS)"
        self.ksk_cut = "(hp_Hlt1TrackMVADecision_TOS)"
        self.id_labels = ["dp","dm"] #D^+, D^-
        self.total_yields = []
        self.indiv_yields = []
        self.asymmetry = []  
        self.chi = {} 
        
        
    def fit(self):   
        self.create_output()
        # for n in tqdm(range(self.nbins)):
        df_kpipi, df_ksk = self.get_weights(self.kpipi_input_files, self.ksk_input_file, self.kpipi_weights_files, self.ksk_weights_files)
        
        for id, id_cut in enumerate(self.particle_id):
            self.datasets(id, id_cut, df_kpipi, df_ksk)
            self.fit_shapes(id) #, sigma_narrow)
            self.histogram(id)
            
        print("---------------------------------------------")
        print("OUTPUT LOG OF THE FIT : ")
        sys.stdout.flush()

        self.simult_constraints(self.ksk_histograms) #or self.kpipi_histograms
        if self.decay == "ksk" :
            self.execute_fit(self.ksk_data, self.ksk_histograms)
        if self.decay == "kpipi" :
            self.execute_fit(self.kpipi_data, self.kpipi_histograms)
        self.plot_results()
        # self.save_asym()
            
        self.ksk_file.Close()
        self.kpipi_file.Close()
        self.json_results()

    def create_output(self):
        # self.asym_file = self.output_files[0]
        self.fit_results_file = self.output_files[0]
        # self.yields_file = self.output_files[2]
        self.asym = [[0, 0]] # asymmetry, uncertainty
        self.sig_yields = [[0, 0]] #nsig_pos, nsig_neg
        self.chi = [0,0]
        self.constr = []
        self.frame = {}

    def get_weights(self, kpipi_input_file, ksk_input_file, kpipi_weights_file, ksk_weights_file):
        # get trees from files and add weights as friend tree
        self.kpipi_file = ROOT.TFile(f"./{self.root_files_directory}/"+kpipi_input_file)
        self.ksk_file = ROOT.TFile(f"./{self.root_files_directory}/"+ksk_input_file)
        kpipi_tree = self.kpipi_file.Get("DecayTree")
        ksk_tree = self.ksk_file.Get("DecayTree")
        kpipi_weights_file = ROOT.TFile(kpipi_weights_file)
        kpipi_weights_tree = kpipi_weights_file.Get(f"weights_kpipi_tree_{self.n}")
        kpipi_tree.AddFriend(kpipi_weights_tree)
        ksk_weights_file = ROOT.TFile(ksk_weights_file)
        ksk_weights_tree = ksk_weights_file.Get(f"weights_ksk_tree_{self.n}")
        ksk_tree.AddFriend(ksk_weights_tree)
        df_kpipi = RDF(kpipi_tree)
        df_ksk = RDF(ksk_tree)
        for k,kin in enumerate(self.kinematic_vars):
            if kin in df_kpipi.GetColumnNames() : df_kpipi = df_kpipi.Redefine(kin, self.formula_vars_tree_kpipi[k])
            else : df_kpipi = df_kpipi.Define(kin, self.formula_vars_tree_kpipi[k])
            if kin in df_ksk.GetColumnNames() : df_ksk = df_ksk.Redefine(kin, self.formula_vars_tree_ksk[k])
            else : df_ksk = df_ksk.Define(kin, self.formula_vars_tree_ksk[k])

        # intialize dicts to save parameters (avoid memory issues)
        self.signal_shapes, self.background_shapes = {}, {} # function to fit to
        self.total_models, self.signal_models, self.background_models = {}, {}, {} #RooAddPdf of the shapes
        self.signal_pdfs, self.background_pdfs = {}, {}
        self.extra_parameters = []
        self.alphas, self.sigmas = {}, {}
        self.sig_p, self.bkg_p, self.sig_p_str, self.bkg_p_str = [],[],[],[]
        self.kpipi_histograms, self.ksk_histograms = {}, {}
        self.kpipi_data, self.ksk_data = {}, {}
        self.total_signal_yield, self.signal_asymmetry = {}, {}
        
        return df_kpipi, df_ksk
        ## make sure all the variables are defined in the dataframes

    def datasets(self, id, id_cut, df_kpipi, df_ksk):
        self.ds_kpipi = df_kpipi.Filter(self.particle_id_kpipi[id]).AsNumpy(columns=self.kinematic_vars)
        self.ds_ksk = df_ksk.Filter(id_cut).AsNumpy(columns=self.kinematic_vars)
        self.ds_kpipi_weights = df_kpipi.Filter(self.particle_id_kpipi[id]).AsNumpy(columns=[f"weights_kpipi_{self.n}"])
        self.ds_ksk_weights = df_ksk.Filter(id_cut).AsNumpy(columns=[f"weights_ksk_{self.n}"])
        
    def fit_shapes(self, id):
        # print("Fitting bin ", self.n, " for decay ", self.decay, " with id ", id)
        # Use JSON parameters for the specific bin
        if self.w != 0 : parameters_file = f"bin{self.n}"
        else : parameters_file = "bin-1"
        self.signal_shapes[id] = CrystalBall_Gaussian(self.id_labels[id], self.decay, parameters_file) 
        self.background_shapes[id] = Exponential(self.id_labels[id], self.decay, parameters_file)
        
        sig_model, sig1, sig2, frac, alphaL, alphaR, sigmaL, sigmaR = self.signal_shapes[id].init_shape(
            self.mass, bin_name=parameters_file
        )

        self.signal_models[id] = sig_model
        self.extra_parameters += [sig1, sig2, frac]
        self.n_sig, self.n_bkg = {},{}
        self.bkg_model = self.background_shapes[id].init_shape(self.mass, bin_name=f"bin{self.n}")

        self.background_models[id] = self.bkg_model
        self.alphas[id] = [alphaL, alphaR]
        self.sigmas[id] = [sigmaL, sigmaR]
        self.fit = {}
        
    def histogram(self, id):
        self.kpipi_histograms[id] = ROOT.TH1F("kpipi_Dp_M_"+str(id), self.var, self.bins,  np.array(self.bin_edges_dp))
        self.ksk_histograms[id] = ROOT.TH1F("ksk_Dp_M_"+str(id), self.var, self.bins,  np.array(self.bin_edges_dp))

        for i in tqdm(range(len(self.ds_kpipi[self.var]))):
            entry = self.ds_kpipi[self.var][i]
            weight = self.ds_kpipi_weights[f'weights_kpipi_{self.n}'][i]
            self.kpipi_histograms[id].Fill(entry, weight)
        for i in tqdm(range(len(self.ds_ksk[self.var]))):
            entry = self.ds_ksk[self.var][i]
            weight = self.ds_ksk_weights[f'weights_ksk_{self.n}'][i]
            self.ksk_histograms[id].Fill(entry, weight)
            
        self.kpipi_histograms[id].SetDirectory(0)
        self.ksk_histograms[id].SetDirectory(0)
        
        self.kpipi_data[id] = ROOT.RooDataHist(f"kpipi_data_{id}", "", ROOT.RooArgList(self.mass),self.kpipi_histograms[id])
        self.ksk_data[id] = ROOT.RooDataHist(f"ksk_data_{id}", "", ROOT.RooArgList(self.mass),self.ksk_histograms[id])
        
    def simult_constraints(self, histograms):
        projected_histogram = histograms[0] #dm
        projected_histogram.Add(histograms[1]) #dp
        projected_histogram.SetDirectory(0)
        self.total_signal_yield = ROOT.RooRealVar("n_signal_total", "", 0.4*projected_histogram.GetEntries(), 0.0, 1.5*projected_histogram.GetEntries())
        self.signal_asymmetry = ROOT.RooRealVar("signal_asymmetry", "", 0.0, -0.3, 0.3);

        self.n_signal_pos = ROOT.RooFormulaVar ( "signal_yield_positive", "signal_yield_positive", 
                                                            "(1+@0)*@1*0.5", 
                                                            ROOT.RooArgList ( self.signal_asymmetry, 
                                                                        self.total_signal_yield) );

        self.n_signal_minus = ROOT.RooFormulaVar ( "signal_yield_negative", "signal_yield_negative",  
                                                            "(1-@0)*@1*0.5", 
                                                            ROOT.RooArgList ( self.signal_asymmetry, 
                                                                        self.total_signal_yield ) );
        
        self.diff_alphaL = ROOT.RooFormulaVar ( "sq_alphaL", "sq_alphaL", 
                                                "(@0-@1)*(@0-@1)",
                                                ROOT.RooArgList ( self.alphas[0][0], self.alphas[1][0] ) );
        self.diff_alphaR = ROOT.RooFormulaVar ( "sq_alphaR", "sq_alphaR",
                                                "(@0-@1)*(@0-@1)",
                                                ROOT.RooArgList ( self.alphas[0][1], self.alphas[1][1] ) );
        self.diff_aL = ROOT.RooGaussian ( "diff_aL", "diff_aL",
                                        self.diff_alphaL, 
                                        ROOT.RooFit.RooConst ( 0.0 ), 
                                        ROOT.RooFit.RooConst ( 0.01 ) );
        self.diff_aR = ROOT.RooGaussian ( "diff_aR", "diff_aR",
                                        self.diff_alphaR, 
                                        ROOT.RooFit.RooConst ( 0.0 ), 
                                        ROOT.RooFit.RooConst ( 0.01 ) );
        self.diff_sigmaL = ROOT.RooFormulaVar ( "sq_sigmaL", "sq_sigmaL",
                                                "(@0-@1)*(@0-@1)",
                                                ROOT.RooArgList ( self.sigmas[0][0], self.sigmas[1][0] ) );
        self.diff_sigmaR = ROOT.RooFormulaVar ( "sq_sigmaR", "sq_sigmaR",
                                                "(@0-@1)*(@0-@1)",
                                                ROOT.RooArgList ( self.sigmas[0][1], self.sigmas[1][1] ) );
        self.diff_sL = ROOT.RooGaussian ( "diff_sL", "diff_sL",
                                        self.diff_sigmaL, 
                                        ROOT.RooFit.RooConst ( 0.0 ), 
                                        ROOT.RooFit.RooConst ( 0.01 ) );
        self.diff_sR = ROOT.RooGaussian ( "diff_sR", "diff_sR",
                                        self.diff_sigmaR, 
                                        ROOT.RooFit.RooConst ( 0.0 ), 
                                        ROOT.RooFit.RooConst ( 0.01 ) );

        self.constr += [self.diff_aL]
        self.constr += [self.diff_aR]
        # self.constr += [self.diff_sL]
        # self.constr += [self.diff_sR]
        
    def execute_fit(self, datasets, histogram):
        ns = [self.n_signal_pos, self.n_signal_minus]
        for id, id_cut in enumerate(self.particle_id):
            self.fit[id] = ExecuteFit(self.var, self.bins, self.mass_range, self.mass, id_cut,self.id_labels[id])
            self.total_models[id], self.signal_pdfs[id], self.background_pdfs[id] = self.fit[id].init_model(self.signal_models[id],self.background_models[id],datasets[id],histogram[id], self.constr, ns[id])
            self.n_sig[id], self.n_bkg[id] = self.fit[id].get_params()
               
        self.category = ROOT.RooCategory("tag", "tag")
        self.category.defineType("positive")
        self.category.defineType("negative")

        self.sim_pdf = ROOT.RooSimultaneous("SimPdf", "simultaneous pdf", self.category)#ROOT.RooArgList(total_models[0], total_models[1]),category )
        self.sim_pdf.addPdf(self.total_models[0], "positive")
        self.sim_pdf.addPdf(self.total_models[1], "negative")

        self.sim_dataset = ROOT.RooDataHist(f"CategorisedDataset_{self.n}",f"CategorisedDataset_{self.n}",
                                       ROOT.RooArgList( self.mass ),
                                       RF.Index(self.category),
                                       RF.Import("positive", datasets[0]),
                                       RF.Import("negative", datasets[1]))
        
        attempt = 0
            
        while attempt < 10:
            self.fit_result = self.sim_pdf.fitTo(self.sim_dataset, RF.Save(), 
                                    RF.Extended(True), 
                                    RF.Minimizer("Minuit"),
                                    RF.Strategy(2),
                                    RF.PrintEvalErrors(-1),
                                    RF.ExternalConstraints(self.constr),
                                    RF.SumW2Error(True)) #, ROOT.RooFit.ExternalConstraints(constraints))

            if self.fit_result.status() == 0 and self.fit_result.covQual() >= 3:
                print(f"Fit succeeded on attempt {attempt + 1} (Status: {self.fit_result.status()}, covQual: {self.fit_result.covQual()})")
                attempt = 10
            
            else : 
                print(f"Fit failed on attempt {attempt + 1} (Status: {self.fit_result.status()}, covQual: {self.fit_result.covQual()})")
                for id in range(len(self.particle_id)):
                    for i in self.signal_shapes[id].get_param():
                        print(f"{i.GetName()} : {i.getVal()}")
                    for i in self.background_shapes[id].get_param():
                        print(f"{i.GetName()} : {i.getVal()}")
                attempt += 1

        print("ASYMMETRY")
        print(self.signal_asymmetry.getVal())
        self.asym = [self.signal_asymmetry.getVal(),self.signal_asymmetry.getError()]
        self.sig_yields = [round(self.n_sig[0].getVal()),round(self.n_sig[1].getVal())]

    def plot_results(self):
        print(f"Plotting results for bin {self.n}")
        for id, id_cut in enumerate(self.particle_id):
            self.cat = self.category.lookupType(id).GetName()
            print("CATEGORY + NAME : ")
            print(self.cat)
            print(self.signal_pdfs[id].GetName())
            self.sig_p.append(self.signal_shapes[id].get_param())
            self.bkg_p.append(self.background_shapes[id].get_param())
            self.sig_p_str.append(self.signal_shapes[id].get_param_str())
            self.bkg_p_str.append(self.background_shapes[id].get_param_str())
            
            self.frame[id] = self.mass.frame()

            self.sim_dataset.plotOn(self.frame[id], 
                            RF.Name(f"data_{self.n}bin_{self.id_labels[id]}"), 
                            RF.Cut("tag==tag::"+self.cat), 
                            RF.MarkerStyle(ROOT.kPlus))
            
            self.sim_pdf.plotOn(self.frame[id], 
                                    RF.LineColor(ROOT.kBlue-2),
                                    RF.Name(f"model_{self.n}bin_{self.id_labels[id]}"),
                                    # RF.Cut("tag==tag::"+self.cat),
                                    RF.Slice(self.category,self.cat), 
                                    # RF.MoveToBack(),
                                    RF.ProjWData (ROOT.RooArgSet(self.category), self.sim_dataset))
            self.sim_pdf.plotOn(self.frame[id],
                                    RF.Components(self.signal_pdfs[id].GetName()),
                                    RF.Name(f"signal_{self.n}bin_{self.id_labels[id]}"),
                                    RF.LineStyle(ROOT.kDashed), 
                                    RF.LineColor(ROOT.kCyan+2),
                                    RF.Slice(self.category,self.cat),
                                    RF.MoveToBack(),
                                    RF.ProjWData(ROOT.RooArgSet(self.category), self.sim_dataset)) #, ROOT.RooFit.Components(signal_pdfs[id]), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kCyan+2),ROOT.RooFit.Name(f"signal_{k}bin_{id_labels[id]}"))
            self.sim_pdf.plotOn(self.frame[id],
                                    RF.Components(self.background_pdfs[id].GetName()),
                                    RF.Name(f"background_{self.n}bin_{self.id_labels[id]}"),
                                    RF.LineStyle(ROOT.kDashed), 
                                    RF.LineColor(ROOT.kRed),
                                    RF.Slice(self.category,self.cat),
                                    RF.MoveToBack(),
                                    RF.ProjWData(ROOT.RooArgSet(self.category), self.sim_dataset))
            self.chi[id] = (self.frame[id].chiSquare(f"model_{self.n}bin_{self.id_labels[id]}", f"data_{self.n}bin_{self.id_labels[id]}", nFitParam = len(self.sig_p[id])+len(self.bkg_p[id])))

            pad1 = ROOT.TPad ("data_pad", "data_pad", 0.01, 0.35, 0.98, 0.99);
            ROOT.SetOwnership(pad1, False)
            pad1.SetBottomMargin(0);
            pad1.SetTopMargin(0.05);
            pad2 = ROOT.TPad("pull_pad", "pull_pad", 0.01, 0.01, 0.98, 0.33);
            ROOT.SetOwnership(pad2, False)
            pad2.SetTopMargin(0.0);
            pad2.SetBottomMargin(0.35);
            self.frame[id].GetXaxis().SetTitle("");
            self.frame[id].GetYaxis().SetTitle("Events");
            self.frame[id].GetYaxis().SetTitleSize(0.08);
            self.frame[id].GetYaxis().SetLabelSize(0.06);
            # pad1.SetLogy()

            self.h_resid_errors = self.frame[id].pullHist (f"data_{self.n}bin_{self.id_labels[id]}", f"model_{self.n}bin_{self.id_labels[id]}", True);
            pull_graph = ROOT.TGraph ( self.h_resid_errors.GetN (), self.h_resid_errors.GetX (), self.h_resid_errors.GetY () );
            pull_graph.GetXaxis ().SetLimits ( pull_graph.GetX ()[0], pull_graph.GetX ()[pull_graph.GetN () - 1] );
            pull_graph.GetYaxis().SetRangeUser(-5,5);
            pull_graph.GetXaxis().SetNdivisions(505);
            pull_graph.GetYaxis().SetNdivisions(505);

            pull_graph.GetYaxis().SetTitle("Pull");
            if id == 0 : pull_graph.GetXaxis().SetTitle( "#it{m} (#it{D}_{#it{s}}^{+}) [MeV/#it{c}]");
            if id == 1 : pull_graph.GetXaxis().SetTitle( "#it{m} (#it{D}_{#it{s}}^{-}) [MeV/#it{c}]");
            # pull_graph.GetXaxis().SetTitle( x_title );
            pull_graph.GetXaxis().SetTitleSize(0.16);
            pull_graph.GetXaxis().SetLabelSize(0.13);
            pull_graph.GetYaxis().SetTitleSize(0.16);
            pull_graph.GetYaxis().SetLabelSize(0.13);
            pull_graph.GetYaxis().SetTitleOffset(0.4);
            pull_graph.GetXaxis().SetTitleOffset(0.9);
            pull_graph.SetFillStyle(1001)
            pull_graph.SetFillColor(ROOT.kBlack)

            dot_functie = ROOT.TF1("plus_twee_sigma", "2", pull_graph.GetX ()[0], pull_graph.GetX ()[pull_graph.GetN () - 1] )
            dot_functie.SetLineColor(ROOT.kRed)
            dot_functie.SetLineStyle( ROOT.kDotted )

            min_twee_sigma = ROOT.TF1("min_twee_sigma", "-2", pull_graph.GetX ()[0], pull_graph.GetX ()[pull_graph.GetN () - 1] )
            min_twee_sigma.SetLineColor(ROOT.kRed)
            min_twee_sigma.SetLineStyle( ROOT.kDotted )
            
            self.canvas.cd(1);
            # canvas.Clear()
            pad1.Draw();
            pad2.Draw()

            pad1.cd ();
            self.frame[id].SetMinimum(0.5);
            self.frame[id].Draw();
            
            pad2.cd ();
            pull_graph.Draw("AB");
            dot_functie.Draw("SAME")
            min_twee_sigma.Draw("SAME")
            pull_graph.Draw("B SAME");

            self.canvas.cd(); 
            if self.decay == "ksk" : 
                text2 = ROOT.TLatex(0.2, 0.85, "LHCb 2024 Block 1")
                text2.SetTextSize(0.05)
                text2.Draw("same")
                text3 = ROOT.TLatex(0.2, 0.78, "#it{D}_{#it{s}}^{\pm} #rightarrow ##it{K}_{#it{S}} #it{K}")
                text3.SetTextSize(0.055)
                text3.Draw("same")
                if self.w == 0 :
                    text = ROOT.TLatex(0.65, 0.83, "#it{K}_{#it{S}} #it{K} unweighted")
                else : 
                    text_1 = ROOT.TLatex(0.65, 0.88, f"Dalitz Bin {self.n}")
                    text_1.SetTextSize(0.05)
                    text_1.Draw("same")
                    text = ROOT.TLatex(0.65, 0.83, f"{self.w} re-weightings")
                text.SetTextSize(0.05)
                text.Draw("same")
                    

            legend = ROOT.TLegend(0.65, 0.6, 0.8, 0.8)
            legend.SetTextSize(0.045)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.AddEntry(self.frame[id].findObject(f"model_{self.n}bin_{self.id_labels[id]}"), "Fit","l")
            legend.AddEntry(self.frame[id].findObject(f"signal_{self.n}bin_{self.id_labels[id]}"), "Signal", "l")#signal_models[cc].GetName(),"l")
            # for i,p in enumerate(self.sig_p[id]):
            #     legend.AddEntry('None', self.sig_p_str[id][i]+": "+str(round(p.getVal(),1)), "")
            legend.AddEntry(self.frame[id].findObject(f"background_{self.n}bin_{self.id_labels[id]}"), "Background", "l")#background_models[cc].GetName(),"l")
            # for i,p in enumerate(self.bkg_p[id]):
            #     legend.AddEntry('None',self.bkg_p_str[id][i]+": "+str(round(p.getVal(),4)),"")
            legend.AddEntry('None', "#chi^{2}: "+str(round(self.chi[id],2)), "p")
            legend.Draw("same")
            self.canvas.Update()
            
            with open(self.fit_results_file, "a") as f:
                f.write(f"Weighting {self.w} for {self.id_labels[id]} : "+"\n")
                f.write(f"Signal yield for {self.id_labels[id]} : " + str(round(self.n_sig[id].getVal()))+" +/- "+str(round(self.n_sig[id].getPropagatedError(self.fit_result)))+"\n")
                f.write(f"Background yield for {self.id_labels[id]} : " + str(round(self.n_bkg[id].getVal()))+" +/- "+str(round(self.n_bkg[id].getPropagatedError(self.fit_result)))+"\n")
                f.write("Asymmetry = "+str(round(self.signal_asymmetry.getVal(),6))+" +/- "+str(round(self.signal_asymmetry.getError(),6))+"\n")
                f.write("Chi2 : "+str(round(self.chi[id],2))+"\n")
                f.write("Fit Status : "+str(self.fit_result.status())+"\n")
                f.write("Cov Qual : "+str(self.fit_result.covQual())+"\n")
                f.write("Parameters : "+"\n")
                for i in self.signal_shapes[id].get_param():
                    f.write(f"{i.GetName()} : {round(i.getVal(),4)} +/- {round(i.getError(),4)}\n")
                for i in self.bkg_p[id]:
                    f.write(f"{i.GetName()} : {round(i.getVal(),4)} +/- {round(i.getError(),4)}\n")
                f.write("--------------------"+"\n")
                
            self.canvas.Print(f'./{self.folder}/{self.plot_file}.pdf')
            self.frame[id].remove(f"data_{self.n}bin_{self.id_labels[id]}")

            
    def json_results(self):
        print("Saving final results ")
        
        indiv_yields = [self.sig_yields[0], self.sig_yields[1]]
        total_sig = np.sum(self.sig_yields)
        self.total_yields.append(total_sig)
        self.asymmetry.append([self.asym[0], self.asym[1]])
        n_err = [self.n_sig[0].getPropagatedError(self.fit_result),self.n_sig[1].getPropagatedError(self.fit_result)]
        asym_err = 2/(self.n_sig[0].getVal()+ self.n_sig[1].getVal())**2 * np.sqrt(self.n_sig[1].getVal()**2 * n_err[0]**2 + self.n_sig[0].getVal()**2 * n_err[1]**2)
        # asymmetry_2.append([round(self.asym[0],5), round(asym_err,5)])
        if self.w == 0 :
            self.results_dict[f"ksk"] = {
                "yields": indiv_yields,
                "Araw": [self.asymmetry[0][0],self.asymmetry[0][1]],
                "chi2": [round(self.chi[0],2), round(self.chi[1],2)],
            }
        else :     
            self.results_dict[f"{self.w}w_bin{self.n}"] = {
                "yields": indiv_yields,
                "Araw": [self.asymmetry[0][0],self.asymmetry[0][1]],
                "chi2": [round(self.chi[0],2), round(self.chi[1],2)],
            }


