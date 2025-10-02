import ROOT
from ROOT import TMath
import numpy as np
import logging
import json
ROOT.gROOT.ProcessLine(".L extras/RooCruijff.C");

"""
functions necessary for fitting
"""

def format_scientific(value):
    """Format a number in scientific notation."""
    return f"{value:.2e}"

def load_fitting_parameters(json_file="fitting_parameters.json", bin_name=None):
    """
    Load fitting parameters from JSON file.
    
    Args:
        json_file (str): Path to the JSON file containing fitting parameters
        bin_name (str): Specific bin to load (e.g., "bin0"). If None, returns all bins.
    
    Returns:
        dict: Dictionary containing the parameters for the specified bin or all bins
    """
    try:
        with open(json_file, 'r') as f:
            params = json.load(f)
        
        if bin_name:
            if bin_name in params:
                return params[bin_name]
            else:
                raise KeyError(f"Bin {bin_name} not found in {json_file}")
        else:
            return params
            
    except FileNotFoundError:
        print(f"Error: {json_file} not found")
        return {}
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in {json_file}")
        return {}

def get_param_value(params_dict, param_name, default_value=None):
    """
    Get parameter value from the parameters dictionary.
    
    Args:
        params_dict (dict): Dictionary containing parameters
        param_name (str): Name of the parameter to retrieve
        default_value: Default value if parameter not found
        
    Returns:
        float: Parameter value
    """
    if param_name in params_dict:
        return float(params_dict[param_name])
    else:
        print(f"Warning: Parameter {param_name} not found, using default value {default_value}")
        return default_value if default_value is not None else 0.0      

class CrystalBall_Gaussian():
    def __init__(self, name, decay, bin_name=None, json_file="fitting_parameters.json"):
        self.name = name
        self.decay = decay
        
        # Load parameters from JSON if bin_name is provided
        if bin_name:
            params = load_fitting_parameters(json_file, bin_name)
            # Get central values from JSON
            mean_val = get_param_value(params, "mean", 1968.5)
            sigmaL_val = get_param_value(params, "sigmaL", 10.5)
            sigmaR_val = get_param_value(params, "sigmaR", 9.8)
            mean_g_val = get_param_value(params, "mean_g", 1968.5)
            sigma_g_val = get_param_value(params, "sigma_g", 6.4)
            self.alphaL_val = get_param_value(params, "alphaL", 1.4)
            self.alphaR_val = get_param_value(params, "alphaR", 2.1)
            self.nL_val = get_param_value(params, "nL", 3.0)
            self.nR_val = get_param_value(params, "nR", 3.0)
            # Store constants info for later use
            self.constants_str = params.get("constants", "")
        else:
            # Default values if no JSON parameters provided
            mean_val = 1968.5
            sigmaL_val = 10.5
            sigmaR_val = 9.8
            mean_g_val = 1968.5
            sigma_g_val = 6.4
            self.alphaL_val = 1.4
            self.alphaR_val = 2.1
            self.nL_val = 3.0
            self.nR_val = 3.0
            self.constants_str = ""
        
        self.mean = ROOT.RooRealVar("mean_"+str(self.name), "Mean Mass", mean_val, 1920, 2025)
        self.sigmaL = ROOT.RooRealVar("sigmaL_"+str(self.name), "Left Tail Width", sigmaL_val, 1, 20)
        self.sigmaR = ROOT.RooRealVar("sigmaR_"+str(self.name), "Right Tail Width", sigmaR_val, 1, 20)
        # self.alphaR = ROOT.RooRealVar("alphaR_"+str(self.name), "Right Tail Parameter", 1.5, 0.5, 5)
        # self.alphaL = ROOT.RooRealVar("alphaL_"+str(self.name), "Left Tail Parameter", 1.5, 0.5, 5)
        # self.nL = ROOT.RooRealVar("nL_"+str(self.name), "Left Tail Exponent", 3.5, 0.5, 10)
        # self.nR = ROOT.RooRealVar("nR_"+str(self.name), "Right Tail Exponent", 5.0, 0.5, 10)
        self.mean_g = ROOT.RooRealVar("mean_g_"+str(self.name), "Mean of Gaussian", mean_g_val, 1920, 2025)
        self.sigma_g = ROOT.RooRealVar("sigma_g_"+str(self.name), "Width of Gaussian", sigma_g_val, 1, 20)
        self.frac = ROOT.RooRealVar("sig_fraction_"+str(self.name), "Fraction", 0.4, 0.0, 1.0)
        
    def init_shape (self,x, alphaL=None, alphaR=None, nL=None, nR=None, bin_name=None, json_file="fitting_parameters.json"):
        print("---------------------------------------------")
        print("Signal : Double Sided Crystal Ball + Gaussian")
        self.x = x
        
        # Load alpha and n parameters from JSON if bin_name is provided
        if bin_name:
            params = load_fitting_parameters(json_file, bin_name)
            # Override the stored values with JSON values
            self.alphaL_val = get_param_value(params, "alphaL", self.alphaL_val)
            self.alphaR_val = get_param_value(params, "alphaR", self.alphaR_val)
            self.nL_val = get_param_value(params, "nL", self.nL_val)
            self.nR_val = get_param_value(params, "nR", self.nR_val)
        
        # Create or use provided RooRealVar objects
        if alphaL is None:
            self.alphaL = ROOT.RooRealVar("alphaL_"+str(self.name), "Left Tail Parameter", self.alphaL_val, 0.5, 8)  # Wider range
        else:
            self.alphaL = alphaL
            
        if alphaR is None:
            self.alphaR = ROOT.RooRealVar("alphaR_"+str(self.name), "Right Tail Parameter", self.alphaR_val, 0.5, 8)  # Wider range
        else:
            self.alphaR = alphaR
            
        if nL is None:
            self.nL = ROOT.RooRealVar("nL_"+str(self.name), "Left Tail Exponent", self.nL_val, 0.5, 15)  # Wider range
        else:
            self.nL = nL
            
        if nR is None:
            self.nR = ROOT.RooRealVar("nR_"+str(self.name), "Right Tail Exponent", self.nR_val, 0.5, 15)  # Wider range
        else:
            self.nR = nR
            
        self.set_constants(bin_name=bin_name, json_file=json_file)
            
        self.cb = ROOT.RooCrystalBall("Crystal Ball_"+str(self.name), "Crystal Ball Signal", self.x, self.mean, self.sigmaL, self.sigmaR, self.alphaL, self.nL, self.alphaR, self.nR)
        self.gauss = ROOT.RooGaussian("Gaussian_"+str(self.name), "Gaussian Signal 1", self.x, self.mean_g, self.sigma_g)
        self.sig_model = ROOT.RooAddPdf("Crystal Ball + Gaussian Signal_"+str(self.name), "Crystal Ball + Gaussian Signal", [self.cb, self.gauss], [self.frac])
        return self.sig_model, self.cb, self.gauss, self.frac, self.alphaL, self.alphaR, self.sigmaL, self.sigmaR

    def set_constants(self, params=None, bin_name=None, json_file="fitting_parameters.json"):
        self.constants = []

        # Combine parameters from JSON (if any) and from explicit input
        param_names = set()

        # From JSON file if bin_name is given
        if bin_name:
            json_params = load_fitting_parameters(json_file, bin_name)
            constants_str = json_params.get("constants", "")
            if constants_str != "":
                json_const_list = [name.strip() for name in constants_str.split(",")]
                param_names.update(json_const_list)

        # From explicit parameter list
        if params:
            param_names.update(params)
            
        for name in param_names:
        # Support names like 'tau_kpkm' when param is called self.tau_kpkm
            try:
                param_obj = getattr(self, name)
                if hasattr(param_obj, "setConstant"):
                    param_obj.setConstant(True)
                    self.constants.append(param_obj.GetName())
                    print(f"Set '{name}' to constant")
            except AttributeError:
                print(f"Warning: parameter '{name}' not found in self")
    
    def return_constants(self):
        return self.constants
        
    def get_param(self):
        return [self.mean_g, self.sigma_g, self.mean, self.sigmaL, self.sigmaR, self.alphaL, self.alphaR, self.nL, self.nR]
    
    def get_param_name(self):
        return [self.mean_g.GetName(), self.sigma_g.GetName(), self.mean.GetName(), self.sigmaL.GetName(), self.sigmaR.GetName(), self.alphaL.GetName(), self.alphaR.GetName(), self.nL.GetName(), self.nR.GetName()]
    
    def get_param_str(self):
        return ["#mu_{G}", "#sigma_{G}", "#mu", "#sigma_{L}", "#sigma_{R}", "#alpha_{L}", "#alpha_{R}", "n_{L}", "n_{R}"]

class Exponential():
    def __init__(self, name, decay, bin_name=None, json_file="fitting_parameters.json"):
        self.name = name
        self.decay = decay
        
        # Load tau parameter from JSON if bin_name is provided
        if bin_name:
            params = load_fitting_parameters(json_file, bin_name)
            tau_val = get_param_value(params, "tau", 0.001)
        else:
            tau_val = 0.001

        self.tau = ROOT.RooRealVar("tau_"+str(self.name), "Background slope", tau_val, -0.1, 0.1)

    def init_shape (self,x, bin_name=None, json_file="fitting_parameters.json"):
        print("Background : Exponential")
        self.x = x
        self.set_constants(bin_name=bin_name, json_file=json_file)
        return ROOT.RooExponential("Exponential Background "+str(self.name), "Exponential Background", self.x, self.tau)#, [self.cnstr_tau]
        
    def set_constants(self, params=None, bin_name=None, json_file="fitting_parameters.json"):
        self.constants = []

        # Combine parameters from JSON (if any) and from explicit input
        param_names = set()

        # From JSON file if bin_name is given
        if bin_name:
            json_params = load_fitting_parameters(json_file, bin_name)
            constants_str = json_params.get("constants", "")
            if constants_str:
                json_const_list = [name.strip() for name in constants_str.split(",")]
                print(json_const_list)
                param_names.update(json_const_list)

        # From explicit parameter list
        if params:
            param_names.update(params)
            
        for name in param_names:
        # Support names like 'tau_kpkm' when param is called self.tau_kpkm
            try:
                param_obj = getattr(self, name)
                if hasattr(param_obj, "setConstant"):
                    param_obj.setConstant(True)
                    self.constants.append(param_obj.GetName())
                    print(f"Set '{name}' to constant")
            except AttributeError:
                print(f"Warning: parameter '{name}' not found in self")

   
    def return_constants(self):
        return self.constants

    def get_param(self):
        return [self.tau]
    
    def get_param_name(self):
        return [self.tau.GetName()]
    
    def get_param_str(self):
        return ["#tau"]  
    
class Bernstein():
    def __init__(self, name, decay):
        self.name = name
        self.decay = decay
        self.a1 = ROOT.RooRealVar("a1_"+str(self.name), "a1", 1, 0, 10)
        self.a2 = ROOT.RooRealVar("a2_"+str(self.name), "a2", 1, 0, 10)
        self.a3 = ROOT.RooRealVar("a3_"+str(self.name), "a3",1, 0, 10)
        self.a4 = ROOT.RooRealVar("a4_"+str(self.name), "a4", 1, 0, 100)

    def init_shape (self,x):
        self.x = x
        # self.__init__()
        # if self.decay == "ksk" :
        return ROOT.RooBernstein("Bernstein_Background_"+str(self.name), "Bernstein_Background", self.x, [self.a1, self.a2, self.a3, self.a4])

    def set_constants(self, params):
        constants = []
        for param in params:
            if param == "a1_"+str(self.name):
                print("setting a1 to constant")
                self.a1.setConstant(True)
                constants.append(param)
            if param == "a2_"+str(self.name):
                print("setting a2 to constant")
                self.a2.setConstant(True)
                constants.append(param)
            if param == "a3_"+str(self.name):
                print("setting a3 to constant")
                self.a3.setConstant(True)
            if param == "a4_"+str(self.name):
                print("setting a4 to constant")
                self.a4.setConstant(True)
            
        return constants

    def get_param(self):
        return [self.a1, self.a2, self.a3, self.a4]
    
    def get_param_name(self):
        return [self.a1.GetName(), self.a2.GetName(), self.a3.GetName(), self.a4.GetName()]

    def get_param_str(self):
        return ["a1", "a2", "a3", "a4"]

class Cruijff:
    def __init__(self):
        self.mean = ROOT.RooRealVar("mean", "Mean Mass", 1968, 1925, 2025)
        self.alphaL = ROOT.RooRealVar("alphaL", "Left Tail Parameter", 1.5, 0.5, 5)
        self.alphaR = ROOT.RooRealVar("alphaR", "Right Tail Parameter", 1.5, 0.5, 5)
        self.sigmaL = ROOT.RooRealVar("sigmaL", "Left Tail Width", 10, 1, 20)
        self.sigmaR = ROOT.RooRealVar("sigmaR", "Right Tail Width", 10, 1, 20)
        self.mean1 = ROOT.RooRealVar("mean1", "Mean of Gaussian 1", 1968, 1925, 2025)
        self.sigma1 = ROOT.RooRealVar("sigma1", "Width of Gaussian 1", 10, 1, 20)
        self.frac = ROOT.RooRealVar("sig_fraction", "Fraction", 0.5, 0, 1)
        # self.frac.setConstant(True)
        
    def init_shape (self, x):
        self.x = x
        self.__init__()
        self.cruijff = ROOT.RooCruijff("Cruijff Signal", "Cruijff Signal", self.x, self.mean, self.sigmaL, self.sigmaR, self.alphaL, self.alphaR)
        self.gauss = ROOT.RooGaussian("Gaussian", "Gaussian Signal 1", self.x, self.mean1, self.sigma1)
        self.sig_model = ROOT.RooAddPdf("Cruijff + Gaussian Signal", "Cruijff + Gaussian Signal", [self.cruijff, self.gauss], [self.frac])
        self.eff_mean = self.frac.getVal() * self.mean.getVal() + (1 - self.frac.getVal()) * self.mean1.getVal()
        self.cr_width = np.sqrt((self.sigmaL.getVal()**2 + self.sigmaR.getVal()**2)/2)
        self.eff_width = np.sqrt(self.frac.getVal() * self.cr_width**2 + (1 - self.frac.getVal()) * self.sigma1.getVal()**2 + self.frac.getVal()*(1-self.frac.getVal())*(self.mean.getVal()-self.mean1.getVal())**2)
        return self.sig_model, self.cruijff, self.gauss, self.eff_mean, self.eff_width, self.frac
    
    def get_param(self):
        return [self.mean, self.sigmaL, self.sigmaR, self.alphaL, self.alphaR, self.mean1, self.sigma1]
    
    def get_param_str(self):
        return ["#mu", "#sigma_{L}", "#sigma_{R}", "#alpha_{L}", "#alpha_{R}", "#mu_{G}", "#sigma_{G}"]
    
class DoubleGaussian():
    def __init__(self,name,decay):
        self.name = name
        self.decay = decay
        self.mean1 = ROOT.RooRealVar("mean1_"+str(self.name), "Mean of Gaussian 1", 1968, 1925, 2025)
        self.sigma1 = ROOT.RooRealVar("sigma1_"+str(self.name), "Width of Gaussian 1", 10, 1, 20)
        self.mean2 = ROOT.RooRealVar("mean2_"+str(self.name), "Mean of Gaussian 2", 1968, 1925, 2025)
        self.sigma2 = ROOT.RooRealVar("sigma2_"+str(self.name), "Width of Gaussian 2", 5, 1, 20)
        self.frac = ROOT.RooRealVar("sig_fraction_"+str(self.name), "Fraction", 0.3, 0, 1)
        self.cnstr_sig1 = ROOT.RooGaussian("cnstr_sig1_"+str(self.name),"",self.sigma1,ROOT.RooFit.RooConst(9),ROOT.RooFit.RooConst(1))
        self.cnstr_sig2 = ROOT.RooGaussian("cnstr_sig2_"+str(self.name),"",self.sigma2,ROOT.RooFit.RooConst(6),ROOT.RooFit.RooConst(1))
        self.constraints = [self.cnstr_sig1,self.cnstr_sig2]
        
    def set_constants(self, params):
        constants = []
        for param in params:
            if param == "mean1_"+str(self.name):
                print("setting mean1 to constant")
                self.mean1.setConstant(True)
                constants.append(param)
        return constants
        
    def init_shape (self, x):
        self.x = x
        self.gauss1 = ROOT.RooGaussian("Gaussian 1_"+str(self.name), "Gaussian Signal 1", self.x, self.mean1, self.sigma1)
        self.gauss2 = ROOT.RooGaussian("Gaussian 2_"+str(self.name), "Gaussian Signal 2", self.x, self.mean2, self.sigma2)
        self.sig_model = ROOT.RooAddPdf("Double Gaussian Signal_"+str(self.name), "Double Gaussian Signal", [self.gauss1, self.gauss2], [self.frac])
        self.eff_mean = self.frac.getVal() * self.mean1.getVal() + (1 - self.frac.getVal()) * self.mean2.getVal()
        self.eff_width = np.sqrt(self.frac.getVal() * self.sigma1.getVal()**2 + (1 - self.frac.getVal()) * self.sigma2.getVal()**2)
        return self.sig_model, self.gauss1, self.gauss2, self.frac, self.constraints
    
    def get_param(self):
        return [self.mean1, self.sigma1, self.mean2, self.sigma2]
    
    def get_param_name(self):
        return [self.mean1.GetName(), self.sigma1.GetName(), self.mean2.GetName(), self.sigma2.GetName(), self.frac.GetName()]
    
    def get_param_str(self):
        return ["#mu_{1}", "#sigma_{1}", "#mu_{2}", "#sigma_{2}"]
               
class ExecuteFit :
    def __init__(self, var, bins, mass_range, x, cut, name="") : 
        self.name = name
        self.var = var
        self.bins = bins
        self.mass_range = mass_range
        self.x = x
        self.cut = cut

    def init_model (self, sig_model, bkg_model, data, hist, constraints, nsig) :
        self.sig_model = sig_model
        self.bkg_model = bkg_model
        self.data = data
        self.hist = hist
        self.constraints = constraints
        if self.hist.Integral() == 0: 
            logging.fatal("Empty histogram")
            exit()
        self.nbkg = ROOT.RooRealVar("nbkg_"+str(self.name), "Background Yield", self.hist.Integral() * 0.1, 1.0, 0.8 * self.hist.Integral())
        self.nsig = nsig
        self.background = ROOT.RooExtendPdf("background_"+str(self.name),"background",self.bkg_model, self.nbkg)
        self.signal = ROOT.RooExtendPdf("signal_"+str(self.name),"signal",self.sig_model, self.nsig)
        self.model = ROOT.RooAddPdf("model_"+str(self.name),"model",[self.signal, self.background])
        return self.model, self.signal, self.background
        
    def fit (self) :
        self.fit_result = self.model.fitTo(self.data, ROOT.RooFit.Save(True), ROOT.RooFit.Extended(), ROOT.RooFit.ExternalConstraints(self.constraints),ROOT.RooFit.Minimizer("Minuit")) # ROOT.RooFit.SumW2Error(True)
        return self.fit_result, self.nsig, self.nbkg
    
    def get_params(self):
        return [self.nsig, self.nbkg]
