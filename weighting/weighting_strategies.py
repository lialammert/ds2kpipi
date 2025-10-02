import math
import numpy as np
import copy

"""
Different weighting strategies for reweighting datasets.
Compute binning for the variables used in the reweighting.
Default number of bins is 30 for both weighting and plotting.
Can be adjusted by passing `nbins_weigh` and `nbins_plot` parameters.
If additional variables are needed, they can be added through the `add_vars` parameter.
Create additional strategies by subclassing `WeightingStrategy` and implementing the `get_variables` method.
Name of the strategy should be unique, with month and year eg. '_may25'
"""

def compute_binning(weighting_variables, nbins_weigh=30, nbins_plot=30, add_vars = None) :
    kinematic_vars = ["hp_P","hp_ETA","hp_PHI","hp_PT","Dp_P","Dp_ETA","Dp_PHI","Dp_PT"] 
    formula_vars_tree_kpipi = ["TMath::Sqrt(Kp_PX*Kp_PX + Kp_PY*Kp_PY + Kp_PZ*Kp_PZ)","0.5*TMath::Log((hp_P + Kp_PZ)/(hp_P - Kp_PZ))", "TMath::ATan2(Kp_PY, Kp_PX)", "TMath::Sqrt(Kp_PX*Kp_PX + Kp_PY*Kp_PY)", "TMath::Sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY + Dp_PZ*Dp_PZ)", "0.5*TMath::Log((Dp_P + Dp_PZ)/(Dp_P - Dp_PZ))", "TMath::ATan2(Dp_PY, Dp_PX)", "TMath::Sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY)"]
    formula_vars_tree_ksk = ["TMath::Sqrt(hp_PX*hp_PX + hp_PY*hp_PY + hp_PZ*hp_PZ)","0.5*TMath::Log((hp_P + hp_PZ)/(hp_P - hp_PZ))", "TMath::ATan2(hp_PY, hp_PX)", "TMath::Sqrt(hp_PX*hp_PX + hp_PY*hp_PY)", "TMath::Sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY + Dp_PZ*Dp_PZ)", "0.5*TMath::Log((Dp_P + Dp_PZ)/(Dp_P - Dp_PZ))", "TMath::ATan2(Dp_PY, Dp_PX)", "TMath::Sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY)"]
    
    eta_range = [2.25, 4.5] 
    p_range = [10e3, 120e3]
    pt_range = [2.5e3, 12e3]
    phi_range = [-math.pi-0.05, math.pi+0.05]
    d_eta_range = [2, 4.5]
    d_p_range = [20e3, 220e3]
    d_pt_range = [2.7e3, 15e3]
    d_phi_range = [-math.pi-0.05, math.pi+0.05]
    #here add ranges for new variables

    variable_ranges = [p_range, eta_range, phi_range, pt_range, d_p_range, d_eta_range, d_phi_range, d_pt_range]
    
    if add_vars is not None:
        kinematic_vars += add_vars[1]
        formula_vars_tree_kpipi += add_vars[0]
        formula_vars_tree_ksk += add_vars[1]
        variable_ranges += []  # Add ranges for new variables

    weight_vars_flat = [x for pair in weighting_variables for x in pair]
    all_variables = copy.deepcopy(kinematic_vars)
    all_variables += [x for x in weight_vars_flat if x not in kinematic_vars]  # Ensure no duplicates

    bin_edges_weigh = [0]*len(all_variables)
    bin_edges_plot = [0]*len(all_variables)
    for ik in range(len(all_variables)) :
        bin_edges_weigh[ik] = np.array(np.linspace(variable_ranges[ik][0], variable_ranges[ik][1], nbins_weigh+1))
        bin_edges_plot[ik] = np.array(np.linspace(variable_ranges[ik][0], variable_ranges[ik][1], nbins_plot+1))

    return kinematic_vars, formula_vars_tree_kpipi, formula_vars_tree_ksk,variable_ranges, bin_edges_weigh, bin_edges_plot

class WeightingStrategy():
    
    def __init__(self, name, max_iter, first_coarse = False, order=None, binning=None):
        self.name = name
        self.max_iter = max_iter
        self.first_coarse = first_coarse  # If True, use coarse binning for the first iteration
        self.order = [0,0,0]
        if binning is not None : 
            self.binning = binning  
        else :
            self.binning = compute_binning(self.get_variables(), nbins_weigh=30, nbins_plot=30, add_vars=None)
                                        #    add_vars=[['Kp_PX', 'Kp_PY', 'Dp_PX', 'Dp_PY'],['hp_PX','hp_PY', 'Dp_PX', 'Dp_PY']])
        if self.first_coarse : self.coarse_binning = compute_binning(self.get_variables(), nbins_weigh=5, nbins_plot=5)

class Vega_may25(WeightingStrategy):
    """ 
        Vega weighting strategy : 
        D (ETA, PT) -> K (ETA, PT) -> K (ETA, PHI)
        from ksk to kpipi
    """
    def get_variables(self):
        return [['Dp_ETA', 'Dp_PT'], ['hp_ETA', 'hp_PT'], ['hp_ETA','hp_PHI']]
        
class Lyra_may25(WeightingStrategy):
    """ 
        Lyra weighting strategy : 
        K(ETA, PT) -> K (ETA, PHI) -> D (ETA, PT)
        from ksk to kpipi
    """
    
    def get_variables(self):
        return [['hp_ETA', 'hp_PT'], ['hp_ETA','hp_PHI'],['Dp_ETA', 'Dp_PT']]
    
class Draco_jun25(WeightingStrategy):
    """ 
        Draco weighting strategy : 
        D (ETA, PT)
        from ksk to kpipi
    """
    def get_variables(self):
        return [['Dp_ETA', 'Dp_PT']]
    
class Kepler_jun25(WeightingStrategy):
    """ 
        Kepler weighting strategy : 
        K (ETA, PT) -> K (ETA, PHI)
        from ksk to kpipi
    """
    
    def get_variables(self):
        return [['hp_ETA', 'hp_PT'], ['hp_ETA', 'hp_PHI']]
    
class Ursa_jun25(WeightingStrategy):
    """ 
        Ursa weighting strategy : 
        D (ETA, PT) 
        from kpipi to ksk
        -> K (ETA, PT) -> K (ETA, PHI)
        from ksk to kpipi
    """
    def __init__(self, name, max_iter, first_coarse = False, order=None, binning=None):
        self.name = name
        self.max_iter = max_iter
        self.order = [0,0,0]
        self.first_coarse = first_coarse  # If True, use coarse binning for the first iteration
        if binning is not None : 
            self.binning = binning  
        else :
            self.binning = compute_binning(self.get_variables(), nbins_weigh=30, nbins_plot=30)
        if self.first_coarse : self.coarse_binning = compute_binning(self.get_variables(), nbins_weigh=5, nbins_plot=5)
    
    def get_variables(self):
        return [['Dp_ETA', 'Dp_PT'], ['hp_ETA', 'hp_PT'], ['hp_ETA','hp_PHI']]
    
class Aries_jul25(WeightingStrategy):
    """ 
        Aries weighting strategy : 
        K (ETA, PX) -> K (ETA, PY) -> D (ETA, PT)
        from kpipi to ksk
    """
    def get_variables(self):
        return [['hp_ETA', 'hp_PX'], ['hp_ETA','hp_PY'],['Dp_ETA', 'Dp_PT']]

STRATEGIES = {
    "Vega_may25": Vega_may25,
    "Lyra_may25": Lyra_may25,
    "Draco_jun25": Draco_jun25,
    "Kepler_jun25": Kepler_jun25,
    "Ursa_jun25": Ursa_jun25,
    "Aries_jul25": Aries_jul25
}