import ROOT

"""
To apply fiducial kinematic cuts on kpipi and/or ksk files.
Change input/output folder names as needed.
Add or change variables being cut on, formulas, ranges.
"""


tos = "kp"

input_folder = "Dalitz_Bins_Kp_TOS"
output_folder = "FCUT_Dalitz_Bins"

for file in ["ksk", "kpipi"]:
    if file == "ksk": 
        print("Applying cuts for KSK")
        ns = 1
    if file == "kpipi":
        print("Applying cuts for Kpipi")
        ns = 20
    for n in range(ns):
        print(f"file number {n}")
        if tos == "kp" : 
            if file == "ksk":
                input_file = f"./{input_folder}/hlt1_kp_sweights.root"
            if file == "kpipi":
                input_file = f"./{input_folder}/dalitz_bin_{n}_hlt1_kp.root"
        tree_name = "DecayTree"
        df = ROOT.RDataFrame(tree_name, input_file)
        df = df.Define("Dp_PT", "TMath::Sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY)")
        df = df.Define("Dp_P", "TMath::Sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY+ Dp_PZ*Dp_PZ)")
        df = df.Define("Dp_ETA", "atanh(Dp_PZ/sqrt(Dp_PX*Dp_PX + Dp_PY*Dp_PY + Dp_PZ*Dp_PZ))")
        df = df.Define("pip_ETA", "atanh(pip_PZ/sqrt(pip_PX*pip_PX + pip_PY*pip_PY + pip_PZ*pip_PZ))")
        df = df.Define("pim_ETA", "atanh(pim_PZ/sqrt(pim_PX*pim_PX + pim_PY*pim_PY + pim_PZ*pim_PZ))")
        df = df.Define("Dp_PHI", "TMath::ATan2(Dp_PY, Dp_PX)")
        if file == "ksk":
            df = df.Define("hp_PT", "TMath::Sqrt(hp_PX*hp_PX + hp_PY*hp_PY)")
            df = df.Define("hp_P", "TMath::Sqrt(hp_PX*hp_PX + hp_PY*hp_PY + hp_PZ*hp_PZ)")
            df = df.Define("hp_ETA", "atanh(hp_PZ/sqrt(hp_PX*hp_PX + hp_PY*hp_PY + hp_PZ*hp_PZ))")
            df = df.Define("hp_PHI", "TMath::ATan2(hp_PY, hp_PX)")
        if file == "kpipi":
            df = df.Define("Kp_PT", "TMath::Sqrt(Kp_PX*Kp_PX + Kp_PY*Kp_PY)")
            df = df.Define("hp_P", "TMath::Sqrt(Kp_PX*Kp_PX + Kp_PY*Kp_PY + Kp_PZ*Kp_PZ)")
            df = df.Define("hp_ETA", "atanh(Kp_PZ/sqrt(Kp_PX*Kp_PX + Kp_PY*Kp_PY + Kp_PZ*Kp_PZ))")
            df = df.Define("hp_PHI", "TMath::ATan2(Kp_PY, Kp_PX)")
            
        # Define cuts
        # mass_cut = "(Dp_M > 1920) && (Dp_M < 2025)"
        # kpipi_cut = "(pippim_Hlt1TwoTrackMVADecision_TOS)"
        # ksk_cut = "(KS0_Hlt1TwoTrackMVADecision_TOS || KS0_Hlt1TwoTrackKsDecision_TOS || KS0_Hlt1TrackMVADecision_TOS)"
        
        if tos == "kp":
            eta = "(hp_ETA > 2.25 && hp_ETA < 4.5) && (pip_ETA > 2.25 && pip_ETA < 4.5) && (pim_ETA > 2.25 && pim_ETA < 4.5) && (Dp_ETA > 2 && Dp_ETA < 4.3) "
            p = "&& (hp_P > 10e3 && hp_P < 120e3) && (Dp_P > 20e3 && Dp_P < 220e3)"
            pt = "&& (hp_PT > 2.8e3 && hp_PT < 12e3) && (Dp_PT > 2.7e3 && Dp_PT < 15e3)"
        if tos == "pippim" :
            eta = "(hp_ETA > 2.25 && hp_ETA < 4.5) && (Dp_ETA > 2.25 && Dp_ETA < 4.35)"
            p = "&& (hp_P > 3e3 && hp_P < 100e3) && (Dp_P > 15e3 && Dp_P < 225e3)"
            pt = "&& (hp_PT > 300 && hp_PT < 5.5e3) && (Dp_PT > 2.5e3 && Dp_PT < 15e3)"
        
        cut = eta + p + pt
        df_filtered = df.Filter(cut)

        # Save to new ROOT file
        if file == "ksk":
            output_file = f"./{output_folder}/fcut_hlt1_kp_sweights.root"
        if file == "kpipi":
            output_file = f"./{output_folder}/fcut_dalitz_bin_{n}_hlt1_kp.root"
        df_filtered.Snapshot(tree_name, output_file)
