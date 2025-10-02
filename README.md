This README explains the use of the source files and scripts developed during my Master's thesis at CERN in the LHCb collaboration. The aim of the work was to contribute to the measurements of Charge-Parity (CP) asymmetries in the charm sector, more specifically in the Cabibbo suppressed Ds -> Kpipi decay. This was done using Run 3 (2024) data from the LHCb experiment at CERN.
This thesis focused on the use of a control channel (Ds -> KSK) to compute nuisance asymmetries using statistical tools for data analysis, mainly pyROOT.


# Computing the raw asymmetry of Ds->KsK with Run3 data

Files used : 
    Block 1 Dalitz binned KTOS kpipi files - `/eos/lhcb/lbdt3/user/njurik/Ds2Kpipi/data/block1/ds2kpipi/Hlt1Split_sWeights/dalitz_bin_*_hlt1_kp.root`
    ksk file - `/eos/lhcb/lbdt3/user/njurik/Ds2Kpipi/data/block1/ds2ksk/SimpleSelection_sWeights/hlt1_kp_sweights.root`

In all source files and scripts used, a small text on the top explains what the classes/functions do.

## Reweighting
The reweighting uses kinematic variables in the KsK and Kpipi datasets to match the sWeighted distributions. 
Source files : `./ds2kpipi/run3/acp_ds2kpipi/lia_ksk/weighting/`
    The different weighting strategies (can adapt number of iterations, variables being reweighted, binning) : `weighting_strategies.py`
    To perform fiducial (kinematic) cuts on the kpipi and ksk files : `fiducial_cuts.py`
    The weighting algorithm is performed in 5 steps : 
        - create the datasets using RDF : `step1_create_datasets.py`
        - find weights per bin by taking the ratio of histograms : `step2_weights_per_bin.py`
        - get weights per value by matching to event : `step3_weights_per_val.py`
        - save weights in ROOT files : `step4_save_weights.py` - one ROOT file per iteration which contains one TTree per Dalitz bin.
        - plot histogram at the end of the iteration for (ETA,PT,P,PHI) of the kaon and the Ds : `step5_plot_histograms.py` (`lhcbStyle.C` file needed, found in `./extras`)
    The weighting is executed from the `EXE_reweighting.py` file (to adapt).

Scripts : `./ds2kpipi/run3/analysis/lia_ksk/`
    The submission file is `sub_jobs_reweighting.py` in which 3 parameters need to be changed : 
        - the weighting strategy chosen: `strategy = "Vega_may25"`
        - the maximum number of iterations performed (corresponds to number of weight files saved) : `max_iter = 8`
        - the name given to the folder : `name = "Vega_sWeights"`
    A bash (.sh) script is then submitted, containing the line : `python EXE_reweighting.py -i "$1" -n "$2" -s "$3" -na "$4" -mi "$5"` and calls the file that executes the weighting algorithm.
    The submission files automatically merges the root files from individual bins after waiting until all the reweighting is finished.

## Mass Fitting
This is a simultaneous fitting procedure for the Ds+ and Ds- mass ditributions. 
Source files : `./ds2kpipi/run3/acp_ds2kpipi/lia_ksk/fitting/`
    The fitting is performed using different PDFs defined in :`fitting.py`  
        default is signal : Crystal Ball + Gaussian // background : Negative Exponential
    The SimultaneousFitter class is in `simulataneous_fitting.py``
    The executing file is `EXE_fitting.py` (to adapt).

Scripts : `./ds2kpipi/run3/analysis/lia_ksk/`
    The submission file is `sub_jobs_fitting.py` and takes parameters : 
        - the weighting strategy chosen: `strategy = "Vega_may25"`
        - the maximum number of iterations performed (corresponds to number of weight files saved) : `max_iter = 8`
        - the name given to the folder : `name = "Vega_sWeights"`
        - the decay being fitted `decay = "ksk"` (ksk is default, doesn't need to be changed)
    A bash (.sh) script is then submitted, containing the line : `python EXE_fitting.py -s "$1" -na "$2" -mi "$3" -d "$4" -b "$5"` and calls the file that executes the simultaneous fitting.

## Ds Production Asymmetry
This allows to compute the production asymmetry of the Ds using previously measured values.
Source file : `./ds2kpipi_cpv/run3/acp_ds2kpipi/lia_ksk/per_bin_prod_asym.py` 
    Computes the normalised entries of the binned Kpipi and weighted KsK samples in the bins of the Ds production asymmetries measurements.
Script : `./ds2kpipi_cpv/run3/analysis/lia_ksk/compute_aprod.py`
    Computes and saves the central values and errors of the Ds production asymmetry for each reweighting.

