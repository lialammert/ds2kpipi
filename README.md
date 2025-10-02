# ds2kpipi
Master Thesis @ CERN - Analysis for CP asymmetry measurement of Ds2Kpipi decay

Measurement of raw asymmetry on KsK sample.
Matching of kinematic distributions between Kpipi and KsK samples.
Mass fitting for the Ds+ and the Ds- of the reweighted KsK sample using a simultaneous fitting procedure.

## reweighting
Submit using `./sub_scripts/sub_jobs_reweighting.py`, calls `./weighting/EXE_reweighting.py`
to modify : fiducial cuts, weighting strategy, binning

## fitting
Submit using `./sub_scripts/sub_jobs/fitting.py`, calls `./fitting/EXE_fitting.py`
to modify : fitting pdfs used

## extras
### Ds production asymmetry
Execute `./scripts/prod_acp.py` which will call `./scripts/per_bin_prod_asym.py` to compute the production asymmetry of the Ds in each bin.
### Difference in uncertainties
Execute `./scripts/uncertainties.py` to compare the uncertainty after weighting on kpipi and ksk with and without considering the Ds production asymmetry. The difference should tend to zero as the sample is reweighted.
### Dalitz plot
Execute `./scripts/dalitz_plot.py` to plot the Dalitz plot of the kpipi sample, with bins or not.


