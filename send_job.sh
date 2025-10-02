#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=24:00:00
#SBATCH --nodelist=lphepc122
#SBATCH --job-name=2d_W
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

## simultaneous fitting
# python /home/lammert/ds2kpipi/fitting/EXE_fitting.py -s "$1" -na "$2" -mi "$3" -d "$4" -b "$5"

## 2d weighting_separated
# python /home/lammert/ds2kpipi/weighting/EXE_reweighting.py -i "$1" -n "$2" -s "$3" -na "$4" -mi "$5"

## Ds prod asym
# python /home/lammert/ds2kpipi/scripts/prod_acp.py


