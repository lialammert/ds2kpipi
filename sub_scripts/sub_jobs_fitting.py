import subprocess
from weighting.weighting_strategies import STRATEGIES

"""
Submission script for mass fitting for different number of weighting iterations.
Strategy (and potentially bin number) should be changed manually in the script.
This script calls a bash script `send_job.sh` that handles the actual job submission to the cluster or wherever.
"""

if __name__ == '__main__':
    strategy = "Vega_may25"
    max_iter = 8
    name = "Vega_sWeights"
    decay = "ksk"
    bins = [0,1,2,3,4,5,6,7,8,9,10,11,15,16]
    weighting_strat = STRATEGIES[strategy](name, max_iter)
    
    for bin in bins :
        print(f"Submitting the fitting for decay {decay}, bin {bin}")
        result = subprocess.run(f"sbatch send_job.sh {strategy} {name} {max_iter} {decay} {bin}", shell=True, capture_output=True, text=True)
        job_id = result.stdout.strip().split()[-1]
        print(f"Submitted job {job_id}")