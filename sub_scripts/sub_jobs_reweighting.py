import os
import subprocess
import time
from weighting.weighting_strategies import STRATEGIES

"""
Submission script for reweighting jobs.
Runs on each bin and create .root weights ntuples for ksk and kpipi data.
When all jobs are completed, hadds the weights together and cleans up temporary directories.
Strategy can be changed as needed.
This script calls a bash script `send_job.sh` that handles the actual job submission to the cluster or wherever.
"""

if __name__ == '__main__':
    kpipi_files = [
    "dalitz_bin_0_hlt1_kp.root",
    "dalitz_bin_1_hlt1_kp.root",
    "dalitz_bin_2_hlt1_kp.root",
    "dalitz_bin_3_hlt1_kp.root",
    "dalitz_bin_4_hlt1_kp.root",
    "dalitz_bin_5_hlt1_kp.root",
    "dalitz_bin_6_hlt1_kp.root",
    "dalitz_bin_7_hlt1_kp.root",
    "dalitz_bin_8_hlt1_kp.root",
    "dalitz_bin_9_hlt1_kp.root",
    "dalitz_bin_10_hlt1_kp.root",
    "dalitz_bin_11_hlt1_kp.root",
    "dalitz_bin_12_hlt1_kp.root",
    "dalitz_bin_13_hlt1_kp.root",
    "dalitz_bin_14_hlt1_kp.root",
    "dalitz_bin_15_hlt1_kp.root",
    "dalitz_bin_16_hlt1_kp.root",
    "dalitz_bin_17_hlt1_kp.root",
    "dalitz_bin_18_hlt1_kp.root",
    "dalitz_bin_19_hlt1_kp.root"
    ]
    
    job_ids = []
    
    strategy = "Vega_may25"
    max_iter = 8
    name = "Vega_sWeights"
    weighting_strat = STRATEGIES[strategy](name, max_iter)  # Change to Vega_may25, Lyra_may25, Draco_jun25, Kepler_jun25 as needed
    folder = weighting_strat.name
    w = weighting_strat.max_iter
    
    for i in range(w+1):
        os.system(f'mkdir ./weights_ksk_{i}')
        os.system(f'mkdir ./weights_kpipi_{i}')

    for i, k in enumerate(kpipi_files):
            
        print(f"Submitting for bin {i}")
        result = subprocess.run(f"sbatch send_job.sh {k} {i} {strategy} {name} {max_iter}", shell=True, capture_output=True, text=True)
        job_id = result.stdout.strip().split()[-1]
        print(f"Submitted job {job_id}")
        job_ids.append(job_id)
        
    while True:
        still_running = False
        for job_id in job_ids:
            check = subprocess.run(f"squeue -j {job_id}", shell=True, capture_output=True, text=True)
            if job_id in check.stdout:
                still_running = True
                break
        if not still_running:
            break
        time.sleep(15)  # Wait before polling again
        
    print("All jobs completed. Assembling weights.")
    for i in range(w+1):
        os.system(f'hadd -f ./{folder}/{i}_iter/weights_ksk_{i}.root ./weights_ksk_{i}/weights_ksk_*.root')
        os.system(f'hadd -f ./{folder}/{i}_iter/weights_kpipi_{i}.root ./weights_kpipi_{i}/weights_kpipi_*.root')
        os.system(f'rm -r ./weights_ksk_{i}')
        os.system(f'rm -r ./weights_kpipi_{i}')
