#!/bin/bash

## Call the observed-CCC simulated annealing script 100 times.  
## Each run of the script yields one SA-iteration per all 45 region pairs. 

## Nhung Hoang, June 2021 

#SBATCH --mail-user=nhung.hoang@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --account=rubinov_lab
#SBATCH --time=6:00:00
#SBATCH --mem=2GB
#SBATCH --array=0-99
#SBATCH --output=/data/rubinov_lab/brain_genomics_project/platypus/outs_ccc32/%a.out

ml GCC/8.2.0 OpenMPI/3.1.4 Python/3.7.2
source /data/rubinov_lab/brain_genomics_project/platypus/venv_analyses/bin/activate
python -u /data/rubinov_lab/brain_genomics_project/platypus/run_ccc_sa_32.py ${SLURM_ARRAY_TASK_ID} 
deactivate 
