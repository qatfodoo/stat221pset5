#!/bin/bash
#SBATCH -J iid-em # name for job array
#SBATCH -o ./dump/q4task.out #Standard output
#SBATCH -e ./dump/q4task.err #Standard error
#SBATCH -p general #Partition
#SBATCH -t 10:00:00 #Running time of 10 hours.
#SBATCH --mem-per-cpu 10000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
#SBATCH --mail-user=kuatefodouop@fas.harvard.edu
#SBATCH --mail-type=ALL    # When to send mail

# argument = job id

Rscript kuatefodouop_pset5_q4_task.R $SLURM_ARRAY_TASK_ID
