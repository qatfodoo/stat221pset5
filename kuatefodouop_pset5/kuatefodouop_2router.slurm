#!/bin/bash
#SBATCH -J iid-em # name for job array
#SBATCH -o ./dump/q9task.out #Standard output
#SBATCH -e ./dump/q9task.err #Standard error
#SBATCH -p general #Partition
#SBATCH -t 10:00:00 #Running time of 10 hours.
#SBATCH --mem-per-cpu 10000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
#SBATCH --mail-user=kuatefodouop@fas.harvard.edu
#SBATCH --mail-type=ALL    # When to send mail

# argument = job id

Rscript kuatefodouop_2router.R $SLURM_ARRAY_TASK_ID
