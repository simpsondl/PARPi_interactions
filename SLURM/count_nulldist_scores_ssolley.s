#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=16:00:00
#SBATCH --mem=8G
#SBATCH --job-name="count_nulldist"
#SBATCH --output=/Genomics/grid/users/ds65/array_logs/%x-%A_%a.out

#######################################
### Submit job with following command:
###
### sbatch --array=M-N --export=NUM=num,REC=record,INTERACTIONSCORES=interactionscores,SEARCH=search count_nulldist_scores.s
###
### where
###### M,N are the job indices to create
###### NUM is the number of records to process in each job
###### REC is the line number of starting gene pair minus 1
###### INTERACTIONSCORES is the file containing scores to get counts for
###### SEARCH is file path for null distribution file to check
###
########################################

module purge

pwd; date

NUM_START=$(( $REC + ($SLURM_ARRAY_TASK_ID - 1) * $NUM + 1 ))
NUM_END=$(( $REC + $SLURM_ARRAY_TASK_ID * $NUM ))

read -d '' cmd << EOF
~/scripts/sort_pval_awk_optimized_ssolley.sh $NUM_START $NUM_END $INTERACTIONSCORES $SEARCH
EOF

echo ${cmd}
eval ${cmd}

date
