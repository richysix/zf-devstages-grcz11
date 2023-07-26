#!/usr/bin/env bash
# rscript-array.sh - Script to run a file of Rscript commands as an array job
# It loads the R module and then runs a single line based on the TASK_ID
# from a supplied file. Filename is the argument to the script.
# The second argument is an optional job name for the success/failure message
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G
#$ -t 1-10

source bash_functions.sh

module load R/4.2.1

if [[ -z $2 ]]; then
    JOB=${JOB_NAME}
else
    JOB=$2
fi

line=`sed "${SGE_TASK_ID}q;d" $1`
eval $line
SUCCESS=$?

verbose=1
error_checking $SUCCESS "job ${JOB}, task ${SGE_TASK_ID} SUCCEEDED." "job ${JOB}, task ${SGE_TASK_ID} FAILED: $SUCCESS"
exit $SUCCESS

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2023. Queen Mary University of London.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007
