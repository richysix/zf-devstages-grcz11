#!/usr/bin/env bash
# indexbam-array.sh - Script to run indexbam as an array job
# see indexbam.sh for details
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -t 1-96

USAGE="indexbam-array.sh"

source bash_functions.sh

./indexbam.sh ${SGE_TASK_ID}
SUCCESS=$?

verbose=1
error_checking $SUCCESS "job indexbam, task ${SGE_TASK_ID} succeeded." "job indexbam, task ${SGE_TASK_ID} failed: $SUCCESS"

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
