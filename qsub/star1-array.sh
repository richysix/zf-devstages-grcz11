#!/usr/bin/env bash
# star1-array.sh - Script to run STAR as an array job
# see star1.sh for details
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=24G
#$ -t 1-96

USAGE="star1-array.sh"

source bash_functions.sh

star1.sh ${SGE_TASK_ID}
SUCCESS=$?

verbose=1
error_checking $SUCCESS "job star1, task ${SGE_TASK_ID} succeeded." "job star1, task ${SGE_TASK_ID} failed: $SUCCESS"

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
