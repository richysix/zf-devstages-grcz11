#!/usr/bin/env bash
# qorts1-array.sh - Script to run qorts1 as an array job
# see qorts1.sh for details
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=240:0:0
#$ -l h_vmem=24G
#$ -t 1-96

USAGE="qorts1-array.sh"

source bash_functions.sh

./qorts1.sh ${SGE_TASK_ID}
SUCCESS=$?

verbose=1
error_checking $SUCCESS "job qorts1, task ${SGE_TASK_ID} succeeded." "job qorts1, task ${SGE_TASK_ID} failed: $SUCCESS"

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
