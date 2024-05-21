#!/usr/bin/env bash
# samtools-array.sh - Script to run samtools as an array job
# see samtools.sh for details
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -t 1-96

USAGE="samtools-array.sh"

source bash_functions.sh

./samtools.sh ${SGE_TASK_ID}
SUCCESS=$?

verbose=1
error_checking $SUCCESS "job samtools, task ${SGE_TASK_ID} succeeded." "job samtools, task ${SGE_TASK_ID} failed: $SUCCESS"

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
