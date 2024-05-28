#!/usr/bin/env bash
# bam2cram-array.sh - Script to run bam2cram as array job
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -t 1-10

source bash_functions.sh

line=`sed "${SGE_TASK_ID}q;d" $1`
bam=`echo $line | awk '{ print $1 }'`
ref=`echo $line | awk '{ print $2 }'`
cram=`echo $line | awk '{ print $3 }'`

sh bam2cram.sh -v $bam $ref $cram

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
