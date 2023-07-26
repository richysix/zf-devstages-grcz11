#!/usr/bin/env bash
# star1.sh - Script to run STAR on a set of FASTQ files
# expects a file name fastq.tsv to exist in the working directory
# This should contains three columns
# sample names
# comma-separated list of read1 fastq filenames
# comma-separated list of read2 fastq filenames
# Input to the script is a line number to index into the fastq.tsv file
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=24G

USAGE="star1.sh LINE_NUM"

line=`sed "${1}q;d" fastq.tsv`
sample=`echo $line | awk '{ print $1 }'`
fastq1=`echo $line | awk '{ print $2 }'`
fastq2=`echo $line | awk '{ print $3 }'`
mkdir $sample

module load STAR
STAR \
--runThreadN 1 \
--genomeDir ../reference/grcz11 \
--readFilesIn $fastq1 $fastq2 \
--readFilesCommand zcat \
--outFileNamePrefix $sample/ \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate

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
