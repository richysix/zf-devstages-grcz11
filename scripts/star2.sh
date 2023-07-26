#!/usr/bin/env bash
# star2.sh - Script to run STAR on a set of FASTQ files
# expects a file name fastq.tsv to exist in the working directory
# This should contains three columns
# sample names
# comma-separated list of read1 fastq filenames
# comma-separated list of read2 fastq filenames
# This script is for the second pass of the two-pass mapping strategy with STAR
# The results of the first pass should be in a directory names star1
# Input to the script is a line number to index into the fastq.tsv file
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=24G

USAGE="star2.sh LINE_NUM"

# expects a file named fastq.tsv
if [[ ! -e fastq.tsv ]]; then
    echo "File fastq.tsv not found!"
    exit 2
fi

line=`sed "${1}q;d" fastq.tsv`
sample=`echo $line | awk '{ print $1 }'`
fastq1=`echo $line | awk '{ print $2 }'`
fastq2=`echo $line | awk '{ print $3 }'`
mkdir -p $sample

module load STAR
STAR \
--runThreadN 1 \
--genomeDir ../reference/grcz11 \
--readFilesIn $fastq1 $fastq2 \
--readFilesCommand zcat \
--outFileNamePrefix $sample/ \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--sjdbFileChrStartEnd `find ../star1 | grep SJ.out.tab$ | sort | tr '\n' ' '`


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
