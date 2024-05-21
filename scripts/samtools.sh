#!/usr/bin/env bash
# samtools.sh - Script to use samtools to make bam files
# that contain only numeric chromosomes
# Input to the script is a line number to index into the samples.tsv file
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G

module purge
module load samtools/1.20

line=`sed "${1}q;d" samples.tsv`
sample=`echo $line | awk '{ print $1 }'`
samtools view -b -o $sample.chr.bam ../star2/$sample/Aligned.sortedByCoord.out.bam `seq 25`
