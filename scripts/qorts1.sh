#!/usr/bin/env bash
# qorts1.sh - Script to run the first stage of QoRTs
# Input to the script is a line number to index into the samples.tsv file
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=240:0:0
#$ -l h_vmem=24G

module purge
module load QoRTs

line=`sed "${1}q;d" samples.tsv`
sample=`echo $line | awk '{ print $1 }'`
condition=`echo $line | awk '{ print $2 }'`
reads=`echo $line | awk '{ print $3 }'`
mkdir $sample
qorts QC \
--stranded \
--chromSizes chr-sizes.tsv \
--genomeFA ../reference/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa \
--seqReadCt $reads \
--title $sample \
$sample.chr.bam \
Danio_rerio.GRCz11.111.chr.gtf \
$sample
