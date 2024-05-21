#!/usr/bin/env bash
# qorts2.sh - Script to run 2nd QoRTs job
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G

source bash_functions.sh

module load QoRTs
Rscript qorts2.R

SUCCESS=$?
verbose=1
error_checking $SUCCESS "job qorts2 succeeded." "job qorts2 failed: $SUCCESS"

