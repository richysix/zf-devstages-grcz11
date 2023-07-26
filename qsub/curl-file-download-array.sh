#!/usr/bin/env bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=100M
#$ -t 1-3

sh ../scripts/curl-file-download.sh "curl.${SGE_TASK_ID}"
