#!/usr/bin/env bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=100M
#$ -o rclone-copy-star2.o
#$ -e rclone-copy-star2.e

module load rclone
rclone copy --filter "+ *cram*" --filter "- *" star2/ sharepoint-qmul:Documents/Projects/zf-stages-grcz11/e109/star2/
rclone copy /data/SBBS-BuschLab/genomes/GRCz11/GRCz11.fa sharepoint-qmul:Documents/Projects/zf-stages-grcz11/e109/
