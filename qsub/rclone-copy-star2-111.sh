#!/usr/bin/env bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -o rclone-copy-star2.o
#$ -e rclone-copy-star2.e

module load rclone
rclone copy --verbose --filter "+ *cram*" --filter "- *" star2/ sharepoint-qmul-buschlab:Projects/zf-stages-grcz11/111/star2/
rclone copy --verbose /data/scratch/bty114/zf-stages-grcz11/111/reference/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa sharepoint-qmul-buschlab:Projects/zf-stages-grcz11/111/

rclone copy --verbose deseq2-all/ sharepoint-qmul-buschlab:Projects/zf-stages-grcz11/111/deseq2-all/
rclone copy --filter "+ SJ.star*" --filter "- *" --verbose star2/ sharepoint-qmul-buschlab:Projects/zf-stages-grcz11/111/star2/
