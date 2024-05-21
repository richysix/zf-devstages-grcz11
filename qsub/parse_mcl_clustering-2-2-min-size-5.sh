#!/usr/bin/env bash
# parse_mcl_clustering-2-2-min-size-5.sh - Run parse_mcl_output.pl on 2.2
# clustering and remove clusters smaller then 5
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G

module load Perl/5.43.0

GENE_NAME_FILE='gene-names.tsv'
uid=$( echo $RANDOM | md5sum | head -c 6 )
head -n1 $GENE_NAME_FILE > gene-info-$uid.tsv
join -t$'\t' <( sort -t$'\t' -k1,1 $GENE_NAME_FILE ) \
<( awk '{print $2}' all-cor-70-knn400.tab | sort -u ) >> gene-info-$uid.tsv

OUTPUT_BASE='all-cor-70-knn400'
INFLATION_SUFFIX='2-2'
perl ~/checkouts/bioinf-gen/parse_mcl_output.pl \
--info_file gene-info-$uid.tsv \
--info_file_node_id_col 1 \
--info_file_node_name_col 2 \
--graphml_file $OUTPUT_BASE.mci.I${INFLATION_SUFFIX}-min5.graphml \
--cluster_size_threshold 5 \
$OUTPUT_BASE.tab $OUTPUT_BASE.mci \
$OUTPUT_BASE.mci.I${INFLATION_SUFFIX} \
$OUTPUT_BASE.mci.I${INFLATION_SUFFIX}-min5.tsv

# remove gene info file
rm gene-info-$uid.tsv

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
