#!/usr/bin/env bash
# run-mcl-clustering.sh - Script to create a clustered network
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=4G
#$ -o mcl-clustering.o
#$ -e mcl-clustering.e

source bash_functions.sh

USAGE="Script to create a clustered correlation network

run-mcl-clustering.sh [options] input_file gene_name_for_id_file

input_file: tab-separated file with 3 columns
1. Gene id for node 1
2. Gene id for node 2
3. Edge weight

gene_name_for_id_file: tab-separated file with 2 columns
1. Gene id
2. Gene Name"

OPTIONS="Options:
    -i    Inflation value [default: 1.4]
    -o    Output basename [default: input file name minus .tsv]
    -t    Edge weight threshold [default: 0.85]
    -c    Cluster size threshold
    -d    print debugging info
    -v    verbose output
    -q    turn verbose output off
    -h    print help info"

# default values
debug=0
verbose=1
INFLATION=1.4
THRESHOLD=0.85
CLUSTER_SIZE_THRESHOLD=1

while getopts ":i:o:t:c:dhqv" opt; do
  case $opt in
    i)  INFLATION=$OPTARG  ;;
    o)  OUTPUT_BASE=$OPTARG  ;;
    t)  THRESHOLD=$OPTARG  ;;
    c)  CLUSTER_SIZE_THRESHOLD=$OPTARG  ;;
    d)  debug=1  ;;
    h)  usage "" ;;
    v)  verbose=1 ;;
    q)  verbose=0 ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
  esac
done
shift "$(($OPTIND -1))"

# unpack arguments
INPUT_FILE=$1
GENE_NAME_FILE=$2
if [[ -z $OUTPUT_BASE ]]; then
    OUTPUT_BASE=$( echo $INPUT_FILE | sed -e 's|\.tsv||')
fi
INFLATION_SUFFIX=$( echo $INFLATION | sed -e 's|\.|-|' )

module load MCL

if [[ ! -e $OUTPUT_BASE.mci || ! -e $OUTPUT_BASE.tab ]]; then
    mcxload --stream-mirror -abc $INPUT_FILE -o $OUTPUT_BASE.mci -write-tab $OUTPUT_BASE.tab -tf "abs(),add(-$THRESHOLD)"
    SUCCESS=$?
    error_checking $SUCCESS "mcxload SUCCEEDED." "mcxload FAILED: $SUCCESS"
fi

# CLUSTER
mcl $OUTPUT_BASE.mci -I $INFLATION -o $OUTPUT_BASE.mci.I${INFLATION_SUFFIX}
SUCCESS=$?
error_checking $SUCCESS "mcl SUCCEEDED." "mcl FAILED: $SUCCESS"

# CREATE GENE INFO FILE FOR NODES IN GRAPH
# make unique ID so that jobs run concurrently don't clash
uid=$( echo $RANDOM | md5sum | head -c 6 )
head -n1 $GENE_NAME_FILE > gene-info-$uid.tsv
join -t$'\t' <( sort -t$'\t' -k1,1 $GENE_NAME_FILE ) \
<( awk '{print $2}' $OUTPUT_BASE.tab | sort -u ) >> gene-info-$uid.tsv

# PARSE THE OUTPUT
module load Perl/5.43.0
perl ~/checkouts/bioinf-gen/parse_mcl_output.pl \
--info_file gene-info-$uid.tsv \
--info_file_node_id_col 1 \
--info_file_node_name_col 2 \
--graphml_file $OUTPUT_BASE.mci.I${INFLATION_SUFFIX}-min${CLUSTER_SIZE_THRESHOLD}.graphml \
--cluster_size_threshold $CLUSTER_SIZE_THRESHOLD \
$OUTPUT_BASE.tab $OUTPUT_BASE.mci \
$OUTPUT_BASE.mci.I${INFLATION_SUFFIX} \
$OUTPUT_BASE.mci.I${INFLATION_SUFFIX}-min${CLUSTER_SIZE_THRESHOLD}.tsv

# remove gene info file
rm gene-info-$uid.tsv

SUCCESS=$?
error_checking $SUCCESS "job mcl SUCCEEDED." "job mcl FAILED: $SUCCESS"
exit $SUCCESS

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
