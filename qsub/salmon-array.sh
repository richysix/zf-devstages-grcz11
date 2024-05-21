#!/usr/bin/env bash
# salmon-array.sh - Script to run salmon as an array job
# expects a file name fastq.tsv to exist in the working directory
# This should contains three columns
# sample names
# comma-separated list of read1 fastq filenames
# comma-separated list of read2 fastq filenames
# Each line fo the fastq.tsv file will be run as a separate job
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=12G
#$ -t 1-96

USAGE="salmon-array.sh [options] "

source bash_functions.sh

OPTIONS="Options:
    -g    genemap, gtf [Danio_rerio.GRCz11.111.gtf]
    -i    index [Danio_rerio.GRCz11.salmon]
    -l    library type [ISR]
    -p    number of threads to use [4]
    -v    verbose output
    -h    print help info"

# default values
GTF=Danio_rerio.GRCz11.111.gtf
INDEX=Danio_rerio.GRCz11.salmon
LIB_TYPE=ISR
threads=4
verbose=0
while getopts ":g:i:l:p:vh" opt; do
  case $opt in
    g) GTF=$OPTARG ;;
    i) INDEX=$OPTARG ;;
    l) LIB_TYPE=$OPTARG ;;
    p) threads=$OPTARG ;;
    v) verbose=1 ;;
    h)  usage "" ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
  esac
done
shift "$(($OPTIND -1))"

OPT=""
if [[ $verbose -gt 0 ]]; then
    OPT="-v"
fi

CMD="./salmon.sh -g $GTF \
-i $INDEX \
-l $LIB_TYPE \
-p $threads \
$OPT ${SGE_TASK_ID}"

if [[ $verbose -gt 0 ]]; then
    echo $CMD
fi

eval $CMD
SUCCESS=$?
verbose=1
error_checking $SUCCESS "job salmon, task ${SGE_TASK_ID} succeeded." "job salmon, task ${SGE_TASK_ID} failed: $SUCCESS"

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2024. Queen Mary University of London.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007
