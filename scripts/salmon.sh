#!/usr/bin/env bash
# salmon.sh - Script to run salmon on a set of FASTQ files
# expects a file name fastq.tsv to exist in the working directory
# This should contains three columns
# sample names
# comma-separated list of read1 fastq filenames
# comma-separated list of read2 fastq filenames
# Input to the script is a line number to index into the fastq.tsv file
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=12G

USAGE="salmon.sh [options] TASK_NUMBER"

source bash_functions.sh

OPTIONS="Options:
    -g    genemap, gtf [Danio_rerio.GRCz11.111.gtf]
    -i    index [Danio_rerio.GRCz11.salmon]
    -l    library type [ISR]
    -p    number of threads to use [4]
    -v    verbose output
    -d    print debugging info
    -h    print help info"

# default values
GTF=Danio_rerio.GRCz11.111.gtf
INDEX=Danio_rerio.GRCz11.salmon
LIB_TYPE=ISR
threads=4
verbose=0
debug=0
while getopts ":g:i:l:p:vdh" opt; do
  case $opt in
    g) GTF=$OPTARG ;;
    i) INDEX=$OPTARG ;;
    l) LIB_TYPE=$OPTARG ;;
    p) threads=$OPTARG ;;
    v) verbose=1 ;;
    d) debug=1 ;;
    h)  usage "" ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
  esac
done
shift "$(($OPTIND -1))"

line=`sed "${1}q;d" fastq.tsv`
sample=`echo $line | awk '{ print $1 }'`
fastq1=`echo $line | awk '{ print $2 }' | sed -e 's|,| |g' `
fastq2=`echo $line | awk '{ print $3 }' | sed -e 's|,| |g' `
mkdir -p $sample

module purge
module load salmon

CMD="salmon quant \
-l $LIB_TYPE \
-i $INDEX  \
-1 $fastq1 \
-2 $fastq2 \
-o $sample \
-p $threads \
-g $GTF \
--seqBias --gcBias --posBias"

if [[ $verbose -gt 0 ]]; then
    echo $CMD
fi

if [[ $debug -gt 0 ]]; then
    echo "gtf file = $GTF
library type = $LIB_TYPE
index = $INDEX
threads = $threads
verbose = $verbose" 1>&2
    exit 1
fi

eval $CMD
SUCCESS=$?
verbose=1
error_checking $SUCCESS "salmon command succeeded." "salmon command failed: $SUCCESS"

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
