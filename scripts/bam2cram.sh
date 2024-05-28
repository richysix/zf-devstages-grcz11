#!/usr/bin/env bash
# bam2cram.sh - Script to convert bam files to cram files
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G

source bash_functions.sh

USAGE="bam2cram.sh [options] bam_file ref_file cram"

OPTIONS="Options:
    -d    print debugging info
    -q    turn verbose output off
    -v    verbose output
    -h    print help info"

# default values
debug=0
verbose=1

while getopts ":dhqv" opt; do
  case $opt in
    d)  debug=1  ;;
    h)  usage "" $USAGE $OPTIONS ;;
    q)  verbose=0 ;;
    v)  verbose=1 ;;
    \?) usage "Invalid option: -$OPTARG" $USAGE $OPTIONS ;;
    :)  usage "Option -$OPTARG requires an argument!" $USAGE $OPTIONS ;;
  esac
done
shift "$(($OPTIND -1))"

# unpack arguments
BAM=$1
REF=$2
CRAM=$3

if [[ -z $REF ]];then
    REF=$SHARED/genomes/GRCz11/GRCz11.fa
fi

if [[ -z $CRAM ]];then
    CRAM=$( echo $BAM | sed -e 's|\.bam$|.cram|' )
fi

if [[ $debug -eq 1 ]]; then
    echo "BAM: $BAM
    REF: $REF
    CRAM: $CRAM" >&2
    exit 1
fi

module load samtools/1.17

samtools view -h -C --reference $REF -o $CRAM --write-index $BAM
SUCCESS=$?

error_checking $SUCCESS "job bam2cram SUCCEEDED." "job bam2cram FAILED: $SUCCESS"

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
