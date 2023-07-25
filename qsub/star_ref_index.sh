#!/usr/bin/env bash
# star-ref-index.sh - Script to create an index file for STAR
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=240:0:0
#$ -l h_vmem=8G

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

USAGE="star_ref_index.sh [options] REF_URL GTF_URL"

OPTIONS="Options:
    -n    Name (grcz11)
    -d    print debugging info
    -v    quite output
    -v    verbose output [default]
    -h    print help info"

# default values
debug=0
verbose=1
NAME=grcz11

while getopts ":n:dhqv" opt; do
  case $opt in
    n) NAME=$OPTARG ;;
    d) debug=1 ;;
    h) echo ""; echo "$USAGE"; echo "$OPTIONS"; exit 1 ;;
    q) verbose=0 ;;
    v) verbose=1 ;;
    \?) echo ""; echo "Invalid option: -$OPTARG" >&2; echo "$USAGE" >&2; echo "$OPTIONS" >&2; exit 1 ;;
    :) echo ""; echo "Option -$OPTARG requires an argument!" >&2; echo "$USAGE" >&2; echo "$OPTIONS" >&2; exit 1 ;;
  esac
done
shift "$(($OPTIND -1))"

# unpack arguments
REF_URL=$1
GTF_URL=$2

REF_GZ_FILE=$( echo $REF_URL | sed -e 's|^.*/||' )
REF_FILE=$( echo $REF_GZ_FILE | sed -e 's|\.gz$||' )

GTF_GZ_FILE=$( echo $GTF_URL | sed -e 's|^.*/||' )
GTF_FILE=$( echo $GTF_GZ_FILE | sed -e 's|\.gz$||' )

if [[ $debug -eq 1 ]]; then
    echo "REF URL: $REF_URL
    ZIPPED FILE: $REF_GZ_FILE
    FILE: $REF_FILE
    GTF URL: $GTF_URL
    ZIPPED FILE: $GTF_GZ_FILE
    FILE: $GTF_FILE
    INDEX NAME: $NAME"
    exit 0
fi

# define function for errors
# arg1 exit code
# arg2 Message for success
# arg3 Message for failure
function error_checking (){
    if [ $1 -eq 0 ]
    then
        if [[ -n $2 && $verbose -eq 1 ]]; then echo $2 ; fi
    else
        date >& 2
        echo $3 >&2
        exit $1
    fi
}

# check if reference file exists. if not download it
if [[ ! -e $REF_FILE ]]; then
    # download
    wget $REF_URL
    error_checking $? "Download of reference file, $REF_URL, succeeded." "Download of reference file, $REF_URL, failed: $?"
    # unzip
    gunzip $REF_FILE
    error_checking $? "Unzipping reference file succeeded." "Unzipping reference file failed: $?"
else
    if [[ $verbose -eq 1 ]]; then echo "Reference file already exists. No need to download" ; fi
fi

if [[ ! -e $GTF_FILE ]]; then
    wget $GTF_URL
    error_checking $? "Download of transcriptome gtf file, $GTF_URL, succeeded." "Download of transcriptome gtf file, $GTF_URL, failed: $?"
    gunzip $GTF_FILE
    error_checking $? "Unzipping transcriptome gtf file succeeded." "Unzipping transcriptome gtf file failed: $?"
else
    if [[ $verbose -eq 1 ]]; then echo "Transcriptome gtf file already exists. No need to download" ; fi
fi

mkdir $NAME
error_checking $? "Created directory, $NAME." "Couldn't create directory, $NAME: $?"

module load STAR
CMD="STAR \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir $NAME \
--genomeFastaFiles $REF_FILE \
--sjdbGTFfile $GTF_FILE \
--sjdbOverhang 74"

eval $CMD
error_checking $? "STAR index command succeeded." "STAR index command failed: $?"
