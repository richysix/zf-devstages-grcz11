#!/usr/bin/env bash
# salmon_index.sh - Script to create an index file for salmon
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=1:0:0
#$ -l h_vmem=4G

USAGE="salmon_index.sh [options] fasta_file index_name"

source bash_functions.sh

OPTIONS="Options:
    -f    decoys file [NULL]
    -k    kmer number [31]
    -p    number of threads to use [4]
    -v    verbose output
    -d    print debugging info
    -h    print help info"

# default values
kmer=31
threads=4
verbose=0
debug=0
while getopts ":f:k:p:vdh" opt; do
  case $opt in
    f) decoy_file=$OPTARG ;;
    k) kmer=$OPTARG ;;
    p) threads=$OPTARG ;;
    v) verbose=1 ;;
    d) debug=1 ;;
    h)  usage "" ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
  esac
done
shift "$(($OPTIND -1))"

if [[ -z $decoy_file ]]; then
  DECOYS=""
else
  DECOYS="-d $decoy_file"
fi

module purge
module load salmon

CMD="salmon index \
--transcripts $1 
$DECOYS \
-k $kmer \
--index $2 \
--threads $threads"

if [[ $verbose -gt 0 ]]; then
    echo $CMD
fi

if [[ $debug -gt 0 ]]; then
    echo "decoy file = $decoy_file
kmer = $kmer
threads = $threads
verbose = $verbose" 1>&2
    exit 1
fi

eval $CMD
error_checking $? "salmon index command succeeded." "salmon index command failed: $?"
