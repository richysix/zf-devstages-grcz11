#!/usr/bin/env bash
# curl-file-download.sh - Script to download files using curl

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

USAGE="curl-file-download.sh [options] input_file

The input file should be a tab-separated file of urls and optionally file names"

OPTIONS="Options:
    -d    print debugging info
    -v    verbose output
    -h    print help info"

# default values
debug=0
verbose=0

while getopts ":dhvq" opt; do
  case $opt in
    d) debug=1 ;;
    h)  echo $USAGE; echo $OPTIONS ;;
    v)  verbose=1 ;;
    q)  verbose=0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 2 ;;
    :)  echo "Option -$OPTARG requires an argument!" >&2; exit 2 ;;
  esac
done
shift "$(($OPTIND -1))"

# unpack arguments
INPUT_FILE=$1

while read line
do
    url=$( echo $line | awk '{ print $1 }' )
    filename=$( echo $line | awk '{ print $2 }' )
    if [[ -z $filename ]]; then
        filename=$( echo $url | sed -e 's|^.*/||')
    fi
    if [[ -e $filename ]]; then continue; fi
    CMD="curl -sS -o $filename.tmp $url && mv $filename.tmp $filename"
    if [[ $debug -eq 1 ]]; then
        echo URL:$url FILE:$filename
        echo "$CMD"
    fi
    eval $CMD
    SUCCESS=$?
    if [[ $SUCCESS -eq 0 ]]; then
      echo "Download of file, $filename, succeeded."
    else
      echo "Download of file, $filename, failed: $?"
      exit $SUCCESS
    fi
    sleep 2
done < $INPUT_FILE
