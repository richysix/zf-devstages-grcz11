# Zebrafish Developmental Stages Data 

Remap the Zebrafish Developmental Stages date published in [White et al 2017](https://doi.org/10.7554/eLife.30860)
to GRCz11.

## Setup

Environment variables
```
gitdir=$HOME/checkouts/zf-devstages-grcz11
basedir=/data/scratch/$USER/zf-stages-grcz11/109/
```

Download bash_functions.sh and place in path
```
cd $HOME/bin # This needs to be one of the directories in the $PATH variable
wget https://raw.githubusercontent.com/richysix/bioinf-gen/master/bash_functions.sh
```

Create directories

```
mkdir -p $basedir $gitdir/scripts $gitdir/qsub
# set up symlinks
cd $basedir
ln -s $gitdir/qsub qsub
ln -s $gitdir/scripts scripts
```

## Index genome with STAR

Download job script from uge repository
```
cd $gitdir/qsub
wget https://raw.githubusercontent.com/richysix/uge-job-scripts/8d9b85ebc9e1ded8c1ba315c5d4b627f0e4a2b9c/star_ref_index.sh
```

Index the e109 genome and annotation
```
cd $basedir
mkdir -p reference/grcz11
cd reference/
qsub $gitdir/qsub/star_ref_index.sh -v -n grcz11 \
ftp://ftp.ensembl.org/pub/release-109/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz \
ftp://ftp.ensembl.org/pub/release-109/gtf/danio_rerio/Danio_rerio.GRCz11.109.gtf.gz
```

## Download FASTQ data

Download project data from ENA:
```
cd $basedir
project=PRJEB12982
fields="secondary_sample_accession,run_accession,fastq_md5,fastq_ftp,sample_alias,sample_description"
curl -SsL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$project&result=read_run&fields=$fields&format=tsv" \
> ena.tsv
```

Download sample info from Figshare
```
wget https://figshare.com/ndownloader/files/21662469
mv 21662469 zfs-rnaseq-samples.tsv
```

Join sample info to ENA info
```
# header
head -1 ena.tsv | sed -e 's|secondary_sample_accession|accession_number|' | \
join -t$'\t' -1 3 -2 2 <( head -1 zfs-rnaseq-samples.tsv ) - > ena-sample-names.tsv
# samples
join -t$'\t' -1 3 -2 2 <( sort -t$'\t' -k3,3 zfs-rnaseq-samples.tsv ) <( sort -t$'\t' -k2,2 ena.tsv ) | \
sort -t$'\t' -k4,4 >> ena-sample-names.tsv
```

Create download files
```
dir=fastq
mkdir $dir
cut -f11 ena-sample-names.tsv | grep -v fastq_ftp | sed -e 's|;|\n|' | \
sed -e 's|^\(.*/\)\(.*\)$|\1\2\t\2|' > $dir/curl
split -l72 -d $dir/curl $dir/curl.
rename $dir/curl.0 $dir/curl. $dir/curl.*
num_tasks=$( find $dir -type f -name "curl.*" | wc -l )
mv $dir/curl.0 $dir/curl.${num_tasks}
rm $dir/curl
cd $dir
```

Download job script and run as array job
(`curl-file-download.sh` is in the `scripts` directory
```
cd $gitdir/qsub
https://raw.githubusercontent.com/richysix/uge-job-scripts/8d9b85ebc9e1ded8c1ba315c5d4b627f0e4a2b9c/curl-file-download-array.sh
qsub -t 1-${num_tasks} $gitdir/qsub/curl-file-download-array.sh
```

Check md5sums
Create a file of checksums and filenames
```
paste <(cut -f10 ../ena-sample-names.tsv | sed -e 's/;/\n/') \
<(cut -f11 ../ena-sample-names.tsv | sed -e 's/;/\n/') | \
grep -v fastq_ftp | sed -e 's|\t.*/|  |' > md5sum.txt
```

Download check-md5sum and run job
```
cd $gitdir/qsub
wget https://raw.githubusercontent.com/richysix/uge-job-scripts/8d9b85ebc9e1ded8c1ba315c5d4b627f0e4a2b9c/check-md5sums.sh
cd $basedir/$dir
qsub ../qsub/check-md5sums.sh
```

Check there are no lines in md5sum.out
Lines are only output if the checksum doesn't match
```
wc -l md5sum.out 
0 md5sum.out
```
