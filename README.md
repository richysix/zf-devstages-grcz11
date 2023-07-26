# Zebrafish Developmental Stages Data 

Remap the Zebrafish Developmental Stages date published in [White et al 2017](https://doi.org/10.7554/eLife.30860)
to GRCz11.

## Setup

Environment variables
```
gitdir=$HOME/checkouts/zf-devstages-grcz11
basedir=/data/scratch/$USER/zf-stages-grcz11/109
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
`curl-file-download.sh` is in the `scripts` directory
```
cd $gitdir/qsub
https://raw.githubusercontent.com/richysix/uge-job-scripts/8d9b85ebc9e1ded8c1ba315c5d4b627f0e4a2b9c/curl-file-download-array.sh
qsub -t 1-${num_tasks} $gitdir/qsub/curl-file-download-array.sh
```

### Check md5sums  
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

## Run STAR

Create file of samples to fastq files
```
cd ..
dir=star1
mkdir -p $dir
grep -v 'sample_accession' ena-sample-names.tsv | perl -ane 'push @{$h{$F[6]}}, $F[8]; END { foreach my $k (keys %h)
{ printf "%s\t%s\t%s\n", $k, (join ",", map { "../fastq/${_}_1.fastq.gz" } sort @{$h{$k}}), (join ",", map { "../fastq/${_}_2.fastq.gz" } sort @{$h{$k}}) } }' | \
sort -V > $dir/fastq.tsv
```

Download star1 job files
```
cd $gitdir/qsub
wget https://raw.githubusercontent.com/richysix/uge-job-scripts/8d9b85ebc9e1ded8c1ba315c5d4b627f0e4a2b9c/star1-array.sh
cd ../scripts/
wget https://raw.githubusercontent.com/richysix/uge-job-scripts/8d9b85ebc9e1ded8c1ba315c5d4b627f0e4a2b9c/star1.sh
```

Run all samples as array job
```
cd $basedir/$dir
# symlink star1 script
ln -s $gitdir/scripts/star1.sh star1.sh
# run job
num_tasks=$( wc -l fastq.tsv | awk '{ print $1 }' )
qsub -t 1-${num_tasks} $gitdir/qsub/star1-array.sh
```

Check whether they succeeded or failed
```
grep -i -E 'succeeded|failed' star1-array.sh.o* | awk '{ print $5 }' | sort | uniq -c
     90 succeeded.
cd ..
```

### Remap using discovered splice junctions

Download star2 job files
```
cd $gitdir/qsub
wget https://raw.githubusercontent.com/richysix/uge-job-scripts/8d9b85ebc9e1ded8c1ba315c5d4b627f0e4a2b9c/star2-array.sh
cd ../scripts/
wget https://raw.githubusercontent.com/richysix/uge-job-scripts/8d9b85ebc9e1ded8c1ba315c5d4b627f0e4a2b9c/star2.sh
```

Copy file of fastq files to star2 directory
```
cd $basedir
dir=star2
mkdir -p $dir
cp star1/fastq.tsv $dir
```

Run STAR as array job
```
cd $dir
# symlink star2 script
ln -s $gitdir/scripts/star2.sh star2.sh

qsub -t 1-${num_tasks} $gitdir/qsub/star2-array.sh
```

Check jobs completed
```
grep -i -E 'succeeded|failed' star2-array.sh.o* | awk '{ print $5 }' | sort | uniq -c
     90 SUCCEEDED.
```
Check error files
```
cat star2-array.sh.e* | wc -l
0
```

## DESeq2

Get Ensembl e109 annotation
```
# download scripts
cd $gitdir/qsub
wget https://raw.githubusercontent.com/richysix/uge-job-scripts/50a25a47cce6654f860ff350dfe6466e9f1c2da0/get_ensembl_gene_annotation.sh 
cd $gitdir/scripts
wget https://raw.githubusercontent.com/iansealy/bio-misc/4e27d60323907d55d37a1ec8f468ca771f542f78/get_ensembl_gene_annotation.pl

mkdir annotation
qsub  qsub/get_ensembl_gene_annotation.sh -e 109 -s 'Danio rerio' \
-o annotation/annotation.txt -f scripts/get_ensembl_gene_annotation.pl
```

### Run DESeq2 to get normalised counts

Run DESeq2 with samples divided into Maternal and Zygotic. Split is irrelevant,
but the DESeq2 script requires two conditions to run.

Create samples file
```
dir=deseq2-all
mkdir $dir
awk -F"\t" '{ print $7 "\t" $6 }' zfs-rnaseq-samples.tsv | grep -v sampleName > $dir/samples.tsv
```

Create file of commands to run with Rscript
```
con=$( cut -f4,6 zfs-rnaseq-samples.tsv | grep -v condition | uniq | sed -e 's|ZFS:0*||' | \
awk '{if($1 <= 11){ print $2 }}' | tr '\n' ',' | sed -e 's|,$||' )
exp=$( cut -f4,6 zfs-rnaseq-samples.tsv | grep -v condition | uniq | sed -e 's|ZFS:0*||' | \
awk '{if($1 > 11){ print $2 }}' | tr '\n' ',' | sed -e 's|,$||' )
echo "Rscript $gitdir/scripts/deseq2.R -s $dir/samples.tsv -e $exp -c $con \
-a $basedir/annotation/annotation.txt -d $basedir/star2 -o $dir" > deseq2.txt
```

Run rscript job script
```
cd $gitdir/qsub
wget https://raw.githubusercontent.com/richysix/uge-job-scripts/5161dae7f8688aca015ea809af0ce231790e9e60/rscript-array.sh
cd $basedir
qsub -t 1 -o $dir/deseq-all.o -e $dir/deseq-all.e $gitdir/qsub/rscript-array.sh deseq2.txt deseq2
```
