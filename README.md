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

### Archive mapped bams

Make file of bam files names with ZFS stage cram file names
```
cut -f4,7 zfs-rnaseq-samples.tsv | grep -v condition | \
sed -e 's|-\([1-5]\)$|-\1\t\1|; s|ZFS:|zfs-|; s|00000||;' | \
awk -F"\t" '{ print "star2/" $2 "/Aligned.sortedByCoord.out.bam\t" "'$SHARED'" "/genomes/GRCz11/GRCz11.fa\tstar2/" $1 "-" $3 ".cram" }' > bam2cram.txt
```

Run bam2cram array job
```
# symlink to bam2cram script
ln -s ~/checkouts/uge-job-scripts/bam2cram.sh bam2cram.sh
qsub -t 1-90 ~/checkouts/uge-job-scripts/bam2cram-array.sh bam2cram.txt

grep -ihE 'SUCCEEDED|FAILED' bam2cram-array.sh.[oe]30790* | awk '{print $3}' | sort | uniq -c
     90 SUCCEEDED.
```

Copy to archive
```
echo '#!/usr/bin/env bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -o rclone-copy-star2.o
#$ -e rclone-copy-star2.e

module load rclone
rclone copy --filter "+ *cram*" --filter "- *" star2/ sharepoint-qmul:Projects/zf-stages-grcz11/e109/star2/
rclone copy /data/SBBS-BuschLab/genomes/GRCz11/GRCz11.fa sharepoint-qmul:Projects/zf-stages-grcz11/e109/' > qsub/rclone-copy-star2.sh

qsub qsub/rclone-copy-star2.sh
```

Check files transferred properly
```
rclone check --filter "+ *cram*" --filter "- *" star2/ sharepoint-qmul:Projects/zf-stages-grcz11/e109/star2/
2023/08/30 12:18:18 NOTICE: OneDrive root 'Projects/zf-stages-grcz11/e109/star2': 0 differences found
2023/08/30 12:18:18 NOTICE: OneDrive root 'Projects/zf-stages-grcz11/e109/star2': 180 matching files

# check file sizes match
rclone ls sharepoint-qmul:Documents/Projects/zf-stages-grcz11/e109/star2/ > transferred-files.txt
ls -l star2/*.cram* | sed -e 's|star2/||' | sort -k9,9 | join -1 9 -2 2 - <( sort -k2,2 transferred-files.txt ) | awk '{ if($6 != $10 ){ print $0 } }'
```

## Make correlation network and cluster with MCL

```
cd scripts
wget https://raw.githubusercontent.com/richysix/bioinf-gen/f13590412d28c443c700f270d90db8923f13f552/create_cor_network.R

cd $basedir/$dir
echo "Rscript $gitdir/scripts/create_cor_network.R samples.tsv all.tsv \
--threshold 0.7 --expansion 2 --inflation 1.4 \
--output all-cor-long-70.tsv --clusters_file all-cor-long-70-clusters.tsv" > cor.txt
qsub -t 1 -l h_vmem=128G -o cor.o -e cor.e $gitdir/qsub/rscript-array.sh cor.txt cor
```

Decide on best value of correlation coeff cut-off to use
```
qlogin -l h_vmem=8G

cd /data/scratch/bty114/zf-stages-grcz11/109/deseq2-all/
module load MCL

threshold=0.7
suffix=$( echo $threshold | perl -lane 'print $F[0] * 100' )
base=all-cor
mcxload --stream-mirror -abc all-cor-long-70.tsv -o $base-$suffix.mci -write-tab $base-$suffix.tab -tf 'abs()'

mcx query -imx $base-$suffix.mci --vary-correlation -vary-threshold 0.7/1/15
[mclIO] reading <all-cor-70.mci>
.......................................
[mclIO] read native interchange 26007x26007 matrix with 62072692 entries
-------------------------------------------------------------------------------
 L       Percentage of nodes in the largest component
 D       Percentage of nodes in components of size at most 3 [-div option]
 R       Percentage of nodes not in L or D: 100 - L -D
 S       Percentage of nodes that are singletons
 E       Fraction of edges retained (input graph has 62072692)
 cce     Expected size of component, nodewise [ sum(sz^2) / sum^2(sz) ]
 EW      Edge weight traits (mean, median and IQR)
 ND      Node degree traits [mean, median and IQR]
 CCF     Clustering coefficient (scale 1-100)
 eff     Induced component efficiency relative to start graph (scale 1-1000)
Total number of nodes: 26007
----------------------------------------------------------------------------------------------
  L   D   R   S     E     cce  EWmean   EWmed   EWiqr  NDmean   NDmed  NDiqr CCF  eff  Cutoff 
----------------------------------------------------------------------------------------------
 99   0   1   0 1.000   25434    0.81    0.81    0.12  2386.8  1569.5 3775.5   -    -     0.70
 96   2   1   2 0.900   24044    0.82    0.82    0.11  2147.8  1335.5 3314.0   -    -     0.72
 93   5   2   4 0.802   22545    0.84    0.83    0.10  1913.6  1100.5 2838.0   -    -     0.74
 90   7   3   6 0.706   21059    0.85    0.84    0.09  1684.7   884.5 2369.5   -    -     0.76
 87   9   3   9 0.613   19910    0.86    0.85    0.08  1462.3   683.5 2037.5   -    -     0.78
 85  12   3  12 0.523   18657    0.87    0.87    0.07  1248.1   503.5 1779.5   -    -     0.80
 81  15   3  14 0.437   17273    0.88    0.88    0.07  1043.8   345.5 1499.5   -    -     0.82
 78  19   3  18 0.356   15968    0.89    0.89    0.06   850.2   219.5 1187.5   -    -     0.84
 75  22   3  21 0.280   14568    0.91    0.90    0.05   668.5   124.5  847.5   -    -     0.86
 70  27   3  27 0.210   12684    0.92    0.92    0.04   501.2    60.5  555.5   -    -     0.88
 64  33   3  33 0.147   10719    0.93    0.93    0.03   351.4    23.5  322.8   -    -     0.90
 58  40   2  40 0.093    8602    0.94    0.94    0.03   222.5     6.5  147.5   -    -     0.92
 48  50   2  50 0.050    5968    0.96    0.96    0.02   118.3     1.5   46.5   -    -     0.94
 34  64   2  63 0.019    3046    0.97    0.97    0.01    45.0     0.5    7.5   -    -     0.96
 10  83   7  82 0.003     284    0.99    0.98    0.01     7.0     0.5    0.5   -    -     0.98
----------------------------------------------------------------------------------------------
```

