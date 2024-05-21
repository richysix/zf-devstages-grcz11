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
# qsub -t 1 -l h_vmem=128G -l h_rt=240:0:0 -o cor.o -e cor.e $gitdir/qsub/rscript-array.sh cor.txt cor
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

Run MCL clustering for a variety of 
```
# run one job per threshold first
# otherwise jobs will clash all trying to create the same file
for threshold in 0.8 0.86 0.9
do
  threshold_suffix=$( echo $threshold | perl -lane 'print $F[0] * 100' )
  if [[ ! -e all-cor-long-${threshold_suffix}.tsv ]]; then
    awk -F"\t" '{if(NR == 1){ print $0 }else{ if($3 > '$threshold' ){ print $0 } }}' all-cor-long-70.tsv > all-cor-long-${threshold_suffix}.tsv
  fi
  for INFLATION in 1.4
  do
    qsub -o mcl-clustering-$threshold_suffix-$INFLATION.o -e mcl-clustering-$threshold_suffix-$INFLATION.e \
    $gitdir/qsub/run-mcl-clustering.sh -i $INFLATION \
    -o all-cor-${threshold_suffix} -t $threshold all-cor-long-${threshold_suffix}.tsv gene-names.tsv
  done
done

# run the rest of the jobs with different inflation parameters
for threshold in 0.8 0.86 0.9
do
  threshold_suffix=$( echo $threshold | perl -lane 'print $F[0] * 100' )
  if [[ ! -e all-cor-long-${threshold_suffix}.tsv ]]; then
    awk -F"\t" '{if(NR == 1){ print $0 }else{ if($3 > '$threshold' ){ print $0 } }}' all-cor-long-70.tsv > all-cor-long-${threshold_suffix}.tsv
  fi
  for INFLATION in 1.8 2.0 2.2 4.0
  do
    qsub -o mcl-clustering-$threshold_suffix-$INFLATION.o -e mcl-clustering-$threshold_suffix-$INFLATION.e \
    $gitdir/qsub/run-mcl-clustering.sh -i $INFLATION \
    -o all-cor-${threshold_suffix} -t $threshold all-cor-long-${threshold_suffix}.tsv gene-names.tsv
  done
done

for threshold in 0.8 0.86 0.9
do
  threshold_suffix=$( echo $threshold | perl -lane 'print $F[0] * 100' )
  echo $threshold
  clm info all-cor-$threshold_suffix.mci all-cor-$threshold_suffix.mci.I[0-9]-[0-9]
done
0.8
efficiency=0.43917 massfrac=0.93145 areafrac=0.08476  source=all-cor-80.mci.I1-4 clusters=8467 max=6848  ctr=2486.1 avg=3.5 min=1 DGI=6848 TWI=3 TWL=2106 sgl=8076 qrt=8312
===
efficiency=0.48883 massfrac=0.88921 areafrac=0.06055  source=all-cor-80.mci.I1-8 clusters=8616 max=6141  ctr=1776.4 avg=3.4 min=1 DGI=6141 TWI=7 TWL=671 sgl=8078 qrt=8408
===
efficiency=0.50045 massfrac=0.87969 areafrac=0.05761  source=all-cor-80.mci.I2-0 clusters=8683 max=5981  ctr=1690.2 avg=3.4 min=1 DGI=5981 TWI=7 TWL=668 sgl=8082 qrt=8456
===
efficiency=0.51833 massfrac=0.86358 areafrac=0.05296  source=all-cor-80.mci.I2-2 clusters=8732 max=5856  ctr=1553.8 avg=3.4 min=1 DGI=5856 TWI=9 TWL=486 sgl=8081 qrt=8483
===
efficiency=0.58919 massfrac=0.77489 areafrac=0.03836  source=all-cor-80.mci.I4-0 clusters=9249 max=5096  ctr=1125.6 avg=3.2 min=1 DGI=5096 TWI=30 TWL=68 sgl=8163 qrt=8885
0.86
efficiency=0.50980 massfrac=0.95017 areafrac=0.05115  source=all-cor-86.mci.I1-4 clusters=12178 max=5346  ctr=1500.6 avg=2.4 min=1 DGI=5346 TWI=10 TWL=95 sgl=11792 qrt=12004
===
efficiency=0.56456 massfrac=0.89814 areafrac=0.03156  source=all-cor-86.mci.I1-8 clusters=12371 max=4709  ctr=926.3 avg=2.4 min=1 DGI=4709 TWI=55 TWL=22 sgl=11792 qrt=12115
===
efficiency=0.57908 massfrac=0.88667 areafrac=0.02919  source=all-cor-86.mci.I2-0 clusters=12461 max=4539  ctr=856.9 avg=2.4 min=1 DGI=4539 TWI=76 TWL=18 sgl=11793 qrt=12177
===
efficiency=0.59041 massfrac=0.87518 areafrac=0.02724  source=all-cor-86.mci.I2-2 clusters=12535 max=4306  ctr=799.6 avg=2.3 min=1 DGI=4306 TWI=98 TWL=16 sgl=11796 qrt=12223
===
efficiency=0.66028 massfrac=0.78957 areafrac=0.01650  source=all-cor-86.mci.I4-0 clusters=13194 max=3506  ctr=484.8 avg=2.2 min=1 DGI=3506 TWI=362 TWL=6 sgl=11877 qrt=12714
0.9
efficiency=0.20067 massfrac=0.90552 areafrac=0.16022  source=all-cor-90.mci.I1-4 clusters=417 max=6384  ctr=2876.5 avg=43.0 min=2 DGI=413 TWI=2 TWL=1888 sgl=0 qrt=230
===
efficiency=0.28103 massfrac=0.83654 areafrac=0.10388  source=all-cor-90.mci.I1-8 clusters=634 max=5070  ctr=1865.4 avg=28.3 min=1 DGI=619 TWI=3 TWL=1206 sgl=1 qrt=363
===
efficiency=0.31514 massfrac=0.80634 areafrac=0.08881  source=all-cor-90.mci.I2-0 clusters=724 max=4700  ctr=1594.8 avg=24.8 min=1 DGI=695 TWI=4 TWL=906 sgl=1 qrt=409
===
efficiency=0.34957 massfrac=0.77595 areafrac=0.07253  source=all-cor-90.mci.I2-2 clusters=816 max=4184  ctr=1302.6 avg=22.0 min=1 DGI=759 TWI=6 TWL=558 sgl=4 qrt=477
===
efficiency=0.46231 massfrac=0.62436 areafrac=0.05033  source=all-cor-90.mci.I4-0 clusters=1639 max=3656  ctr=904.3 avg=11.0 min=1 DGI=1283 TWI=20 TWL=84 sgl=115 qrt=1096

# check jobs ran successfully
grep -lE SUCCEEDED mcl-clustering-[89]*o | wc -l
15
```

Also, try pruning edges with k-nearest neighbours
```
mcx query -imx $base-$suffix.mci -vary-knn 100/500/100
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
k-NN     The knn parameter
Total number of nodes: 26007
----------------------------------------------------------------------------------------------
  L   D   R   S     E     cce  EWmean   EWmed   EWiqr  NDmean   NDmed  NDiqr CCF  eff    k-NN 
----------------------------------------------------------------------------------------------
 96   3   1   3 0.076   23898    0.90    0.91    0.10   182.1   143.5  307.5   -    -      500
 95   4   1   3 0.059   23504    0.90    0.91    0.09   141.1   107.5  237.5   -    -      400
 94   5   1   5 0.043   22819    0.91    0.92    0.09   101.9    74.5  170.5   -    -      300
 92   7   1   6 0.027   21777    0.91    0.93    0.08    64.7    44.5  106.5   -    -      200
 86  12   1  11 0.013   19369    0.92    0.94    0.07    30.2    20.5   48.5   -    -      100
----------------------------------------------------------------------------------------------
```

Try k-NN 400. It keeps ~6% edges and has median node degree ~100

```
qlogin -l h_vmem=8G

cd /data/scratch/bty114/zf-stages-grcz11/109/deseq2-all
module load MCL

threshold=0.7
suffix=$( echo $threshold | perl -lane 'print $F[0] * 100' )
base=all-cor

mcx query -imx $base-$suffix.mci -vary-knn 100/500/100
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
k-NN     The knn parameter
Total number of nodes: 26007
----------------------------------------------------------------------------------------------
  L   D   R   S     E     cce  EWmean   EWmed   EWiqr  NDmean   NDmed  NDiqr CCF  eff    k-NN 
----------------------------------------------------------------------------------------------
 96   3   1   3 0.076   23898    0.90    0.91    0.10   182.1   143.5  307.5   -    -      500
 95   4   1   3 0.059   23504    0.90    0.91    0.09   141.1   107.5  237.5   -    -      400
 94   5   1   5 0.043   22819    0.91    0.92    0.09   101.9    74.5  170.5   -    -      300
 92   7   1   6 0.027   21777    0.91    0.93    0.08    64.7    44.5  106.5   -    -      200
 86  12   1  11 0.013   19369    0.92    0.94    0.07    30.2    20.5   48.5   -    -      100
---------------------------------------------------------------------------------------------

knn_threshold=400
mcx alter -imx $base-$suffix.mci -tf "abs(),#knn($knn_threshold)" -o $base-$suffix-knn${knn_threshold}.mci

for INFLATION in 1.4 1.6 2.0 4.0
do
INFLATION_SUFFIX=$( echo $INFLATION | perl -lane 'print $F[0] * 10' )
mcl $base-$suffix-knn${knn_threshold}.mci -I ${INFLATION} -o $base-$suffix-knn${knn_threshold}.mci.I${INFLATION_SUFFIX}
done

clm info all-cor-70-knn400.mci all-cor-70-knn400.mci.I[0-9][0-9]
efficiency=0.15629 massfrac=0.83000 areafrac=0.07443  source=all-cor-70-knn400.mci.I14 clusters=1384 max=4108  ctr=1936.5 avg=18.8 min=1 DGI=1330 TWI=4 TWL=2371 sgl=907 qrt=1127
===
efficiency=0.19235 massfrac=0.77197 areafrac=0.04415  source=all-cor-70-knn400.mci.I16 clusters=1619 max=2576  ctr=1149.2 avg=16.1 min=1 DGI=1385 TWI=7 TWL=1152 sgl=909 qrt=1254
===
efficiency=0.27028 massfrac=0.66856 areafrac=0.01725  source=all-cor-70-knn400.mci.I20 clusters=2261 max=1287  ctr=449.6 avg=11.5 min=1 DGI=1157 TWI=18 TWL=516 sgl=1026 qrt=1797
===
efficiency=0.26274 massfrac=0.34114 areafrac=0.00248  source=all-cor-70-knn400.mci.I40 clusters=7212 max=411  ctr=65.5 avg=3.6 min=1 DGI=411 TWI=255 TWL=14 sgl=4384 qrt=6518

```

Forgot to subtract the correlation threshold. Redo
```
qlogin -l h_vmem=8G

cd /data/scratch/bty114/zf-stages-grcz11/109/deseq2-all/
module load MCL

knn_threshold=400
threshold=0.7
threshold_suffix=$( echo $threshold | perl -lane 'print $F[0] * 100' )

OUTPUT_BASE=all-cor-${threshold_suffix}-knn${knn_threshold}
mcxload --stream-mirror -abc all-cor-long-${threshold_suffix}.tsv \
-o $OUTPUT_BASE.mci -write-tab $OUTPUT_BASE.tab -tf "abs(),#knn($knn_threshold),add(-$threshold)"
# exit
^D

# run clustering 
for INFLATION in 1.4 1.6 1.8 2.0 2.2 4.0
do
  qsub -o mcl-clustering-$threshold_suffix-knn${knn_threshold}-$INFLATION.o \
  -e mcl-clustering-$threshold_suffix-knn${knn_threshold}-$INFLATION.e \
  $gitdir/qsub/run-mcl-clustering.sh -i $INFLATION \
  -o all-cor-${threshold_suffix}-knn${knn_threshold} -t $threshold all-cor-long-${threshold_suffix}.tsv gene-names.tsv
done

# number of nodes in clusters of size > 99
for INFLATION in 1.4 1.6 1.8 2.0 2.2 4.0
do
  INFLATION_SUFFIX=$( echo $INFLATION | sed -e 's|\.|-|' )
  cut -f3 all-cor-70-knn400.mci.I${INFLATION_SUFFIX}.tsv | uniq -c | sort -grk1,2 | \
  awk 'BEGIN{sum = 0} {if($1 > 99){ sum = sum + $1 }} END{ print sum }' 
done
20141
19000
18117
16737
15196
3939

# number of nodes in clusters of size > 4
for INFLATION in 1.4 1.6 1.8 2.0 2.2 4.0
do
  INFLATION_SUFFIX=$( echo $INFLATION | sed -e 's|\.|-|' )
  cut -f3 all-cor-70-knn400.mci.I${INFLATION_SUFFIX}.tsv | uniq -c | sort -grk1,2 | \
  awk 'BEGIN{sum = 0} {if($1 > 4){ sum = sum + $1 }} END{ print sum }' 
done
22910
22507
22176
21442
20509
12296
```

```
clm info all-cor-70-knn400.mci all-cor-70-knn400.mci.I[0-9]-[0-9]
efficiency=0.29987 massfrac=0.83982 areafrac=0.05446  source=all-cor-70-knn400.mci.I1-4 clusters=6139 max=3596  ctr=1597.8 avg=4.8 min=1 DGI=3596 TWI=5 TWL=1433 sgl=5547 qrt=5882
===
efficiency=0.32930 massfrac=0.79819 areafrac=0.03286  source=all-cor-70-knn400.mci.I1-6 clusters=6354 max=2264  ctr=964.5 avg=4.6 min=1 DGI=2264 TWI=9 TWL=1004 sgl=5548 qrt=6029
===
efficiency=0.36753 massfrac=0.74634 areafrac=0.01759  source=all-cor-70-knn400.mci.I1-8 clusters=6547 max=1428  ctr=516.6 avg=4.5 min=1 DGI=1428 TWI=16 TWL=467 sgl=5558 qrt=6155
===
efficiency=0.39474 massfrac=0.70847 areafrac=0.01267  source=all-cor-70-knn400.mci.I2-0 clusters=6911 max=1049  ctr=372.6 avg=4.2 min=1 DGI=1049 TWI=23 TWL=244 sgl=5639 qrt=6490
===
efficiency=0.41563 massfrac=0.67331 areafrac=0.00968  source=all-cor-70-knn400.mci.I2-2 clusters=7480 max=886  ctr=284.7 avg=3.9 min=1 DGI=886 TWI=33 TWL=151 sgl=5939 qrt=7050
===
efficiency=0.42246 massfrac=0.43981 areafrac=0.00147  source=all-cor-70-knn400.mci.I4-0 clusters=14156 max=452  ctr=44.0 avg=2.1 min=1 DGI=452 TWI=1317 TWL=3 sgl=11130 qrt=13544
```

Rerun parse_mcl_output.pl on 2.2 clustering and remove clusters smaller then 5
```
qsub ../qsub/parse_mcl_clustering-2-2-min-size-5.sh
```