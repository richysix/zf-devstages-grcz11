# Zebrafish Developmental Stages Data 

Remap the Zebrafish Developmental Stages date published in [White et al 2017](https://doi.org/10.7554/eLife.30860)
to GRCz11.

## Setup

```
gitdir=$HOME/checkouts/zf-devstages-grcz11
basedir=/data/scratch/$USER/zf-stages-grcz11/109/
```

### Create directories

```
mkdir -p $basedir $gitdir/scripts $gitdir/qsub
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
