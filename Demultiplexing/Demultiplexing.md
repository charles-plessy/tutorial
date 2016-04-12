Demultiplexing with TagDust 2
=============================

In this example, we will install TagDust 2, download nanoCAGE data from DDBJ,
and demultiplex it.  These commands have been tested on a virtual machine
running Linux.

## If $HOME/bin is not available, create it.

```
mkdir $HOME/bin
```

Log out, and log in again.  On Linux systems such as Debian, $HOME/bin is now in the $PATH.

## If TagDust 2 is not available, install it.

Here we are downloading version 2.31.

```
wget http://downloads.sourceforge.net/project/tagdust/tagdust-2.31.tar.gz
tar xvfz tagdust-2.31.tar.gz 
cd tagdust-2.31
./configure
make
cp src/tagdust $HOME/bin
cd
```

## Download the sequences from DDBJ.

These are Read 1 and Read 2 from a HiSeq lane where nanoCAGE libraries have
been multiplexed using a barcode that spans the first 6 bases of Read1.
Illumina's TruSeq indexing system has not been used here, which is why we need
to demultiplex by hand.  The reason for using barcodes in Read 1 instead of
TruSeq indexes is to allow for multiplexing the library construction at an
early stage.

{{{
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004180/DRX044600/DRR049557_1.fastq.bz2
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004180/DRX044600/DRR049557_2.fastq.bz2
}}}

## Demultiplex.

### Create an _architecture_ file for TagDust.

```
cat > DRR049557.arch <<__END__
tagdust -1 B:ACATGA,ATCATA,CACGTG,CGATGA,GAGATA,GCTCTC,GTATGA,TCGATA,AGTAGC,ATCGCA,CACTCT,CTGAGC,GAGCGT,GCTGCA,TATAGC,CACGAT,CTGACG -2 F:NNNNNNNN -3 S:TATAGGG -4 R:N
tagdust -1 R:N
__END__
```

### Run TagDust in paired-end mode with the architecture file.

```
tagdust -t 4 -ref tagdust.fa -arch DRR049557.arch -o DRR049557 DRR049557_1.fastq.bz2 DRR049557_2.fastq.bz2
```

`-t 4` allocates 4 cores to the task.  You can adapt this number to the machine where you run that command.

## Rename the samples according to their barcodes.

Here is the list associating barcodes to sample names.  It was constructed
using the information submitted to DDBJ
(<https://trace.ddbj.nig.ac.jp/DRASearch/experiment?acc=DRX044600>).

```
cat > DRR049557.samples.txt <<__END__
ACATGA HaCaT_r1
ATCATA HaCaT_r2
CACGTG HaCaT_r3
CGATGA c33a_r1
GAGATA c33a_r2
GCTCTC c33a_r3
GTATGA Hela_r1
TCGATA Hela_r2
AGTAGC Hela_r3
ATCGCA Caski_r1
CACTCT Caski_r2
CTGAGC Caski_r3
GAGCGT SiHa_r1
GCTGCA SiHa_r2
TATAGC SiHa_r3
CACGAT W12E_r1
CTGACG W12E_r2
__END__
```

Read the list and rename with a loop.

```
cat DRR049557.samples.txt |
  while read from to
  do
    prename -f "s/$from/$to/" *$from*
  done
```
