Demultiplexing with TagDust 2
=============================

In this example, we will install [TagDust 2](https://sourceforge.net/projects/tagdust/),
download [nanoCAGE](http://population-transcriptomics.org/nanoCAGE/) data from DDBJ,
and demultiplex it.  These commands have been tested on a virtual machine
running Linux.

## If $HOME/bin is not available, create it.

```
mkdir -p $HOME/bin
```

Log out, and log in again.  On Linux systems such as Debian, `$HOME/bin` is now in the `$PATH`.

## If TagDust 2 is not available, install it.

Here we are downloading version 2.33.  This is the first version where the option `-show_finger_seq`
is available.

```
wget http://downloads.sourceforge.net/project/tagdust/tagdust-2.33.tar.gz
tar xvfz tagdust-2.33.tar.gz 
cd tagdust-2.33
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

```
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004180/DRX044600/DRR049557_1.fastq.bz2
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004180/DRX044600/DRR049557_2.fastq.bz2
```

## Demultiplex.

The sequence of Read 1 should match with `BBBBBBNNNNNNNNTATAGGG`, where `B`
is a barcode base, `N` is a UMI base, and `TATAGGG` is a linker sequence that
is the same in every read.  The rest of the bases are cDNA sequences to be aligned
on the genome.

### Create an _architecture_ file for TagDust.

The _architecture_ files declare the structure of the reads in a machine-readable
format, where `B`, `F`, `S` and `R` bloks declare _barcodes_, _fingerprints_ (UMIs),
_spacers_ (linkers), and the _rna_ sequence to be aligned.  More details can be
found in the user manual (the `doc/User-Manual.pdf` file in TagDust's source
archive) or in [Lassmann T., 2015](10.1186/s12859-015-0454-y).

```
cat > DRR049557.arch <<__END__
tagdust -1 B:ACATGA,ATCATA,CACGTG,CGATGA,GAGATA,GCTCTC,GTATGA,TCGATA,AGTAGC,ATCGCA,CACTCT,CTGAGC,GAGCGT,GCTGCA,TATAGC,CACGAT,CTGACG -2 F:NNNNNNNN -3 S:TATAGGG -4 R:N
tagdust -1 R:N
__END__
```

### Run TagDust in paired-end mode with the architecture file.

```
tagdust -t 4 -show_finger_seq -ref tagdust.fa -arch DRR049557.arch -o DRR049557 DRR049557_1.fastq.bz2 DRR049557_2.fastq.bz2
```

`-t 4` allocates 4 cores to the task.  You can adapt this number to the machine where you run that command.

The file [`tagdust.fa`](./tagdust.fa) contains sequences of nanoCAGE linkers and of Nextera sequencing
linkers.

As per Illumina's [ Illumina Adapter Sequences Document
](http://support.illumina.com/downloads/illumina-customer-sequence-letter.html),
here is the copyright notice for the sequences that have a name starting with
"Nextera".

> Oligonucleotide sequences © 2016 Illumina, Inc. All rights reserved.
> Derivative works created by Illumina customers are authorized for use with
> Illumina instruments and products only. All other uses are strictly
> prohibited.

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

## Result

The following files should be produced:

```
DRR049557_BC_Caski_r1_READ1.fq
DRR049557_BC_Caski_r1_READ2.fq
DRR049557_BC_Caski_r2_READ1.fq
DRR049557_BC_Caski_r2_READ2.fq
DRR049557_BC_Caski_r3_READ1.fq
DRR049557_BC_Caski_r3_READ2.fq
DRR049557_BC_HaCaT_r1_READ1.fq
DRR049557_BC_HaCaT_r1_READ2.fq
DRR049557_BC_HaCaT_r2_READ1.fq
DRR049557_BC_HaCaT_r2_READ2.fq
DRR049557_BC_HaCaT_r3_READ1.fq
DRR049557_BC_HaCaT_r3_READ2.fq
DRR049557_BC_Hela_r1_READ1.fq
DRR049557_BC_Hela_r1_READ2.fq
DRR049557_BC_Hela_r2_READ1.fq
DRR049557_BC_Hela_r2_READ2.fq
DRR049557_BC_Hela_r3_READ1.fq
DRR049557_BC_Hela_r3_READ2.fq
DRR049557_BC_SiHa_r1_READ1.fq
DRR049557_BC_SiHa_r1_READ2.fq
DRR049557_BC_SiHa_r2_READ1.fq
DRR049557_BC_SiHa_r2_READ2.fq
DRR049557_BC_SiHa_r3_READ1.fq
DRR049557_BC_SiHa_r3_READ2.fq
DRR049557_BC_W12E_r1_READ1.fq
DRR049557_BC_W12E_r1_READ2.fq
DRR049557_BC_W12E_r2_READ1.fq
DRR049557_BC_W12E_r2_READ2.fq
DRR049557_BC_c33a_r1_READ1.fq
DRR049557_BC_c33a_r1_READ2.fq
DRR049557_BC_c33a_r2_READ1.fq
DRR049557_BC_c33a_r2_READ2.fq
DRR049557_BC_c33a_r3_READ1.fq
DRR049557_BC_c33a_r3_READ2.fq
DRR049557_un_READ1.fq
DRR049557_un_READ2.fq
```

The `tagdust` command also produces a summary file, `DRR049557_logfile.txt`.
In our test run, its content was:

```
[2016-07-07 09:29:55]	Tagdust 2.33, Copyright (C) 2013 Timo Lassmann <timolassmann@gmail.com>
[2016-07-07 09:29:55]	cmd: tagdust -t 4 -show_finger_seq -ref tagdust.fa -arch DRR049557.arch -o DRR049557 DRR049557_1.fastq.bz2 DRR049557_2.fastq.bz2 
[2016-07-07 09:29:55]	Start Run
--------------------------------------------------
[2016-07-07 09:29:55]	Looking at file:DRR049557_1.fastq.bz2
[2016-07-07 09:29:55]	Searching for best architecture in file 'DRR049557.arch'
[2016-07-07 09:30:26]	Using: -1 B:ACATGA,ATCATA,CACGTG,CGATGA,GAGATA,GCTCTC,GTATGA,TCGATA,AGTAGC,ATCGCA,CACTCT,CTGAGC,GAGCGT,GCTGCA,TATAGC,CACGAT,CTGACG -2 F:NNNNNNNN -3 S:TATAGGG -4 R:N 
[2016-07-07 09:30:26]	1.00 Confidence.
[2016-07-07 09:30:26]	Looking at file:DRR049557_2.fastq.bz2
[2016-07-07 09:30:26]	Searching for best architecture in file 'DRR049557.arch'
[2016-07-07 09:30:53]	Using: -1 R:N 
[2016-07-07 09:30:53]	1.00 Confidence.
[2016-07-07 09:31:10]	Determining threshold for read0.
[2016-07-07 09:31:20]	Long sequence found. Need to realloc model...
[2016-07-07 09:34:31]	Selected Threshold:: 2.658036
[2016-07-07 09:34:31]	Determining threshold for read1.
[2016-07-07 09:34:41]	Long sequence found. Need to realloc model...
[2016-07-07 09:34:42]	Selected Threshold:: 3.051751
[2016-07-07 09:35:01]	Detected casava 1.8 format.
[2016-07-08 00:24:06]	Done.

[2016-07-08 00:24:06]	DRR049557_1.fastq.bz2	Input file 0.
[2016-07-08 00:24:06]	DRR049557_2.fastq.bz2	Input file 1.
[2016-07-08 00:24:06]	217752577	total input reads
[2016-07-08 00:24:06]	3.05	selected threshold
[2016-07-08 00:24:06]	187155993	successfully extracted
[2016-07-08 00:24:06]	85.9%	extracted
[2016-07-08 00:24:06]	12041093	problems with architecture
[2016-07-08 00:24:06]	1251417	barcode / UMI not found
[2016-07-08 00:24:06]	0	too short
[2016-07-08 00:24:06]	10318	low complexity
[2016-07-08 00:24:06]	17293756	match artifacts:
[2016-07-08 00:24:06]	1750534	RT_(without_random_bases)
[2016-07-08 00:24:06]	15346	empty_(TS_linker_+_RT_reverse-complemented)
[2016-07-08 00:24:06]	470963	Nextera_501
[2016-07-08 00:24:06]	108	Nextera_502
[2016-07-08 00:24:06]	1349	Nextera_503
[2016-07-08 00:24:06]	391	Nextera_504
[2016-07-08 00:24:06]	14	Nextera_505
[2016-07-08 00:24:06]	1	Nextera_506
[2016-07-08 00:24:06]	73	Nextera_507
[2016-07-08 00:24:06]	3	Nextera_508
[2016-07-08 00:24:06]	11706684	Nextera_701
[2016-07-08 00:24:06]	374315	Nextera_702
[2016-07-08 00:24:06]	2972225	Nextera_703
[2016-07-08 00:24:06]	16	Nextera_704
[2016-07-08 00:24:06]	302	Nextera_705
[2016-07-08 00:24:06]	151	Nextera_706
[2016-07-08 00:24:06]	416	Nextera_707
[2016-07-08 00:24:06]	28	Nextera_708
[2016-07-08 00:24:06]	137	Nextera_709
[2016-07-08 00:24:06]	16	Nextera_710
[2016-07-08 00:24:06]	445	Nextera_711
[2016-07-08 00:24:06]	169	Nextera_712
[2016-07-08 00:24:06]	60	Nextera_501_Reversed:
[2016-07-08 00:24:06]	2	Nextera_502_Reversed:
[2016-07-08 00:24:06]	2	Nextera_503_Reversed:
[2016-07-08 00:24:06]	3	Nextera_504_Reversed:
[2016-07-08 00:24:06]	1	Nextera_505_Reversed:
[2016-07-08 00:24:06]	2	Nextera_703_Reversed:
```

The resulting FASTQ reads look like this:

```
@DRR049557.7 HWI-ST549:177:C6Y04ACXX:2:1101:1316:1986 length=51;FP:GTCAGGGG;RQ:30.19
GGGGAATCAGGGTTCGATTCCGGAGAGGGN
+
JJDDDDDDDDDDDDDDDEDDEDDDDDDDD#
```

The sequence of the fingerprint is indicated by the `FP:` tag in the semicolon-separated part of the read name.
