Demultiplexing nanoCAGE reads
=============================

_Multiplexing_ is done by adding artificial sequences to molecules, to encode their
sample of origin, so that multiple samples can be sequenced together (to save costs).
_Demultiplexing_ is to retreive this information from the sequence reads, and recode
it as an addition to the read name or by grouping the reads in separate files.  Some
multiplexing strategies are supported natively by some sequencers, but a lot of
custom designs also exist.

In this example, we will demultiplex [nanoCAGE](https://population-transcriptomics.org/nanoCAGE/)
data from DDBJ using install [TagDust 2](https://sourceforge.net/projects/tagdust/),
after removing PCR duplicates with [clumpify](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify-guide/).
These commands have been tested on a virtual machine running Linux.

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
## If [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) is not available, install it.

```
wget --trust-server-names https://sourceforge.net/projects/bbmap/files/latest/download
tar xvfz tar xvfz BBMap_*.tar.gz
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

## Deduplicate

```
./bbmap/clumpify.sh in=DRR049557_1.fastq.bz2 in2=DRR049557_2.fastq.bz2 out=DRR049557_1.demul.fastq.gz out2=DRR049557_2.demul.fastq.gz shortname=t dedupe=t dupesubs=2 addcount=t
```

output:

```
java version "1.7.0_131"
OpenJDK Runtime Environment (IcedTea 2.6.9) (7u131-2.6.9-2~deb8u1)
OpenJDK 64-Bit Server VM (build 24.131-b00, mixed mode)
java -ea -Xmx67638m -Xms67638m -cp /home/plessy/CAGEscan_3.0_realign2/bbmap/current/ clump.Clumpify in=DRR049557_1.fastq.bz2 in2=DRR049557_2.fastq.bz2 out=DRR049557_1.demul.fastq.gz out2=DRR049557_2.demul.fastq.gz shortname=t dedupe=t dupesubs=2 addcount=t
Executing clump.Clumpify [in=DRR049557_1.fastq.bz2, in2=DRR049557_2.fastq.bz2, out=DRR049557_1.demul.fastq.gz, out2=DRR049557_2.demul.fastq.gz, shortname=t, dedupe=t, dupesubs=2, addcount=t]


Clumpify version 37.55
Read Estimate:          628733779
Memory Estimate:        300563 MB
Memory Available:       53113 MB
Set groups to 53
Executing clump.KmerSplit [in1=DRR049557_1.fastq.bz2, in2=DRR049557_2.fastq.bz2, out=DRR049557_1.demul_clumpify_p1_temp%_7c966ffa9f3ece8a.fastq.gz, out2=null, groups=53, ecco=false, addname=f, shortname=t, unpair=false, repair=f, namesort=f, ow=true, dedupe=t]

Reset INTERLEAVED to false because paired input files were specified.
Set INTERLEAVED to false
Input is being processed as paired
Writing interleaved.
Made a comparator with k=31, seed=1, border=1, hashes=4
Time:                         	1760.816 seconds.
Reads Processed:        435m 	247.33k reads/sec
Bases Processed:      22210m 	12.61m bases/sec
Executing clump.KmerSort3 [in1=DRR049557_1.demul_clumpify_p1_temp%_7c966ffa9f3ece8a.fastq.gz, in2=null, out=DRR049557_1.demul.fastq.gz, out2=DRR049557_2.demul.fastq.gz, groups=53, ecco=f, addname=false, shortname=f, unpair=f, repair=false, namesort=false, ow=true, dedupe=t]

Making comparator.
Made a comparator with k=31, seed=1, border=1, hashes=4
Making 2 fetch threads.
Starting threads.
Fetching reads.
Fetched 3734838 reads: 	12.909 seconds.
Making clumps.
Clump time: 	0.161 seconds.
Deduping.
Dedupe time: 	0.814 seconds.
Writing.
Fetching reads.
Fetched 3930647 reads: 	0.000 seconds.
Making clumps.
Clump time: 	0.239 seconds.
Deduping.
Dedupe time: 	6.698 seconds.
Writing.
Fetching reads.
Fetched 3766454 reads: 	14.516 seconds.
Making clumps.
Clump time: 	0.173 seconds.
Deduping.
Dedupe time: 	0.729 seconds.
Writing.
Fetching reads.
Fetched 4870043 reads: 	2.077 seconds.
Making clumps.
Clump time: 	0.376 seconds.
Deduping.
Dedupe time: 	4.397 seconds.
Writing.
Fetching reads.
Fetched 4267639 reads: 	9.465 seconds.
Making clumps.
Clump time: 	0.260 seconds.
Deduping.
Dedupe time: 	2.107 seconds.
Writing.
Fetching reads.
Fetched 3712221 reads: 	0.000 seconds.
Making clumps.
Clump time: 	0.188 seconds.
Deduping.
Dedupe time: 	12.646 seconds.
Writing.
Fetching reads.
Fetched 4078710 reads: 	13.333 seconds.
Making clumps.
Clump time: 	0.157 seconds.
Deduping.
Dedupe time: 	0.936 seconds.
Writing.
Fetching reads.
Fetched 4074164 reads: 	1.285 seconds.
Making clumps.
Clump time: 	0.210 seconds.
Deduping.
Dedupe time: 	1.336 seconds.
Writing.
Fetching reads.
Fetched 3648178 reads: 	15.841 seconds.
Making clumps.
Clump time: 	0.130 seconds.
Deduping.
Dedupe time: 	0.415 seconds.
Writing.
Fetching reads.
Fetched 4019769 reads: 	2.093 seconds.
Making clumps.
Clump time: 	0.180 seconds.
Deduping.
Dedupe time: 	0.613 seconds.
Writing.
Fetching reads.
Fetched 3814138 reads: 	17.609 seconds.
Making clumps.
Clump time: 	0.136 seconds.
Deduping.
Dedupe time: 	0.450 seconds.
Writing.
Fetching reads.
Fetched 4220795 reads: 	0.000 seconds.
Making clumps.
Clump time: 	0.163 seconds.
Deduping.
Dedupe time: 	0.979 seconds.
Writing.
Fetching reads.
Fetched 4358672 reads: 	19.616 seconds.
Making clumps.
Clump time: 	0.145 seconds.
Deduping.
Dedupe time: 	0.806 seconds.
Writing.
Fetching reads.
Fetched 4557703 reads: 	0.000 seconds.
Making clumps.
Clump time: 	0.249 seconds.
Deduping.
Dedupe time: 	2.897 seconds.
Writing.
Fetching reads.
Fetched 3741505 reads: 	15.113 seconds.
Making clumps.
Clump time: 	0.145 seconds.
Deduping.
Dedupe time: 	0.702 seconds.
Writing.
Fetching reads.
Fetched 4289768 reads: 	2.488 seconds.
Making clumps.
Clump time: 	0.274 seconds.
Deduping.
Dedupe time: 	1.044 seconds.
Writing.
Fetching reads.
Fetched 3916702 reads: 	13.397 seconds.
Making clumps.
Clump time: 	0.145 seconds.
Deduping.
Dedupe time: 	0.499 seconds.
Writing.
Fetching reads.
Fetched 3754220 reads: 	0.954 seconds.
Making clumps.
Clump time: 	0.237 seconds.
Deduping.
Dedupe time: 	2.925 seconds.
Writing.
Fetching reads.
Fetched 3636764 reads: 	12.835 seconds.
Making clumps.
Clump time: 	2.806 seconds.
Deduping.
Dedupe time: 	0.333 seconds.
Writing.
Fetching reads.
Fetched 4031662 reads: 	1.662 seconds.
Making clumps.
Clump time: 	0.270 seconds.
Deduping.
Dedupe time: 	0.825 seconds.
Writing.
Fetching reads.
Fetched 3667871 reads: 	14.132 seconds.
Making clumps.
Clump time: 	0.146 seconds.
Deduping.
Dedupe time: 	0.401 seconds.
Writing.
Fetching reads.
Fetched 3912289 reads: 	0.000 seconds.
Making clumps.
Clump time: 	0.204 seconds.
Deduping.
Dedupe time: 	0.845 seconds.
Writing.
Fetching reads.
Fetched 3640157 reads: 	17.439 seconds.
Making clumps.
Clump time: 	0.153 seconds.
Deduping.
Dedupe time: 	0.395 seconds.
Writing.
Fetching reads.
Fetched 4083911 reads: 	1.254 seconds.
Making clumps.
Clump time: 	0.260 seconds.
Deduping.
Dedupe time: 	1.121 seconds.
Writing.
Fetching reads.
Fetched 4221907 reads: 	16.722 seconds.
Making clumps.
Clump time: 	0.170 seconds.
Deduping.
Dedupe time: 	4.061 seconds.
Writing.
Fetching reads.
Fetched 4273035 reads: 	0.000 seconds.
Making clumps.
Clump time: 	0.198 seconds.
Deduping.
Dedupe time: 	1.055 seconds.
Writing.
Fetching reads.
Fetched 3674670 reads: 	15.864 seconds.
Making clumps.
Clump time: 	0.138 seconds.
Deduping.
Dedupe time: 	0.450 seconds.
Writing.
Fetching reads.
Fetched 3951075 reads: 	0.000 seconds.
Making clumps.
Clump time: 	0.197 seconds.
Deduping.
Dedupe time: 	1.089 seconds.
Writing.
Fetching reads.
Fetched 3989937 reads: 	18.870 seconds.
Making clumps.
Clump time: 	0.176 seconds.
Deduping.
Dedupe time: 	0.972 seconds.
Writing.
Fetching reads.
Fetched 5098672 reads: 	1.877 seconds.
Making clumps.
Clump time: 	15.550 seconds.
Deduping.
Dedupe time: 	2.143 seconds.
Writing.
Fetching reads.
Fetched 4018745 reads: 	11.249 seconds.
Making clumps.
Clump time: 	0.185 seconds.
Deduping.
Dedupe time: 	0.578 seconds.
Writing.
Fetching reads.
Fetched 6944917 reads: 	13.619 seconds.
Making clumps.
Clump time: 	0.646 seconds.
Deduping.
Dedupe time: 	6.677 seconds.
Writing.
Fetching reads.
Fetched 3973142 reads: 	0.000 seconds.
Making clumps.
Clump time: 	0.318 seconds.
Deduping.
Dedupe time: 	1.203 seconds.
Writing.
Fetching reads.
Fetched 3859868 reads: 	6.863 seconds.
Making clumps.
Clump time: 	0.124 seconds.
Deduping.
Dedupe time: 	0.457 seconds.
Writing.
Fetching reads.
Fetched 3910877 reads: 	11.144 seconds.
Making clumps.
Clump time: 	0.122 seconds.
Deduping.
Dedupe time: 	0.545 seconds.
Writing.
Fetching reads.
Fetched 4076702 reads: 	4.728 seconds.
Making clumps.
Clump time: 	0.168 seconds.
Deduping.
Dedupe time: 	0.623 seconds.
Writing.
Fetching reads.
Fetched 3959347 reads: 	15.109 seconds.
Making clumps.
Clump time: 	0.170 seconds.
Deduping.
Dedupe time: 	0.514 seconds.
Writing.
Fetching reads.
Fetched 5220396 reads: 	0.000 seconds.
Making clumps.
Clump time: 	0.329 seconds.
Deduping.
Dedupe time: 	5.364 seconds.
Writing.
Fetching reads.
Fetched 3856875 reads: 	13.424 seconds.
Making clumps.
Clump time: 	0.145 seconds.
Deduping.
Dedupe time: 	0.471 seconds.
Writing.
Fetching reads.
Fetched 4914934 reads: 	3.800 seconds.
Making clumps.
Clump time: 	0.373 seconds.
Deduping.
Dedupe time: 	5.099 seconds.
Writing.
Fetching reads.
Fetched 3878886 reads: 	8.116 seconds.
Making clumps.
Clump time: 	0.138 seconds.
Deduping.
Dedupe time: 	0.532 seconds.
Writing.
Fetching reads.
Fetched 3881790 reads: 	3.793 seconds.
Making clumps.
Clump time: 	0.149 seconds.
Deduping.
Dedupe time: 	0.545 seconds.
Writing.
Fetching reads.
Fetched 3518762 reads: 	11.311 seconds.
Making clumps.
Clump time: 	0.153 seconds.
Deduping.
Dedupe time: 	0.377 seconds.
Writing.
Fetching reads.
Fetched 3673340 reads: 	3.251 seconds.
Making clumps.
Clump time: 	0.150 seconds.
Deduping.
Dedupe time: 	0.364 seconds.
Writing.
Fetching reads.
Fetched 4302883 reads: 	12.704 seconds.
Making clumps.
Clump time: 	2.050 seconds.
Deduping.
Dedupe time: 	0.647 seconds.
Writing.
Fetching reads.
Fetched 3934267 reads: 	0.720 seconds.
Making clumps.
Clump time: 	0.189 seconds.
Deduping.
Dedupe time: 	0.691 seconds.
Writing.
Fetching reads.
Fetched 3814885 reads: 	14.683 seconds.
Making clumps.
Clump time: 	0.130 seconds.
Deduping.
Dedupe time: 	0.477 seconds.
Writing.
Fetching reads.
Fetched 4202796 reads: 	1.591 seconds.
Making clumps.
Clump time: 	0.207 seconds.
Deduping.
Dedupe time: 	1.284 seconds.
Writing.
Fetching reads.
Fetched 3615950 reads: 	12.948 seconds.
Making clumps.
Clump time: 	0.167 seconds.
Deduping.
Dedupe time: 	0.443 seconds.
Writing.
Fetching reads.
Fetched 4315804 reads: 	3.579 seconds.
Making clumps.
Clump time: 	0.190 seconds.
Deduping.
Dedupe time: 	3.146 seconds.
Writing.
Fetching reads.
Fetched 3672194 reads: 	12.653 seconds.
Making clumps.
Clump time: 	0.129 seconds.
Deduping.
Dedupe time: 	0.401 seconds.
Writing.
Fetching reads.
No more reads to fetch.
Adding poison.
Fetched 5742151 reads: 	2.546 seconds.
Making clumps.
Clump time: 	0.324 seconds.
Deduping.
Dedupe time: 	3.361 seconds.
Writing.
Fetching reads.
Encountered poison; count=1
A fetch thread finished.
No more reads to fetch.
Adding poison.
Fetched 3524950 reads: 	6.771 seconds.
Making clumps.
Clump time: 	0.131 seconds.
Deduping.
Dedupe time: 	0.374 seconds.
Writing.
Closing fetch threads.
A fetch thread finished.
Closed fetch threads.
Waiting for writing to complete.
Write time: 	2.160 seconds.
Done!
Time:                         	539.753 seconds.
Reads Processed:        435m 	806.86k reads/sec
Bases Processed:      22210m 	41.15m bases/sec

Reads In:          435505154
Clumps Formed:      17115599
Duplicates Found:  238273388
Total time: 	2304.207 seconds.
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

`-t 4` allocates 4 cores to the task.  You can adapt this number to the machine
where you run that command.  Note that even when using multiple cores, TagDust may
take hours before finishing.  Look at the timestamps in the [example output](#example-log-file)
below for an example.

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
### Example log file

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
