


Simple interrogation of FANTOM5's SDRF files
============================================

Some of the FANTOM5 metadata is distributed in the _Sample and Data
Relationship Format_ ([SDRF](https://en.wikipedia.org/wiki/SDRF)) format along
with the alignment files in <http://fantom.gsc.riken.jp/5/datafiles/>.  Each
line describes one FANTOM5 library, except the first line that is a header.

Here is one example downloaded from the FANTOM5 website.


```sh
curl http://fantom.gsc.riken.jp/5/datafiles/phase1.2/basic/human.cell_line.LQhCAGE/00_human.cell_line.LQhCAGE.hg19.assay_sdrf.txt
```

```
Extract Name	Comment [rna_tube]	Comment [sample_name]	Comment [organism]	Protocol REF	Parameter [rna_extraction]	Extract Name	Comment [rna_id]	Comment [Material Type]	Comment [comment_on_rna]	Protocol REF	Comment [lsid]	Parameter [library_protocol]	Library Name	Comment [library_id]	Protocol REF	Parameter [sequence_protocol]	Parameter [machine_name]	Parameter [run_name]	Parameter [flowcell.channel]	File Name	Comment [sequence_raw_file]	Protocol REF	Parameter [sex]	File Name
10789-110H6	10789-110H6	acute myeloid leukemia (FAB M0) cell line:Kasumi-3	Homo sapiens	standard whole cell	standard whole cell	10789-110H6.rna	10789-110H6.rna	total RNA	comment:Sample category:"FANTOM5 cell line"; company:"JAPAN HEALTH SCIENCES FOUNDATION - Health Science Research Resources Bank"; cell catalog:"jcrb1004 ";; rna_sample_type:total RNA; extract_method:standard whole cell; RIN:NULL; A260A280:1.9; A260A230:1.01	OP-HELISCOPE-CAGE-v5.0	1006	OP-HELISCOPE-CAGE-v5.0	CNhs13241	CNhs13241	OP-HELISCOPE-sequencing-v1.0	OP-HELISCOPE-sequencing-v1.0	NA	NA	fc2.ch06	/sequencedata/heliscope/srf/2011_05_17_R2_30Q_1100FOV	/sequencedata/heliscope/srf/2011_05_17_R2_30Q_1100FOV	delve_align	male	acute%20myeloid%20leukemia%20%28FAB%20M0%29%20cell%20line%3aKasumi-3.CNhs13241.10789-110H6.hg19.nobarcode.bam
10722-110A2	10722-110A2	leiomyosarcoma cell line:Hs 5.T	Homo sapiens	standard whole cell	standard whole cell	10722-110A2.rna	10722-110A2.rna	total RNA	comment:Sample category:"FANTOM5 cell line"; company:"ATCC"; cell catalog:"CRL-7822";; rna_sample_type:total RNA; extract_method:standard whole cell; RIN:NULL; A260A280:1.95; A260A230:1.39	OP-HELISCOPE-CAGE-v5.0	976	OP-HELISCOPE-CAGE-v5.0	CNhs12192	CNhs12192	OP-HELISCOPE-sequencing-v1.0	OP-HELISCOPE-sequencing-v1.0	NA	NA	fc1.ch21	/sequencedata/heliscope/srf/2011_01_28_R2_30Q_1100FOV	/sequencedata/heliscope/srf/2011_01_28_R2_30Q_1100FOV	delve_align	female	leiomyosarcoma%20cell%20line%3aHs%205%2eT.CNhs12192.10722-110A2.hg19.nobarcode.bam
10696-109G3	10696-109G3	mesodermal tumor cell line:HIRS-BM	Homo sapiens	standard whole cell	standard whole cell	10696-109G3.rna	10696-109G3.rna	total RNA	comment:Sample category:"FANTOM5 cell line"; company:"RIKEN Bioresource centre"; cell catalog:"RCB0978 ";; rna_sample_type:total RNA; extract_method:standard whole cell; RIN:NULL; A260A280:2.1; A260A230:1.49	OP-HELISCOPE-CAGE-v5.0	976	OP-HELISCOPE-CAGE-v5.0	CNhs12191	CNhs12191	OP-HELISCOPE-sequencing-v1.0	OP-HELISCOPE-sequencing-v1.0	NA	NA	fc1.ch24	/sequencedata/heliscope/srf/2011_01_28_R2_30Q_1100FOV	/sequencedata/heliscope/srf/2011_01_28_R2_30Q_1100FOV	delve_align	female	mesodermal%20tumor%20cell%20line%3aHIRS-BM.CNhs12191.10696-109G3.hg19.nobarcode.bam
10730-110B1	10730-110B1	non-small cell lung cancer cell line:NCI-H1385	Homo sapiens	standard whole cell	standard whole cell	10730-110B1.rna	10730-110B1.rna	total RNA	comment:Sample category:"FANTOM5 cell line"; company:"ATCC"; cell catalog:"CRL-5867";; rna_sample_type:total RNA; extract_method:standard whole cell; RIN:NULL; A260A280:2; A260A230:1.02	OP-HELISCOPE-CAGE-v5.0	976	OP-HELISCOPE-CAGE-v5.0	CNhs12193	CNhs12193	OP-HELISCOPE-sequencing-v1.0	OP-HELISCOPE-sequencing-v1.0	NA	NA	fc1.ch18	/sequencedata/heliscope/srf/2011_01_28_R2_30Q_1100FOV	/sequencedata/heliscope/srf/2011_01_28_R2_30Q_1100FOV	delve_align	female	non-small%20cell%20lung%20cancer%20cell%20line%3aNCI-H1385.CNhs12193.10730-110B1.hg19.nobarcode.bam
```


As you see there are long lines with multiple fields.  We will use the
3<sup>rd</sup> (`Comment [sample_name]`) and the 14<sup>th</sup> (`Library
Name`).  Here is their contents in the above example.


```sh
curl http://fantom.gsc.riken.jp/5/datafiles/phase1.2/basic/human.cell_line.LQhCAGE/00_human.cell_line.LQhCAGE.hg19.assay_sdrf.txt | cut -f 3,14
```

```
Comment [sample_name]	Library Name
acute myeloid leukemia (FAB M0) cell line:Kasumi-3	CNhs13241
leiomyosarcoma cell line:Hs 5.T	CNhs12192
mesodermal tumor cell line:HIRS-BM	CNhs12191
non-small cell lung cancer cell line:NCI-H1385	CNhs12193
```


Simple shell functions to return or query library IDs
-----------------------------------------------------

The file `sdrf.sh` distributed here defines two simple shell functions,
`SDRFlib` and `SDRFdesc`.  These functions expect the FANTOM5 SDRF files to be
in the current directory or in a directory indicated by the environment
variable `F5`.  In the examples later, we load these functions with the command
`source sdrf.sh`.


```sh
cat sdrf.sh
```

```
# Return all libraries.

function SDRFlib {
for SDRF in $(find ${F5-.} -name '*sdrf.txt')
do
  grep "$1" $SDRF | cut -f14
done
}

# Returns the description of a library.

function SDRFdesc {
for SDRF in $(find ${F5-.} -name '*sdrf.txt')
do
  grep "$1" $SDRF | cut -f3
done
}
```


Download the SDRF files
-----------------------

The FANTOM5 data is organised in multiple subdirectories, each containing a
single SDRF file.  Here, we will download all the SDRF files in a single
directory.


```sh
PHASE='phase1.2'
echo 'mget */*sdrf.txt' |
  lftp http://fantom.gsc.riken.jp/5/datafiles/$PHASE/basic
wc -l *sdrf.txt
```

```
##       5 00_human.cell_line.LQhCAGE.hg19.assay_sdrf.txt
##     258 00_human.cell_line.hCAGE.hg19.assay_sdrf.txt
##      13 00_human.fractionation.hCAGE.hg19.assay_sdrf.txt
##      51 00_human.primary_cell.LQhCAGE.hg19.assay_sdrf.txt
##     489 00_human.primary_cell.hCAGE.hg19.assay_sdrf.txt
##      30 00_human.timecourse.hCAGE.hg19.assay_sdrf.txt
##     153 00_human.tissue.hCAGE.hg19.assay_sdrf.txt
##       2 00_mouse.cell_line.hCAGE.mm9.assay_sdrf.txt
##      27 00_mouse.primary_cell.LQhCAGE.mm9.assay_sdrf.txt
##      85 00_mouse.primary_cell.hCAGE.mm9.assay_sdrf.txt
##      29 00_mouse.qualitycontrol.hCAGE.mm9.assay_sdrf.txt
##      21 00_mouse.timecourse.hCAGE.mm9.assay_sdrf.txt
##       2 00_mouse.tissue.LQhCAGE.mm9.assay_sdrf.txt
##     236 00_mouse.tissue.hCAGE.mm9.assay_sdrf.txt
##    1401 total
```


_Note: the SDRF files have one header line, so the total number of libraries is
not the total number of lines._

Library names (column 14) start with `CNhs` for HeliScopeCAGE libraries.  Phase
1.2 only contains HeliScopeCAGE libraries.


```sh
cut -f 14 *sdrf.txt |
  grep -v 'Library Name' |   # Remove header lines.
  cut -c 1-4 |
  sort | uniq -c | sort -n
```

```
##    1387 CNhs
```


Examples of use
---------------

Find all libraries with _pancreas_ in their description.


```bash
source sdrf.sh

SDRFlib pancreas
```

```
## CNhs11756
## CNhs10486
## CNhs11012
## CNhs11042
## CNhs11003
## CNhs10599
## CNhs10580
## CNhs11105
## CNhs11138
## CNhs11139
## CNhs11136
## CNhs11094
## CNhs11182
## CNhs11731
## CNhs11732
## CNhs11733
## CNhs11814
```


Find the library name of all the libraries with _pancreas_ in their
description.


```bash
source sdrf.sh

for lib in $(SDRFlib pancreas);
do
  echo -ne "$lib\t"
  SDRFdesc $lib
done
```

```
## CNhs11756	pancreas, adult, donor1
## CNhs10486	pancreas, adult
## CNhs11012	pancreas, embryo E14
## CNhs11042	pancreas, embryo E15
## CNhs11003	pancreas, embryo E16
## CNhs10599	pancreas, embryo E17
## CNhs10580	pancreas, embryo E18
## CNhs11105	pancreas, neonate N00
## CNhs11138	pancreas, neonate N01
## CNhs11139	pancreas, neonate N02
## CNhs11136	pancreas, neonate N16
## CNhs11094	pancreas, neonate N25
## CNhs11182	pancreas, neonate N30
## CNhs11731	embryonic pancreas cell line:1B2C6
## CNhs11732	embryonic pancreas cell line:1C3D3
## CNhs11733	embryonic pancreas cell line:1C3IKEI
## CNhs11814	embryonic pancreas cell line:2C6
```


Conclusion
----------

We will use these functions to get list the name of all human and mouse FANTOM5
libraries in the next tutorial.

Note also that there is more metadata and more expressive ways to interrogate
it, see <http://fantom.gsc.riken.jp/views/> for details.
