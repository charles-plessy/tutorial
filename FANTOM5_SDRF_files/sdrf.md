


Simplistic interrogation of FANTOM5's SDRF files
================================================

Simple shell functions to return or query library IDs
-----------------------------------------------------

should also work with $F5 pointed to the full mirror ?


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


_Note: the SDRF files have one header line, so the total number of libraries is not the total number of lines._

Library names are in column 14 and start with `CNhs` for HeliScopeCAGE
libraries.  Phase 1.2 only contains HeliScopeCAGE libraries.


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

Find example that also return sRNA libraries


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


