---
title: "FANTOM5 DPI peak shape"
author: "Charles Plessy"
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---

FANTOM5 promoter peak shape
===========================

The FANTOM5 expression atlas paper (Forrest et al., 2014) describes promoter
models produced by Decomposition-based Peak Identification (DPI) on CAGE data.
The DPI peaks have different widths according to the underlying CAGE signal;
the goal of this page is to use this width as a proxy to classify the promoters
in _sharp_ or _broad_ categories.

There are multiple limitations to this approach.  In particular, promoter shape
can change across tissues, and also across time in the same tissue (for instance
during development).  Indeed, the promoter shape analysis in Forrest 2014 (ext.
data fig. 2) was performed with a more robust approach, using the 0.1–0.9
interquantile width of the CAGE expression data accross all samples.  Moreover,
genes may have multiple promoters of different shapes.  Nevertheless, there are
cases, such as when handlings large groups of genes, where it may be useful to
estimate if these genes tend to have sharper or broader promoters.  This is the
aim of the approach presented here.


Width per DPI peak
------------------

First, run the example in `../FANTOM5_DPI_BED_file_with_annotation/`.

Then, load the result files.


```r
(human <- rtracklayer::import.bed("../FANTOM5_DPI_BED_file_with_annotation/hg19.cage_peak_phase1and2combined_anncoord.bed"))
```

```
## GRanges object with 201802 ranges and 2 metadata columns:
##            seqnames               ranges strand |
##               <Rle>            <IRanges>  <Rle> |
##        [1]     chr1     [564572, 564600]      + |
##        [2]     chr1     [564640, 564649]      + |
##        [3]     chr1     [565267, 565278]      + |
##        [4]     chr1     [565479, 565483]      + |
##        [5]     chr1     [565510, 565541]      + |
##        ...      ...                  ...    ... .
##   [201798]     chrY [28817277, 28817283]      - |
##   [201799]     chrY [58856052, 58856076]      + |
##   [201800]     chrY [58884632, 58884640]      + |
##   [201801]     chrY [58903886, 58903887]      + |
##   [201802]     chrY [59027989, 59028000]      + |
##                                   name     score
##                            <character> <numeric>
##        [1]                 p1@MTND1P23      2398
##        [2]                 p3@MTND1P23       220
##        [3]                 p3@MTND2P28       535
##        [4]                 p4@MTND2P28       106
##        [5]                 p1@MTND2P28      3594
##        ...                         ...       ...
##   [201798] p@chrY:28817276..28817283,-        39
##   [201799] p@chrY:58856051..58856076,+       422
##   [201800] p@chrY:58884631..58884640,+       129
##   [201801] p@chrY:58903885..58903887,+       131
##   [201802] p@chrY:59027988..59028000,+       190
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```

```r
(mouse <- rtracklayer::import.bed("../FANTOM5_DPI_BED_file_with_annotation/mm9.cage_peak_phase1and2combined_anncoord.bed"))
```

```
## GRanges object with 158966 ranges and 2 metadata columns:
##            seqnames             ranges strand |                      name
##               <Rle>          <IRanges>  <Rle> |               <character>
##        [1]     chr1 [3309586, 3309588]      - | p@chr1:3309585..3309588,-
##        [2]     chr1 [3367868, 3367870]      - | p@chr1:3367867..3367870,-
##        [3]     chr1 [3479231, 3479234]      - | p@chr1:3479230..3479234,-
##        [4]     chr1 [3644977, 3644980]      - | p@chr1:3644976..3644980,-
##        [5]     chr1 [3657916, 3657919]      - | p@chr1:3657915..3657919,-
##        ...      ...                ...    ... .                       ...
##   [158962]     chrY [2858252, 2858257]      + | p@chrY:2858251..2858257,+
##   [158963]     chrY [2858345, 2858362]      + | p@chrY:2858344..2858362,+
##   [158964]     chrY [2877093, 2877118]      + | p@chrY:2877092..2877118,+
##   [158965]     chrY [2890363, 2890366]      + | p@chrY:2890362..2890366,+
##   [158966]     chrY [2898757, 2898767]      - | p@chrY:2898756..2898767,-
##                score
##            <numeric>
##        [1]       153
##        [2]       569
##        [3]       163
##        [4]       493
##        [5]       904
##        ...       ...
##   [158962]       157
##   [158963]       277
##   [158964]        62
##   [158965]       101
##   [158966]       120
##   -------
##   seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Then, just use Bioconductor to calculate the width


```r
human$width <- GenomicRanges::width(human)
mouse$width <- GenomicRanges::width(mouse)
```


Promoter shape per gene
------------------------

Here we attribute to a given gene symbol the properties if its promoter with the
highest expression, which by definition has a name starting with `p1`.



```r
human <- human[substr(human$name, 1, 3) == "p1@"]
human$name  <- factor(sub("p1@", "", human$name))
human$name  <- factor(sub(",.*", "", human$name))
human
```

```
## GRanges object with 25444 ranges and 3 metadata columns:
##           seqnames               ranges strand |            name     score
##              <Rle>            <IRanges>  <Rle> |        <factor> <numeric>
##       [1]     chr1     [564572, 564600]      + |        MTND1P23      2398
##       [2]     chr1     [565510, 565541]      + |        MTND2P28      3594
##       [3]     chr1     [566902, 566923]      + |      uc001aaz.2      1281
##       [4]     chr1     [568913, 568941]      + | ENST00000467115     12638
##       [5]     chr1     [569902, 569922]      + | ENST00000416718      7865
##       ...      ...                  ...    ... .             ...       ...
##   [25440]     chrY [21906595, 21906622]      - |           KDM5D     33978
##   [25441]     chrY [22737620, 22737669]      + |          EIF1AY    127236
##   [25442]     chrY [22918022, 22918033]      + |          RPS4Y2       284
##   [25443]     chrY [23613728, 23613735]      - |         CYorf17        39
##   [25444]     chrY [27869489, 27869491]      + | ENST00000516617      1080
##               width
##           <integer>
##       [1]        29
##       [2]        32
##       [3]        22
##       [4]        29
##       [5]        21
##       ...       ...
##   [25440]        28
##   [25441]        50
##   [25442]        12
##   [25443]         8
##   [25444]         3
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```

```r
mouse <- mouse[substr(mouse$name, 1, 3) == "p1@"]
mouse$name  <- factor(sub("p1@", "", mouse$name))
mouse$name  <- factor(sub(",.*", "", mouse$name))
mouse
```

```
## GRanges object with 22636 ranges and 3 metadata columns:
##           seqnames             ranges strand |          name     score
##              <Rle>          <IRanges>  <Rle> |      <factor> <numeric>
##       [1]     chr1 [3661753, 3661814]      - |          Xkr4      4638
##       [2]     chr1 [4350338, 4350393]      - |           Rp1      4733
##       [3]     chr1 [4483680, 4483695]      - |         Sox17     13842
##       [4]     chr1 [4775753, 4775825]      - |        Mrpl15     83671
##       [5]     chr1 [4775854, 4775883]      + | A930006A01Rik       115
##       ...      ...                ...    ... .           ...       ...
##   [22632]     chrY [ 347031,  347056]      + |       Eif2s3y     72031
##   [22633]     chrY [ 582129,  582195]      - |           Uty     39950
##   [22634]     chrY [ 623007,  623058]      - |         Ddx3y     43152
##   [22635]     chrY [ 796215,  796225]      - |         Usp9y       179
##   [22636]     chrY [1426342, 1426344]      - |          Zfy2       143
##               width
##           <integer>
##       [1]        62
##       [2]        56
##       [3]        16
##       [4]        73
##       [5]        30
##       ...       ...
##   [22632]        26
##   [22633]        67
##   [22634]        52
##   [22635]        11
##   [22636]         3
##   -------
##   seqinfo: 22 sequences from an unspecified genome; no seqlengths
```


Save lists
----------


```r
write.table( file = "human_top_promoter_width.txt"
           , GenomicRanges::mcols(human)
           , row.names = FALSE, quote = FALSE, sep = "\t")

write.table( file = "mouse_top_promoter_width.txt"
           , GenomicRanges::mcols(mouse)
           , row.names = FALSE, quote = FALSE, sep = "\t")
```


Session information
-------------------


```r
sessionInfo()
```

```
## R version 3.3.3 (2017-03-06)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 9 (stretch)
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.9                knitr_1.20                
##  [3] XVector_0.14.0             magrittr_1.5              
##  [5] GenomicAlignments_1.10.0   GenomicRanges_1.26.2      
##  [7] BiocGenerics_0.20.0        IRanges_2.8.1             
##  [9] BiocParallel_1.8.1         lattice_0.20-34           
## [11] stringr_1.3.1              GenomeInfoDb_1.10.3       
## [13] tools_3.3.3                grid_3.3.3                
## [15] SummarizedExperiment_1.4.0 parallel_3.3.3            
## [17] Biobase_2.34.0             htmltools_0.3.5           
## [19] yaml_2.1.14                rprojroot_1.3-2           
## [21] digest_0.6.11              Matrix_1.2-7.1            
## [23] rtracklayer_1.34.1         S4Vectors_0.12.1          
## [25] bitops_1.0-6               RCurl_1.95-4.8            
## [27] evaluate_0.10              rmarkdown_1.10            
## [29] stringi_1.2.4              backports_1.1.2           
## [31] Biostrings_2.42.1          Rsamtools_1.26.1          
## [33] stats4_3.3.3               XML_3.98-1.5
```
