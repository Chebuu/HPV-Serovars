Simulating PCR Diagnostics to Discriminate HPV Serotypes by DNA Sequence
================

`---`

``` yaml
title: Title of submission to PLOS journal
author:
  - name: Alice Anonymous
    email: alice@example.com
    affiliation: Some Institute of Technology
    corresponding: alice@example.com
  - name: Bob Security
    email: bob@example.com
    affiliation: Another University
address:
  - code: Some Institute of Technology
    address: Department, Street, City, State, Zip
  - code: Another University
    address: Department, Street, City, State, Zip
abstract: |
  Lorem ipsum dolor sit amet, consectetur adipiscing elit.
  
author_summary: |
  Lorem ipsum dolor sit amet, consectetur adipiscing elit.
bibliography: citations.bib
includes:
  in_header: 
    - preamble.tex
    - common.tex
output: rticles::plos_article
```

`---`

    This PNAS journal template is provided to help you write your work in the
    correct journal format. Instructions for use are provided below.
    
    Note: please start your introduction WITHOUT including the word
    "Introduction" as a section heading (except for math articles in the
    Physical Sciences section); this heading is implied in the first
    paragraphs.

#### Installing OligAarrayAux

Download
[OligoArrayAux](http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux)
and see the following vignette for installation instructions.

``` r
library(devtools)

devtools::install_github('chebuu/Design-Group-Specific-Primers')

library(Design_Group_Specific_Primers)

vignette('Installing-OligoArrayAux', package = 'Design_Group_Specific_Primers')
```

# Introduction

The first infectious cause of cancer was identified in 1911 by Peyton
Rous who demonstrated the transmissibility of a tumor in fowl through
injection of sub 0.2 micron filtrate of the tumor Rous (1911).

``` r
if (!requireNamespace("BiocManager", quietly = T))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")
```

``` r
library(dplyr)
library(tidyr)
library(stringr)

library(Biostrings)
library(DECIPHER)
```

``` r
browseVignettes("DECIPHER")
```

# Data

I included GenBank queries for all the HPV DNA in this repo. If you want
to run them just paste in the search bar at:
<https://www.ncbi.nlm.nih.gov/nuccore>

## All HPV

<https://www.ncbi.nlm.nih.gov/nuccore>

``` sql
 "Human papillomavirus"[Primary Organism] 
 AND viruses[filter] 
 NOT Polyamides[All Fields] 
 NOT Method[All Fields] 
 NOT Patent[All Fields]
```

## Complete Genomes

<https://www.ncbi.nlm.nih.gov/nuccore>

``` sql
"Human papillomavirus"[Primary Organism] 
AND "complete genome"[All Fields]
NOT Isolate[Title] 
```

## Oral Isolates

<https://www.ncbi.nlm.nih.gov/nuccore>

``` sql
https://www.ncbi.nlm.nih.gov/nuccore
     "Human papillomavirus"[Primary Organism] 
 AND viruses[filter] 
 NOT Polyamides[All Fields] 
 NOT Method[All Fields] 
 NOT Patent[All Fields] 
 AND Oral[All Fields] 
```

### L1 gene

A collection of *l1* DNA sequences (avg. \~600bp) from various oral
isolates of HPV (Extracted from results of GenBank query above) are
aligned and loaded in SQLite (RAM) along with annotations from GenBank
that describe the serotype of each isolate.

``` r
doAlignment <- function(xset) {
    useqs <- unique(xset)
    uidxs <- match(xset, useqs)
    aseqs <- AlignSeqs(useqs, verbose=F) 
    aseqs[uidxs]
}

dbConn <- dbConnect(SQLite(), ':memory:')

Oral.L1.seqs <- readDNAStringSet('HPV.oral.L1.fasta') %>% doAlignment
Oral.L1.vars <- read.csv('HPV.oral.L1.csv', stringsAsFactors = F)

head(Oral.L1.seqs)
```

    ##   A DNAStringSet instance of length 6
    ##     width seq                                               names               
    ## [1]   607 -----------------------...----------------------- gi|944543704|gb|K...
    ## [2]   607 -----------------------...----------------------- gi|944543703|gb|K...
    ## [3]   607 -----------------------...----------------------- gi|944543701|gb|K...
    ## [4]   607 -----------------------...----------------------- gi|944543699|gb|K...
    ## [5]   607 -----------------------...----------------------- gi|944543697|gb|K...
    ## [6]   607 -----------------------...----------------------- gi|944543695|gb|K...

``` r
Seqs2DB(Oral.L1.seqs, 'XStringSet', dbConn, '')
```

    ## Adding 31 sequences to the database.
    ## 
    ## 31 total sequences in table Seqs.
    ## Time difference of 0.11 secs

``` r
Add2DB(Oral.L1.vars %>% mutate(identifier = SVAR), dbConn)
```

    ## Expression:
    ## alter table Seqs add column GI INTEGER
    ## 
    ## Expression:
    ## update Seqs set GI = :GI where row_names = :row_names
    ## 
    ## Expression:
    ## alter table Seqs add column SVAR INTEGER
    ## 
    ## Expression:
    ## update Seqs set SVAR = :SVAR where row_names = :row_names
    ## 
    ## Expression:
    ## alter table Seqs add column ISO INTEGER
    ## 
    ## Expression:
    ## update Seqs set ISO = :ISO where row_names = :row_names
    ## 
    ## Expression:
    ## alter table Seqs add column GENE INTEGER
    ## 
    ## Expression:
    ## update Seqs set GENE = :GENE where row_names = :row_names
    ## 
    ## Expression:
    ## update Seqs set identifier = :identifier where row_names = :row_names
    ## 
    ## Added to table Seqs:  "GI" and "SVAR" and "ISO" and "GENE" and "identifier".
    ## 
    ## Time difference of 0.03 secs

``` r
dbGetQuery(dbConn, "select * from Seqs") %>% head(4)
```

    ##   row_names identifier
    ## 1         1      HPV18
    ## 2         2      HPV18
    ## 3         3      HPV16
    ## 4         4      HPV18
    ##                                                                                      description
    ## 1 gi|944543704|gb|KT365847.1| Human papillomavirus isolate HPV18-14 L1 protein gene, partial cds
    ## 2 gi|944543703|gb|KT365846.1| Human papillomavirus isolate HPV18-13 L1 protein gene, partial cds
    ## 3  gi|944543701|gb|KT365845.1| Human papillomavirus isolate HPV16-5 L1 protein gene, partial cds
    ## 4 gi|944543699|gb|KT365844.1| Human papillomavirus isolate HPV18-12 L1 protein gene, partial cds
    ##          GI  SVAR ISO GENE
    ## 1 944543704 HPV18  14   L1
    ## 2 944543703 HPV18  13   L1
    ## 3 944543701 HPV16   5   L1
    ## 4 944543699 HPV18  12   L1

Design HPV16-specific F/W primers using `DECIPHER` utilities.

``` r
tiles.L1 <- TileSeqs(
  dbFile = dbConn,
  minLength = 18,
  maxLength = 29,
  minCoverage = 0.8
)
```

    ## ================================================================================
    ## 
    ## Time difference of 7.29 secs

``` r
print(
  tiles.L1[,c(6,11)] %>% sample_n(6)
)
```

    ##   misprime                   target_site
    ## 1     TRUE AAGCAGACACAGTTATGTATTTTGGGCTG
    ## 2    FALSE GGTGCTACTGGTCATCCATATTTTAATAG
    ## 3    FALSE CCCATTGGTGAGCATTGGGCCAAGGGCAC
    ## 4    FALSE       TCATGCTGGTCAGCCTGGTGAGT
    ## 5    FALSE AAATGCTAGTGCTTATGCAGCAAATGCAG
    ## 6    FALSE GAAAGACTCCAACGACGCAGAGAAACACA

``` r
oligos.L1 <- DesignPrimers(
  tiles = tiles.L1,
  identifier = 'HPV16',
  worstScore = -1E3,
  maxPermutations = 5,
  minGroupCoverage = 0.85,
  minCoverage = 0.85,
  minLength = 20,
  maxLength = 28
)
```

    ## 
    ## HPV16 (7 candidate primers):
    ## ================================================================================
    ## 
    ## Time difference of 5.9 secs

``` r
head(oligos.L1)
```

    ##     identifier start_forward start_reverse start_aligned_forward
    ## 201      HPV16           249           275                   248
    ## 231      HPV16           284           302                   278
    ## 232      HPV16           287           304                   279
    ## 236      HPV16           288           306                   283
    ## 237      HPV16           288           306                   284
    ## 238      HPV16           288           307                   285
    ##     start_aligned_reverse permutations_forward permutations_reverse
    ## 201                   276                    5                    5
    ## 231                   306                    5                    5
    ## 232                   307                    5                    5
    ## 236                   311                    5                    5
    ## 237                   312                    5                    5
    ## 238                   313                    5                    5
    ##     score_forward score_reverse          forward_primer.1
    ## 201  -6.27037....  -0.00038.... TTGGTTTTCCTGACACCTCATTTTA
    ## 231  -0.00333....  -0.00082....      ACACAGCGGCTGGTTTGGGC
    ## 232  -0.00932....  -0.00024....      CACAGCGGCTGGTTTGGGCC
    ## 236  -0.20587....  -0.00023....      GCGGCTGGTTTGGGCCTGTG
    ## 237  -0.01501....  -0.00047....      CGGCTGGTTTGGGCCTGTGT
    ## 238  -6.07531....  -0.00115....      GGCTGGTTTGGGCCTGTGTA
    ##                forward_primer.2         forward_primer.3
    ## 201 AAAGTGAAGTTCCACTGGATATTTGTA TGCCTTAGTGGACCCTACAGTATA
    ## 231     GAAAAGGAAAGACTAGTGTGGGC     CGCGAGCGATTAGTGTGGAA
    ## 232       AAGGAAAGACTAGTGTGGGCC    CGCGAGCGATTAGTGTGGAAA
    ## 236        AAGACTAGTGTGGGCCTGTG GCGAGCGATTAGTGTGGAAATTTA
    ## 237        AGACTAGTGTGGGCCTGTGC     AGACTGGTTTGGGCCTGTAG
    ## 238        GACTAGTGTGGGCCTGTGCA     GACTGGTTTGGGCCTGTAGA
    ##                 forward_primer.4           forward_primer.5
    ## 201 AATTTGCTTTAGCAGATATGTCAGTCTA  TCGCCTTGGTAGATATGAATGTCTA
    ## 231         CATGAGCGTTTAGTGTGGCG       AAGGAAAGACTGGTTTGGGC
    ## 232         ATGAGCGTTTAGTGTGGCGT       AGGAAAGACTGGTTTGGGCC
    ## 236         AAGACTGGTTTGGGCCTGTA      AGCGTTTAGTGTGGCGTTTAC
    ## 237        GCGTTTAGTGTGGCGTTTACG  GCGAGCGATTAGTGTGGAAATTTAC
    ## 238       GCGTTTAGTGTGGCGTTTACGT GCGAGCGATTAGTGTGGAAATTTACA
    ##            reverse_primer.1             reverse_primer.2
    ## 201 GAGGTGTCAGGAAAACCAAACTT TACAAATATCCAGTGGAACTTCACTTTT
    ## 231    AGCCGCTGTGTATCTGGATT    ACACTAGTCTTTCCTTTTCTGGGTT
    ## 232    CAGCCGCTGTGTATCTGGAT    CACACTAGTCTTTCCTTTTCTGGGT
    ## 236    AAACCAGCCGCTGTGTATCT     GCCCACACTAGTCTTTCCTTTTCT
    ## 237    CAAACCAGCCGCTGTGTATC      GCCCACACTAGTCTTTCCTTTTC
    ## 238    CCAAACCAGCCGCTGTGTAT      GGCCCACACTAGTCTTTCCTTTT
    ##             reverse_primer.3             reverse_primer.4
    ## 201 TGTAGGGTCCACTAAGGCAAATTT AGACTGACATATCTGCTAAAGCAAATTT
    ## 231     AATCGCTCGCGGTCTGGGTT         AAACGCTCATGGTCAGAGTT
    ## 232     TAATCGCTCGCGGTCTGGGT      CACTAAACGCTCATGGTCAGAGT
    ## 236     ACACTAATCGCTCGCGGTCT     GCCCAAACCAGTCTTTCCTTTTCA
    ## 237  GCCCAAACCAGTCTTTCCTTTTC        CCACACTAAACGCTCATGGTC
    ## 238   GCCCAAACCAGTCTTTCCTTTT        GCCACACTAAACGCTCATGGT
    ##                reverse_primer.5 forward_efficiency.1 forward_efficiency.2
    ## 201 GACATTCATATCTACCAAGGCGAATCT            0.8339538            0.8284785
    ## 231   AAACCAGTCTTTCCTTTTCAGGATT            0.9976087            0.8209337
    ## 232  CCAAACCAGTCTTTCCTTTTCAGGAT            0.9987960            0.8996088
    ## 236      CCACACTAAACGCTCATGGTCA            0.9987474            0.9127847
    ## 237        CACACTAATCGCTCGCGGTC            0.9957214            0.9361697
    ## 238        CCACACTAATCGCTCGCGGT            0.9859330            0.9767059
    ##     forward_efficiency.3 forward_efficiency.4 forward_efficiency.5
    ## 201            0.8905282            0.8600553            0.8531185
    ## 231            0.8767872            0.8416041            0.8159230
    ## 232            0.9008893            0.8601513            0.8485258
    ## 236            0.9127556            0.9438210            0.8215438
    ## 237            0.8015350            0.8500537            0.9127556
    ## 238            0.9105089            0.9458488            0.9484353
    ##     reverse_efficiency.1 reverse_efficiency.2 reverse_efficiency.3
    ## 201            0.9283584            0.8692461            0.9066815
    ## 231            0.8941688            0.8974904            0.9989823
    ## 232            0.8768771            0.8541983            0.9957801
    ## 236            0.9203204            0.8333086            0.9243034
    ## 237            0.9137364            0.8209337            0.9233699
    ## 238            0.9500386            0.9443670            0.8935778
    ##     reverse_efficiency.4 reverse_efficiency.5 forward_coverage.1
    ## 201            0.8127476            0.9368589          0.4444444
    ## 231            0.8755932            0.8085530          0.4444444
    ## 232            0.8794141            0.8211595          0.4444444
    ## 236            0.9405149            0.8775017          0.4444444
    ## 237            0.8451844            0.9707554          0.4444444
    ## 238            0.9113278            0.9687265          0.4444444
    ##     forward_coverage.2 forward_coverage.3 forward_coverage.4 forward_coverage.5
    ## 201          0.1111111          0.1111111          0.1111111          0.1111111
    ## 231          0.1111111          0.1111111          0.1111111          0.1111111
    ## 232          0.1111111          0.1111111          0.1111111          0.1111111
    ## 236          0.1111111          0.1111111          0.1111111          0.1111111
    ## 237          0.1111111          0.1111111          0.1111111          0.1111111
    ## 238          0.1111111          0.1111111          0.1111111          0.1111111
    ##     reverse_coverage.1 reverse_coverage.2 reverse_coverage.3 reverse_coverage.4
    ## 201          0.4444444          0.1111111          0.1111111          0.1111111
    ## 231          0.4444444          0.1111111          0.1111111          0.1111111
    ## 232          0.4444444          0.1111111          0.1111111          0.1111111
    ## 236          0.4444444          0.1111111          0.1111111          0.1111111
    ## 237          0.4444444          0.1111111          0.1111111          0.1111111
    ## 238          0.4444444          0.1111111          0.1111111          0.1111111
    ##     reverse_coverage.5                              mismatches_forward
    ## 201          0.1111111                                                
    ## 231          0.1111111 HPV18 (0.0986%) HPV11 (0.0393%) HPV17 (0.192%) 
    ## 232          0.1111111                  HPV18 (0.437%) HPV17 (0.495%) 
    ## 236          0.1111111                   HPV18 (20.5%) HPV17 (0.114%) 
    ## 237          0.1111111                   HPV18 (1.26%) HPV17 (0.238%) 
    ## 238          0.1111111                                                
    ##                  mismatches_reverse
    ## 201                HPV17 (0.0296%) 
    ## 231 HPV4 (0.0656%) HPV17 (0.0125%) 
    ## 232                 HPV4 (0.0189%) 
    ## 236                HPV17 (0.0207%) 
    ## 237                HPV17 (0.0426%) 
    ## 238                  HPV17 (0.11%)

Plot hybridization curves for the reverse compliment of each primer pair
(5) in each primer set (7). Honestly, I don’t really understand the
dimension of the dataframe output by `DECIPHER::DesignPimers`. If I
remember correctly, it’s a 7x17 matrix where each cell holds a list of 5
oligos. Unfortunately, I can’t print the matrix to stdout, which makes
debugging slow, so I’ll just deal with this issue later I guess.

``` r
plotMeltCurves <- function(temps, effs) {
  plot(
    temps, effs[,1], ylim=c(0,1),
    ylab="Hybridization Efficiency",
    xlab=expression(paste("Temperature (", degree, "C)", sep="")),
    type="l", lwd=2, col="Blue", main="Denaturation Plot"
  )
  lines(temps, effs[,2], col="Red", lwd=2)
  abline(h=0.5, lty=2, lwd=2, col="Orange")
  abline(v=64, lty=2, lwd=2, col="Green")
  legend(
    "topright",
    legend = c(
      "Forward Primer",
      "Reverse Primer",
      "50% Efficiency",
      "Annealing Temperature"
    ),
    col = c("Blue", "Red", "Orange", "Green"),
    lwd = c(2, 2, 2, 2), lty = c(1, 1, 2, 2)
  )
}

meltCurves <- 
  function(primers, target=reverseComplement(DNAStringSet(primers)), temps=60:75, P=4e-7, ions=.225, doPlot=TRUE, ...) {
    fxn <- function(temp) CalculateEfficiencyPCR(primers, target, temp, P=P, ions=ions, ...)
    effs <- matrix(unlist(lapply(temps, fxn)), ncol=2, byrow=TRUE)
    if (doPlot) plotMeltCurves(temps, effs)
    list(temps = temps, effs = effs)
  }

nsets <- nrow(oligos.L1)
npair <- ncol(oligos.L1)

for (i in 1:nsets)
  for (j in 1:npair)
    tryCatch(
      {
        c(
          oligos.L1$forward_primer[i,j],
          oligos.L1$reverse_primer[i,j]
        ) %>% meltCurves
      },
      error = function(e) NULL
        # # ¿ --- It's supposed to be 7o x 17p --- ?
        # warning(sprintf('%s [iteration i=%s j=%s]', e, i, j))
    )
```

![](HPV-Serovars_files/figure-gfm/melt_curves-1.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-2.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-3.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-4.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-5.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-6.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-7.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-8.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-9.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-10.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-11.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-12.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-13.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-14.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-15.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-16.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-17.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-18.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-19.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-20.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-21.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-22.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-23.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-24.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-25.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-26.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-27.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-28.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-29.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-30.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-31.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-32.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-33.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-34.png)<!-- -->![](HPV-Serovars_files/figure-gfm/melt_curves-35.png)<!-- -->

Some more primer designs (RFLP, sequencing, etc.):

``` r
TYPE <- 'sequence'
MIN_LENGTH <- 15
MAX_LENGTH <- 25
MIN_SIZE <- 60
MAX_SIZE <- 100
RESOLUTION <- 3
LEVELS <- 2
ENZYMES <- NULL

DesignSignatures(
  dbConn,
  type=TYPE,
  minLength=MIN_LENGTH,
  maxLength=MAX_LENGTH,
  minProductSize=MIN_SIZE,
  maxProductSize=MAX_SIZE,
  resolution=RESOLUTION,
  levels=LEVELS,
  enzymes=ENZYMES
)
```

    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.41 secs
    ## 
    ## Designing primer sequences based on the group 'HPV18':
    ## ================================================================================
    ## 
    ## Time difference of 59.77 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 13.69 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 2.8 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.55 secs

    ##                forward_primer            reverse_primer score coverage products
    ## 1      CCAACGACGCAGAGAAACACAA  ATGTCTTGCAATGTTGCCTTAGGT     0      0.2        0
    ## 2     CCAACGACGCAGAGAAACACAAG  ATGTCTTGCAATGTTGCCTTAGGT     0      0.2        0
    ## 3       GCTTTGAGGATCCAACACGGC    GCAGTGAAGTGTTCAGTTCCGT     0      0.2        0
    ## 4       GCTTTGAGGATCCAACACGGC    GTGTTCAGTTCCGTGCACAGAT     0      0.2        0
    ## 5       GCTTTGAGGATCCAACACGGC   GCAGTGAAGTGTTCAGTTCCGTG     0      0.2        0
    ## 6       GCTTTGAGGATCCAACACGGC     AGTGTTCAGTTCCGTGCACAG     0      0.2        0
    ## 7       GCTTTGAGGATCCAACACGGC    AGTGAAGTGTTCAGTTCCGTGC     0      0.2        0
    ## 8       GCTTTGAGGATCCAACACGGC     GTGTTCAGTTCCGTGCACAGA     0      0.2        0
    ## 9       GCTTTGAGGATCCAACACGGC    GTGAAGTGTTCAGTTCCGTGCA     0      0.2        0
    ## 10      GCTTTGAGGATCCAACACGGC    GAAGTGTTCAGTTCCGTGCACA     0      0.2        0
    ## 11          TCGTGCTGCAACCGAGC  ATGTCTTGCAATGTTGCCTTAGGT     0      0.2        0
    ## 12          TCGTGCTGCAACCGAGC GTTGCCTTAGGTCCATGCATACTTA     0      0.2        0
    ## 13      CAGGAAMGACTCCAACGACGC  ATGTCTTGCAATGTTGCCTTAGGT     0      0.2        0
    ## 14      CAGGAAMGACTCCAACGACGC GTTGCCTTAGGTCCATGCATACTTA     0      0.2        0
    ## 15        GGAAMGACTCCAACGACGC  ATGTCTTGCAATGTTGCCTTAGGT     0      0.2        0
    ## 16        GGAAMGACTCCAACGACGC GTTGCCTTAGGTCCATGCATACTTA     0      0.2        0
    ## 17         CCCCCCCGCCAACTACTA    GTGCAGCATCCTTTTGACAGGT     0      0.2        0
    ## 18         CCCCCCCGCCAACTACTA  CCTTATTTTCAGCYGGTGCAGCAT     0      0.2        0
    ## 19         CCCCCCCGCCAACTACTA        TTTCAGCYGGTGCAGCAT     0      0.2        0
    ## 20         CCCCCCCGCCAACTACTA GTGCAGCATCCTTTTGACAGGTAAT     0      0.2        0
    ## 21         CCCCCCCGCCAACTACTA CCTTATTTTCAGCYGGTGCAGCATC     0      0.2        0
    ## 22         CCCCCCCGCCAACTACTA   TCCTTATTTTCAGCYGGTGCAGC     0      0.2        0
    ## 23         CCCCCCCGCCAACTACTA    CCTTATTTTCAGCYGGTGCAGC     0      0.2        0
    ## 24         CCCCCCCGCCAACTACTA   GTGCAGCATCCTTTTGACAGGTA     0      0.2        0
    ## 25         CCCCCCCGCCAACTACTA   CCTTATTTTCAGCYGGTGCAGCA     0      0.2        0
    ## 26         CCCCCCCGCCAACTACTA         TTTCAGCYGGTGCAGCA     0      0.2        0
    ## 27         CCCCCCCGCCAACTACTA  GTGCAGCATCCTTTTGACAGGTAA     0      0.2        0
    ## 28   TCTCCTGTACCTGGGCAATATGAT TCYTCAACATGTCTGCTATACTGCT     0      0.2        0
    ## 29   TCTCCTGTACCTGGGCAATATGAT  CCTCAACATGTCTGCTATACTGCC     0      0.2        0
    ## 30  GCATGCTGCATGCCATAAATGTATA TCYAATGTGTCTCCATACACAGAGT     0      0.2        0
    ## 31      TGCCGCCACGTCTAATGTTTC GGCACAGCCCAAAATACATAACTGT     0      0.2        0
    ## 32      TGCCGCCACGTCTAATGTTTC  GGGCACAGCCCAAAATACATAACT     0      0.2        0
    ## 33      TGCCGCCACGTCTAATGTTTC GGGCACAGCCCAAAATACATAACTG     0      0.2        0
    ## 34      TGCCGCCACGTCTAATGTTTC       CAATAGCAGGGGCACAGCC     0      0.2        0
    ## 35      TGCCGCCACGTCTAATGTTTC  GGGGCACAGCCCAAAATACATAAC     0      0.2        0
    ## 36      TGCCGCCACGTCTAATGTTTC          GCAGGGKYACAGCCCA     0      0.2        0
    ## 37  CAACGACGCAGAGAAACACAAGTAT  ATGTCTTGCAATGTTGCCTTAGGT     0      0.2        0
    ## 38  GATGTGAGAAACRCACCACAATACT      GCTTGTAGGGTCGCCGTGTT     0      0.2        0
    ## 39  GATGTGAGAAACRCACCACAATACT       GCTTGTAGGGTCGCCGTGT     0      0.2        0
    ## 40  GATGTGAGAAACRCACCACAATACT  CACAGATCAGGTAGCTTGTAGGGT     0      0.2        0
    ## 41  GATGTGAGAAACRCACCACAATACT     AGTTCCGTGCACAGATCAGGT     0      0.2        0
    ## 42  GATGTGAGAAACRCACCACAATACT    GTGTTCAGTTCCGTGCACAGAT     0      0.2        0
    ## 43  GATGTGAGAAACRCACCACAATACT      CTTGTAGGGTCGCCGTGTTG     0      0.2        0
    ## 44  GATGTGAGAAACRCACCACAATACT        GCTTGTAGGGTCGCCGTG     0      0.2        0
    ## 45  GATGTGAGAAACRCACCACAATACT        GTAGGGTCGCCGTGTTGG     0      0.2        0
    ## 46  GATGTGAGAAACRCACCACAATACT  GCACAGATCAGGTAGCTTGTAGGG     0      0.2        0
    ## 47  GATGTGAGAAACRCACCACAATACT     CAGTTCCGTGCACAGATCAGG     0      0.2        0
    ## 48  GATGTGAGAAACRCACCACAATACT   AGATCAGGTAGCTTGTAGGGTCG     0      0.2        0
    ## 49  GATGTGAGAAACRCACCACAATACT       GTAGCTTGTAGGGTCGCCG     0      0.2        0
    ## 50  GATGTGAGAAACRCACCACAATACT    GTTCCGTGCACAGATCAGGTAG     0      0.2        0
    ## 51  GATGTGAGAAACRCACCACAATACT   GTTCAGTTCCGTGCACAGATCAG     0      0.2        0
    ## 52  GATGTGAGAAACRCACCACAATACT  ACAGATCAGGTAGCTTGTAGGGTC     0      0.2        0
    ## 53  GATGTGAGAAACRCACCACAATACT    TGTTCAGTTCCGTGCACAGATC     0      0.2        0
    ## 54  GATGTGAGAAACRCACCACAATACT     TCAGGTAGCTTGTAGGGTCGC     0      0.2        0
    ## 55  GATGTGAGAAACRCACCACAATACT      AGGTAGCTTGTAGGGTCGCC     0      0.2        0
    ## 56  GATGTGAGAAACRCACCACAATACT    AGTTCCGTGCACAGATCAGGTA     0      0.2        0
    ## 57  GATGTGAGAAACRCACCACAATACT     GTGTTCAGTTCCGTGCACAGA     0      0.2        0
    ## 58  GATGTGAGAAACRCACCACAATACT    GTTCAGTTCCGTGCACAGATCA     0      0.2        0
    ## 59   TCTGAGGACGTTAGGGACAATGTG GGCACAGCCCAAAATACATAACTGT     0      0.2        0
    ## 60   TCTGAGGACGTTAGGGACAATGTG  GGGCACAGCCCAAAATACATAACT     0      0.2        0
    ## 61   TCTGAGGACGTTAGGGACAATGTG GGGCACAGCCCAAAATACATAACTG     0      0.2        0
    ## 62   TCTGAGGACGTTAGGGACAATGTG       CAATAGCAGGGGCACAGCC     0      0.2        0
    ## 63   TCTGAGGACGTTAGGGACAATGTG  GGGGCACAGCCCAAAATACATAAC     0      0.2        0
    ## 64   TCTGAGGACGTTAGGGACAATGTG          GCAGGGKYACAGCCCA     0      0.2        0
    ## 65    CAACGACGCAGAGAAACACAAGT  ATGTCTTGCAATGTTGCCTTAGGT     0      0.2        0
    ## 66   GCAGGTACTATGGGTGACACTGTG GCACGCATACCTGTGCCTTTAATAT     0      0.2        0
    ## 67   GCAGGTACTATGGGTGACACTGTG   GCACGCATACCTGTGCCTTTAAT     0      0.2        0
    ## 68   GCAGGTACTATGGGTGACACTGTG    GCACGCATACCTGTGCCTTTAA     0      0.2        0
    ## 69    CGGCTGGTTTTATGTACAGGCTA   CCATATCCGACCCTGTGTCTGTT     0      0.2        0
    ## 70    CGGCTGGTTTTATGTACAGGCTA    CCATATCCGACCCTGTGTCTGT     0      0.2        0
    ## 71    CGGCTGGTTTTATGTACAGGCTA    ACCATATCCGACCCTGTGTCTG     0      0.2        0
    ## 72    CGGCTGGTTTTATGTACAGGCTA  CGACCCTGTGTCTGTTGCATTTTC     0      0.2        0
    ## 73  ACTCTGTGTATGGAGACACATTGGA    GCACCGCAGGCACCTTATTAAT     0      0.2        0
    ## 74  ACTCTGTGTATGGAGACACATTGGA     GCACCGCAGGCACCTTATTAA     0      0.2        0
    ## 75  ACTCTGTGTATGGAGACACATTGGA GCACCGCAGGCACCTTATTAATAAA     0      0.2        0
    ## 76    AGATGTGAGAAACRCACCACAAT      GCTTGTAGGGTCGCCGTGTT     0      0.2        0
    ## 77    AGATGTGAGAAACRCACCACAAT       GCTTGTAGGGTCGCCGTGT     0      0.2        0
    ## 78    AGATGTGAGAAACRCACCACAAT  CACAGATCAGGTAGCTTGTAGGGT     0      0.2        0
    ## 79    AGATGTGAGAAACRCACCACAAT     AGTTCCGTGCACAGATCAGGT     0      0.2        0
    ## 80    AGATGTGAGAAACRCACCACAAT    GTGTTCAGTTCCGTGCACAGAT     0      0.2        0
    ## 81    AGATGTGAGAAACRCACCACAAT      CTTGTAGGGTCGCCGTGTTG     0      0.2        0
    ## 82    AGATGTGAGAAACRCACCACAAT        GCTTGTAGGGTCGCCGTG     0      0.2        0
    ## 83    AGATGTGAGAAACRCACCACAAT        GTAGGGTCGCCGTGTTGG     0      0.2        0
    ## 84    AGATGTGAGAAACRCACCACAAT  GCACAGATCAGGTAGCTTGTAGGG     0      0.2        0
    ## 85    AGATGTGAGAAACRCACCACAAT     CAGTTCCGTGCACAGATCAGG     0      0.2        0
    ## 86    AGATGTGAGAAACRCACCACAAT   AGATCAGGTAGCTTGTAGGGTCG     0      0.2        0
    ## 87    AGATGTGAGAAACRCACCACAAT       GTAGCTTGTAGGGTCGCCG     0      0.2        0
    ## 88    AGATGTGAGAAACRCACCACAAT    GTTCCGTGCACAGATCAGGTAG     0      0.2        0
    ## 89    AGATGTGAGAAACRCACCACAAT   GTTCAGTTCCGTGCACAGATCAG     0      0.2        0
    ## 90    AGATGTGAGAAACRCACCACAAT     AGTGTTCAGTTCCGTGCACAG     0      0.2        0
    ## 91    AGATGTGAGAAACRCACCACAAT  ACAGATCAGGTAGCTTGTAGGGTC     0      0.2        0
    ## 92    AGATGTGAGAAACRCACCACAAT    TGTTCAGTTCCGTGCACAGATC     0      0.2        0
    ## 93    AGATGTGAGAAACRCACCACAAT    AGTGAAGTGTTCAGTTCCGTGC     0      0.2        0
    ## 94    AGATGTGAGAAACRCACCACAAT     TCAGGTAGCTTGTAGGGTCGC     0      0.2        0
    ## 95    AGATGTGAGAAACRCACCACAAT      AGGTAGCTTGTAGGGTCGCC     0      0.2        0
    ## 96    AGATGTGAGAAACRCACCACAAT    TGAAGTGTTCAGTTCCGTGCAC     0      0.2        0
    ## 97    AGATGTGAGAAACRCACCACAAT    AGTTCCGTGCACAGATCAGGTA     0      0.2        0
    ## 98    AGATGTGAGAAACRCACCACAAT     GTGTTCAGTTCCGTGCACAGA     0      0.2        0
    ## 99    AGATGTGAGAAACRCACCACAAT    GTTCAGTTCCGTGCACAGATCA     0      0.2        0
    ## 100   AGATGTGAGAAACRCACCACAAT    GTGAAGTGTTCAGTTCCGTGCA     0      0.2        0
    ##     similar_signatures        missing_signatures
    ## 1                      HPV16, HPV11, HPV4, HPV17
    ## 2                      HPV16, HPV11, HPV4, HPV17
    ## 3                      HPV16, HPV11, HPV4, HPV17
    ## 4                      HPV16, HPV11, HPV4, HPV17
    ## 5                      HPV16, HPV11, HPV4, HPV17
    ## 6                      HPV16, HPV11, HPV4, HPV17
    ## 7                      HPV16, HPV11, HPV4, HPV17
    ## 8                      HPV16, HPV11, HPV4, HPV17
    ## 9                      HPV16, HPV11, HPV4, HPV17
    ## 10                     HPV16, HPV11, HPV4, HPV17
    ## 11                     HPV16, HPV11, HPV4, HPV17
    ## 12                     HPV16, HPV11, HPV4, HPV17
    ## 13                     HPV16, HPV11, HPV4, HPV17
    ## 14                     HPV16, HPV11, HPV4, HPV17
    ## 15                     HPV16, HPV11, HPV4, HPV17
    ## 16                     HPV16, HPV11, HPV4, HPV17
    ## 17                     HPV16, HPV11, HPV4, HPV17
    ## 18                     HPV16, HPV11, HPV4, HPV17
    ## 19                     HPV16, HPV11, HPV4, HPV17
    ## 20                     HPV16, HPV11, HPV4, HPV17
    ## 21                     HPV16, HPV11, HPV4, HPV17
    ## 22                     HPV16, HPV11, HPV4, HPV17
    ## 23                     HPV16, HPV11, HPV4, HPV17
    ## 24                     HPV16, HPV11, HPV4, HPV17
    ## 25                     HPV16, HPV11, HPV4, HPV17
    ## 26                     HPV16, HPV11, HPV4, HPV17
    ## 27                     HPV16, HPV11, HPV4, HPV17
    ## 28                     HPV16, HPV11, HPV4, HPV17
    ## 29                     HPV16, HPV11, HPV4, HPV17
    ## 30                     HPV16, HPV11, HPV4, HPV17
    ## 31                     HPV16, HPV11, HPV4, HPV17
    ## 32                     HPV16, HPV11, HPV4, HPV17
    ## 33                     HPV16, HPV11, HPV4, HPV17
    ## 34                     HPV16, HPV11, HPV4, HPV17
    ## 35                     HPV16, HPV11, HPV4, HPV17
    ## 36                     HPV16, HPV11, HPV4, HPV17
    ## 37                     HPV16, HPV11, HPV4, HPV17
    ## 38                     HPV16, HPV11, HPV4, HPV17
    ## 39                     HPV16, HPV11, HPV4, HPV17
    ## 40                     HPV16, HPV11, HPV4, HPV17
    ## 41                     HPV16, HPV11, HPV4, HPV17
    ## 42                     HPV16, HPV11, HPV4, HPV17
    ## 43                     HPV16, HPV11, HPV4, HPV17
    ## 44                     HPV16, HPV11, HPV4, HPV17
    ## 45                     HPV16, HPV11, HPV4, HPV17
    ## 46                     HPV16, HPV11, HPV4, HPV17
    ## 47                     HPV16, HPV11, HPV4, HPV17
    ## 48                     HPV16, HPV11, HPV4, HPV17
    ## 49                     HPV16, HPV11, HPV4, HPV17
    ## 50                     HPV16, HPV11, HPV4, HPV17
    ## 51                     HPV16, HPV11, HPV4, HPV17
    ## 52                     HPV16, HPV11, HPV4, HPV17
    ## 53                     HPV16, HPV11, HPV4, HPV17
    ## 54                     HPV16, HPV11, HPV4, HPV17
    ## 55                     HPV16, HPV11, HPV4, HPV17
    ## 56                     HPV16, HPV11, HPV4, HPV17
    ## 57                     HPV16, HPV11, HPV4, HPV17
    ## 58                     HPV16, HPV11, HPV4, HPV17
    ## 59                     HPV16, HPV11, HPV4, HPV17
    ## 60                     HPV16, HPV11, HPV4, HPV17
    ## 61                     HPV16, HPV11, HPV4, HPV17
    ## 62                     HPV16, HPV11, HPV4, HPV17
    ## 63                     HPV16, HPV11, HPV4, HPV17
    ## 64                     HPV16, HPV11, HPV4, HPV17
    ## 65                     HPV16, HPV11, HPV4, HPV17
    ## 66                     HPV16, HPV11, HPV4, HPV17
    ## 67                     HPV16, HPV11, HPV4, HPV17
    ## 68                     HPV16, HPV11, HPV4, HPV17
    ## 69                     HPV16, HPV11, HPV4, HPV17
    ## 70                     HPV16, HPV11, HPV4, HPV17
    ## 71                     HPV16, HPV11, HPV4, HPV17
    ## 72                     HPV16, HPV11, HPV4, HPV17
    ## 73                     HPV16, HPV11, HPV4, HPV17
    ## 74                     HPV16, HPV11, HPV4, HPV17
    ## 75                     HPV16, HPV11, HPV4, HPV17
    ## 76                     HPV16, HPV11, HPV4, HPV17
    ## 77                     HPV16, HPV11, HPV4, HPV17
    ## 78                     HPV16, HPV11, HPV4, HPV17
    ## 79                     HPV16, HPV11, HPV4, HPV17
    ## 80                     HPV16, HPV11, HPV4, HPV17
    ## 81                     HPV16, HPV11, HPV4, HPV17
    ## 82                     HPV16, HPV11, HPV4, HPV17
    ## 83                     HPV16, HPV11, HPV4, HPV17
    ## 84                     HPV16, HPV11, HPV4, HPV17
    ## 85                     HPV16, HPV11, HPV4, HPV17
    ## 86                     HPV16, HPV11, HPV4, HPV17
    ## 87                     HPV16, HPV11, HPV4, HPV17
    ## 88                     HPV16, HPV11, HPV4, HPV17
    ## 89                     HPV16, HPV11, HPV4, HPV17
    ## 90                     HPV16, HPV11, HPV4, HPV17
    ## 91                     HPV16, HPV11, HPV4, HPV17
    ## 92                     HPV16, HPV11, HPV4, HPV17
    ## 93                     HPV16, HPV11, HPV4, HPV17
    ## 94                     HPV16, HPV11, HPV4, HPV17
    ## 95                     HPV16, HPV11, HPV4, HPV17
    ## 96                     HPV16, HPV11, HPV4, HPV17
    ## 97                     HPV16, HPV11, HPV4, HPV17
    ## 98                     HPV16, HPV11, HPV4, HPV17
    ## 99                     HPV16, HPV11, HPV4, HPV17
    ## 100                    HPV16, HPV11, HPV4, HPV17

``` r
# TYPE <- 'length'
# LEVELS <- 2
# MIN_SIZE <- 200
# MAX_SIZE <- 1400
# RESOLUTION <- c(
#   seq(200, 700, 3),
#   seq(705, 1000, 5),
#   seq(1010, 1400, 10)
# )
# 
# TYPE <- 'melt'
# MIN_SIZE <- 55
# MAX_SIZE <- 400
# 
# data(RESTRICTION_ENZYMES) 
# LEVELS <- 10
# RESOLUTION <- seq(75, 100, 0.25)
# ENZYMES <- RESTRICTION_ENZYMES["MslI"]
# 
# primers <- DesignSignatures(
#   dbConn,
#   type=TYPE,
#   minProductSize=MIN_SIZE,
#   maxProductSize=MAX_SIZE,
#   resolution=RESOLUTION,
#   levels=LEVELS,
#   enzymes=ENZYMES
# )
```

## Anal Isolates

<https://www.ncbi.nlm.nih.gov/nuccore>

``` sql
 "Human papillomavirus"[Primary Organism]
 AND viruses[filter] 
 NOT Polyamides[All Fields] 
 NOT Method[All Fields] 
 NOT Patent[All Fields]
 AND Anal[All Fields]
```

## Cutaneous Isolates

<https://www.ncbi.nlm.nih.gov/nuccore>

``` sql
"Human papillomavirus"[Primary Organism]
NOT Polyamides[All Fields] 
NOT Method[All Fields] 
NOT Patent[All Fields]
AND Anal[All Fields]
AND "Complete Genome"[All Fields]
```

``` r
tryCatch({dbDisconnect(dbconn)}, error=warning)
tryCatch({dbDisconnect(dbConn)}, error=warning)
```

<div id="refs" class="references">

<div id="ref-Rous_1911">

Rous, Peyton. 1911. “A Sarcoma of the Fowl Transmissible by an Agent
Separable from the Tumor Cells.” *The Journal of Experimental Medicine*
13 (4): 397–411. <https://doi.org/10.1084/jem.13.4.397>.

</div>

</div>
