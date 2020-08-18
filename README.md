Designing PCR Diagnostics that Discriminate HPV Serotypes by DNA
Sequence
================
chebuu

<!-- --- -->

<!-- ... -->

<!-- bibliography: citations.bib -->

<!-- includes: -->

<!--   in_header:  -->

<!--     - preamble.tex -->

<!--     - common.tex -->

<!-- ... -->

<!-- --- -->

# Introduction

The first infectious cause of cancer was identified in 1911 by Peyton
Rous who demonstrated the transmissibility of a tumor in fowl through
injection of sub 0.2 micron filtrate of the tumor .

#### Installing this package

``` r
devtools::install_github('Chebuu/HPV-Serovars')
library(HPVSerovars)
```

(Not run) This vignette is packaged and accessible via:
`browseVignettes('HPVSerovars')`.

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

### Install Bioconductor packages

<https://www.bioconductor.org/>

``` r
# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = T))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")

# Open Biostrings vignettes
library(Biostrings)
browseVignettes("Biostrings")

# Open DECIPHER docs and vignettes
library(DECIPHER)
help("DECIPHER")
browseVignettes("DECIPHER")
```

### Load packages

``` r
library(dplyr)
library(tidyr)
library(stringr)

library(Biostrings)
library(DECIPHER)

library(HPVSerovars)
```

#### Make an alignment utility:

``` r
doAlignment <- function(xset) {
  useqs <- unique(xset)
  uidxs <- match(xset, useqs)
  aseqs <- AlignSeqs(useqs, verbose=F) 
  aseqs[uidxs]
}
```

# Wiki Pheno Data

I copied the tables from [Wiki:
Human\_papillomavirus\_infection](https://en.wikipedia.org/wiki/Human_papillomavirus_infection)
and included them in the R package assoc. w/ this document
`github.com/Chebuu/HPV-Serovars`.

``` r
data(HPV.Sequelae)
data(HPV.GenitalCancerRisk)
```

The Wiki data contains phenotypic categorizations of several HPV
strains. It’s useful for testing.

``` r
HPV <- c(HPV.Sequelae, HPV.GenitalCancerRisk)
p.1 <- unique(unlist(HPV))

knitr::kable(
  sprintf('HPV-%s', sort(p.1)), col.names = '',
  caption = sprintf(
    'Unique HPV serotypes (N=%s) (Total=%s): ', 
    length(p.1), length(HPV)
  )
)
```

|        |
| :----- |
| HPV-1  |
| HPV-2  |
| HPV-3  |
| HPV-4  |
| HPV-6  |
| HPV-7  |
| HPV-10 |
| HPV-11 |
| HPV-13 |
| HPV-16 |
| HPV-18 |
| HPV-22 |
| HPV-26 |
| HPV-28 |
| HPV-31 |
| HPV-32 |
| HPV-33 |
| HPV-35 |
| HPV-39 |
| HPV-42 |
| HPV-44 |
| HPV-45 |
| HPV-51 |
| HPV-52 |
| HPV-53 |
| HPV-56 |
| HPV-58 |
| HPV-59 |
| HPV-60 |
| HPV-63 |
| HPV-66 |
| HPV-68 |
| HPV-73 |
| HPV-82 |

Unique HPV serotypes (N=34) (Total=13):

``` r
hist(unlist(HPV), max(unlist(HPV)), col='blue'); head(HPV, 4)
```

![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpcsFCFZ/preview-772336fa426.dir/Manuscript_files/figure-gfm/wikiprint-1.png)<!-- -->

    ## $`Common warts`
    ## [1]  2  7 22
    ## 
    ## $`Plantar warts`
    ## [1]  1  2  4 63
    ## 
    ## $`Flat warts`
    ## [1]  3 10 28
    ## 
    ## $`Anogenital warts`
    ## [1]  6 11 42 44

# NCBI Data

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
"Human papillomavirus"[Primary Organism] 
AND viruses[filter] 
NOT Polyamides[All Fields] 
NOT Method[All Fields] 
NOT Patent[All Fields] 
AND Oral[All Fields] 
```

### L1 gene

A collection of *l1* DNA sequences (avg. \~600bp) from various oral
isolates of HPV (extracted from results of GenBank query above) are
aligned and loaded in SQLite (RAM) along with annotations from GenBank
that describe the serotype of each isolate.

#### Design HPV16-specific F/R primers using `DECIPHER` utilities.

``` r
data(Oral.L1.seqs)
data(Oral.L1.vars)

dbConn <- dbConnect(RSQLite::SQLite(), ':memory:')

Seqs2DB(Oral.L1.seqs, 'XStringSet', dbConn, '')
```

    ## Adding 31 sequences to the database.
    ## 
    ## 31 total sequences in table Seqs.
    ## Time difference of 0.19 secs

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

    ## Warning: Factors converted to character

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

    ## Warning: Factors converted to character

    ## Expression:
    ## update Seqs set identifier = :identifier where row_names = :row_names

    ## Warning: Factors converted to character

    ## Added to table Seqs:  "GI" and "SVAR" and "ISO" and "GENE" and "identifier".
    ## 
    ## Time difference of 0.06 secs

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

``` r
tiles.L1 <- TileSeqs(
  dbFile = dbConn,
  minLength = 18,
  maxLength = 29,
  minCoverage = 0.8
)
```

    ## Warning in TileSeqs(dbFile = dbConn, minLength = 18, maxLength = 29, minCoverage
    ## = 0.8): Skipped due to multiple width sequences: HPV18

    ## Warning in TileSeqs(dbFile = dbConn, minLength = 18, maxLength = 29, minCoverage
    ## = 0.8): Skipped due to multiple width sequences: HPV16

    ## ================================================================

    ## Warning in TileSeqs(dbFile = dbConn, minLength = 18, maxLength = 29, minCoverage
    ## = 0.8): Skipped due to multiple width sequences: HPV17

    ## 
    ## 
    ## Time difference of 2.76 secs

``` r
print(
  tiles.L1[,c(6,11)] %>% sample_n(6)
)
```

    ##   misprime                   target_site
    ## 1    FALSE TCCATATTTTAATAGATACTATGATACTG
    ## 2    FALSE GTGTTGCCAGATCCTAACAAGTTTGCATT
    ## 3    FALSE CTGTTTGACCCCACTACACAGCGTTTAGT
    ## 4    FALSE CAATTGAATTAATGCATACAGTGATACAG
    ## 5    FALSE AGACTCCTTGCTGTGGGACATCCATATTA
    ## 6    FALSE TGGAGGTAGGCAGGGGTCAACCTTTAGGC

``` r
oligos.L1 <- DesignPrimers(
  tiles = tiles.L1,
  identifier = 'HPV4',
  worstScore = -1E3,
  maxPermutations = 5,
  minGroupCoverage = 0.85,
  minCoverage = 0.85,
  minLength = 20,
  maxLength = 28
)
```

    ## 
    ## HPV4 (324 candidate primers):
    ## ================================================================================
    ## 
    ## Time difference of 5.9 secs

``` r
head(oligos.L1)
```

    ##   identifier start_forward start_reverse start_aligned_forward
    ## 3       HPV4            37            52                    30
    ## 4       HPV4            37            53                    31
    ## 5       HPV4            38            54                    32
    ## 6       HPV4            38            55                    33
    ## 7       HPV4            39            55                    34
    ## 8       HPV4            40            55                    35
    ##   start_aligned_reverse permutations_forward permutations_reverse score_forward
    ## 3                    58                    1                    1             0
    ## 4                    59                    1                    1             0
    ## 5                    60                    1                    1             0
    ## 6                    61                    1                    1             0
    ## 7                    62                    1                    1             0
    ## 8                    63                    1                    1             0
    ##   score_reverse         forward_primer.1 forward_primer.2 forward_primer.3
    ## 3             0   TGGTCCCAAAGGTTTCAGCAAA             <NA>             <NA>
    ## 4             0  TGGTCCCAAAGGTTTCAGCAAAT             <NA>             <NA>
    ## 5             0  GGTCCCAAAGGTTTCAGCAAATC             <NA>             <NA>
    ## 6             0 GGTCCCAAAGGTTTCAGCAAATCA             <NA>             <NA>
    ## 7             0 GTCCCAAAGGTTTCAGCAAATCAG             <NA>             <NA>
    ## 8             0 TCCCAAAGGTTTCAGCAAATCAGT             <NA>             <NA>
    ##   forward_primer.4 forward_primer.5        reverse_primer.1 reverse_primer.2
    ## 3             <NA>             <NA> GAAACCTTTGGGACCACAATTTT             <NA>
    ## 4             <NA>             <NA> TGAAACCTTTGGGACCACAATTT             <NA>
    ## 5             <NA>             <NA> CTGAAACCTTTGGGACCACAATT             <NA>
    ## 6             <NA>             <NA> GCTGAAACCTTTGGGACCACAAT             <NA>
    ## 7             <NA>             <NA>  GCTGAAACCTTTGGGACCACAA             <NA>
    ## 8             <NA>             <NA>   GCTGAAACCTTTGGGACCACA             <NA>
    ##   reverse_primer.3 reverse_primer.4 reverse_primer.5 forward_efficiency.1
    ## 3             <NA>             <NA>             <NA>            0.8695543
    ## 4             <NA>             <NA>             <NA>            0.8997148
    ## 5             <NA>             <NA>             <NA>            0.9171716
    ## 6             <NA>             <NA>             <NA>            0.9356067
    ## 7             <NA>             <NA>             <NA>            0.8303409
    ## 8             <NA>             <NA>             <NA>            0.9161453
    ##   forward_efficiency.2 forward_efficiency.3 forward_efficiency.4
    ## 3                   NA                   NA                   NA
    ## 4                   NA                   NA                   NA
    ## 5                   NA                   NA                   NA
    ## 6                   NA                   NA                   NA
    ## 7                   NA                   NA                   NA
    ## 8                   NA                   NA                   NA
    ##   forward_efficiency.5 reverse_efficiency.1 reverse_efficiency.2
    ## 3                   NA            0.8433773                   NA
    ## 4                   NA            0.8369640                   NA
    ## 5                   NA            0.8697235                   NA
    ## 6                   NA            0.9333216                   NA
    ## 7                   NA            0.9122828                   NA
    ## 8                   NA            0.8831295                   NA
    ##   reverse_efficiency.3 reverse_efficiency.4 reverse_efficiency.5
    ## 3                   NA                   NA                   NA
    ## 4                   NA                   NA                   NA
    ## 5                   NA                   NA                   NA
    ## 6                   NA                   NA                   NA
    ## 7                   NA                   NA                   NA
    ## 8                   NA                   NA                   NA
    ##   forward_coverage.1 forward_coverage.2 forward_coverage.3 forward_coverage.4
    ## 3                  1                 NA                 NA                 NA
    ## 4                  1                 NA                 NA                 NA
    ## 5                  1                 NA                 NA                 NA
    ## 6                  1                 NA                 NA                 NA
    ## 7                  1                 NA                 NA                 NA
    ## 8                  1                 NA                 NA                 NA
    ##   forward_coverage.5 reverse_coverage.1 reverse_coverage.2 reverse_coverage.3
    ## 3                 NA                  1                 NA                 NA
    ## 4                 NA                  1                 NA                 NA
    ## 5                 NA                  1                 NA                 NA
    ## 6                 NA                  1                 NA                 NA
    ## 7                 NA                  1                 NA                 NA
    ## 8                 NA                  1                 NA                 NA
    ##   reverse_coverage.4 reverse_coverage.5 mismatches_forward mismatches_reverse
    ## 3                 NA                 NA                                      
    ## 4                 NA                 NA                                      
    ## 5                 NA                 NA                                      
    ## 6                 NA                 NA                                      
    ## 7                 NA                 NA                                      
    ## 8                 NA                 NA

Honestly, I don’t really understand the dimensions of the dataframe
output by `DECIPHER::DesignPimers`. If I remember correctly, it’s a 7x17
matrix where each cell holds a list of 5 oligos. Unfortunately, I can’t
print the matrix to stdout, which makes debugging slow, so I’ll just
deal with this issue later I guess…

Plot hybridization curves for the reverse compliment of each primer pair
(5) in each primer set (7).

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

More primer designs (RFLP, sequencing, etc.):

``` r
TYPE <- 'sequence'
MIN_LENGTH <- 15
MAX_LENGTH <- 25
MIN_SIZE <- 60
MAX_SIZE <- 100
LEVELS <- 2
RESOLUTION <- 3
ENZYMES <- NULL

DesignSignatures(
  dbConn,
  type = TYPE,
  minLength = MIN_LENGTH,
  maxLength = MAX_LENGTH,
  minProductSize = MIN_SIZE,
  maxProductSize = MAX_SIZE,
  resolution = RESOLUTION,
  levels = LEVELS,
  enzymes = ENZYMES
)
```

    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.86 secs
    ## 
    ## Designing primer sequences based on the group 'HPV18':
    ## ================================================================================
    ## 
    ## Time difference of 85.35 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 17.73 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 3.51 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.01 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 2.01 secs

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
data(Oral.L1.seqs)
data(Oral.L1.vars)

dbConn <- dbConnect(SQLite(), ':memory:')

print(Oral.L1.seqs)

Seqs2DB(Oral.L1.seqs, 'XStringSet', dbConn, '')
Add2DB(Oral.L1.vars %>% mutate(identifier = SVAR), dbConn)

dbGetQuery(dbConn, "select * from Seqs") %>% head(4)

data(RESTRICTION_ENZYMES)

FOCUS_ID <- 'HPV11'

TYPE <- 'melt'
LEVELS <- 4
MIN_LENGTH <- 15
MAX_LENGTH <- 25
MIN_SIZE <- 60
MAX_SIZE <- 200
RESOLUTION <- seq(75, 300, 15)

RFLP.Oal.L1 <- lapply(
  rep(sample(1:3, 3), 3),  # Random RE combinations
  function(i)
    DesignSignatures(
      dbConn,
      type = TYPE,
      identifier = '',
      focusID = FOCUS_ID,
      enzymes = RESTRICTION_ENZYMES[i],
      minProductSize = MIN_SIZE,
      maxProductSize = MAX_SIZE,
      resolution = RESOLUTION, 
      levels = LEVELS
    )
)

head(RFLP.Oal.L1)

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
```

#### RFLP optimization by grid search

``` r
FOCUS_ID <- 'HPV18'

params <- expand.grid(
  enzymes = lapply(1:3, function(x) 
    RESTRICTION_ENZYMES[
      sample(1:length(RESTRICTION_ENZYMES), 3)
    ]
  ),
  minProductSize = seq(80, 150, length.out = 3),
  maxProductSize = seq(500, 5000, length.out = 3),
  resolution = seq(1, 2, length.out = 3), 
  levels = seq(5, 20, length.out = 3)
) %>% mutate(type = rep('melt', 243))

gridsearch <- apply(params, 1, function(...) 
    list(
      arg = as.list(...), 
      out = do.call(
        DesignSignatures, 
        append(
          list(
            dbFile = dbConn, 
            tblName = 'Seqs', 
            identifier = '', 
            focusID = FOCUS_ID
          ), 
          as.list(...)
        )
      )
    )
) %>% bind_rows

objective <- function(x) {
  print(x)
}
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

# PaVE Data

Reference genomes for all serovars listed in `data(HPV.Sequelae)`and
`data(HPV.GenitalCancerRisk)`were selected from among the results of the
following [PaVEPaVE](https://pave.niaid.nih.gov) query:

    https://pave.niaid.nih.gov/#search/search_database/kw
    ?dbNamespace=Genomes
    &includeNR=true
    &refCloneOnly=false
    &sort=Locus_ID
    &sortType=true
    &page=600&start=1
    &text=Human&showTable=1
    &

``` r
# data(PaVE.Wikiv)
data(HPV.REF.PaVE.seqs)
data(HPV.REF.PaVE.vars)

testVars <- lapply(grep('[1-9]*', HPV.REF.PaVE.vars, value = T), function(x) {
  sapply(HPV.Sequelae, function(y) any(x %in% y))
})

dbConn <- dbConnect(SQLite(), ':memory:')

DECIPHER::Seqs2DB(HPV.REF.PaVE.seqs, 'XStringSet', dbConn, '')
```

    ## Adding 218 sequences to the database.
    ## 
    ## 218 total sequences in table Seqs.
    ## Time difference of 0.11 secs

``` r
DECIPHER::Add2DB(
  data.frame(
    # identifier = testVars
    identifier = lapply(HPV.REF.PaVE.vars, function(x) grep('[1-9]', x, perl = T))
  ), 
  dbConn, verbose = F
)
```

#### Genital cancer risk groups (PaVE REF genomes)

HPV reference genomes from PaVE categorized by risk of genital (but not
anal) cancer.

``` r
library(HPVSerovars)

data(gA) # High risk
data(gB) # Med risk
data(gC) # Some risk

dbDisconnect(dbConn)

dbConn <- dbConnect(RSQLite::SQLite(), ':memory:')

Seqs2DB(
  c(gA,gB,gC), 
  type = 'XStringSet', 
  dbFile = dbConn, 
  identifier = ''
)
```

    ## Adding 18 sequences to the database.
    ## 
    ## 18 total sequences in table Seqs.
    ## Time difference of 0.06 secs

``` r
Add2DB(
  data.frame(
    identifier = c(
      rep('A', length(gA)), 
      rep('B', length(gB)), 
      rep('C', length(gC))
    )
  ), 
  dbConn
)
```

    ## Expression:
    ## update Seqs set identifier = :identifier where row_names = :row_names

    ## Warning: Factors converted to character

    ## Added to table Seqs:  "identifier".
    ## 
    ## Time difference of 0.01 secs

``` r
gABC.synt <- FindSynteny(dbConn)
```

    ## ================================================================================
    ## 
    ## Time difference of 0.84 secs

``` r
summary(gABC.synt)
```

    ##  A.Length  A.Class  A.Mode  B.Length  B.Class  B.Mode 
    ##    4      -none-   numeric  16376    -none-   numeric 
    ##  320      -none-   numeric      8    -none-   numeric 
    ##  240      -none-   numeric    480    -none-   numeric 
    ##  C.Length  C.Class  C.Mode 
    ##  12936    -none-   numeric 
    ##  26560    -none-   numeric 
    ##      6    -none-   numeric

``` r
pairs(gABC.synt, boxBlocks=TRUE)
```

![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpcsFCFZ/preview-772336fa426.dir/Manuscript_files/figure-gfm/risk_groups-1.png)<!-- -->

``` r
unlist(AlignSynteny(gABC.synt, dbConn)) %>% head(4)
```

    ## ================================================================================
    ## 
    ## Time difference of 18.3 secs

    ## $`A/B`
    ## DNAStringSetList of length 32
    ## [["1"]] A=AGGGCGTAACCGAAATCGGTTGAACCGAAACCGGTTAGTATAAAAGCAGACATTTTATGCACCAAAA...
    ## [["2"]] A=ATGGCGCGCTTTGACGATCCAAAGCAACGACCCTACAAGCTACCAGATTTGTGCACAGAATTGAATA...
    ## [["3"]] A=ATGGCGCGCTTTGAGGATCCAACACGGCGACCCTACAAGCTACCTGATCTGTGCACGGAACTGAACA...
    ## [["4"]] A=ACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGCAACAGTTACTGCGACGTGAGGTAT...
    ## [["5"]] A=AAAAAAGGGTGTAACCGAAAACGGTTGCAACCAAAAACGGTGCATATAAAA-GCTTTGTGGAAAAGT...
    ## [["6"]] A=AAAAAAGGGAGTAACCGAAAACGGTCGGGACCGAAAACGGTGTATATAAAA-GATGTGAGAAACACA...
    ## [["7"]] A=AAAAAAGTAGGGAGTGACCGAAAGTGGT-GAACCGAAAACGGTTGGTATATAA-------AGCACAT...
    ## [["8"]] A=ATTTACAGATTTAACAATAGTATATAGGGACGACACACCACACGGAGTGTGTACAAAATGTTTAAGA...
    ## [["9"]] A=AAAAAAGGGAGTAACCGAAAACGGTCGGGACCGAAAACGGTGTATATAAAAG---ATGTGAGAAACA...
    ## [["10"]] A=ATGAACTAAGATTGAATTGTGTCTACTGCAAAGGTCAGTTAACAGAAACAGAGGTATTAGATTTTG...
    ## ...
    ## <22 more elements>
    ## 
    ## $`A/C`
    ## DNAStringSetList of length 24
    ## [["1"]] A=AACGACCCTACAAGCTACCAGATTTGTGCACAGAATTGAATACATCACTACAAGACGTATCTATTGC...
    ## [["2"]] A=CGACCCTACAAGCTACCTGATCTGTGCACGGAACTGAACACTTCACTGCAAGACATAGAAATAACCT...
    ## [["3"]] A=AACCGAAACCGGTTAGTATAAAAGCAGACATTTTATGCACCAAAAGAGAACTGCAATGTTTCAGGAC...
    ## [["4"]] A=TAACAATTATACTACATAAAAAAGGGTGTAACCGAAAACGGTTGCAACCAAAAACGGTGCATATAAA...
    ## [["5"]] A=AAAAGGGAGTAACCGAAAACGGTCGGGACCGAAAACGGTGTATATAAAAG----ATGTGAGAAACAC...
    ## [["6"]] A=AGGGAGTAACCGAAA-ACGGTCGGGACCGAAAACGGTGTATATAAAAG--ATGTGAGAAACACACCA...
    ## [["7"]] A=AGGGCGTAACCGAAATCGGTTGA--ACCGAAACCGGTTAGTATAAAAGCAGACATTTTATGCACCAA...
    ## [["8"]] A=ACCAAAAACGGTGCATATAAAAGCTTTGTGGA------AAAGTGCATTACAGGATGGCGCGCT--TT...
    ## [["9"]] A=AACACTGCAAGAAATTGTATTGCATTTGGAACCTCAGAATGAATTAGATCCTGTTGACCTGTTGTGT...
    ## [["10"]] A=GAACTTCCAGACAGCGGGTATGGCAATACTGAAG------TGGAAACGCAGCAGATGGTACAGGTA...
    ## ...
    ## <14 more elements>
    ## 
    ## $`B/C`
    ## DNAStringSetList of length 48
    ## [["1"]] B=GAAAGTTTCAATCATACTTTTATATATTGGGAGTGACCGAAAAGGGTTTAAGACCGAAAACGGTACA...
    ## [["2"]] B=ATGGCGCGATTTCACAATCCTGCAGAACGGCCATACAAATTGCCAGACCTGTGCACAACGCTGGACA...
    ## [["3"]] B=ACAATTATCTTGTAAAAACTAGGGTGTAACCGAAAAGGGTTATGACCGAAAACGGTGCATATAAAAG...
    ## [["4"]] B=CAATCATACTTTTATATATTGGGAGTGACCGAAAAGGGTTTAAGACCGAAAACGGTACATATAAAAG...
    ## [["5"]] B=AGGGAAAGACCACGAACGCTGCATGAATTATGTGAAGCTTTGAACGTTTCTATGCACAATATACAGG...
    ## [["6"]] B=AACGACCATACAAACTGCCTGATTTGAGCACAACATTGAATATTCCTCTGCATGATATTCGCATCAA...
    ## [["7"]] B=AGGGTGTAACCGAAAACGGTCT-GACCGAAACCGGTGCATATATAA--AGCAGACATTTTTTGGTAG...
    ## [["8"]] B=AGGGTGTAACCGAAAGCGGTT-CAACCGAAAACGGTGCATATATAAAGCAAACATTTTGCAGTAAGG...
    ## [["9"]] B=AATTAAGATATTATAGAGATTCCGTGTATGGAGAAACATTAGAGGCTGAAACCAAGACACCGTTACA...
    ## [["10"]] B=AGGGAGTGACCGAAAACGGTCG-TACCGAAAACGGTTGCCATAAAAGCAGAAGTGCACAAAAAAGC...
    ## ...
    ## <38 more elements>

### Using OpenPrimer

#### Installing OpenPrimer Dependencies

From the [openPimer
vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/openPrimeR/inst/doc/openPrimeR_vignette.html):

openPrimeR requires external programs for some features, particularly
for computing the physicochemical properties of primers. Please make
sure you have the following tools installed on your system such that
they are in your system’s path:

[MELTING](http://www.ebi.ac.uk/biomodels/tools/melting/) (\>= 5.1.1):
For melting temperature computations.
[ViennaRNA](http://www.tbi.univie.ac.at/RNA/) (\>= 2.2.4): For secondary
structure prediction.
[OligoArrayAux](http://unafold.rna.albany.edu/OligoArrayAux.php) (\>=
3.8): For primer efficiency computations as performed by DECIPHER.
[MAFFT](http://mafft.cbrc.jp/alignment/software/) (\>= 7.305): For
computing multiple sequence alignments. [Pandoc](http://pandoc.org/)
(\>= 1.19.1): For creating PDF reports.

If you would like to be able to access the immunoglobulin repository
IMGT from the openPrimeR Shiny app, you should additionally fulfill the
following dependencies:

[PhantomJS](http://phantomjs.org/) (\>= 2.1): For headless website
calls. [Python](http://www.python.org/) (\>=2.7.9) and the
[selenium](http://selenium-python.readthedocs.io/) (\>=3.0.1) module:
For data extraction scripts. openPrimeR will automatically check for all
dependencies and inform you about any missing dependencies when the
package is attached.

Note that the tool is still functional if there are missing external
programs. However, we recommend that all dependencies are fulfilled to
guarantee the best user experience.

``` r
library(openPrimeR)
```
