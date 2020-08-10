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
```

    ## Skipping install of 'HPVSerovars' from a github remote, the SHA1 (17418f04) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
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
package is attached:

``` r
# library(openPrimeR)
```

Note that the tool is still functional if there are missing external
programs. However, we recommend that all dependencies are fulfilled to
guarantee the best user experience.

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

#### Make alignment utils:

``` r
doAlignment <- function(xset) {
  useqs <- unique(xset)
  uidxs <- match(xset, useqs)
  aseqs <- AlignSeqs(useqs, verbose=F) 
  aseqs[uidxs]
}
```

#### Quick-ref

Some references to remind me of packages, methods, etc.:

    # CRAN
    # https://rdrr.io/cran/dpcR/f/vignettes/overview.Rmd
    
    # BIOC
    # https://bioconductor.org/packages/release/bioc/html/twoddpcr.html
    # https://bioconductor.org/packages/release/bioc/html/ddCt.html
    # https://bioconductor.org/packages/3.5/bioc/html/mBPCR.html
    
    # Other
    # https://doi.org/10.1007/s13402-017-0331-y

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

![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/wikiprint-1.png)<!-- -->

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
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(stringr)

library(Biostrings)
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 3.6.3

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
library(DECIPHER)
```

    ## Loading required package: RSQLite

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
    ## Time difference of 0.31 secs

``` r
DECIPHER::Add2DB(
  data.frame(
    # identifier = testVars
    identifier = lapply(HPV.REF.PaVE.vars, function(x) grep('[1-9]', x, perl = T))
  ), 
  dbConn, verbose = F
)

# # DEL
# library(Biostrings)
# PaVE.Wikiv<-Biostrings::readDNAStringSet('../inst/extdata/HPV.REFS.PaVE.wikivars.fasta')
# usethis::use_data(PaVE.Wikiv)
#
# https://www.youtube.com/watch?v=zu8Z9-kebKQ&list=PLjiXAZO27elDHGlQwfd06r7coiFtpPkvI&index=5
#
# # END
```

# NCBI Data

``` r
library(dplyr)
library(tidyr)
library(stringr)

library(Biostrings)
library(DECIPHER)

library(HPVSerovars)
```

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

``` r
# library(HPVSerovars)
# 
# data(Oral.L1.seqs)
# data(Oral.L1.vars)

load('../data/Oral.L1.seqs.rda')
load('../data/Oral.L1.vars.rda')

dbConn <- dbConnect(SQLite(), ':memory:')

# Oral.L1.seqs <- readDNAStringSet('../inst/extdata/HPV.oral.L1.fasta') %>% doAlignment
# Oral.L1.vars <- read.csv('../inst/extdata/HPV.oral.L1.csv', stringsAsFactors = F)

Seqs2DB(Oral.L1.seqs, 'XStringSet', dbConn, '')
```

    ## Adding 31 sequences to the database.
    ## 
    ## 31 total sequences in table Seqs.
    ## Time difference of 0.02 secs

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
    ## Time difference of 0.04 secs

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

Design HPV16-specific F/R primers using `DECIPHER` utilities.

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
    ## Time difference of 2.92 secs

``` r
print(
  tiles.L1[,c(6,11)] %>% sample_n(6)
)
```

    ##   misprime                   target_site
    ## 1    FALSE TCCTAAGCAGGGTGACAATGTGGATAGCA
    ## 2    FALSE CTAAGCAGGTACAAATGCTTATTGTGGGC
    ## 3    FALSE AGATTAGTTTGGCGACTGGTAGGACTTCA
    ## 4    FALSE GCCAGCAGTTCTAGACTCCTTGCTGTGGG
    ## 5    FALSE GGGCATCCATTGCTAAACAAATATGATGA
    ## 6    FALSE GCAGTTCTAGACTCCTTGCTGTGGGACAT

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
    ## Time difference of 7.1 secs

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

![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-1.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-2.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-3.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-4.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-5.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-6.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-7.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-8.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-9.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-10.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-11.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-12.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-13.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-14.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-15.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-16.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-17.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-18.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-19.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-20.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-21.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-22.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-23.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-24.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-25.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-26.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-27.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-28.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-29.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-30.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-31.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-32.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-33.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-34.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-35.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-36.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-37.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-38.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-39.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-40.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-41.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-42.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-43.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-44.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-45.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-46.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-47.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-48.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-49.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-50.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-51.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-52.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-53.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-54.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-55.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-56.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-57.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-58.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-59.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-60.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-61.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-62.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-63.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-64.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-65.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-66.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-67.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-68.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-69.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-70.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-71.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-72.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-73.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-74.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-75.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-76.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-77.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-78.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-79.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-80.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-81.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-82.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-83.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-84.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-85.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-86.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-87.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-88.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-89.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-90.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-91.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-92.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-93.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-94.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-95.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-96.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-97.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-98.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-99.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-100.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-101.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-102.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-103.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-104.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-105.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-106.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-107.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-108.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-109.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-110.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-111.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-112.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-113.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-114.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-115.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-116.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-117.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-118.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-119.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-120.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-121.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-122.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-123.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-124.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-125.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-126.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-127.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-128.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-129.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-130.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-131.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-132.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-133.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-134.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-135.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-136.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-137.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-138.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-139.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-140.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-141.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-142.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-143.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-144.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-145.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-146.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-147.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-148.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-149.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-150.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-151.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-152.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-153.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-154.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-155.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-156.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-157.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-158.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-159.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-160.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-161.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-162.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-163.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-164.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-165.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-166.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-167.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-168.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-169.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-170.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-171.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-172.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-173.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-174.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-175.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-176.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-177.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-178.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-179.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-180.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-181.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-182.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-183.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-184.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-185.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-186.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-187.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-188.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-189.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-190.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-191.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-192.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-193.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-194.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-195.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-196.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-197.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-198.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-199.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-200.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-201.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-202.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-203.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-204.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-205.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-206.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-207.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-208.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-209.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-210.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-211.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-212.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-213.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-214.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-215.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-216.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-217.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-218.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-219.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-220.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-221.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-222.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-223.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-224.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-225.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-226.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-227.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-228.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-229.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-230.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-231.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-232.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-233.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-234.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-235.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-236.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-237.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-238.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-239.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-240.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-241.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-242.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-243.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-244.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-245.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-246.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-247.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-248.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-249.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-250.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-251.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-252.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-253.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-254.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-255.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-256.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-257.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-258.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-259.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-260.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-261.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-262.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-263.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-264.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-265.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-266.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-267.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-268.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-269.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-270.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-271.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-272.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-273.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-274.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-275.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-276.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-277.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-278.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-279.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-280.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-281.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-282.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-283.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-284.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-285.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-286.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-287.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-288.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-289.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-290.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-291.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-292.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-293.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-294.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-295.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-296.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-297.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-298.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-299.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-300.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-301.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-302.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-303.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-304.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-305.png)<!-- -->![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpGWPYOU/preview-15858733632bf.dir/Manuscript_files/figure-gfm/melt_curves-306.png)<!-- -->

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
    ## Time difference of 0.48 secs
    ## 
    ## Designing primer sequences based on the group 'HPV18':
    ## ================================================================================
    ## 
    ## Time difference of 67.94 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 14.97 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 3.04 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.01 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.71 secs

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
```

    ##   A DNAStringSet instance of length 31
    ##      width seq                                              names               
    ##  [1]   313 GGTATACATGGTACACATTATTA...ATATAGAGTATTTAGGGTGCAG gi|944543704|gb|K...
    ##  [2]   312 TATACCGCATGCTGCATGCCATA...CACAAGTATAATATTAAGTATG gi|944543703|gb|K...
    ##  [3]   316 ATTAGAGTTAATAAACACAGTTA...TTTAGCCAGTTCAAATTATTTT gi|944543701|gb|K...
    ##  [4]   321 TATCACAGGGCGATTGCCCCCCT...ACAGGTATGCGTGCTTCACCTG gi|944543699|gb|K...
    ##  [5]   313 CGGATGAATATGTTGCACGCACA...CTTTATTAAATAAATTGGATGA gi|944543697|gb|K...
    ##  ...   ... ...
    ## [27]   437 TTGAAGTACGGAATGGGCAGGAC...GTTAATAAACTCCATTATCGAA gi|440573440|gb|K...
    ## [28]   437 TTAATGTATATGAAGGTAACGAG...GTTAGTTAACAGCACTATTCAG gi|440573438|gb|K...
    ## [29]   440 ATGATGTCAGATCTAGTGATGGC...ATTAAAAAATACAGTAATTGAA gi|440573436|gb|K...
    ## [30]   437 TTGATGTTAGAGACACTGTGGAT...ATTAGTAAATACTGTAATTGAA gi|440573434|gb|K...
    ## [31]   440 TTGAAGTGCGTAATGGTTTGGGT...GTTAATAAATACAGTAATTGAA gi|440573432|gb|K...

``` r
Seqs2DB(Oral.L1.seqs, 'XStringSet', dbConn, '')
```

    ## Adding 31 sequences to the database.
    ## 
    ## 31 total sequences in table Seqs.
    ## Time difference of 0.03 secs

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
    ## Time difference of 0.02 secs

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
data(RESTRICTION_ENZYMES)

FOCUS_ID <- 'HPV11'

TYPE <- 'melt'
LEVELS <- 4
MIN_LENGTH <- 15
MAX_LENGTH <- 25
MIN_SIZE <- 60
MAX_SIZE <- 200
RESOLUTION <- seq(75, 300, 15)

lapply(
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
) %>% print
```

    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.39 secs
    ## 
    ## Designing primer sequences based on the group 'HPV11':
    ## ================================================================================
    ## 
    ## Time difference of 18.67 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 3.85 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 2.67 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.01 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.71 secs
    ## 
    ## Finding the best restriction enzyme:
    ## ================================================================================
    ## 
    ## Time difference of 2.39 secs
    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.36 secs
    ## 
    ## Designing primer sequences based on the group 'HPV11':
    ## ================================================================================
    ## 
    ## Time difference of 18.14 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 4.1 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 2.31 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.02 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.69 secs
    ## 
    ## Finding the best restriction enzyme:
    ## ================================================================================
    ## 
    ## Time difference of 2.35 secs
    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.36 secs
    ## 
    ## Designing primer sequences based on the group 'HPV11':
    ## ================================================================================
    ## 
    ## Time difference of 18.32 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 4.72 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 2.84 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.02 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.86 secs
    ## 
    ## Finding the best restriction enzyme:
    ## ================================================================================
    ## 
    ## Time difference of 2.15 secs
    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.4 secs
    ## 
    ## Designing primer sequences based on the group 'HPV11':
    ## ================================================================================
    ## 
    ## Time difference of 21.91 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 4.9 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 3.5 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.01 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.98 secs
    ## 
    ## Finding the best restriction enzyme:
    ## ================================================================================
    ## 
    ## Time difference of 2.97 secs
    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.58 secs
    ## 
    ## Designing primer sequences based on the group 'HPV11':
    ## ================================================================================
    ## 
    ## Time difference of 22.65 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 4.64 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 3.17 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.01 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.98 secs
    ## 
    ## Finding the best restriction enzyme:
    ## ================================================================================
    ## 
    ## Time difference of 2.36 secs
    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.4 secs
    ## 
    ## Designing primer sequences based on the group 'HPV11':
    ## ================================================================================
    ## 
    ## Time difference of 22.19 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 3.7 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 2.35 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.01 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.7 secs
    ## 
    ## Finding the best restriction enzyme:
    ## ================================================================================
    ## 
    ## Time difference of 2.28 secs
    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.34 secs
    ## 
    ## Designing primer sequences based on the group 'HPV11':
    ## ================================================================================
    ## 
    ## Time difference of 17.81 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 4.07 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 2.61 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.01 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.68 secs
    ## 
    ## Finding the best restriction enzyme:
    ## ================================================================================
    ## 
    ## Time difference of 1.94 secs
    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.35 secs
    ## 
    ## Designing primer sequences based on the group 'HPV11':
    ## ================================================================================
    ## 
    ## Time difference of 17.57 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 3.96 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 2.73 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.01 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.66 secs
    ## 
    ## Finding the best restriction enzyme:
    ## ================================================================================
    ## 
    ## Time difference of 1.93 secs
    ## Tallying 8-mers for 5 groups:
    ## ================================================================================
    ## 
    ## Time difference of 0.36 secs
    ## 
    ## Designing primer sequences based on the group 'HPV11':
    ## ================================================================================
    ## 
    ## Time difference of 17.87 secs
    ## 
    ## Selecting the most common primer sequences:
    ## ================================================================================
    ## 
    ## Time difference of 3.94 secs
    ## 
    ## Determining PCR products from each group:
    ## ================================================================================
    ## 
    ## Time difference of 2.66 secs
    ## 
    ## Scoring primer pair combinations:
    ## ================================================================================
    ## 
    ## Time difference of 0.01 secs
    ## 
    ## Choosing optimal forward and reverse pairs:
    ## ================================================================================
    ## 
    ## Time difference of 1.71 secs
    ## 
    ## Finding the best restriction enzyme:
    ## ================================================================================
    ## 
    ## Time difference of 1.66 secs
    ## [[1]]
    ##            forward_primer             reverse_primer score coverage products
    ## 1   GTTTGACCCCACTACACAGCG       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 2    TGACCCCACTACACAGCGTT       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 3   GTTTGACCCCACTACACAGCG   CACTAACACCAACGCCTAAAGGTT     0      0.2        0
    ## 4   GTTTGACCCCACTACACAGCG    CACTAACACCAACGCCTAAAGGT     0      0.2        0
    ## 5   GTTTGACCCCACTACACAGCG      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 6   GTTTGACCCCACTACACAGCG       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 7   GTTTGACCCCACTACACAGCG        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 8   GTTTGACCCCACTACACAGCG    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 9   GTTTGACCCCACTACACAGCG   ACTAACACCAACGCCTAAAGGTTG     0      0.2        0
    ## 10  GTTTGACCCCACTACACAGCG       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 11  GTTTGACCCCACTACACAGCG    CCACTAACRCCAACRCCTAAAGG     0      0.2        0
    ## 12  GTTTGACCCCACTACACAGCG      GGATGCCCACTAACACCAACG     0      0.2        0
    ## 13  GTTTGACCCCACTACACAGCG    CCCACTAACACCAACGCCTAAAG     0      0.2        0
    ## 14  GTTTGACCCCACTACACAGCG       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 15  GTTTGACCCCACTACACAGCG       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 16  GTTTGACCCCACTACACAGCG        TGCCCACTAACACCAACGC     0      0.2        0
    ## 17  GTTTGACCCCACTACACAGCG         TGACCCCTGCCTACCTCC     0      0.2        0
    ## 18  GTTTGACCCCACTACACAGCG      CAACGCCTAAAGGTTGACCCC     0      0.2        0
    ## 19  GTTTGACCCCACTACACAGCG      CCAACGCCTAAAGGTTGACCC     0      0.2        0
    ## 20  GTTTGACCCCACTACACAGCG        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 21  GTTTGACCCCACTACACAGCG      ACCAACGCCTAAAGGTTGACC     0      0.2        0
    ## 22  GTTTGACCCCACTACACAGCG     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 23  GTTTGACCCCACTACACAGCG    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 24  GTTTGACCCCACTACACAGCG     ACACCAACGCCTAAAGGTTGAC     0      0.2        0
    ## 25  GTTTGACCCCACTACACAGCG     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 26  GTTTGACCCCACTACACAGCG  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 27  GTTTGACCCCACTACACAGCG        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 28  GTTTGACCCCACTACACAGCG     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 29  GTTTGACCCCACTACACAGCG     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 30  GTTTGACCCCACTACACAGCG      GCCCACTAACACCAACGCCTA     0      0.2        0
    ## 31  GTTTGACCCCACTACACAGCG   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 32  GTTTGACCCCACTACACAGCG   CTAACACCAACGCCTAAAGGTTGA     0      0.2        0
    ## 33  GTTTGACCCCACTACACAGCG         CACCYCTRCCTACCTCCA     0      0.2        0
    ## 34  GTTTGACCCCACTACACAGCG     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 35  GTTTGACCCCACTACACAGCG  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 36  GTTTGACCCCACTACACAGCG     GCCCACTAACACCAACGCCTAA     0      0.2        0
    ## 37  GTTTGACCCCACTACACAGCG  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 38  GTTTGACCCCACTACACAGCG        GACCCCTGCCTACCTCCAA     0      0.2        0
    ## 39  GTTTGACCCCACTACACAGCG    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 40  GTTTGACCCCACTACACAGCG     CCCACTAACACCAACGCCTAAA     0      0.2        0
    ## 41   TGACCCCACTACACAGCGTT      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 42   TGACCCCACTACACAGCGTT       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 43   TGACCCCACTACACAGCGTT        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 44   TGACCCCACTACACAGCGTT    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 45   TGACCCCACTACACAGCGTT       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 46   TGACCCCACTACACAGCGTT       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 47   TGACCCCACTACACAGCGTT       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 48   TGACCCCACTACACAGCGTT        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 49   TGACCCCACTACACAGCGTT     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 50   TGACCCCACTACACAGCGTT    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 51   TGACCCCACTACACAGCGTT     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 52   TGACCCCACTACACAGCGTT  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 53   TGACCCCACTACACAGCGTT     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 54   TGACCCCACTACACAGCGTT     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 55   TGACCCCACTACACAGCGTT   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 56   TGACCCCACTACACAGCGTT     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 57   TGACCCCACTACACAGCGTT  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 58   TGACCCCACTACACAGCGTT  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 59   TGACCCCACTACACAGCGTT    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 60  TCCTTGCTGTGGGACATCCAT      CGCCCATACTAAACGCTGTGT     0      0.2        0
    ## 61  TCCTTGCTGTGGGACATCCAT      CACGCCCATACTAAACGCTGT     0      0.2        0
    ## 62  TCCTTGCTGTGGGACATCCAT        TGCCTACCTCCAACCCTGT     0      0.2        0
    ## 63  TCCTTGCTGTGGGACATCCAT  GGGGTCAAACAGAGATGAATCAGGT     0      0.2        0
    ## 64  TCCTTGCTGTGGGACATCCAT   CGCCCATACTAAACGCTGTGTAGT     0      0.2        0
    ## 65  TCCTTGCTGTGGGACATCCAT       GCACGCCCATACTAAACGCT     0      0.2        0
    ## 66  TCCTTGCTGTGGGACATCCAT        CCTGCCTACCTCCAACCCT     0      0.2        0
    ## 67  TCCTTGCTGTGGGACATCCAT        CCCTGTGCACGCCCATACT     0      0.2        0
    ## 68  TCCTTGCTGTGGGACATCCAT   TGTGTAGTGGGGTCAAACAGAGAT     0      0.2        0
    ## 69  TCCTTGCTGTGGGACATCCAT TGTAGTGGGGTCAAACAGAGATGAAT     0      0.2        0
    ## 70  TCCTTGCTGTGGGACATCCAT      ACGCCCATACTAAACGCTGTG     0      0.2        0
    ## 71  TCCTTGCTGTGGGACATCCAT       TGCCTACCTCCAACCCTGTG     0      0.2        0
    ## 72  TCCTTGCTGTGGGACATCCAT   GCCCATACTAAACGCTGTGTAGTG     0      0.2        0
    ## 73  TCCTTGCTGTGGGACATCCAT      GCACGCCCATACTAAACGCTG     0      0.2        0
    ## 74  TCCTTGCTGTGGGACATCCAT       CCTGCCTACCTCCAACCCTG     0      0.2        0
    ## 75  TCCTTGCTGTGGGACATCCAT  TGTGTAGTGGGGTCAAACAGAGATG     0      0.2        0
    ## 76  TCCTTGCTGTGGGACATCCAT   CATACTAAACGCTGTGTAGTGGGG     0      0.2        0
    ## 77  TCCTTGCTGTGGGACATCCAT   GGGGTCAAACAGAGATGAATCAGG     0      0.2        0
    ## 78  TCCTTGCTGTGGGACATCCAT         CCTCCAACCCTGTGCACG     0      0.2        0
    ## 79  TCCTTGCTGTGGGACATCCAT    CGCCCATACTAAACGCTGTGTAG     0      0.2        0
    ## 80  TCCTTGCTGTGGGACATCCAT   GCTGTGTAGTGGGGTCAAACAGAG     0      0.2        0
    ## 81  TCCTTGCTGTGGGACATCCAT  GTGGGGTCAAACAGAGATGAATCAG     0      0.2        0
    ## 82  TCCTTGCTGTGGGACATCCAT     GCTGTGTAGTGGGGTCAAACAG     0      0.2        0
    ## 83  TCCTTGCTGTGGGACATCCAT     CTAAACGCTGTGTAGTGGGGTC     0      0.2        0
    ## 84  TCCTTGCTGTGGGACATCCAT       GCCTACCTCCAACCCTGTGC     0      0.2        0
    ## 85  TCCTTGCTGTGGGACATCCAT          TCCAACCCTGTGCACGC     0      0.2        0
    ## 86  TCCTTGCTGTGGGACATCCAT          CCAACCCTGTGCACGCC     0      0.2        0
    ## 87  TCCTTGCTGTGGGACATCCAT        CCCTGCCTACCTCCAACCC     0      0.2        0
    ## 88  TCCTTGCTGTGGGACATCCAT        CCCCTGCCTACCTCCAACC     0      0.2        0
    ## 89  TCCTTGCTGTGGGACATCCAT         CCCTGTGCACGCCCATAC     0      0.2        0
    ## 90  TCCTTGCTGTGGGACATCCAT      CCTACCTCCAACCCTGTGCAC     0      0.2        0
    ## 91  TCCTTGCTGTGGGACATCCAT        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 92  TCCTTGCTGTGGGACATCCAT     CCTGTGCACGCCCATACTAAAC     0      0.2        0
    ## 93  TCCTTGCTGTGGGACATCCAT      CGCTGTGTAGTGGGGTCAAAC     0      0.2        0
    ## 94  TCCTTGCTGTGGGACATCCAT     CGCCCATACTAAACGCTGTGTA     0      0.2        0
    ## 95  TCCTTGCTGTGGGACATCCAT GGGGTCAAACAGAGATGAATCAGGTA     0      0.2        0
    ## 96  TCCTTGCTGTGGGACATCCAT       CCCTGTGCACGCCCATACTA     0      0.2        0
    ## 97  TCCTTGCTGTGGGACATCCAT          CCCTGTGCACGCCCATA     0      0.2        0
    ## 98  TCCTTGCTGTGGGACATCCAT   TGTAGTGGGGTCAAACAGAGATGA     0      0.2        0
    ## 99  TCCTTGCTGTGGGACATCCAT    TGTGTAGTGGGGTCAAACAGAGA     0      0.2        0
    ## 100 TCCTTGCTGTGGGACATCCAT    GCTGTGTAGTGGGGTCAAACAGA     0      0.2        0
    ##     similar_signatures        missing_signatures enzyme digest_score fragments
    ## 1                      HPV18, HPV16, HPV4, HPV17 Acc65I            0         2
    ## 2                      HPV18, HPV16, HPV4, HPV17 Acc65I            0         2
    ## 3                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 4                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 5                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 6                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 7                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 8                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 9                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 10                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 11                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 12                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 13                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 14                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 15                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 16                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 17                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 18                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 19                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 20                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 21                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 22                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 23                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 24                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 25                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 26                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 27                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 28                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 29                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 30                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 31                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 32                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 33                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 34                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 35                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 36                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 37                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 38                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 39                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 40                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 41                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 42                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 43                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 44                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 45                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 46                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 47                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 48                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 49                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 50                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 51                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 52                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 53                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 54                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 55                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 56                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 57                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 58                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 59                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 60                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 61                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 62                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 63                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 64                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 65                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 66                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 67                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 68                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 69                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 70                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 71                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 72                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 73                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 74                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 75                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 76                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 77                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 78                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 79                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 80                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 81                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 82                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 83                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 84                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 85                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 86                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 87                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 88                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 89                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 90                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 91                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 92                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 93                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 94                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 95                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 96                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 97                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 98                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 99                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 100                    HPV18, HPV16, HPV4, HPV17                   0         0
    ## 
    ## [[2]]
    ##            forward_primer             reverse_primer score coverage products
    ## 1   GTTTGACCCCACTACACAGCG   CACTAACACCAACGCCTAAAGGTT     0      0.2        0
    ## 2   GTTTGACCCCACTACACAGCG    CACTAACACCAACGCCTAAAGGT     0      0.2        0
    ## 3   GTTTGACCCCACTACACAGCG      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 4   GTTTGACCCCACTACACAGCG       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 5   GTTTGACCCCACTACACAGCG        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 6   GTTTGACCCCACTACACAGCG    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 7   GTTTGACCCCACTACACAGCG   ACTAACACCAACGCCTAAAGGTTG     0      0.2        0
    ## 8   GTTTGACCCCACTACACAGCG       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 9   GTTTGACCCCACTACACAGCG    CCACTAACRCCAACRCCTAAAGG     0      0.2        0
    ## 10  GTTTGACCCCACTACACAGCG      GGATGCCCACTAACACCAACG     0      0.2        0
    ## 11  GTTTGACCCCACTACACAGCG    CCCACTAACACCAACGCCTAAAG     0      0.2        0
    ## 12  GTTTGACCCCACTACACAGCG       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 13  GTTTGACCCCACTACACAGCG       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 14  GTTTGACCCCACTACACAGCG        TGCCCACTAACACCAACGC     0      0.2        0
    ## 15  GTTTGACCCCACTACACAGCG         TGACCCCTGCCTACCTCC     0      0.2        0
    ## 16  GTTTGACCCCACTACACAGCG       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 17  GTTTGACCCCACTACACAGCG      CAACGCCTAAAGGTTGACCCC     0      0.2        0
    ## 18  GTTTGACCCCACTACACAGCG      CCAACGCCTAAAGGTTGACCC     0      0.2        0
    ## 19  GTTTGACCCCACTACACAGCG        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 20  GTTTGACCCCACTACACAGCG      ACCAACGCCTAAAGGTTGACC     0      0.2        0
    ## 21  GTTTGACCCCACTACACAGCG     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 22  GTTTGACCCCACTACACAGCG    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 23  GTTTGACCCCACTACACAGCG     ACACCAACGCCTAAAGGTTGAC     0      0.2        0
    ## 24  GTTTGACCCCACTACACAGCG     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 25  GTTTGACCCCACTACACAGCG  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 26  GTTTGACCCCACTACACAGCG        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 27  GTTTGACCCCACTACACAGCG     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 28  GTTTGACCCCACTACACAGCG     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 29  GTTTGACCCCACTACACAGCG      GCCCACTAACACCAACGCCTA     0      0.2        0
    ## 30  GTTTGACCCCACTACACAGCG   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 31  GTTTGACCCCACTACACAGCG   CTAACACCAACGCCTAAAGGTTGA     0      0.2        0
    ## 32  GTTTGACCCCACTACACAGCG         CACCYCTRCCTACCTCCA     0      0.2        0
    ## 33  GTTTGACCCCACTACACAGCG     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 34  GTTTGACCCCACTACACAGCG  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 35  GTTTGACCCCACTACACAGCG     GCCCACTAACACCAACGCCTAA     0      0.2        0
    ## 36  GTTTGACCCCACTACACAGCG  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 37  GTTTGACCCCACTACACAGCG        GACCCCTGCCTACCTCCAA     0      0.2        0
    ## 38  GTTTGACCCCACTACACAGCG    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 39  GTTTGACCCCACTACACAGCG     CCCACTAACACCAACGCCTAAA     0      0.2        0
    ## 40   TGACCCCACTACACAGCGTT      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 41   TGACCCCACTACACAGCGTT       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 42   TGACCCCACTACACAGCGTT        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 43   TGACCCCACTACACAGCGTT    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 44   TGACCCCACTACACAGCGTT       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 45   TGACCCCACTACACAGCGTT       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 46   TGACCCCACTACACAGCGTT       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 47   TGACCCCACTACACAGCGTT       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 48   TGACCCCACTACACAGCGTT        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 49   TGACCCCACTACACAGCGTT     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 50   TGACCCCACTACACAGCGTT    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 51   TGACCCCACTACACAGCGTT     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 52   TGACCCCACTACACAGCGTT  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 53   TGACCCCACTACACAGCGTT     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 54   TGACCCCACTACACAGCGTT     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 55   TGACCCCACTACACAGCGTT   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 56   TGACCCCACTACACAGCGTT     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 57   TGACCCCACTACACAGCGTT  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 58   TGACCCCACTACACAGCGTT  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 59   TGACCCCACTACACAGCGTT    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 60  TCCTTGCTGTGGGACATCCAT      CGCCCATACTAAACGCTGTGT     0      0.2        0
    ## 61  TCCTTGCTGTGGGACATCCAT      CACGCCCATACTAAACGCTGT     0      0.2        0
    ## 62  TCCTTGCTGTGGGACATCCAT        TGCCTACCTCCAACCCTGT     0      0.2        0
    ## 63  TCCTTGCTGTGGGACATCCAT  GGGGTCAAACAGAGATGAATCAGGT     0      0.2        0
    ## 64  TCCTTGCTGTGGGACATCCAT   CGCCCATACTAAACGCTGTGTAGT     0      0.2        0
    ## 65  TCCTTGCTGTGGGACATCCAT       GCACGCCCATACTAAACGCT     0      0.2        0
    ## 66  TCCTTGCTGTGGGACATCCAT        CCTGCCTACCTCCAACCCT     0      0.2        0
    ## 67  TCCTTGCTGTGGGACATCCAT        CCCTGTGCACGCCCATACT     0      0.2        0
    ## 68  TCCTTGCTGTGGGACATCCAT   TGTGTAGTGGGGTCAAACAGAGAT     0      0.2        0
    ## 69  TCCTTGCTGTGGGACATCCAT TGTAGTGGGGTCAAACAGAGATGAAT     0      0.2        0
    ## 70  TCCTTGCTGTGGGACATCCAT      ACGCCCATACTAAACGCTGTG     0      0.2        0
    ## 71  TCCTTGCTGTGGGACATCCAT       TGCCTACCTCCAACCCTGTG     0      0.2        0
    ## 72  TCCTTGCTGTGGGACATCCAT   GCCCATACTAAACGCTGTGTAGTG     0      0.2        0
    ## 73  TCCTTGCTGTGGGACATCCAT      GCACGCCCATACTAAACGCTG     0      0.2        0
    ## 74  TCCTTGCTGTGGGACATCCAT       CCTGCCTACCTCCAACCCTG     0      0.2        0
    ## 75  TCCTTGCTGTGGGACATCCAT  TGTGTAGTGGGGTCAAACAGAGATG     0      0.2        0
    ## 76  TCCTTGCTGTGGGACATCCAT   CATACTAAACGCTGTGTAGTGGGG     0      0.2        0
    ## 77  TCCTTGCTGTGGGACATCCAT   GGGGTCAAACAGAGATGAATCAGG     0      0.2        0
    ## 78  TCCTTGCTGTGGGACATCCAT         CCTCCAACCCTGTGCACG     0      0.2        0
    ## 79  TCCTTGCTGTGGGACATCCAT    CGCCCATACTAAACGCTGTGTAG     0      0.2        0
    ## 80  TCCTTGCTGTGGGACATCCAT   GCTGTGTAGTGGGGTCAAACAGAG     0      0.2        0
    ## 81  TCCTTGCTGTGGGACATCCAT  GTGGGGTCAAACAGAGATGAATCAG     0      0.2        0
    ## 82  TCCTTGCTGTGGGACATCCAT     GCTGTGTAGTGGGGTCAAACAG     0      0.2        0
    ## 83  TCCTTGCTGTGGGACATCCAT     CTAAACGCTGTGTAGTGGGGTC     0      0.2        0
    ## 84  TCCTTGCTGTGGGACATCCAT       GCCTACCTCCAACCCTGTGC     0      0.2        0
    ## 85  TCCTTGCTGTGGGACATCCAT          TCCAACCCTGTGCACGC     0      0.2        0
    ## 86  TCCTTGCTGTGGGACATCCAT          CCAACCCTGTGCACGCC     0      0.2        0
    ## 87  TCCTTGCTGTGGGACATCCAT        CCCTGCCTACCTCCAACCC     0      0.2        0
    ## 88  TCCTTGCTGTGGGACATCCAT        CCCCTGCCTACCTCCAACC     0      0.2        0
    ## 89  TCCTTGCTGTGGGACATCCAT         CCCTGTGCACGCCCATAC     0      0.2        0
    ## 90  TCCTTGCTGTGGGACATCCAT      CCTACCTCCAACCCTGTGCAC     0      0.2        0
    ## 91  TCCTTGCTGTGGGACATCCAT        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 92  TCCTTGCTGTGGGACATCCAT     CCTGTGCACGCCCATACTAAAC     0      0.2        0
    ## 93  TCCTTGCTGTGGGACATCCAT      CGCTGTGTAGTGGGGTCAAAC     0      0.2        0
    ## 94  TCCTTGCTGTGGGACATCCAT     CGCCCATACTAAACGCTGTGTA     0      0.2        0
    ## 95  TCCTTGCTGTGGGACATCCAT GGGGTCAAACAGAGATGAATCAGGTA     0      0.2        0
    ## 96  TCCTTGCTGTGGGACATCCAT       CCCTGTGCACGCCCATACTA     0      0.2        0
    ## 97  TCCTTGCTGTGGGACATCCAT          CCCTGTGCACGCCCATA     0      0.2        0
    ## 98  TCCTTGCTGTGGGACATCCAT   TGTAGTGGGGTCAAACAGAGATGA     0      0.2        0
    ## 99  TCCTTGCTGTGGGACATCCAT    TGTGTAGTGGGGTCAAACAGAGA     0      0.2        0
    ## 100 TCCTTGCTGTGGGACATCCAT    GCTGTGTAGTGGGGTCAAACAGA     0      0.2        0
    ##     similar_signatures        missing_signatures enzyme digest_score fragments
    ## 1                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 2                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 3                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 4                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 5                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 6                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 7                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 8                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 9                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 10                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 11                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 12                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 13                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 14                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 15                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 16                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 17                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 18                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 19                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 20                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 21                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 22                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 23                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 24                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 25                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 26                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 27                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 28                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 29                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 30                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 31                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 32                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 33                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 34                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 35                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 36                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 37                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 38                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 39                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 40                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 41                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 42                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 43                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 44                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 45                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 46                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 47                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 48                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 49                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 50                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 51                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 52                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 53                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 54                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 55                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 56                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 57                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 58                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 59                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 60                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 61                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 62                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 63                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 64                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 65                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 66                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 67                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 68                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 69                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 70                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 71                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 72                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 73                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 74                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 75                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 76                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 77                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 78                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 79                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 80                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 81                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 82                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 83                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 84                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 85                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 86                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 87                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 88                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 89                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 90                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 91                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 92                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 93                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 94                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 95                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 96                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 97                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 98                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 99                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 100                    HPV18, HPV16, HPV4, HPV17                   0         0
    ## 
    ## [[3]]
    ##            forward_primer             reverse_primer score coverage products
    ## 1   GTTTGACCCCACTACACAGCG   CACTAACACCAACGCCTAAAGGTT     0      0.2        0
    ## 2   GTTTGACCCCACTACACAGCG    CACTAACACCAACGCCTAAAGGT     0      0.2        0
    ## 3   GTTTGACCCCACTACACAGCG      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 4   GTTTGACCCCACTACACAGCG       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 5   GTTTGACCCCACTACACAGCG        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 6   GTTTGACCCCACTACACAGCG    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 7   GTTTGACCCCACTACACAGCG   ACTAACACCAACGCCTAAAGGTTG     0      0.2        0
    ## 8   GTTTGACCCCACTACACAGCG       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 9   GTTTGACCCCACTACACAGCG    CCACTAACRCCAACRCCTAAAGG     0      0.2        0
    ## 10  GTTTGACCCCACTACACAGCG      GGATGCCCACTAACACCAACG     0      0.2        0
    ## 11  GTTTGACCCCACTACACAGCG    CCCACTAACACCAACGCCTAAAG     0      0.2        0
    ## 12  GTTTGACCCCACTACACAGCG       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 13  GTTTGACCCCACTACACAGCG       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 14  GTTTGACCCCACTACACAGCG        TGCCCACTAACACCAACGC     0      0.2        0
    ## 15  GTTTGACCCCACTACACAGCG         TGACCCCTGCCTACCTCC     0      0.2        0
    ## 16  GTTTGACCCCACTACACAGCG       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 17  GTTTGACCCCACTACACAGCG      CAACGCCTAAAGGTTGACCCC     0      0.2        0
    ## 18  GTTTGACCCCACTACACAGCG      CCAACGCCTAAAGGTTGACCC     0      0.2        0
    ## 19  GTTTGACCCCACTACACAGCG        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 20  GTTTGACCCCACTACACAGCG      ACCAACGCCTAAAGGTTGACC     0      0.2        0
    ## 21  GTTTGACCCCACTACACAGCG     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 22  GTTTGACCCCACTACACAGCG    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 23  GTTTGACCCCACTACACAGCG     ACACCAACGCCTAAAGGTTGAC     0      0.2        0
    ## 24  GTTTGACCCCACTACACAGCG     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 25  GTTTGACCCCACTACACAGCG  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 26  GTTTGACCCCACTACACAGCG        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 27  GTTTGACCCCACTACACAGCG     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 28  GTTTGACCCCACTACACAGCG     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 29  GTTTGACCCCACTACACAGCG      GCCCACTAACACCAACGCCTA     0      0.2        0
    ## 30  GTTTGACCCCACTACACAGCG   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 31  GTTTGACCCCACTACACAGCG   CTAACACCAACGCCTAAAGGTTGA     0      0.2        0
    ## 32  GTTTGACCCCACTACACAGCG         CACCYCTRCCTACCTCCA     0      0.2        0
    ## 33  GTTTGACCCCACTACACAGCG     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 34  GTTTGACCCCACTACACAGCG  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 35  GTTTGACCCCACTACACAGCG     GCCCACTAACACCAACGCCTAA     0      0.2        0
    ## 36  GTTTGACCCCACTACACAGCG  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 37  GTTTGACCCCACTACACAGCG        GACCCCTGCCTACCTCCAA     0      0.2        0
    ## 38  GTTTGACCCCACTACACAGCG    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 39  GTTTGACCCCACTACACAGCG     CCCACTAACACCAACGCCTAAA     0      0.2        0
    ## 40   TGACCCCACTACACAGCGTT      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 41   TGACCCCACTACACAGCGTT       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 42   TGACCCCACTACACAGCGTT        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 43   TGACCCCACTACACAGCGTT    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 44   TGACCCCACTACACAGCGTT       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 45   TGACCCCACTACACAGCGTT       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 46   TGACCCCACTACACAGCGTT       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 47   TGACCCCACTACACAGCGTT       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 48   TGACCCCACTACACAGCGTT        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 49   TGACCCCACTACACAGCGTT     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 50   TGACCCCACTACACAGCGTT    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 51   TGACCCCACTACACAGCGTT     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 52   TGACCCCACTACACAGCGTT  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 53   TGACCCCACTACACAGCGTT     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 54   TGACCCCACTACACAGCGTT     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 55   TGACCCCACTACACAGCGTT   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 56   TGACCCCACTACACAGCGTT     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 57   TGACCCCACTACACAGCGTT  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 58   TGACCCCACTACACAGCGTT  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 59   TGACCCCACTACACAGCGTT    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 60  TCCTTGCTGTGGGACATCCAT      CGCCCATACTAAACGCTGTGT     0      0.2        0
    ## 61  TCCTTGCTGTGGGACATCCAT      CACGCCCATACTAAACGCTGT     0      0.2        0
    ## 62  TCCTTGCTGTGGGACATCCAT        TGCCTACCTCCAACCCTGT     0      0.2        0
    ## 63  TCCTTGCTGTGGGACATCCAT  GGGGTCAAACAGAGATGAATCAGGT     0      0.2        0
    ## 64  TCCTTGCTGTGGGACATCCAT   CGCCCATACTAAACGCTGTGTAGT     0      0.2        0
    ## 65  TCCTTGCTGTGGGACATCCAT       GCACGCCCATACTAAACGCT     0      0.2        0
    ## 66  TCCTTGCTGTGGGACATCCAT        CCTGCCTACCTCCAACCCT     0      0.2        0
    ## 67  TCCTTGCTGTGGGACATCCAT        CCCTGTGCACGCCCATACT     0      0.2        0
    ## 68  TCCTTGCTGTGGGACATCCAT   TGTGTAGTGGGGTCAAACAGAGAT     0      0.2        0
    ## 69  TCCTTGCTGTGGGACATCCAT TGTAGTGGGGTCAAACAGAGATGAAT     0      0.2        0
    ## 70  TCCTTGCTGTGGGACATCCAT      ACGCCCATACTAAACGCTGTG     0      0.2        0
    ## 71  TCCTTGCTGTGGGACATCCAT       TGCCTACCTCCAACCCTGTG     0      0.2        0
    ## 72  TCCTTGCTGTGGGACATCCAT   GCCCATACTAAACGCTGTGTAGTG     0      0.2        0
    ## 73  TCCTTGCTGTGGGACATCCAT      GCACGCCCATACTAAACGCTG     0      0.2        0
    ## 74  TCCTTGCTGTGGGACATCCAT       CCTGCCTACCTCCAACCCTG     0      0.2        0
    ## 75  TCCTTGCTGTGGGACATCCAT  TGTGTAGTGGGGTCAAACAGAGATG     0      0.2        0
    ## 76  TCCTTGCTGTGGGACATCCAT   CATACTAAACGCTGTGTAGTGGGG     0      0.2        0
    ## 77  TCCTTGCTGTGGGACATCCAT   GGGGTCAAACAGAGATGAATCAGG     0      0.2        0
    ## 78  TCCTTGCTGTGGGACATCCAT         CCTCCAACCCTGTGCACG     0      0.2        0
    ## 79  TCCTTGCTGTGGGACATCCAT    CGCCCATACTAAACGCTGTGTAG     0      0.2        0
    ## 80  TCCTTGCTGTGGGACATCCAT   GCTGTGTAGTGGGGTCAAACAGAG     0      0.2        0
    ## 81  TCCTTGCTGTGGGACATCCAT  GTGGGGTCAAACAGAGATGAATCAG     0      0.2        0
    ## 82  TCCTTGCTGTGGGACATCCAT     GCTGTGTAGTGGGGTCAAACAG     0      0.2        0
    ## 83  TCCTTGCTGTGGGACATCCAT     CTAAACGCTGTGTAGTGGGGTC     0      0.2        0
    ## 84  TCCTTGCTGTGGGACATCCAT       GCCTACCTCCAACCCTGTGC     0      0.2        0
    ## 85  TCCTTGCTGTGGGACATCCAT          TCCAACCCTGTGCACGC     0      0.2        0
    ## 86  TCCTTGCTGTGGGACATCCAT          CCAACCCTGTGCACGCC     0      0.2        0
    ## 87  TCCTTGCTGTGGGACATCCAT        CCCTGCCTACCTCCAACCC     0      0.2        0
    ## 88  TCCTTGCTGTGGGACATCCAT        CCCCTGCCTACCTCCAACC     0      0.2        0
    ## 89  TCCTTGCTGTGGGACATCCAT         CCCTGTGCACGCCCATAC     0      0.2        0
    ## 90  TCCTTGCTGTGGGACATCCAT      CCTACCTCCAACCCTGTGCAC     0      0.2        0
    ## 91  TCCTTGCTGTGGGACATCCAT        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 92  TCCTTGCTGTGGGACATCCAT     CCTGTGCACGCCCATACTAAAC     0      0.2        0
    ## 93  TCCTTGCTGTGGGACATCCAT      CGCTGTGTAGTGGGGTCAAAC     0      0.2        0
    ## 94  TCCTTGCTGTGGGACATCCAT     CGCCCATACTAAACGCTGTGTA     0      0.2        0
    ## 95  TCCTTGCTGTGGGACATCCAT GGGGTCAAACAGAGATGAATCAGGTA     0      0.2        0
    ## 96  TCCTTGCTGTGGGACATCCAT       CCCTGTGCACGCCCATACTA     0      0.2        0
    ## 97  TCCTTGCTGTGGGACATCCAT          CCCTGTGCACGCCCATA     0      0.2        0
    ## 98  TCCTTGCTGTGGGACATCCAT   TGTAGTGGGGTCAAACAGAGATGA     0      0.2        0
    ## 99  TCCTTGCTGTGGGACATCCAT    TGTGTAGTGGGGTCAAACAGAGA     0      0.2        0
    ## 100 TCCTTGCTGTGGGACATCCAT    GCTGTGTAGTGGGGTCAAACAGA     0      0.2        0
    ##     similar_signatures        missing_signatures enzyme digest_score fragments
    ## 1                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 2                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 3                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 4                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 5                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 6                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 7                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 8                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 9                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 10                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 11                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 12                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 13                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 14                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 15                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 16                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 17                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 18                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 19                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 20                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 21                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 22                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 23                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 24                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 25                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 26                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 27                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 28                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 29                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 30                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 31                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 32                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 33                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 34                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 35                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 36                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 37                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 38                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 39                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 40                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 41                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 42                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 43                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 44                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 45                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 46                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 47                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 48                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 49                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 50                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 51                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 52                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 53                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 54                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 55                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 56                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 57                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 58                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 59                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 60                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 61                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 62                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 63                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 64                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 65                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 66                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 67                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 68                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 69                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 70                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 71                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 72                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 73                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 74                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 75                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 76                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 77                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 78                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 79                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 80                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 81                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 82                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 83                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 84                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 85                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 86                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 87                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 88                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 89                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 90                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 91                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 92                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 93                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 94                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 95                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 96                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 97                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 98                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 99                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 100                    HPV18, HPV16, HPV4, HPV17                   0         0
    ## 
    ## [[4]]
    ##            forward_primer             reverse_primer score coverage products
    ## 1   GTTTGACCCCACTACACAGCG       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 2    TGACCCCACTACACAGCGTT       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 3   GTTTGACCCCACTACACAGCG   CACTAACACCAACGCCTAAAGGTT     0      0.2        0
    ## 4   GTTTGACCCCACTACACAGCG    CACTAACACCAACGCCTAAAGGT     0      0.2        0
    ## 5   GTTTGACCCCACTACACAGCG      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 6   GTTTGACCCCACTACACAGCG       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 7   GTTTGACCCCACTACACAGCG        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 8   GTTTGACCCCACTACACAGCG    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 9   GTTTGACCCCACTACACAGCG   ACTAACACCAACGCCTAAAGGTTG     0      0.2        0
    ## 10  GTTTGACCCCACTACACAGCG       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 11  GTTTGACCCCACTACACAGCG    CCACTAACRCCAACRCCTAAAGG     0      0.2        0
    ## 12  GTTTGACCCCACTACACAGCG      GGATGCCCACTAACACCAACG     0      0.2        0
    ## 13  GTTTGACCCCACTACACAGCG    CCCACTAACACCAACGCCTAAAG     0      0.2        0
    ## 14  GTTTGACCCCACTACACAGCG       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 15  GTTTGACCCCACTACACAGCG       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 16  GTTTGACCCCACTACACAGCG        TGCCCACTAACACCAACGC     0      0.2        0
    ## 17  GTTTGACCCCACTACACAGCG         TGACCCCTGCCTACCTCC     0      0.2        0
    ## 18  GTTTGACCCCACTACACAGCG      CAACGCCTAAAGGTTGACCCC     0      0.2        0
    ## 19  GTTTGACCCCACTACACAGCG      CCAACGCCTAAAGGTTGACCC     0      0.2        0
    ## 20  GTTTGACCCCACTACACAGCG        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 21  GTTTGACCCCACTACACAGCG      ACCAACGCCTAAAGGTTGACC     0      0.2        0
    ## 22  GTTTGACCCCACTACACAGCG     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 23  GTTTGACCCCACTACACAGCG    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 24  GTTTGACCCCACTACACAGCG     ACACCAACGCCTAAAGGTTGAC     0      0.2        0
    ## 25  GTTTGACCCCACTACACAGCG     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 26  GTTTGACCCCACTACACAGCG  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 27  GTTTGACCCCACTACACAGCG        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 28  GTTTGACCCCACTACACAGCG     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 29  GTTTGACCCCACTACACAGCG     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 30  GTTTGACCCCACTACACAGCG      GCCCACTAACACCAACGCCTA     0      0.2        0
    ## 31  GTTTGACCCCACTACACAGCG   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 32  GTTTGACCCCACTACACAGCG   CTAACACCAACGCCTAAAGGTTGA     0      0.2        0
    ## 33  GTTTGACCCCACTACACAGCG         CACCYCTRCCTACCTCCA     0      0.2        0
    ## 34  GTTTGACCCCACTACACAGCG     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 35  GTTTGACCCCACTACACAGCG  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 36  GTTTGACCCCACTACACAGCG     GCCCACTAACACCAACGCCTAA     0      0.2        0
    ## 37  GTTTGACCCCACTACACAGCG  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 38  GTTTGACCCCACTACACAGCG        GACCCCTGCCTACCTCCAA     0      0.2        0
    ## 39  GTTTGACCCCACTACACAGCG    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 40  GTTTGACCCCACTACACAGCG     CCCACTAACACCAACGCCTAAA     0      0.2        0
    ## 41   TGACCCCACTACACAGCGTT      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 42   TGACCCCACTACACAGCGTT       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 43   TGACCCCACTACACAGCGTT        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 44   TGACCCCACTACACAGCGTT    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 45   TGACCCCACTACACAGCGTT       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 46   TGACCCCACTACACAGCGTT       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 47   TGACCCCACTACACAGCGTT       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 48   TGACCCCACTACACAGCGTT        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 49   TGACCCCACTACACAGCGTT     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 50   TGACCCCACTACACAGCGTT    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 51   TGACCCCACTACACAGCGTT     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 52   TGACCCCACTACACAGCGTT  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 53   TGACCCCACTACACAGCGTT     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 54   TGACCCCACTACACAGCGTT     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 55   TGACCCCACTACACAGCGTT   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 56   TGACCCCACTACACAGCGTT     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 57   TGACCCCACTACACAGCGTT  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 58   TGACCCCACTACACAGCGTT  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 59   TGACCCCACTACACAGCGTT    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 60  TCCTTGCTGTGGGACATCCAT      CGCCCATACTAAACGCTGTGT     0      0.2        0
    ## 61  TCCTTGCTGTGGGACATCCAT      CACGCCCATACTAAACGCTGT     0      0.2        0
    ## 62  TCCTTGCTGTGGGACATCCAT        TGCCTACCTCCAACCCTGT     0      0.2        0
    ## 63  TCCTTGCTGTGGGACATCCAT  GGGGTCAAACAGAGATGAATCAGGT     0      0.2        0
    ## 64  TCCTTGCTGTGGGACATCCAT   CGCCCATACTAAACGCTGTGTAGT     0      0.2        0
    ## 65  TCCTTGCTGTGGGACATCCAT       GCACGCCCATACTAAACGCT     0      0.2        0
    ## 66  TCCTTGCTGTGGGACATCCAT        CCTGCCTACCTCCAACCCT     0      0.2        0
    ## 67  TCCTTGCTGTGGGACATCCAT        CCCTGTGCACGCCCATACT     0      0.2        0
    ## 68  TCCTTGCTGTGGGACATCCAT   TGTGTAGTGGGGTCAAACAGAGAT     0      0.2        0
    ## 69  TCCTTGCTGTGGGACATCCAT TGTAGTGGGGTCAAACAGAGATGAAT     0      0.2        0
    ## 70  TCCTTGCTGTGGGACATCCAT      ACGCCCATACTAAACGCTGTG     0      0.2        0
    ## 71  TCCTTGCTGTGGGACATCCAT       TGCCTACCTCCAACCCTGTG     0      0.2        0
    ## 72  TCCTTGCTGTGGGACATCCAT   GCCCATACTAAACGCTGTGTAGTG     0      0.2        0
    ## 73  TCCTTGCTGTGGGACATCCAT      GCACGCCCATACTAAACGCTG     0      0.2        0
    ## 74  TCCTTGCTGTGGGACATCCAT       CCTGCCTACCTCCAACCCTG     0      0.2        0
    ## 75  TCCTTGCTGTGGGACATCCAT  TGTGTAGTGGGGTCAAACAGAGATG     0      0.2        0
    ## 76  TCCTTGCTGTGGGACATCCAT   CATACTAAACGCTGTGTAGTGGGG     0      0.2        0
    ## 77  TCCTTGCTGTGGGACATCCAT   GGGGTCAAACAGAGATGAATCAGG     0      0.2        0
    ## 78  TCCTTGCTGTGGGACATCCAT         CCTCCAACCCTGTGCACG     0      0.2        0
    ## 79  TCCTTGCTGTGGGACATCCAT    CGCCCATACTAAACGCTGTGTAG     0      0.2        0
    ## 80  TCCTTGCTGTGGGACATCCAT   GCTGTGTAGTGGGGTCAAACAGAG     0      0.2        0
    ## 81  TCCTTGCTGTGGGACATCCAT  GTGGGGTCAAACAGAGATGAATCAG     0      0.2        0
    ## 82  TCCTTGCTGTGGGACATCCAT     GCTGTGTAGTGGGGTCAAACAG     0      0.2        0
    ## 83  TCCTTGCTGTGGGACATCCAT     CTAAACGCTGTGTAGTGGGGTC     0      0.2        0
    ## 84  TCCTTGCTGTGGGACATCCAT       GCCTACCTCCAACCCTGTGC     0      0.2        0
    ## 85  TCCTTGCTGTGGGACATCCAT          TCCAACCCTGTGCACGC     0      0.2        0
    ## 86  TCCTTGCTGTGGGACATCCAT          CCAACCCTGTGCACGCC     0      0.2        0
    ## 87  TCCTTGCTGTGGGACATCCAT        CCCTGCCTACCTCCAACCC     0      0.2        0
    ## 88  TCCTTGCTGTGGGACATCCAT        CCCCTGCCTACCTCCAACC     0      0.2        0
    ## 89  TCCTTGCTGTGGGACATCCAT         CCCTGTGCACGCCCATAC     0      0.2        0
    ## 90  TCCTTGCTGTGGGACATCCAT      CCTACCTCCAACCCTGTGCAC     0      0.2        0
    ## 91  TCCTTGCTGTGGGACATCCAT        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 92  TCCTTGCTGTGGGACATCCAT     CCTGTGCACGCCCATACTAAAC     0      0.2        0
    ## 93  TCCTTGCTGTGGGACATCCAT      CGCTGTGTAGTGGGGTCAAAC     0      0.2        0
    ## 94  TCCTTGCTGTGGGACATCCAT     CGCCCATACTAAACGCTGTGTA     0      0.2        0
    ## 95  TCCTTGCTGTGGGACATCCAT GGGGTCAAACAGAGATGAATCAGGTA     0      0.2        0
    ## 96  TCCTTGCTGTGGGACATCCAT       CCCTGTGCACGCCCATACTA     0      0.2        0
    ## 97  TCCTTGCTGTGGGACATCCAT          CCCTGTGCACGCCCATA     0      0.2        0
    ## 98  TCCTTGCTGTGGGACATCCAT   TGTAGTGGGGTCAAACAGAGATGA     0      0.2        0
    ## 99  TCCTTGCTGTGGGACATCCAT    TGTGTAGTGGGGTCAAACAGAGA     0      0.2        0
    ## 100 TCCTTGCTGTGGGACATCCAT    GCTGTGTAGTGGGGTCAAACAGA     0      0.2        0
    ##     similar_signatures        missing_signatures enzyme digest_score fragments
    ## 1                      HPV18, HPV16, HPV4, HPV17 Acc65I            0         2
    ## 2                      HPV18, HPV16, HPV4, HPV17 Acc65I            0         2
    ## 3                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 4                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 5                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 6                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 7                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 8                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 9                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 10                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 11                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 12                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 13                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 14                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 15                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 16                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 17                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 18                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 19                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 20                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 21                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 22                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 23                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 24                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 25                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 26                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 27                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 28                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 29                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 30                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 31                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 32                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 33                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 34                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 35                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 36                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 37                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 38                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 39                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 40                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 41                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 42                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 43                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 44                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 45                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 46                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 47                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 48                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 49                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 50                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 51                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 52                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 53                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 54                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 55                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 56                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 57                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 58                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 59                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 60                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 61                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 62                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 63                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 64                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 65                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 66                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 67                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 68                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 69                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 70                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 71                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 72                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 73                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 74                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 75                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 76                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 77                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 78                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 79                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 80                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 81                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 82                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 83                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 84                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 85                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 86                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 87                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 88                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 89                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 90                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 91                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 92                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 93                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 94                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 95                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 96                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 97                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 98                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 99                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 100                    HPV18, HPV16, HPV4, HPV17                   0         0
    ## 
    ## [[5]]
    ##            forward_primer             reverse_primer score coverage products
    ## 1   GTTTGACCCCACTACACAGCG   CACTAACACCAACGCCTAAAGGTT     0      0.2        0
    ## 2   GTTTGACCCCACTACACAGCG    CACTAACACCAACGCCTAAAGGT     0      0.2        0
    ## 3   GTTTGACCCCACTACACAGCG      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 4   GTTTGACCCCACTACACAGCG       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 5   GTTTGACCCCACTACACAGCG        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 6   GTTTGACCCCACTACACAGCG    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 7   GTTTGACCCCACTACACAGCG   ACTAACACCAACGCCTAAAGGTTG     0      0.2        0
    ## 8   GTTTGACCCCACTACACAGCG       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 9   GTTTGACCCCACTACACAGCG    CCACTAACRCCAACRCCTAAAGG     0      0.2        0
    ## 10  GTTTGACCCCACTACACAGCG      GGATGCCCACTAACACCAACG     0      0.2        0
    ## 11  GTTTGACCCCACTACACAGCG    CCCACTAACACCAACGCCTAAAG     0      0.2        0
    ## 12  GTTTGACCCCACTACACAGCG       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 13  GTTTGACCCCACTACACAGCG       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 14  GTTTGACCCCACTACACAGCG        TGCCCACTAACACCAACGC     0      0.2        0
    ## 15  GTTTGACCCCACTACACAGCG         TGACCCCTGCCTACCTCC     0      0.2        0
    ## 16  GTTTGACCCCACTACACAGCG       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 17  GTTTGACCCCACTACACAGCG      CAACGCCTAAAGGTTGACCCC     0      0.2        0
    ## 18  GTTTGACCCCACTACACAGCG      CCAACGCCTAAAGGTTGACCC     0      0.2        0
    ## 19  GTTTGACCCCACTACACAGCG        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 20  GTTTGACCCCACTACACAGCG      ACCAACGCCTAAAGGTTGACC     0      0.2        0
    ## 21  GTTTGACCCCACTACACAGCG     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 22  GTTTGACCCCACTACACAGCG    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 23  GTTTGACCCCACTACACAGCG     ACACCAACGCCTAAAGGTTGAC     0      0.2        0
    ## 24  GTTTGACCCCACTACACAGCG     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 25  GTTTGACCCCACTACACAGCG  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 26  GTTTGACCCCACTACACAGCG        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 27  GTTTGACCCCACTACACAGCG     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 28  GTTTGACCCCACTACACAGCG     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 29  GTTTGACCCCACTACACAGCG      GCCCACTAACACCAACGCCTA     0      0.2        0
    ## 30  GTTTGACCCCACTACACAGCG   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 31  GTTTGACCCCACTACACAGCG   CTAACACCAACGCCTAAAGGTTGA     0      0.2        0
    ## 32  GTTTGACCCCACTACACAGCG         CACCYCTRCCTACCTCCA     0      0.2        0
    ## 33  GTTTGACCCCACTACACAGCG     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 34  GTTTGACCCCACTACACAGCG  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 35  GTTTGACCCCACTACACAGCG     GCCCACTAACACCAACGCCTAA     0      0.2        0
    ## 36  GTTTGACCCCACTACACAGCG  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 37  GTTTGACCCCACTACACAGCG        GACCCCTGCCTACCTCCAA     0      0.2        0
    ## 38  GTTTGACCCCACTACACAGCG    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 39  GTTTGACCCCACTACACAGCG     CCCACTAACACCAACGCCTAAA     0      0.2        0
    ## 40   TGACCCCACTACACAGCGTT      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 41   TGACCCCACTACACAGCGTT       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 42   TGACCCCACTACACAGCGTT        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 43   TGACCCCACTACACAGCGTT    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 44   TGACCCCACTACACAGCGTT       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 45   TGACCCCACTACACAGCGTT       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 46   TGACCCCACTACACAGCGTT       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 47   TGACCCCACTACACAGCGTT       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 48   TGACCCCACTACACAGCGTT        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 49   TGACCCCACTACACAGCGTT     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 50   TGACCCCACTACACAGCGTT    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 51   TGACCCCACTACACAGCGTT     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 52   TGACCCCACTACACAGCGTT  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 53   TGACCCCACTACACAGCGTT     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 54   TGACCCCACTACACAGCGTT     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 55   TGACCCCACTACACAGCGTT   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 56   TGACCCCACTACACAGCGTT     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 57   TGACCCCACTACACAGCGTT  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 58   TGACCCCACTACACAGCGTT  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 59   TGACCCCACTACACAGCGTT    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 60  TCCTTGCTGTGGGACATCCAT      CGCCCATACTAAACGCTGTGT     0      0.2        0
    ## 61  TCCTTGCTGTGGGACATCCAT      CACGCCCATACTAAACGCTGT     0      0.2        0
    ## 62  TCCTTGCTGTGGGACATCCAT        TGCCTACCTCCAACCCTGT     0      0.2        0
    ## 63  TCCTTGCTGTGGGACATCCAT  GGGGTCAAACAGAGATGAATCAGGT     0      0.2        0
    ## 64  TCCTTGCTGTGGGACATCCAT   CGCCCATACTAAACGCTGTGTAGT     0      0.2        0
    ## 65  TCCTTGCTGTGGGACATCCAT       GCACGCCCATACTAAACGCT     0      0.2        0
    ## 66  TCCTTGCTGTGGGACATCCAT        CCTGCCTACCTCCAACCCT     0      0.2        0
    ## 67  TCCTTGCTGTGGGACATCCAT        CCCTGTGCACGCCCATACT     0      0.2        0
    ## 68  TCCTTGCTGTGGGACATCCAT   TGTGTAGTGGGGTCAAACAGAGAT     0      0.2        0
    ## 69  TCCTTGCTGTGGGACATCCAT TGTAGTGGGGTCAAACAGAGATGAAT     0      0.2        0
    ## 70  TCCTTGCTGTGGGACATCCAT      ACGCCCATACTAAACGCTGTG     0      0.2        0
    ## 71  TCCTTGCTGTGGGACATCCAT       TGCCTACCTCCAACCCTGTG     0      0.2        0
    ## 72  TCCTTGCTGTGGGACATCCAT   GCCCATACTAAACGCTGTGTAGTG     0      0.2        0
    ## 73  TCCTTGCTGTGGGACATCCAT      GCACGCCCATACTAAACGCTG     0      0.2        0
    ## 74  TCCTTGCTGTGGGACATCCAT       CCTGCCTACCTCCAACCCTG     0      0.2        0
    ## 75  TCCTTGCTGTGGGACATCCAT  TGTGTAGTGGGGTCAAACAGAGATG     0      0.2        0
    ## 76  TCCTTGCTGTGGGACATCCAT   CATACTAAACGCTGTGTAGTGGGG     0      0.2        0
    ## 77  TCCTTGCTGTGGGACATCCAT   GGGGTCAAACAGAGATGAATCAGG     0      0.2        0
    ## 78  TCCTTGCTGTGGGACATCCAT         CCTCCAACCCTGTGCACG     0      0.2        0
    ## 79  TCCTTGCTGTGGGACATCCAT    CGCCCATACTAAACGCTGTGTAG     0      0.2        0
    ## 80  TCCTTGCTGTGGGACATCCAT   GCTGTGTAGTGGGGTCAAACAGAG     0      0.2        0
    ## 81  TCCTTGCTGTGGGACATCCAT  GTGGGGTCAAACAGAGATGAATCAG     0      0.2        0
    ## 82  TCCTTGCTGTGGGACATCCAT     GCTGTGTAGTGGGGTCAAACAG     0      0.2        0
    ## 83  TCCTTGCTGTGGGACATCCAT     CTAAACGCTGTGTAGTGGGGTC     0      0.2        0
    ## 84  TCCTTGCTGTGGGACATCCAT       GCCTACCTCCAACCCTGTGC     0      0.2        0
    ## 85  TCCTTGCTGTGGGACATCCAT          TCCAACCCTGTGCACGC     0      0.2        0
    ## 86  TCCTTGCTGTGGGACATCCAT          CCAACCCTGTGCACGCC     0      0.2        0
    ## 87  TCCTTGCTGTGGGACATCCAT        CCCTGCCTACCTCCAACCC     0      0.2        0
    ## 88  TCCTTGCTGTGGGACATCCAT        CCCCTGCCTACCTCCAACC     0      0.2        0
    ## 89  TCCTTGCTGTGGGACATCCAT         CCCTGTGCACGCCCATAC     0      0.2        0
    ## 90  TCCTTGCTGTGGGACATCCAT      CCTACCTCCAACCCTGTGCAC     0      0.2        0
    ## 91  TCCTTGCTGTGGGACATCCAT        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 92  TCCTTGCTGTGGGACATCCAT     CCTGTGCACGCCCATACTAAAC     0      0.2        0
    ## 93  TCCTTGCTGTGGGACATCCAT      CGCTGTGTAGTGGGGTCAAAC     0      0.2        0
    ## 94  TCCTTGCTGTGGGACATCCAT     CGCCCATACTAAACGCTGTGTA     0      0.2        0
    ## 95  TCCTTGCTGTGGGACATCCAT GGGGTCAAACAGAGATGAATCAGGTA     0      0.2        0
    ## 96  TCCTTGCTGTGGGACATCCAT       CCCTGTGCACGCCCATACTA     0      0.2        0
    ## 97  TCCTTGCTGTGGGACATCCAT          CCCTGTGCACGCCCATA     0      0.2        0
    ## 98  TCCTTGCTGTGGGACATCCAT   TGTAGTGGGGTCAAACAGAGATGA     0      0.2        0
    ## 99  TCCTTGCTGTGGGACATCCAT    TGTGTAGTGGGGTCAAACAGAGA     0      0.2        0
    ## 100 TCCTTGCTGTGGGACATCCAT    GCTGTGTAGTGGGGTCAAACAGA     0      0.2        0
    ##     similar_signatures        missing_signatures enzyme digest_score fragments
    ## 1                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 2                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 3                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 4                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 5                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 6                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 7                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 8                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 9                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 10                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 11                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 12                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 13                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 14                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 15                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 16                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 17                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 18                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 19                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 20                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 21                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 22                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 23                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 24                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 25                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 26                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 27                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 28                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 29                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 30                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 31                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 32                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 33                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 34                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 35                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 36                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 37                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 38                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 39                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 40                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 41                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 42                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 43                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 44                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 45                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 46                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 47                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 48                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 49                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 50                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 51                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 52                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 53                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 54                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 55                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 56                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 57                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 58                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 59                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 60                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 61                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 62                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 63                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 64                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 65                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 66                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 67                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 68                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 69                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 70                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 71                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 72                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 73                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 74                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 75                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 76                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 77                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 78                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 79                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 80                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 81                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 82                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 83                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 84                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 85                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 86                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 87                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 88                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 89                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 90                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 91                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 92                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 93                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 94                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 95                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 96                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 97                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 98                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 99                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 100                    HPV18, HPV16, HPV4, HPV17                   0         0
    ## 
    ## [[6]]
    ##            forward_primer             reverse_primer score coverage products
    ## 1   GTTTGACCCCACTACACAGCG   CACTAACACCAACGCCTAAAGGTT     0      0.2        0
    ## 2   GTTTGACCCCACTACACAGCG    CACTAACACCAACGCCTAAAGGT     0      0.2        0
    ## 3   GTTTGACCCCACTACACAGCG      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 4   GTTTGACCCCACTACACAGCG       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 5   GTTTGACCCCACTACACAGCG        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 6   GTTTGACCCCACTACACAGCG    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 7   GTTTGACCCCACTACACAGCG   ACTAACACCAACGCCTAAAGGTTG     0      0.2        0
    ## 8   GTTTGACCCCACTACACAGCG       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 9   GTTTGACCCCACTACACAGCG    CCACTAACRCCAACRCCTAAAGG     0      0.2        0
    ## 10  GTTTGACCCCACTACACAGCG      GGATGCCCACTAACACCAACG     0      0.2        0
    ## 11  GTTTGACCCCACTACACAGCG    CCCACTAACACCAACGCCTAAAG     0      0.2        0
    ## 12  GTTTGACCCCACTACACAGCG       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 13  GTTTGACCCCACTACACAGCG       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 14  GTTTGACCCCACTACACAGCG        TGCCCACTAACACCAACGC     0      0.2        0
    ## 15  GTTTGACCCCACTACACAGCG         TGACCCCTGCCTACCTCC     0      0.2        0
    ## 16  GTTTGACCCCACTACACAGCG       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 17  GTTTGACCCCACTACACAGCG      CAACGCCTAAAGGTTGACCCC     0      0.2        0
    ## 18  GTTTGACCCCACTACACAGCG      CCAACGCCTAAAGGTTGACCC     0      0.2        0
    ## 19  GTTTGACCCCACTACACAGCG        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 20  GTTTGACCCCACTACACAGCG      ACCAACGCCTAAAGGTTGACC     0      0.2        0
    ## 21  GTTTGACCCCACTACACAGCG     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 22  GTTTGACCCCACTACACAGCG    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 23  GTTTGACCCCACTACACAGCG     ACACCAACGCCTAAAGGTTGAC     0      0.2        0
    ## 24  GTTTGACCCCACTACACAGCG     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 25  GTTTGACCCCACTACACAGCG  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 26  GTTTGACCCCACTACACAGCG        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 27  GTTTGACCCCACTACACAGCG     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 28  GTTTGACCCCACTACACAGCG     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 29  GTTTGACCCCACTACACAGCG      GCCCACTAACACCAACGCCTA     0      0.2        0
    ## 30  GTTTGACCCCACTACACAGCG   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 31  GTTTGACCCCACTACACAGCG   CTAACACCAACGCCTAAAGGTTGA     0      0.2        0
    ## 32  GTTTGACCCCACTACACAGCG         CACCYCTRCCTACCTCCA     0      0.2        0
    ## 33  GTTTGACCCCACTACACAGCG     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 34  GTTTGACCCCACTACACAGCG  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 35  GTTTGACCCCACTACACAGCG     GCCCACTAACACCAACGCCTAA     0      0.2        0
    ## 36  GTTTGACCCCACTACACAGCG  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 37  GTTTGACCCCACTACACAGCG        GACCCCTGCCTACCTCCAA     0      0.2        0
    ## 38  GTTTGACCCCACTACACAGCG    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 39  GTTTGACCCCACTACACAGCG     CCCACTAACACCAACGCCTAAA     0      0.2        0
    ## 40   TGACCCCACTACACAGCGTT      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 41   TGACCCCACTACACAGCGTT       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 42   TGACCCCACTACACAGCGTT        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 43   TGACCCCACTACACAGCGTT    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 44   TGACCCCACTACACAGCGTT       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 45   TGACCCCACTACACAGCGTT       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 46   TGACCCCACTACACAGCGTT       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 47   TGACCCCACTACACAGCGTT       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 48   TGACCCCACTACACAGCGTT        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 49   TGACCCCACTACACAGCGTT     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 50   TGACCCCACTACACAGCGTT    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 51   TGACCCCACTACACAGCGTT     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 52   TGACCCCACTACACAGCGTT  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 53   TGACCCCACTACACAGCGTT     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 54   TGACCCCACTACACAGCGTT     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 55   TGACCCCACTACACAGCGTT   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 56   TGACCCCACTACACAGCGTT     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 57   TGACCCCACTACACAGCGTT  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 58   TGACCCCACTACACAGCGTT  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 59   TGACCCCACTACACAGCGTT    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 60  TCCTTGCTGTGGGACATCCAT      CGCCCATACTAAACGCTGTGT     0      0.2        0
    ## 61  TCCTTGCTGTGGGACATCCAT      CACGCCCATACTAAACGCTGT     0      0.2        0
    ## 62  TCCTTGCTGTGGGACATCCAT        TGCCTACCTCCAACCCTGT     0      0.2        0
    ## 63  TCCTTGCTGTGGGACATCCAT  GGGGTCAAACAGAGATGAATCAGGT     0      0.2        0
    ## 64  TCCTTGCTGTGGGACATCCAT   CGCCCATACTAAACGCTGTGTAGT     0      0.2        0
    ## 65  TCCTTGCTGTGGGACATCCAT       GCACGCCCATACTAAACGCT     0      0.2        0
    ## 66  TCCTTGCTGTGGGACATCCAT        CCTGCCTACCTCCAACCCT     0      0.2        0
    ## 67  TCCTTGCTGTGGGACATCCAT        CCCTGTGCACGCCCATACT     0      0.2        0
    ## 68  TCCTTGCTGTGGGACATCCAT   TGTGTAGTGGGGTCAAACAGAGAT     0      0.2        0
    ## 69  TCCTTGCTGTGGGACATCCAT TGTAGTGGGGTCAAACAGAGATGAAT     0      0.2        0
    ## 70  TCCTTGCTGTGGGACATCCAT      ACGCCCATACTAAACGCTGTG     0      0.2        0
    ## 71  TCCTTGCTGTGGGACATCCAT       TGCCTACCTCCAACCCTGTG     0      0.2        0
    ## 72  TCCTTGCTGTGGGACATCCAT   GCCCATACTAAACGCTGTGTAGTG     0      0.2        0
    ## 73  TCCTTGCTGTGGGACATCCAT      GCACGCCCATACTAAACGCTG     0      0.2        0
    ## 74  TCCTTGCTGTGGGACATCCAT       CCTGCCTACCTCCAACCCTG     0      0.2        0
    ## 75  TCCTTGCTGTGGGACATCCAT  TGTGTAGTGGGGTCAAACAGAGATG     0      0.2        0
    ## 76  TCCTTGCTGTGGGACATCCAT   CATACTAAACGCTGTGTAGTGGGG     0      0.2        0
    ## 77  TCCTTGCTGTGGGACATCCAT   GGGGTCAAACAGAGATGAATCAGG     0      0.2        0
    ## 78  TCCTTGCTGTGGGACATCCAT         CCTCCAACCCTGTGCACG     0      0.2        0
    ## 79  TCCTTGCTGTGGGACATCCAT    CGCCCATACTAAACGCTGTGTAG     0      0.2        0
    ## 80  TCCTTGCTGTGGGACATCCAT   GCTGTGTAGTGGGGTCAAACAGAG     0      0.2        0
    ## 81  TCCTTGCTGTGGGACATCCAT  GTGGGGTCAAACAGAGATGAATCAG     0      0.2        0
    ## 82  TCCTTGCTGTGGGACATCCAT     GCTGTGTAGTGGGGTCAAACAG     0      0.2        0
    ## 83  TCCTTGCTGTGGGACATCCAT     CTAAACGCTGTGTAGTGGGGTC     0      0.2        0
    ## 84  TCCTTGCTGTGGGACATCCAT       GCCTACCTCCAACCCTGTGC     0      0.2        0
    ## 85  TCCTTGCTGTGGGACATCCAT          TCCAACCCTGTGCACGC     0      0.2        0
    ## 86  TCCTTGCTGTGGGACATCCAT          CCAACCCTGTGCACGCC     0      0.2        0
    ## 87  TCCTTGCTGTGGGACATCCAT        CCCTGCCTACCTCCAACCC     0      0.2        0
    ## 88  TCCTTGCTGTGGGACATCCAT        CCCCTGCCTACCTCCAACC     0      0.2        0
    ## 89  TCCTTGCTGTGGGACATCCAT         CCCTGTGCACGCCCATAC     0      0.2        0
    ## 90  TCCTTGCTGTGGGACATCCAT      CCTACCTCCAACCCTGTGCAC     0      0.2        0
    ## 91  TCCTTGCTGTGGGACATCCAT        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 92  TCCTTGCTGTGGGACATCCAT     CCTGTGCACGCCCATACTAAAC     0      0.2        0
    ## 93  TCCTTGCTGTGGGACATCCAT      CGCTGTGTAGTGGGGTCAAAC     0      0.2        0
    ## 94  TCCTTGCTGTGGGACATCCAT     CGCCCATACTAAACGCTGTGTA     0      0.2        0
    ## 95  TCCTTGCTGTGGGACATCCAT GGGGTCAAACAGAGATGAATCAGGTA     0      0.2        0
    ## 96  TCCTTGCTGTGGGACATCCAT       CCCTGTGCACGCCCATACTA     0      0.2        0
    ## 97  TCCTTGCTGTGGGACATCCAT          CCCTGTGCACGCCCATA     0      0.2        0
    ## 98  TCCTTGCTGTGGGACATCCAT   TGTAGTGGGGTCAAACAGAGATGA     0      0.2        0
    ## 99  TCCTTGCTGTGGGACATCCAT    TGTGTAGTGGGGTCAAACAGAGA     0      0.2        0
    ## 100 TCCTTGCTGTGGGACATCCAT    GCTGTGTAGTGGGGTCAAACAGA     0      0.2        0
    ##     similar_signatures        missing_signatures enzyme digest_score fragments
    ## 1                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 2                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 3                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 4                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 5                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 6                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 7                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 8                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 9                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 10                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 11                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 12                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 13                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 14                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 15                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 16                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 17                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 18                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 19                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 20                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 21                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 22                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 23                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 24                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 25                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 26                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 27                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 28                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 29                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 30                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 31                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 32                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 33                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 34                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 35                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 36                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 37                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 38                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 39                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 40                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 41                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 42                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 43                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 44                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 45                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 46                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 47                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 48                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 49                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 50                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 51                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 52                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 53                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 54                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 55                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 56                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 57                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 58                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 59                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 60                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 61                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 62                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 63                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 64                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 65                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 66                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 67                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 68                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 69                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 70                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 71                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 72                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 73                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 74                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 75                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 76                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 77                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 78                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 79                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 80                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 81                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 82                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 83                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 84                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 85                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 86                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 87                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 88                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 89                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 90                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 91                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 92                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 93                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 94                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 95                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 96                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 97                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 98                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 99                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 100                    HPV18, HPV16, HPV4, HPV17                   0         0
    ## 
    ## [[7]]
    ##            forward_primer             reverse_primer score coverage products
    ## 1   GTTTGACCCCACTACACAGCG       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 2    TGACCCCACTACACAGCGTT       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 3   GTTTGACCCCACTACACAGCG   CACTAACACCAACGCCTAAAGGTT     0      0.2        0
    ## 4   GTTTGACCCCACTACACAGCG    CACTAACACCAACGCCTAAAGGT     0      0.2        0
    ## 5   GTTTGACCCCACTACACAGCG      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 6   GTTTGACCCCACTACACAGCG       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 7   GTTTGACCCCACTACACAGCG        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 8   GTTTGACCCCACTACACAGCG    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 9   GTTTGACCCCACTACACAGCG   ACTAACACCAACGCCTAAAGGTTG     0      0.2        0
    ## 10  GTTTGACCCCACTACACAGCG       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 11  GTTTGACCCCACTACACAGCG    CCACTAACRCCAACRCCTAAAGG     0      0.2        0
    ## 12  GTTTGACCCCACTACACAGCG      GGATGCCCACTAACACCAACG     0      0.2        0
    ## 13  GTTTGACCCCACTACACAGCG    CCCACTAACACCAACGCCTAAAG     0      0.2        0
    ## 14  GTTTGACCCCACTACACAGCG       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 15  GTTTGACCCCACTACACAGCG       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 16  GTTTGACCCCACTACACAGCG        TGCCCACTAACACCAACGC     0      0.2        0
    ## 17  GTTTGACCCCACTACACAGCG         TGACCCCTGCCTACCTCC     0      0.2        0
    ## 18  GTTTGACCCCACTACACAGCG      CAACGCCTAAAGGTTGACCCC     0      0.2        0
    ## 19  GTTTGACCCCACTACACAGCG      CCAACGCCTAAAGGTTGACCC     0      0.2        0
    ## 20  GTTTGACCCCACTACACAGCG        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 21  GTTTGACCCCACTACACAGCG      ACCAACGCCTAAAGGTTGACC     0      0.2        0
    ## 22  GTTTGACCCCACTACACAGCG     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 23  GTTTGACCCCACTACACAGCG    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 24  GTTTGACCCCACTACACAGCG     ACACCAACGCCTAAAGGTTGAC     0      0.2        0
    ## 25  GTTTGACCCCACTACACAGCG     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 26  GTTTGACCCCACTACACAGCG  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 27  GTTTGACCCCACTACACAGCG        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 28  GTTTGACCCCACTACACAGCG     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 29  GTTTGACCCCACTACACAGCG     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 30  GTTTGACCCCACTACACAGCG      GCCCACTAACACCAACGCCTA     0      0.2        0
    ## 31  GTTTGACCCCACTACACAGCG   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 32  GTTTGACCCCACTACACAGCG   CTAACACCAACGCCTAAAGGTTGA     0      0.2        0
    ## 33  GTTTGACCCCACTACACAGCG         CACCYCTRCCTACCTCCA     0      0.2        0
    ## 34  GTTTGACCCCACTACACAGCG     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 35  GTTTGACCCCACTACACAGCG  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 36  GTTTGACCCCACTACACAGCG     GCCCACTAACACCAACGCCTAA     0      0.2        0
    ## 37  GTTTGACCCCACTACACAGCG  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 38  GTTTGACCCCACTACACAGCG        GACCCCTGCCTACCTCCAA     0      0.2        0
    ## 39  GTTTGACCCCACTACACAGCG    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 40  GTTTGACCCCACTACACAGCG     CCCACTAACACCAACGCCTAAA     0      0.2        0
    ## 41   TGACCCCACTACACAGCGTT      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 42   TGACCCCACTACACAGCGTT       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 43   TGACCCCACTACACAGCGTT        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 44   TGACCCCACTACACAGCGTT    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 45   TGACCCCACTACACAGCGTT       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 46   TGACCCCACTACACAGCGTT       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 47   TGACCCCACTACACAGCGTT       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 48   TGACCCCACTACACAGCGTT        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 49   TGACCCCACTACACAGCGTT     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 50   TGACCCCACTACACAGCGTT    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 51   TGACCCCACTACACAGCGTT     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 52   TGACCCCACTACACAGCGTT  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 53   TGACCCCACTACACAGCGTT     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 54   TGACCCCACTACACAGCGTT     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 55   TGACCCCACTACACAGCGTT   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 56   TGACCCCACTACACAGCGTT     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 57   TGACCCCACTACACAGCGTT  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 58   TGACCCCACTACACAGCGTT  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 59   TGACCCCACTACACAGCGTT    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 60  TCCTTGCTGTGGGACATCCAT      CGCCCATACTAAACGCTGTGT     0      0.2        0
    ## 61  TCCTTGCTGTGGGACATCCAT      CACGCCCATACTAAACGCTGT     0      0.2        0
    ## 62  TCCTTGCTGTGGGACATCCAT        TGCCTACCTCCAACCCTGT     0      0.2        0
    ## 63  TCCTTGCTGTGGGACATCCAT  GGGGTCAAACAGAGATGAATCAGGT     0      0.2        0
    ## 64  TCCTTGCTGTGGGACATCCAT   CGCCCATACTAAACGCTGTGTAGT     0      0.2        0
    ## 65  TCCTTGCTGTGGGACATCCAT       GCACGCCCATACTAAACGCT     0      0.2        0
    ## 66  TCCTTGCTGTGGGACATCCAT        CCTGCCTACCTCCAACCCT     0      0.2        0
    ## 67  TCCTTGCTGTGGGACATCCAT        CCCTGTGCACGCCCATACT     0      0.2        0
    ## 68  TCCTTGCTGTGGGACATCCAT   TGTGTAGTGGGGTCAAACAGAGAT     0      0.2        0
    ## 69  TCCTTGCTGTGGGACATCCAT TGTAGTGGGGTCAAACAGAGATGAAT     0      0.2        0
    ## 70  TCCTTGCTGTGGGACATCCAT      ACGCCCATACTAAACGCTGTG     0      0.2        0
    ## 71  TCCTTGCTGTGGGACATCCAT       TGCCTACCTCCAACCCTGTG     0      0.2        0
    ## 72  TCCTTGCTGTGGGACATCCAT   GCCCATACTAAACGCTGTGTAGTG     0      0.2        0
    ## 73  TCCTTGCTGTGGGACATCCAT      GCACGCCCATACTAAACGCTG     0      0.2        0
    ## 74  TCCTTGCTGTGGGACATCCAT       CCTGCCTACCTCCAACCCTG     0      0.2        0
    ## 75  TCCTTGCTGTGGGACATCCAT  TGTGTAGTGGGGTCAAACAGAGATG     0      0.2        0
    ## 76  TCCTTGCTGTGGGACATCCAT   CATACTAAACGCTGTGTAGTGGGG     0      0.2        0
    ## 77  TCCTTGCTGTGGGACATCCAT   GGGGTCAAACAGAGATGAATCAGG     0      0.2        0
    ## 78  TCCTTGCTGTGGGACATCCAT         CCTCCAACCCTGTGCACG     0      0.2        0
    ## 79  TCCTTGCTGTGGGACATCCAT    CGCCCATACTAAACGCTGTGTAG     0      0.2        0
    ## 80  TCCTTGCTGTGGGACATCCAT   GCTGTGTAGTGGGGTCAAACAGAG     0      0.2        0
    ## 81  TCCTTGCTGTGGGACATCCAT  GTGGGGTCAAACAGAGATGAATCAG     0      0.2        0
    ## 82  TCCTTGCTGTGGGACATCCAT     GCTGTGTAGTGGGGTCAAACAG     0      0.2        0
    ## 83  TCCTTGCTGTGGGACATCCAT     CTAAACGCTGTGTAGTGGGGTC     0      0.2        0
    ## 84  TCCTTGCTGTGGGACATCCAT       GCCTACCTCCAACCCTGTGC     0      0.2        0
    ## 85  TCCTTGCTGTGGGACATCCAT          TCCAACCCTGTGCACGC     0      0.2        0
    ## 86  TCCTTGCTGTGGGACATCCAT          CCAACCCTGTGCACGCC     0      0.2        0
    ## 87  TCCTTGCTGTGGGACATCCAT        CCCTGCCTACCTCCAACCC     0      0.2        0
    ## 88  TCCTTGCTGTGGGACATCCAT        CCCCTGCCTACCTCCAACC     0      0.2        0
    ## 89  TCCTTGCTGTGGGACATCCAT         CCCTGTGCACGCCCATAC     0      0.2        0
    ## 90  TCCTTGCTGTGGGACATCCAT      CCTACCTCCAACCCTGTGCAC     0      0.2        0
    ## 91  TCCTTGCTGTGGGACATCCAT        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 92  TCCTTGCTGTGGGACATCCAT     CCTGTGCACGCCCATACTAAAC     0      0.2        0
    ## 93  TCCTTGCTGTGGGACATCCAT      CGCTGTGTAGTGGGGTCAAAC     0      0.2        0
    ## 94  TCCTTGCTGTGGGACATCCAT     CGCCCATACTAAACGCTGTGTA     0      0.2        0
    ## 95  TCCTTGCTGTGGGACATCCAT GGGGTCAAACAGAGATGAATCAGGTA     0      0.2        0
    ## 96  TCCTTGCTGTGGGACATCCAT       CCCTGTGCACGCCCATACTA     0      0.2        0
    ## 97  TCCTTGCTGTGGGACATCCAT          CCCTGTGCACGCCCATA     0      0.2        0
    ## 98  TCCTTGCTGTGGGACATCCAT   TGTAGTGGGGTCAAACAGAGATGA     0      0.2        0
    ## 99  TCCTTGCTGTGGGACATCCAT    TGTGTAGTGGGGTCAAACAGAGA     0      0.2        0
    ## 100 TCCTTGCTGTGGGACATCCAT    GCTGTGTAGTGGGGTCAAACAGA     0      0.2        0
    ##     similar_signatures        missing_signatures enzyme digest_score fragments
    ## 1                      HPV18, HPV16, HPV4, HPV17 Acc65I            0         2
    ## 2                      HPV18, HPV16, HPV4, HPV17 Acc65I            0         2
    ## 3                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 4                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 5                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 6                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 7                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 8                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 9                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 10                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 11                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 12                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 13                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 14                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 15                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 16                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 17                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 18                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 19                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 20                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 21                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 22                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 23                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 24                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 25                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 26                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 27                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 28                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 29                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 30                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 31                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 32                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 33                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 34                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 35                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 36                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 37                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 38                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 39                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 40                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 41                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 42                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 43                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 44                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 45                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 46                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 47                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 48                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 49                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 50                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 51                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 52                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 53                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 54                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 55                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 56                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 57                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 58                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 59                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 60                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 61                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 62                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 63                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 64                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 65                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 66                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 67                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 68                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 69                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 70                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 71                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 72                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 73                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 74                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 75                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 76                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 77                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 78                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 79                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 80                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 81                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 82                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 83                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 84                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 85                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 86                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 87                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 88                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 89                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 90                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 91                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 92                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 93                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 94                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 95                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 96                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 97                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 98                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 99                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 100                    HPV18, HPV16, HPV4, HPV17                   0         0
    ## 
    ## [[8]]
    ##            forward_primer             reverse_primer score coverage products
    ## 1   GTTTGACCCCACTACACAGCG   CACTAACACCAACGCCTAAAGGTT     0      0.2        0
    ## 2   GTTTGACCCCACTACACAGCG    CACTAACACCAACGCCTAAAGGT     0      0.2        0
    ## 3   GTTTGACCCCACTACACAGCG      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 4   GTTTGACCCCACTACACAGCG       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 5   GTTTGACCCCACTACACAGCG        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 6   GTTTGACCCCACTACACAGCG    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 7   GTTTGACCCCACTACACAGCG   ACTAACACCAACGCCTAAAGGTTG     0      0.2        0
    ## 8   GTTTGACCCCACTACACAGCG       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 9   GTTTGACCCCACTACACAGCG    CCACTAACRCCAACRCCTAAAGG     0      0.2        0
    ## 10  GTTTGACCCCACTACACAGCG      GGATGCCCACTAACACCAACG     0      0.2        0
    ## 11  GTTTGACCCCACTACACAGCG    CCCACTAACACCAACGCCTAAAG     0      0.2        0
    ## 12  GTTTGACCCCACTACACAGCG       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 13  GTTTGACCCCACTACACAGCG       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 14  GTTTGACCCCACTACACAGCG        TGCCCACTAACACCAACGC     0      0.2        0
    ## 15  GTTTGACCCCACTACACAGCG         TGACCCCTGCCTACCTCC     0      0.2        0
    ## 16  GTTTGACCCCACTACACAGCG       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 17  GTTTGACCCCACTACACAGCG      CAACGCCTAAAGGTTGACCCC     0      0.2        0
    ## 18  GTTTGACCCCACTACACAGCG      CCAACGCCTAAAGGTTGACCC     0      0.2        0
    ## 19  GTTTGACCCCACTACACAGCG        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 20  GTTTGACCCCACTACACAGCG      ACCAACGCCTAAAGGTTGACC     0      0.2        0
    ## 21  GTTTGACCCCACTACACAGCG     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 22  GTTTGACCCCACTACACAGCG    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 23  GTTTGACCCCACTACACAGCG     ACACCAACGCCTAAAGGTTGAC     0      0.2        0
    ## 24  GTTTGACCCCACTACACAGCG     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 25  GTTTGACCCCACTACACAGCG  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 26  GTTTGACCCCACTACACAGCG        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 27  GTTTGACCCCACTACACAGCG     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 28  GTTTGACCCCACTACACAGCG     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 29  GTTTGACCCCACTACACAGCG      GCCCACTAACACCAACGCCTA     0      0.2        0
    ## 30  GTTTGACCCCACTACACAGCG   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 31  GTTTGACCCCACTACACAGCG   CTAACACCAACGCCTAAAGGTTGA     0      0.2        0
    ## 32  GTTTGACCCCACTACACAGCG         CACCYCTRCCTACCTCCA     0      0.2        0
    ## 33  GTTTGACCCCACTACACAGCG     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 34  GTTTGACCCCACTACACAGCG  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 35  GTTTGACCCCACTACACAGCG     GCCCACTAACACCAACGCCTAA     0      0.2        0
    ## 36  GTTTGACCCCACTACACAGCG  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 37  GTTTGACCCCACTACACAGCG        GACCCCTGCCTACCTCCAA     0      0.2        0
    ## 38  GTTTGACCCCACTACACAGCG    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 39  GTTTGACCCCACTACACAGCG     CCCACTAACACCAACGCCTAAA     0      0.2        0
    ## 40   TGACCCCACTACACAGCGTT      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 41   TGACCCCACTACACAGCGTT       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 42   TGACCCCACTACACAGCGTT        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 43   TGACCCCACTACACAGCGTT    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 44   TGACCCCACTACACAGCGTT       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 45   TGACCCCACTACACAGCGTT       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 46   TGACCCCACTACACAGCGTT       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 47   TGACCCCACTACACAGCGTT       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 48   TGACCCCACTACACAGCGTT        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 49   TGACCCCACTACACAGCGTT     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 50   TGACCCCACTACACAGCGTT    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 51   TGACCCCACTACACAGCGTT     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 52   TGACCCCACTACACAGCGTT  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 53   TGACCCCACTACACAGCGTT     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 54   TGACCCCACTACACAGCGTT     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 55   TGACCCCACTACACAGCGTT   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 56   TGACCCCACTACACAGCGTT     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 57   TGACCCCACTACACAGCGTT  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 58   TGACCCCACTACACAGCGTT  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 59   TGACCCCACTACACAGCGTT    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 60  TCCTTGCTGTGGGACATCCAT      CGCCCATACTAAACGCTGTGT     0      0.2        0
    ## 61  TCCTTGCTGTGGGACATCCAT      CACGCCCATACTAAACGCTGT     0      0.2        0
    ## 62  TCCTTGCTGTGGGACATCCAT        TGCCTACCTCCAACCCTGT     0      0.2        0
    ## 63  TCCTTGCTGTGGGACATCCAT  GGGGTCAAACAGAGATGAATCAGGT     0      0.2        0
    ## 64  TCCTTGCTGTGGGACATCCAT   CGCCCATACTAAACGCTGTGTAGT     0      0.2        0
    ## 65  TCCTTGCTGTGGGACATCCAT       GCACGCCCATACTAAACGCT     0      0.2        0
    ## 66  TCCTTGCTGTGGGACATCCAT        CCTGCCTACCTCCAACCCT     0      0.2        0
    ## 67  TCCTTGCTGTGGGACATCCAT        CCCTGTGCACGCCCATACT     0      0.2        0
    ## 68  TCCTTGCTGTGGGACATCCAT   TGTGTAGTGGGGTCAAACAGAGAT     0      0.2        0
    ## 69  TCCTTGCTGTGGGACATCCAT TGTAGTGGGGTCAAACAGAGATGAAT     0      0.2        0
    ## 70  TCCTTGCTGTGGGACATCCAT      ACGCCCATACTAAACGCTGTG     0      0.2        0
    ## 71  TCCTTGCTGTGGGACATCCAT       TGCCTACCTCCAACCCTGTG     0      0.2        0
    ## 72  TCCTTGCTGTGGGACATCCAT   GCCCATACTAAACGCTGTGTAGTG     0      0.2        0
    ## 73  TCCTTGCTGTGGGACATCCAT      GCACGCCCATACTAAACGCTG     0      0.2        0
    ## 74  TCCTTGCTGTGGGACATCCAT       CCTGCCTACCTCCAACCCTG     0      0.2        0
    ## 75  TCCTTGCTGTGGGACATCCAT  TGTGTAGTGGGGTCAAACAGAGATG     0      0.2        0
    ## 76  TCCTTGCTGTGGGACATCCAT   CATACTAAACGCTGTGTAGTGGGG     0      0.2        0
    ## 77  TCCTTGCTGTGGGACATCCAT   GGGGTCAAACAGAGATGAATCAGG     0      0.2        0
    ## 78  TCCTTGCTGTGGGACATCCAT         CCTCCAACCCTGTGCACG     0      0.2        0
    ## 79  TCCTTGCTGTGGGACATCCAT    CGCCCATACTAAACGCTGTGTAG     0      0.2        0
    ## 80  TCCTTGCTGTGGGACATCCAT   GCTGTGTAGTGGGGTCAAACAGAG     0      0.2        0
    ## 81  TCCTTGCTGTGGGACATCCAT  GTGGGGTCAAACAGAGATGAATCAG     0      0.2        0
    ## 82  TCCTTGCTGTGGGACATCCAT     GCTGTGTAGTGGGGTCAAACAG     0      0.2        0
    ## 83  TCCTTGCTGTGGGACATCCAT     CTAAACGCTGTGTAGTGGGGTC     0      0.2        0
    ## 84  TCCTTGCTGTGGGACATCCAT       GCCTACCTCCAACCCTGTGC     0      0.2        0
    ## 85  TCCTTGCTGTGGGACATCCAT          TCCAACCCTGTGCACGC     0      0.2        0
    ## 86  TCCTTGCTGTGGGACATCCAT          CCAACCCTGTGCACGCC     0      0.2        0
    ## 87  TCCTTGCTGTGGGACATCCAT        CCCTGCCTACCTCCAACCC     0      0.2        0
    ## 88  TCCTTGCTGTGGGACATCCAT        CCCCTGCCTACCTCCAACC     0      0.2        0
    ## 89  TCCTTGCTGTGGGACATCCAT         CCCTGTGCACGCCCATAC     0      0.2        0
    ## 90  TCCTTGCTGTGGGACATCCAT      CCTACCTCCAACCCTGTGCAC     0      0.2        0
    ## 91  TCCTTGCTGTGGGACATCCAT        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 92  TCCTTGCTGTGGGACATCCAT     CCTGTGCACGCCCATACTAAAC     0      0.2        0
    ## 93  TCCTTGCTGTGGGACATCCAT      CGCTGTGTAGTGGGGTCAAAC     0      0.2        0
    ## 94  TCCTTGCTGTGGGACATCCAT     CGCCCATACTAAACGCTGTGTA     0      0.2        0
    ## 95  TCCTTGCTGTGGGACATCCAT GGGGTCAAACAGAGATGAATCAGGTA     0      0.2        0
    ## 96  TCCTTGCTGTGGGACATCCAT       CCCTGTGCACGCCCATACTA     0      0.2        0
    ## 97  TCCTTGCTGTGGGACATCCAT          CCCTGTGCACGCCCATA     0      0.2        0
    ## 98  TCCTTGCTGTGGGACATCCAT   TGTAGTGGGGTCAAACAGAGATGA     0      0.2        0
    ## 99  TCCTTGCTGTGGGACATCCAT    TGTGTAGTGGGGTCAAACAGAGA     0      0.2        0
    ## 100 TCCTTGCTGTGGGACATCCAT    GCTGTGTAGTGGGGTCAAACAGA     0      0.2        0
    ##     similar_signatures        missing_signatures enzyme digest_score fragments
    ## 1                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 2                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 3                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 4                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 5                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 6                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 7                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 8                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 9                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 10                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 11                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 12                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 13                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 14                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 15                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 16                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 17                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 18                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 19                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 20                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 21                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 22                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 23                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 24                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 25                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 26                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 27                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 28                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 29                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 30                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 31                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 32                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 33                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 34                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 35                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 36                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 37                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 38                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 39                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 40                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 41                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 42                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 43                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 44                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 45                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 46                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 47                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 48                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 49                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 50                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 51                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 52                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 53                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 54                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 55                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 56                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 57                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 58                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 59                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 60                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 61                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 62                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 63                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 64                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 65                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 66                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 67                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 68                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 69                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 70                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 71                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 72                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 73                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 74                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 75                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 76                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 77                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 78                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 79                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 80                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 81                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 82                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 83                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 84                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 85                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 86                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 87                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 88                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 89                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 90                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 91                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 92                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 93                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 94                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 95                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 96                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 97                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 98                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 99                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 100                    HPV18, HPV16, HPV4, HPV17                   0         0
    ## 
    ## [[9]]
    ##            forward_primer             reverse_primer score coverage products
    ## 1   GTTTGACCCCACTACACAGCG   CACTAACACCAACGCCTAAAGGTT     0      0.2        0
    ## 2   GTTTGACCCCACTACACAGCG    CACTAACACCAACGCCTAAAGGT     0      0.2        0
    ## 3   GTTTGACCCCACTACACAGCG      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 4   GTTTGACCCCACTACACAGCG       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 5   GTTTGACCCCACTACACAGCG        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 6   GTTTGACCCCACTACACAGCG    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 7   GTTTGACCCCACTACACAGCG   ACTAACACCAACGCCTAAAGGTTG     0      0.2        0
    ## 8   GTTTGACCCCACTACACAGCG       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 9   GTTTGACCCCACTACACAGCG    CCACTAACRCCAACRCCTAAAGG     0      0.2        0
    ## 10  GTTTGACCCCACTACACAGCG      GGATGCCCACTAACACCAACG     0      0.2        0
    ## 11  GTTTGACCCCACTACACAGCG    CCCACTAACACCAACGCCTAAAG     0      0.2        0
    ## 12  GTTTGACCCCACTACACAGCG       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 13  GTTTGACCCCACTACACAGCG       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 14  GTTTGACCCCACTACACAGCG        TGCCCACTAACACCAACGC     0      0.2        0
    ## 15  GTTTGACCCCACTACACAGCG         TGACCCCTGCCTACCTCC     0      0.2        0
    ## 16  GTTTGACCCCACTACACAGCG       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 17  GTTTGACCCCACTACACAGCG      CAACGCCTAAAGGTTGACCCC     0      0.2        0
    ## 18  GTTTGACCCCACTACACAGCG      CCAACGCCTAAAGGTTGACCC     0      0.2        0
    ## 19  GTTTGACCCCACTACACAGCG        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 20  GTTTGACCCCACTACACAGCG      ACCAACGCCTAAAGGTTGACC     0      0.2        0
    ## 21  GTTTGACCCCACTACACAGCG     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 22  GTTTGACCCCACTACACAGCG    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 23  GTTTGACCCCACTACACAGCG     ACACCAACGCCTAAAGGTTGAC     0      0.2        0
    ## 24  GTTTGACCCCACTACACAGCG     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 25  GTTTGACCCCACTACACAGCG  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 26  GTTTGACCCCACTACACAGCG        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 27  GTTTGACCCCACTACACAGCG     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 28  GTTTGACCCCACTACACAGCG     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 29  GTTTGACCCCACTACACAGCG      GCCCACTAACACCAACGCCTA     0      0.2        0
    ## 30  GTTTGACCCCACTACACAGCG   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 31  GTTTGACCCCACTACACAGCG   CTAACACCAACGCCTAAAGGTTGA     0      0.2        0
    ## 32  GTTTGACCCCACTACACAGCG         CACCYCTRCCTACCTCCA     0      0.2        0
    ## 33  GTTTGACCCCACTACACAGCG     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 34  GTTTGACCCCACTACACAGCG  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 35  GTTTGACCCCACTACACAGCG     GCCCACTAACACCAACGCCTAA     0      0.2        0
    ## 36  GTTTGACCCCACTACACAGCG  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 37  GTTTGACCCCACTACACAGCG        GACCCCTGCCTACCTCCAA     0      0.2        0
    ## 38  GTTTGACCCCACTACACAGCG    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 39  GTTTGACCCCACTACACAGCG     CCCACTAACACCAACGCCTAAA     0      0.2        0
    ## 40   TGACCCCACTACACAGCGTT      CCTAAAGGTTGACCCCTGCCT     0      0.2        0
    ## 41   TGACCCCACTACACAGCGTT       ACGCCTAAAGGTTGACCCCT     0      0.2        0
    ## 42   TGACCCCACTACACAGCGTT        GGCTGACCMCKGCCTACCT     0      0.2        0
    ## 43   TGACCCCACTACACAGCGTT    TTGTTTAGCAATGGATGCCCACT     0      0.2        0
    ## 44   TGACCCCACTACACAGCGTT       CGCCTAAAGGTTGACCCCTG     0      0.2        0
    ## 45   TGACCCCACTACACAGCGTT       GGCTGACCMCKGCCTACCTC     0      0.2        0
    ## 46   TGACCCCACTACACAGCGTT       GCCTAAAGGTTGACCCCTGC     0      0.2        0
    ## 47   TGACCCCACTACACAGCGTT       CCTAAAGGTTGACCCCTGCC     0      0.2        0
    ## 48   TGACCCCACTACACAGCGTT        AGGTTGACCCCTGCCTACC     0      0.2        0
    ## 49   TGACCCCACTACACAGCGTT     GCAATGGATGCCCACTAACACC     0      0.2        0
    ## 50   TGACCCCACTACACAGCGTT    CCTAAAGGTTGACCCCTGCCTAC     0      0.2        0
    ## 51   TGACCCCACTACACAGCGTT     AGCAATGGATGCCCACTAACAC     0      0.2        0
    ## 52   TGACCCCACTACACAGCGTT  TGTTTAGCAATGGATGCCCACTAAC     0      0.2        0
    ## 53   TGACCCCACTACACAGCGTT     ATGGATGCCCACTAACACCAAC     0      0.2        0
    ## 54   TGACCCCACTACACAGCGTT     CCTAAAGGTTGACCCCTGCCTA     0      0.2        0
    ## 55   TGACCCCACTACACAGCGTT   TTGTTTAGCAATGGATGCCCACTA     0      0.2        0
    ## 56   TGACCCCACTACACAGCGTT     CAATGGATGCCCACTAACACCA     0      0.2        0
    ## 57   TGACCCCACTACACAGCGTT  GTTTAGCAATGGATGCCCACTAACA     0      0.2        0
    ## 58   TGACCCCACTACACAGCGTT  TTGTTTAGCAATGGATGCCCACTAA     0      0.2        0
    ## 59   TGACCCCACTACACAGCGTT    CAATGGATGCCCACTAACACCAA     0      0.2        0
    ## 60  TCCTTGCTGTGGGACATCCAT      CGCCCATACTAAACGCTGTGT     0      0.2        0
    ## 61  TCCTTGCTGTGGGACATCCAT      CACGCCCATACTAAACGCTGT     0      0.2        0
    ## 62  TCCTTGCTGTGGGACATCCAT        TGCCTACCTCCAACCCTGT     0      0.2        0
    ## 63  TCCTTGCTGTGGGACATCCAT  GGGGTCAAACAGAGATGAATCAGGT     0      0.2        0
    ## 64  TCCTTGCTGTGGGACATCCAT   CGCCCATACTAAACGCTGTGTAGT     0      0.2        0
    ## 65  TCCTTGCTGTGGGACATCCAT       GCACGCCCATACTAAACGCT     0      0.2        0
    ## 66  TCCTTGCTGTGGGACATCCAT        CCTGCCTACCTCCAACCCT     0      0.2        0
    ## 67  TCCTTGCTGTGGGACATCCAT        CCCTGTGCACGCCCATACT     0      0.2        0
    ## 68  TCCTTGCTGTGGGACATCCAT   TGTGTAGTGGGGTCAAACAGAGAT     0      0.2        0
    ## 69  TCCTTGCTGTGGGACATCCAT TGTAGTGGGGTCAAACAGAGATGAAT     0      0.2        0
    ## 70  TCCTTGCTGTGGGACATCCAT      ACGCCCATACTAAACGCTGTG     0      0.2        0
    ## 71  TCCTTGCTGTGGGACATCCAT       TGCCTACCTCCAACCCTGTG     0      0.2        0
    ## 72  TCCTTGCTGTGGGACATCCAT   GCCCATACTAAACGCTGTGTAGTG     0      0.2        0
    ## 73  TCCTTGCTGTGGGACATCCAT      GCACGCCCATACTAAACGCTG     0      0.2        0
    ## 74  TCCTTGCTGTGGGACATCCAT       CCTGCCTACCTCCAACCCTG     0      0.2        0
    ## 75  TCCTTGCTGTGGGACATCCAT  TGTGTAGTGGGGTCAAACAGAGATG     0      0.2        0
    ## 76  TCCTTGCTGTGGGACATCCAT   CATACTAAACGCTGTGTAGTGGGG     0      0.2        0
    ## 77  TCCTTGCTGTGGGACATCCAT   GGGGTCAAACAGAGATGAATCAGG     0      0.2        0
    ## 78  TCCTTGCTGTGGGACATCCAT         CCTCCAACCCTGTGCACG     0      0.2        0
    ## 79  TCCTTGCTGTGGGACATCCAT    CGCCCATACTAAACGCTGTGTAG     0      0.2        0
    ## 80  TCCTTGCTGTGGGACATCCAT   GCTGTGTAGTGGGGTCAAACAGAG     0      0.2        0
    ## 81  TCCTTGCTGTGGGACATCCAT  GTGGGGTCAAACAGAGATGAATCAG     0      0.2        0
    ## 82  TCCTTGCTGTGGGACATCCAT     GCTGTGTAGTGGGGTCAAACAG     0      0.2        0
    ## 83  TCCTTGCTGTGGGACATCCAT     CTAAACGCTGTGTAGTGGGGTC     0      0.2        0
    ## 84  TCCTTGCTGTGGGACATCCAT       GCCTACCTCCAACCCTGTGC     0      0.2        0
    ## 85  TCCTTGCTGTGGGACATCCAT          TCCAACCCTGTGCACGC     0      0.2        0
    ## 86  TCCTTGCTGTGGGACATCCAT          CCAACCCTGTGCACGCC     0      0.2        0
    ## 87  TCCTTGCTGTGGGACATCCAT        CCCTGCCTACCTCCAACCC     0      0.2        0
    ## 88  TCCTTGCTGTGGGACATCCAT        CCCCTGCCTACCTCCAACC     0      0.2        0
    ## 89  TCCTTGCTGTGGGACATCCAT         CCCTGTGCACGCCCATAC     0      0.2        0
    ## 90  TCCTTGCTGTGGGACATCCAT      CCTACCTCCAACCCTGTGCAC     0      0.2        0
    ## 91  TCCTTGCTGTGGGACATCCAT        ACCCCTGCCTACCTCCAAC     0      0.2        0
    ## 92  TCCTTGCTGTGGGACATCCAT     CCTGTGCACGCCCATACTAAAC     0      0.2        0
    ## 93  TCCTTGCTGTGGGACATCCAT      CGCTGTGTAGTGGGGTCAAAC     0      0.2        0
    ## 94  TCCTTGCTGTGGGACATCCAT     CGCCCATACTAAACGCTGTGTA     0      0.2        0
    ## 95  TCCTTGCTGTGGGACATCCAT GGGGTCAAACAGAGATGAATCAGGTA     0      0.2        0
    ## 96  TCCTTGCTGTGGGACATCCAT       CCCTGTGCACGCCCATACTA     0      0.2        0
    ## 97  TCCTTGCTGTGGGACATCCAT          CCCTGTGCACGCCCATA     0      0.2        0
    ## 98  TCCTTGCTGTGGGACATCCAT   TGTAGTGGGGTCAAACAGAGATGA     0      0.2        0
    ## 99  TCCTTGCTGTGGGACATCCAT    TGTGTAGTGGGGTCAAACAGAGA     0      0.2        0
    ## 100 TCCTTGCTGTGGGACATCCAT    GCTGTGTAGTGGGGTCAAACAGA     0      0.2        0
    ##     similar_signatures        missing_signatures enzyme digest_score fragments
    ## 1                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 2                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 3                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 4                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 5                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 6                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 7                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 8                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 9                      HPV18, HPV16, HPV4, HPV17                   0         0
    ## 10                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 11                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 12                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 13                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 14                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 15                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 16                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 17                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 18                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 19                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 20                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 21                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 22                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 23                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 24                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 25                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 26                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 27                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 28                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 29                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 30                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 31                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 32                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 33                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 34                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 35                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 36                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 37                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 38                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 39                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 40                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 41                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 42                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 43                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 44                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 45                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 46                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 47                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 48                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 49                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 50                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 51                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 52                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 53                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 54                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 55                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 56                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 57                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 58                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 59                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 60                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 61                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 62                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 63                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 64                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 65                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 66                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 67                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 68                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 69                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 70                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 71                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 72                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 73                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 74                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 75                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 76                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 77                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 78                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 79                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 80                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 81                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 82                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 83                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 84                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 85                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 86                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 87                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 88                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 89                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 90                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 91                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 92                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 93                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 94                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 95                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 96                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 97                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 98                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 99                     HPV18, HPV16, HPV4, HPV17                   0         0
    ## 100                    HPV18, HPV16, HPV4, HPV17                   0         0

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

``` r
# tryCatch({dbDisconnect(dbconn)}, error=warning)
# tryCatch({dbDisconnect(dbConn)}, error=warning)
```
