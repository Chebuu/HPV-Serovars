---
title: Designing PCR Diagnostics that Discriminate HPV Serotypes by DNA Sequence
author: chebuu
output: github_document
---


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
The first infectious cause of cancer was identified in 1911 by Peyton Rous who demonstrated the transmissibility of a tumor in fowl through injection of sub 0.2 micron filtrate of the tumor \cite{Rous_1911}.

#### Installing this package

```r
devtools::install_github('Chebuu/HPV-Serovars')
```

```
## Skipping install of 'HPVSerovars' from a github remote, the SHA1 (6e73ff62) has not changed since last install.
##   Use `force = TRUE` to force installation
```

```r
library(HPVSerovars)
```

(Not run) This vignette is packaged and accessible via: `browseVignettes('HPVSerovars')`.

#### Installing OligAarrayAux

Download [OligoArrayAux](http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux) and see the following vignette for installation instructions.


```r
library(devtools)

devtools::install_github('chebuu/Design-Group-Specific-Primers')

library(Design_Group_Specific_Primers)

vignette('Installing-OligoArrayAux', package = 'Design_Group_Specific_Primers')
```

#### Installing OpenPrimer Dependencies

From the [openPimer vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/openPrimeR/inst/doc/openPrimeR_vignette.html):

openPrimeR requires external programs for some features, particularly for computing the physicochemical properties of primers. Please make sure you have the following tools installed on your system such that they are in your system’s path:

[MELTING](http://www.ebi.ac.uk/biomodels/tools/melting/) (>= 5.1.1): For melting temperature computations.
[ViennaRNA](http://www.tbi.univie.ac.at/RNA/) (>= 2.2.4): For secondary structure prediction.
[OligoArrayAux](http://unafold.rna.albany.edu/OligoArrayAux.php) (>= 3.8): For primer efficiency computations as performed by DECIPHER.
[MAFFT](http://mafft.cbrc.jp/alignment/software/) (>= 7.305): For computing multiple sequence alignments.
[Pandoc](http://pandoc.org/) (>= 1.19.1): For creating PDF reports.

If you would like to be able to access the immunoglobulin repository IMGT from the openPrimeR Shiny app, you should additionally fulfill the following dependencies:

[PhantomJS](http://phantomjs.org/) (>= 2.1): For headless website calls.
[Python](http://www.python.org/) (>=2.7.9) and the [selenium](http://selenium-python.readthedocs.io/) (>=3.0.1) module: For data extraction scripts.
openPrimeR will automatically check for all dependencies and inform you about any missing dependencies when the package is attached:


```r
# library(openPrimeR)
```

Note that the tool is still functional if there are missing external programs. However, we recommend that all dependencies are fulfilled to guarantee the best user experience.

### Install Bioconductor packages
[https://www.bioconductor.org/](https://www.bioconductor.org/)

```r
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


```r
doAlignment <- function(xset) {
  useqs <- unique(xset)
  uidxs <- match(xset, useqs)
  aseqs <- AlignSeqs(useqs, verbose=F) 
  aseqs[uidxs]
}
```

#### Quick-ref
Some references to remind me of packages, methods, etc.:
```{}
# CRAN
# https://rdrr.io/cran/dpcR/f/vignettes/overview.Rmd

# BIOC
# https://bioconductor.org/packages/release/bioc/html/twoddpcr.html
# https://bioconductor.org/packages/release/bioc/html/ddCt.html
# https://bioconductor.org/packages/3.5/bioc/html/mBPCR.html

# Other
# https://doi.org/10.1007/s13402-017-0331-y
```

# Wiki Pheno Data
I copied the tables from [Wiki: Human_papillomavirus_infection ](https://en.wikipedia.org/wiki/Human_papillomavirus_infection) and included them in the R package assoc. w/ this document `github.com/Chebuu/HPV-Serovars`.

```r
data(HPV.Sequelae)
data(HPV.GenitalCancerRisk)
```

The Wiki data contains phenotypic categorizations of several HPV strains. It's useful for testing.


```r
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



Table: Unique HPV serotypes (N=34) (Total=13): 

|       |
|:------|
|HPV-1  |
|HPV-2  |
|HPV-3  |
|HPV-4  |
|HPV-6  |
|HPV-7  |
|HPV-10 |
|HPV-11 |
|HPV-13 |
|HPV-16 |
|HPV-18 |
|HPV-22 |
|HPV-26 |
|HPV-28 |
|HPV-31 |
|HPV-32 |
|HPV-33 |
|HPV-35 |
|HPV-39 |
|HPV-42 |
|HPV-44 |
|HPV-45 |
|HPV-51 |
|HPV-52 |
|HPV-53 |
|HPV-56 |
|HPV-58 |
|HPV-59 |
|HPV-60 |
|HPV-63 |
|HPV-66 |
|HPV-68 |
|HPV-73 |
|HPV-82 |

```r
hist(unlist(HPV), max(unlist(HPV)), col='blue'); head(HPV, 4)
```

![](/private/var/folders/3p/n0cd926d4tq39fmwy8sqmbdw0000gn/T/RtmpA93Z54/preview-1331dc9c1c8d.dir/Manuscript_files/figure-gfm/wikiprint-1.png)<!-- -->

```
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
```


# PaVE Data

Reference genomes for all serovars listed in `data(HPV.Sequelae)`and `data(HPV.GenitalCancerRisk)`were selected from among the results of the following [PaVEPaVE](https://pave.niaid.nih.gov) query:

```{} 
https://pave.niaid.nih.gov/#search/search_database/kw
?dbNamespace=Genomes
&includeNR=true
&refCloneOnly=false
&sort=Locus_ID
&sortType=true
&page=600&start=1
&text=Human&showTable=1
&
```





























