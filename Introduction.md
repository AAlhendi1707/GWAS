# Introduction


This web tutorial is derived from 'A guide to genome-wide association analysis and post-analytic interrogation' (Statistics in Medicine, in review). The tutorial presents fundamental concepts and specific software tools for implementing a complete genome wide association (GWA) analysis, as well as post-analytic visualization and interrogation of potentially novel findings. In this tutorial we use complete GWA data on 1401 individuals from [the PennCATH study of coronary artery disease (CAD).](http://www.ncbi.nlm.nih.gov/pubmed/21239051)

In the steps to follow we begin by demonstrating a method for downloading necessary R packages and setting global parameters as a means for saving progress while working through a GWA analysis. Next, we include quality control steps for both SNP and sample level filtering. The third section is split into principal component calculation for population stratification in statistical modeling, as well as imputation of non-typed SNPs using 1000 Genomes reference genotype data. We then demonstrate strategies to carry out the GWA analysis on the typed data using basic linear modeling functionality in R, as well as imputed data using functionality contained within the `snpStats` package. Finally, we demonstrate means for post-analytic interrogation, including methods for evaluating the performance of statistical models, as well as visualization of the global and subsetted GWAS output.

## Installing necessary packages

```r
# Run this once interactively to download and install BioConductor packages and other packages.

source ("http://bioconductor.org/biocLite.R ")
list.of.packages <- c("snpStats", "SNPRelate","rtracklayer", "biomaRt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


list.of.packages <- c('plyr', 'LDheatmap', 'doParallel', 'ggplot2', 'coin' ,'igraph', 'devtools', 'downloader')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# GenABEL has moved to CRAN archive. The below command for local installation from CRAN archive.
install.packages("https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.7-6.tar.gz", 
type = "source", repos = NULL)


library(devtools)
install_url("http://cran.r-project.org/src/contrib/Archive/postgwas/postgwas_1.11.tar.gz")

```


## Configuring global parameters

Customize and Run [globals.R](R/globals.R)

```r
source("globals.R")

# Downloading support files
# Download and unzip data needed for this tutorial

library(downloader)

download(urlSupport, zipSupport.fn)
unzip(zipSupport.fn, exdir = data.dir)

```

