# RNAseqQC

<!-- badges: start -->
<!-- badges: end -->

The goal of RNAseqQC is to aid quality control of RNAseq data by providing a 
collection of data visualization functions. It allows identification of samples
with unwanted biological or technical effects and to explore differential testing results.

## Installation

You can install the released version of RNAseqQC from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RNAseqQC")
```

## Example

This is a basic example in which we make a library complexity plot and then compare some samples to the median reference of their respective group:

``` r
library(RNAseqQC)
library("DESeq2")

dds <- makeExampleDESeqDataSet(n=10000, m=30)
plot_library_complexity(dds)
vsd <- vst(dds)
dds$condition
plot_sample_MAs(vsd, group = "condition")[c(1,2,16,17)]
```
