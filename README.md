
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scDED
Single-cell dubious embedding detector (scDED): a statistical method for detecting dubious non-linear embeddings.
This package is used for determining the reliability of non-linear dimension reduction embeddings.
It provides functions to detect dubious cells and trustworthy
cells under tSNE and UMAP. Furthermore, based on the number of dubious
cells, functions in this package will provide the best perplexity
parameter under tSNE and best n.neighbors parameter under UMAP.

## Installation

You can install the released version of dubiousdetector from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dubiousdetector")
```

## Example

This is a basic example which shows you how to find the best parameter.
We use pmbc data as a demo:

``` r
suppressPackageStartupMessages(library(dubiousdetector))
data(pbmc.data)
```

### Choose the suitable dimension (num\_pc)

``` r
chooseK(pbmc.data)
```

<img src="man/figures/chooseK_ex.png" width="100%" /> 

### Example for umap

``` r
umap_example <- umap_tsne_process(pbmc.data , num_pc = 10, use_method = "umap",visualization = TRUE)
```

### Example for tsne

``` r
tsne_example <- umap_tsne_process(pbmc.data, num_pc = 10, use_method = "tsne",visualization = TRUE)
```
