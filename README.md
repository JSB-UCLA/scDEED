
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dubiouspack

This package provides functions to detect dubious cells and trustworthy
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
# install.packages("dubiousdetector")
suppressPackageStartupMessages(library(dubiousdetector))
data(pbmc.data)
```

### Choose the suitable dimension (num\_pc)

``` r
chooseK(pbmc.data)
```

<img src="man/figures/README-chooseK-1.png" width="100%" />

### Example for umap

``` r
umap_example <- umap_tsne_process(pbmc.data , num_pc = 10, use_method = "umap",visualization = TRUE)
```

### Example for tsne

``` r
tsne_example <- umap_tsne_process(pbmc.data, num_pc = 10, use_method = "tsne",visualization = TRUE)
```
