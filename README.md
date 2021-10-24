
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scDEED (single-cell dubious embeddings detector): a statistical method for detecting dubious non-linear embeddings
- This package is used to determine the reliability of non-linear dimension reduction embeddings. It provides functions to detect dubious cells and trustworthy cells in tSNE and UMAP embeddings. Furthermore, by minimizing the number of dubious cells, functions in this package find the best perplexity parameter of tSNE and the best n.neighbors parameter of UMAP.

- Choose the suitable dimension for PCA (num_pc)

- Input count matrix should contain cells as columns and genes as rows

## Installation
You can install the released version of scDEED from GitHub with:

``` r
library(devtools)
devtools::install_github("JSB-UCLA/scDEED")
```

## Example

This is a basic example showing how to find the best parameter.
We use pmbc data as a demo:

``` r
suppressPackageStartupMessages(library(scDEED))
data(pbmc.data)
```

### Choose the suitable dimension for PCA (num\_pc)

``` r
chooseK(pbmc.data)
```
ChooseK plot:

<img src="man/figures/chooseK_ex.png" width="100%" /> 

### Example for umap

``` r
umap_example <- umap_tsne_process(pbmc.data , num_pc = 10, use_method = "umap",visualization = TRUE)
```

``` r
head(umap_example$`number of dubious cells corresponding to n.neighbors list`)
```
|   |  n.neighbors |  number.of.dubious.cells |
| - | ------------ | ------------------------ |
| 1 | 	5           | 	42		                    |
| 2 | 	6           | 	25		                    |
| 3 | 	7           | 	57		                    |
| 4 | 	8           | 	20		                    |
| 5 | 	9           | 	6		                     |
| 6 | 	10          | 	32		                    |

``` r
umap_example$`best n.neighbors`
```
9

``` r
umap_example$`UMAP plot with dubious cells`
```

UMAP Plot corresponding to the best n.neighbers, highlighting the dubious cells:
<img src="man/figures/umap_dubious.png" width="100%" /> 

``` r
umap_example$`UMAP plot with trustworthy cells`
```
UMAP Plot corresponding to the best n.neighbers, highlighting the trustworthy cells:
<img src="man/figures/umap_trustworthy.png" width="100%" /> 

``` r
umap_example$`plot. # of dubious embeddings vs parameters`
```
Plot of number of dubious embeddings vs parameters for UMAP:
<img src="man/figures/umap dub em vs parameter.png" width="100%" /> 

### Example for tsne

``` r
tsne_example <- umap_tsne_process(pbmc.data, num_pc = 10, use_method = "tsne",visualization = TRUE)
```

``` r
head(tsne_example$`number of dubious cells corresponding to perplexity list`)
```

|   |  perplexity |  number.of.dubious.cells |
| - | ----------- | ------------------------ |
| 1 | 	20         | 	30                      |
| 2 | 	50         | 	50                      |
| 3 | 	80         | 	20                      |
| 4 | 	110        | 	14                      |
| 5 | 	140        | 	23                      |
| 6 | 	170        | 	23                      |


``` r
tsne_example$`best perplexity`
```
110

``` r
tsne_example$`tSNE plot with dubious cells`
```

tSNE Plot corresponding to the best perplexity, highlighting the dubious cells:
<img src="man/figures/tsne_dubious.png" width="100%" /> 

``` r
tsne_example$`tSNE plot with trustworthy cells`
```
tSNE Plot corresponding to the best perplexity, highlighting the trustworthy cells:
<img src="man/figures/tsne_trustworthy.png" width="100%" /> 

``` r
tsne_example$`plot. # of dubious embeddings vs parameters`
```
Plot of number of dubious embeddings vs parameters for tSNE:
<img src="man/figures/tsne dub em vs parameter.png" width="100%" /> 