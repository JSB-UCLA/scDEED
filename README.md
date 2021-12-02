
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scDEED (single-cell dubious embeddings detector): a statistical method for detecting dubious non-linear embeddings
- This package is used to determine the reliability of non-linear dimension reduction embeddings. It provides functions to detect dubious cells and trustworthy cells in tSNE and UMAP embeddings. Furthermore, by minimizing the number of dubious cells, functions in this package find the best perplexity parameter of tSNE and the best n.neighbors/min.dist parameter of UMAP.

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

<img src="man/figures/10000_cell_pca.png" width="100%" /> 

### Example for umap

``` r
umap_example <- umap_tsne_process(pbmc.data , num_pc = 10, use_method = "umap",visualization = TRUE)
```

``` r
head(umap_example$`UMAP plot with dubious cells - best pair of n.neighbors and min.dist`)
```
|   | n.neighbors | min.dist | number of dubious cells |
|---|-------------|----------|-------------------------|
| 1 | 5           | 0.1      | 2                       |
| 2 | 5           | 0.3      | 7                       |
| 3 | 5           | 0.5      | 17                      |
| 4 | 5           | 0.7      | 4                       |
| 5 | 5           | 0.9      | 4                       |
| 6 | 6           | 0.1      | 5                       |

``` r
umap_example$`best n.neighbors`
```
9.0, 0.1


Comparative UMAP plots of the randomly selected 10000 cells from Hydra dataset under the n.neighbors 50, min.dist 0.7 and the n.neighbors 5, min.dist 0.5
optimized by scDEED:

Before optimization:
<img src="man/figures/10000_umap_original.png" width="100%" />

After optimization:
<img src="man/figures/10000_umap_optimized.png" width="100%" />


``` r
umap_example$`plot. # of dubious embeddings vs pair of n.neighbors and min.dist`
```
Plot of number of dubious embeddings vs pair of n.neighbors and min.dist for UMAP:
<img src="man/figures/nm_vs_dub.png" width="100%" /> 



### Example for tsne

``` r
tsne_example <- umap_tsne_process(pbmc.data, num_pc = 10, use_method = "tsne",visualization = TRUE)
```

``` r
head(tsne_example$`number of dubious cells corresponding to perplexity list`)
```

|   |  perplexity |  number of dubious cells |
| - | ----------- | ------------------------ |
| 1 | 	20         | 	14                      |
| 2 | 	50         | 	9                      |
| 3 | 	80         | 	4                      |
| 4 | 	110        | 	5                      |
| 5 | 	140        | 	5                      |
| 6 | 	170        | 	5                      |


``` r
tsne_example$`best perplexity`
```
80



Comparative tSNE plots of the randomly selected 10000 cells from Hydra dataset under the perplexity 20 and the perplexity 50
optimized by scDEED:

Before optimization:
<img src="man/figures/10000_tsne_original.png" width="100%" /> 

After optimization:
<img src="man/figures/10000_tsne_optimized.png" width="100%" /> 

``` r
tsne_example$`plot. # of dubious embeddings vs parameters`
```
Plot of number of dubious embeddings vs parameters for tSNE:
<img src="man/figures/tsne_per_dub.png" width="100%" /> 
