\name{umap_tsne_process}
\alias{umap_tsne_process}

\title{
Dubious cells detector under tSNE and UAMP
}

\usage{
umap_tsne_process(pbmc, num_pc,
n_neighbors = c(seq(from=5,to=30,by=1),35,40,45,50),
similarity_percent = 0.5, visualization = FALSE,
use_method = "umap",
perplexity = c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50)),
perplexity_score = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pbmc}{
a Seurat object or a count matrix with rows as the genes
}
  \item{num_pc}{
number of principal components
}
  \item{n_neighbors}{
 a list of choices for n.neighbors parameter in UMAP
}
  \item{similarity_percent}{
When measuring the similarity between the UMAP distances and the PCA distances, we look at this similarity_percent proportion of the smallest distances in UMAP space. We do the same thing for tSNE. Default value 50\%
}
  \item{visualization}{
a boolean value specifies whether we will display the plots with dubious/trustworthy cells
}
  \item{use_method}{
'tsne' or 'umap'
}
  \item{perplexity}{
a list of choices for perplexity parameters in tSNE
}
  \item{perplexity_score}{
a value for perplexity parameter as intialization in tSNE, default 30
}
}

\value{
A list contains i) a vector of dubious numbers corresponding to the list of n.neighbors parameters in UMAP or the list of perplexity parameters in tSNE ii) an n.neighber  or a perplexity corresponding to the lowest dubious number.

If the input visulization is set to be TRUE, i) a vector containing indices for the dubious cells will be returned ii) a plot displaying all dubious cells, and a plot displaying all trustworthy cells will be returned, both corresponding to the parameter resulting in the lowest number of dubious cells.
}

\description{
If the user chooses method umap, this function provides the dubious numbers for all provided n.neighber choices and the n.neighber parameter corresponding to the lowest number of dubious cells.
If the user chooses method tsne, this function provides the dubious numbers for all provided perplexity choices and the perplexity parameter corresponding to the lowest number of dubious cells.
}





