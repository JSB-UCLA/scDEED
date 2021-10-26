\name{umap_tsne_process}
\alias{umap_tsne_process}

\title{
Dubious cells detector under tSNE and UMAP
}

\usage{
umap_tsne_process(pbmc, num_pc, n_neighbors = c(seq(from=5,to=30,by=1),35,40,45,50), min.dist = seq(0.1,0.9, by = 0.2), similarity_percent = 0.5, visualization = FALSE, use_method = "umap", perplexity = c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50)), perplexity_score = 30, optimize_neib = TRUE ,optimize_min = TRUE)
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
  \item{min.dist}{
  a list of choices for min.dist parameter in UMAP
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
a list of choices for perplexity parameter in tSNE
}
  \item{perplexity_score}{
a value for perplexity parameter as intialization in tSNE, default 30
}
  \item{optimize_neib}{
a boolean value to specify if n.neighbors will be optimized, default TRUE
  }
  \item{optimize_min}{
a boolean value to specify if min.dist will be optimized, default TRUE
  }
}

\value{
A list contains i) a dataframe containing n.neighbors/min.dist parameters and their corresponding number of dubious cells in UMAP or the dataframe containing perplexity parameters and their corresponding number of dubious cells in tSNE ii) an n.neighber/min.dist or a perplexity corresponding to the lowest dubious number. iii) a list of dubious cells corresponding to best n.neighbors/min.dist or a list of dubious cells corresponding to best perplexity iv) a list of trustworthy cells corresponding to best n.neighbors/min.dist or a list of trustworthy cells corresponding to best perplexity

If the input visulization is set to be TRUE, two more plots will be added, i) a plot displaying all dubious cells corresponding to the parameter(s) resulting in the lowest number of dubious cells; ii) a plot displaying all trustworthy cells corresponding to the parameter(s) resulting in the lowest number of dubious cells. iii) a plot to show the number of dubious embeddings vs parameters
}

\description{
If the user chooses method umap, this function provides the number of dubious cells for all provided n.neighber/min.dist choices, and also the n.neighber/min.dist parameter corresponding to the lowest number of dubious cells.
If the user chooses method tsne, this function provides the number of dubious cells for all provided perplexity choices, and also the perplexity parameter corresponding to the lowest number of dubious cells. 
}





