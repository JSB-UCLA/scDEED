\name{chooseK}
\alias{chooseK}
\title{
Number of PC choosing technique
}
\usage{
chooseK(pbmc)
}

\arguments{
  \item{pbmc}{a Seurat object or a count matrix with rows as the genes}
}
\value{
an elbow plot for users to choose suitable dimensionality, in terms of numbers of Principal Components (num_pc), as the input of umap_tsne_process.
}
\description{
Generate an elbow plot(a ranking of principle components based on the percentage of variance explained by each one) for users to choose num_pc as the input of umap_tsne_process.
}


