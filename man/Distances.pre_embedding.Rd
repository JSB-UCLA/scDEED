\name{Distances.pre_embedding}
\alias{Distances.pre_embedding}

\title{
Calculation of pre-embedding distances
}

\usage{
Distances.pre_embedding(pbmc, pbmc_permuted, K, pre_embedding = 'pca')
}

\arguments{
  \item{pbmc}{
The original Seurat object 
}
  \item{pbmc_permuted}{
The permuted Seurat object, i.e. the output of \code{Permuted}. 
}

  \item{K}{
The number of PCs to use
}
  \item{pre_embedding}{
Which embedding set to use for distance calculation, default = 'pca'. If the user would like to use an alternate method, like ICA, they should perform that method for the original and permuted data, then specify the slot name here
}
}

\value{
Returns a list with two distance matrices, (1) 'pre_embedding_distances': a matrix showing the distances between cells in the original pre-embedding space (2) 'pre_embedding_distances_permuted': a matrix showing the distances between cells in the permuted pre-embedding space

}

\description{
This function calculates the distances between cells in the pre-embedding space for the original and permuted data 

}





