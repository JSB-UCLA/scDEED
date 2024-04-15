\name{Distances.tSNE}
\alias{Distances.tSNE}

\title{
Calculation of cell-cell distances after t-SNE
}

\usage{
Distances.tSNE(pbmc,pbmc.permuted, K, perplexity_score = 40, pre_embedding = 'pca', check_duplicates = T, rerun = T)
}

\arguments{
  \item{pbmc}{
The original Seurat object 
}
  \item{pbmc_permuted}{
The permuted Seurat object, i.e. the output of \texttt{Permuted}. 
}

  \item{K}{
The number of PCs to use
}
  \item{perplexity_score}{
The perplexity hyperparameter for t-SNE, default = 40.
}
  \item{pre_embedding}{
Which embedding use as input for t-SNE, default = 'pca'. If the user would like to use an alternate method, like ICA, they should perform that method for the original and permuted data, then specify the slot name here
}
  \item{check_duplicates}{
This is an argument to Seurat::RunTSNE. Default = T. If there are duplicates in the data, t-SNE will not proceed. If the user believes there are true biological duplicates in the data, they may change this setting to F.
}
 \item{rerun}{
This is a time-saving argument (default = T). If the user has already performed dimension reduction and would only like to check the results of that dimension reduction, then they can use rerun=F so scDEED does not re-run the embedding method on the data. In most cases, rerun=T because if you are optimizing hyperparameters, the function will need to rerun the embedding method. 
  }
}

\value{
Returns a list with two distance matrices, (1) 'reduced_dim_distances': a matrix showing the distances between cells in the original post-embedding space (2) 'reduced_dim_distances_permuted': a matrix showing the distances between cells in the permuted post-embedding space

}

\description{
This function calculates the distances between cells in the post-embedding space for the original and permuted data 

}





