\name{Cell.Similarity}
\alias{Cell.Similarity}

\title{
Calculation of cell similarity scores for the original and permuted data
}

\usage{
Cell.Similarity(pre_embedding_distances,pre_embedding_distances_permuted, reduced_dim_distances,reduced_dim_distances_permuted,similarity_percent = .50)
}

\arguments{
  \item{pre_embedding_distances}{
A matrix containing the cell-cell distances calculated in the pre-embedding space for the original data
}
  \item{pre_embedding_distances_permuted}{
A matrix containing the cell-cell distances calculated in the pre-embedding space for the permuted data
}
  \item{reduced_dim_distances}{
A matrix containing the cell-cell distances calculated in the post-embedding space for the original data
}
  \item{reduced_dim_distances_permuted}{
A matrix containing the cell-cell distances calculated in the post-embedding space for the permuted data
}

   \item{similarity_percent}{
The percentage of cells to consider in the similarity score calculations (default = 0.5). scDEED uses the nearest floor(number of cells * similarity_percent) neighbors in the similarity percent calculations. Intuitively, a higher similarity score considers more cells as neighbors (emphasis on global preservation) while a lower similarity score considers less cells (emphasis on local preservation) 
}

}

\value{
Returns a list with two vectors, (1) 'rho_original': the similarity scores calculated on the original data (2) 'rho_permuted': the similarity scores calculated on the permuted data

}

\description{
This function calculates the similarity scores for all cells on the original and permuted data. The permuted data will be used to estimate the null distribution. 

}





