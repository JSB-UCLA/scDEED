\name{Permuted}
\alias{Permuted}

\title{
Permutation of data 
}

\usage{
Permuted(pbmc, default_assay = 'RNA', slot = 'scale.data', K)
}

\arguments{
  \item{pbmc}{
a Seurat object to be permuted
}
  \item{default_assay}{
Which assay to use for permutation, default='RNA', but depending on the dataset, users may also want to use 'integrated'
}
  \item{slot}{
Which slot to use for permutation, default = 'scale.data'. If the user would like to use an alternate method of pre-processing/normalization, the user should add the matrix to the Seurat object and specify the slot name here
}
  \item{K}{
The number of PCs to use
}
}

\value{
Returns a new Seurat object in which each row has been permuted and PCA has been calculated on the permuted data matrix. 

} 

\description{
This function permutes each row of the data matrix and calculates PCA on the permuted data matrix. By permuting each row of the data matrix, the marginal gene distributions are maintained. 

}





