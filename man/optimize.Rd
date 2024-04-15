\name{optimize}
\alias{optimize}

\title{
Dubious cells detector under tSNE and UMAP for a single hyperparameter setting
}

\usage{
optimize(input_data, input_data.permuted, pre_embedding, reduction.method, K,
                    n, m, perplexity, results.PCA, similarity_percent, dubious_cutoff,
                    trustworthy_cutoff, check_duplicates = T, rerun = T)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_data}{
a Seurat object
}
  \item{input_data.permuted}{
a Seurat object containing permuted data
}
  \item{pre_embedding}{
the slot to use as input for t-SNE and UMAP. If users would like to use a different pre-embedding space, they can add this to the Seurat object and specify the name here. 
}
  \item{reduction.method}{
Which dimension reduction method to use; currently the package is only set up for 'tsne' or 'umap'
}
  \item{K}{
number of principal components
}

  \item{n}{
 input for the n.neighbors parameter in UMAP
}
  \item{m}{
  input for the min.dist parameter in UMAP
  }
  \item{perplexity}{
  input for the perplexity parameter in tSNE
}
  \item{results.PCA}{
  A named list containing the cell-cell distance matrices calculated in the pre-embedding space for the original (name = "pre_embedding_distances") and permuted data (name = "pre_embedding_distances_permuted").  
}

  \item{similarity_percent}{
The percentage of cells to consider in the similarity score calculations (default = 0.5). scDEED uses the nearest floor(number of cells * similarity_percent) neighbors in the similarity percent calculations. Intuitively, a higher similarity score considers more cells as neighbors (emphasis on global preservation) while a lower similarity score considers less cells (emphasis on local preservation) 
}

  \item{dubious_cutoff}{The cutoff for dubious cells (default = 0.05). Cells with scores worse (lower) than the dubious_cutoff percentile of null scores will be considered dubious. A lower dubious_cutoff means that to be considered dubious, cells will have to have lower scores. A higher dubious_cutoff means that cells can score higher and still be considered dubious. It is similar to significance level in hypothesis testing.   }  \item{trustworthy_cutoff}{The cutoff for trustworthy cells (default = 0.95). Cells with scores better (higher) than the trustworthy_cutoff percentile of null scores will be considered trustworthy. A lower trustworthy_cutoff means that to be considered trustworthy, cells will not have to score as high. A higher trustworthy_cutoff means that cells will need to score higher in order to be considered trustworthy. It is similar to significance level in hypothesis testing.   }
  
  \item{check_duplicates}{
This is an argument to \code{Seurat::RunTSNE}. Default = T. If there are duplicates in the data, t-SNE will not proceed. If the user believes there are true biological duplicates in the data, they may change this setting to F.
  }
  \item{rerun}{
This is a time-saving argument (default = T). If the user has already performed dimension reduction and would only like to check the results of that dimension reduction, then they can use rerun=F so scDEED does not re-run the embedding method on the data. In most cases, rerun=T because if you are optimizing hyperparameters, the function will need to rerun the embedding method. 
  }
}

\value{
a vector of 4 items containing (1): number of dubious cells (2): the indices for the dubious cells, separated by commas (3): the indicies for the trustworthy cells, separated by commas (4): the indicies for the intermediate cells, separated by commas

If one of the categories (dubious, trustworthy, or intermediate) is empty, the entry is 'none'

}

\description{
A wrapper function for the scDEED method. It is similar to the function \code{scDEED}, but there are no defaults, and it can only handle one hyperparameter input. Additionally, \code{results.PCA} must be provided (in \code{scDEED}, this will be calculated internally). You can obtain the same results using \code{scDEED} at a single hyperparameter setting. 
}





