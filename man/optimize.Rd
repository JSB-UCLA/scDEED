\name{optimize}
\alias{optimize}

\title{
Dubious cells detector under tSNE and UMAP
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
When measuring the similarity between the UMAP distances and the PCA distances, we look at this similarity_percent proportion of the smallest distances in UMAP space. We do the same thing for tSNE. Default value 0.5
}

  \item{dubious_cutoff}{
the cutoff for dubious cell classification; if cells have scores are below the dubious_cutoff percentile of null scores, they will be considered dubious. Default setting is 0.05. 
  }
  \item{trustworthy_cutoff}{
the cutoff for trustworthy cell classification; if cells have scores are above the trustworthy_cutoff percentile of null scores, they will be considered trustworthy. Default setting is 0.95. 
  }
  
  \item{check_duplicates}{
This is an argument to Seurat::RunTSNE. Default = T. If there are duplicates in the data, t-SNE will not proceed. If the user believes there are true biological duplicates in the data, they may change this setting to F.
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
A wrapper function for the scDEED method. It is similar to the function \texttt{scDEED}, but there are no defaults, and it can only handle one hyperparameter input. Additionally, \texttt{results.PCA} must be provided (in \texttt{scDEED}, this will be calculated internally). You can obtain the same results using \texttt{scDEED} at a single hyperparameter setting. 
}





