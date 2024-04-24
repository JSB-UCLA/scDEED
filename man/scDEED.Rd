\name{scDEED}
\alias{scDEED}

\title{
Dubious cells detector under tSNE and UMAP
}

\usage{
scDEED(input_data, K, n_neighbors = c(5, 20, 30, 40,50), min.dist = c(0.1, 0.4), similarity_percent = 0.5,reduction.method,
                  perplexity = c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50)), pre_embedding = 'pca', slot = 'scale.data', 
                  dubious_cutoff = 0.05, trustworthy_cutoff = 0.95, permuted = NA, check_duplicates = T, rerun = T, default_assay = 'RNA')
  

}

\arguments{
  \item{input_data}{
a Seurat object
}
  \item{K}{
number of principal components
}
  \item{reduction.method}{
Which dimension reduction method to use; currently the package is only set up for 'tsne' or 'umap'
}
  \item{n_neighbors}{
 a vector containing the settings to test for n.neighbors parameter in UMAP. Default setting is c(5, 20, 30, 40,50). To marginally optimize the n_neighbors hyperparameter, the \code{min.dist} setting should be a single value. 
}
  \item{min.dist}{
  a vector containing the settings to test for min.dist parameter in UMAP. Default setting is c(0.1, 0.4). To marginally optimize the min.dist hyperparameter, the \code{n_neighbors} setting should be a single value. 
  }
  \item{similarity_percent}{
The percentage of cells to consider in the similarity score calculations (default = 0.5). scDEED uses the nearest floor(number of cells * similarity_percent) neighbors in the similarity percent calculations. Intuitively, a higher similarity score considers more cells as neighbors (emphasis on global preservation) while a lower similarity score considers less cells (emphasis on local preservation) 
}
  \item{perplexity}{
a vector of choices for perplexity parameter in tSNE. Default setting is c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50))
}
  \item{pre_embedding}{
the slot to use as input for t-SNE and UMAP. Default setting is 'pca'. If users would like to use a different pre-embedding space, like ICA, they should perform that method for the original and permuted data, then specify the slot name here
}
  \item{slot}{
Which slot to use for permutation, default = 'scale.data'. If the user would like to use an alternate method of pre-processing/normalization, the user should add the matrix to the Seurat object and specify the slot name here
  }
  \item{dubious_cutoff}{
The cutoff for dubious cells (default = 0.05). Cells with scores worse (lower) than the dubious_cutoff percentile of null scores will be considered dubious. A lower dubious_cutoff means that to be considered dubious, cells will have to have lower scores. A higher dubious_cutoff means that cells can score higher and still be considered dubious. It is similar to significance level in hypothesis testing. 
  }
  \item{trustworthy_cutoff}{
The cutoff for trustworthy cells (default = 0.95). Cells with scores better (higher) than the trustworthy_cutoff percentile of null scores will be considered trustworthy. A lower trustworthy_cutoff means that to be considered trustworthy, cells will not have to score as high. A higher trustworthy_cutoff means that cells will need to score higher in order to be considered trustworthy. It is similar to significance level in hypothesis testing. 
  }
  \item{permuted}{
The user may provide their own permuted data here. In general, users will not need to provide this object. Users may want to use this if they are going to use a different pre-embedding space and need a different algorithm, i.e. ICA or something similar, to be applied to the permuted data matrix. Since we do not have this functionality in the scDEED package, users would need to permute the data and apply ICA themselves, then provide this permuted object here. For more details, please see our tutorial. 
  }
  \item{check_duplicates}{
This is an argument to \code{Seurat::RunTSNE}. Default = T. If there are duplicates in the data, t-SNE will not proceed. If the user believes there are true biological duplicates in the data, they may change this setting to F.
  }
  \item{rerun}{
This is a time-saving argument (default = T). If the user has already performed dimension reduction and would only like to check the results of that dimension reduction, then they can use rerun=F so scDEED does not re-run the embedding method on the data. In most cases, rerun=T because if you are optimizing hyperparameters, the function will need to rerun the embedding method. 
  }
\item{default_assay}{
The assay to use for optimization; default is 'RNA'
  }

}

\value{

A list contains i) num_dubious: a data.frame containing n.neighbors/min.dist/perplexity hyperparameters and their corresponding number of dubious cells in UMAP/t-SNE respectively ii) full_results: a data.frame containing the number of dubious cells at each hyperparameter setting, as well as the cell classifications in the final 3 columns. The cell classifications are separated by commas. 

}

\description{
If the user chooses method 'umap', this function provides the number of dubious cells for all provided n.neighber/min.dist choices.
If the user chooses method 'tsne', this function provides the number of dubious cells for all provided perplexity choices.
}





