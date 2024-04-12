\name{scDEED}
\alias{scDEED}

\title{
Dubious cells detector under tSNE and UMAP
}

\usage{
scDEED(input_data, K, n_neighbors = c(5, 20, 30, 40,50), min.dist = c(0.1, 0.4), similarity_percent = 0.5,reduction.method,
                  perplexity = c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50)), pre_embedding = 'pca', slot = 'scale.data', 
                  dubious_cutoff = 0.05, trustworthy_cutoff = 0.95, permuted = NA, check_duplicates = T, rerun = T)

}
%- maybe also 'usage' for other objects documented here.
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
 a list of choices for n.neighbors parameter in UMAP. Default setting is c(5, 20, 30, 40,50)
}
  \item{min.dist}{
  a list of choices for min.dist parameter in UMAP. Default setting is c(0.1, 0.4)
  }
  \item{similarity_percent}{
When measuring the similarity between the UMAP distances and the PCA distances, we look at this similarity_percent proportion of the smallest distances in UMAP space. We do the same thing for tSNE. Default value 0.5
}
  \item{perplexity}{
a list of choices for perplexity parameter in tSNE. Default setting is c(seq(from=20,to=410,by=30),seq(from=450,to=800,by=50))
}
  \item{pre_embedding}{
the slot to use as input for t-SNE and UMAP. Default setting is 'pca'. If users would like to use a different pre-embedding space, they can add this to the Seurat object and specify the name here. 
}
  \item{slot}{
The input to PCA, Default setting is 'scale.data'. If users would like a different pre-processing pipeline (other than Seurat) they can follow their desired pipeline, then add this as another layer to the Seurat object and specify the name here. 
  }
  \item{dubious_cutoff}{
the cutoff for dubious cell classification; if cells have scores are below the dubious_cutoff percentile of null scores, they will be considered dubious. Default setting is 0.05. 
  }
  \item{trustworthy_cutoff}{
the cutoff for trustworthy cell classification; if cells have scores are above the trustworthy_cutoff percentile of null scores, they will be considered trustworthy. Default setting is 0.95. 
  }
  \item{permuted}{
The user may provide their own permuted data here. In general, users will not need to provide this object. Users may want to use this if they are going to use a different pre-embedding space and need a different algorithm, i.e. ICA or something similar, to be applied to the permuted data matrix. Since we do not have this functionality in the scDEED package, users would need to permute the data and apply ICA themselves, then provide this permuted object to the function. For more details, please see our tutorial. 
  }
  \item{check_duplicates}{
This is an argument to Seurat::RunTSNE. Default = T. If there are duplicates in the data, t-SNE will not proceed. If the user believes there are true biological duplicates in the data, they may change this setting to F.
  }
  \item{rerun}{
This is a time-saving argument (default = T). If the user has already performed dimension reduction and would only like to check the results of that dimension reduction, then they can use rerun=F so scDEED does not re-run the embedding method on the data. In most cases, rerun=T because if you are optimizing hyperparameters, the function will need to rerun the embedding method. 
  }
}

\value{
num_dubious: This is a data.frame that shows the hyperparameter settings and the number of dubious cells at those hyperparameter settings

full_results: This is a data.frame that shows the hyperparamater settings as well as the dubious/trustworthy/intermediate cells. The cell classifications at each hyperparameter setting are available in the last 3 columns, separated by commas.

A list contains i) num_dubious: a data.frame containing n.neighbors/min.dist/perplexity hyperparameters and their corresponding number of dubious cells in UMAP/t-SNE respectively ii) full_results: a data.frame containing the number of dubious cells at each hyperparameter setting, as well as the cell classifications in the final 3 columns. The cell classifications are separated by commas. 

}

\description{
If the user chooses method umap, this function provides the number of dubious cells for all provided n.neighber/min.dist choices.
If the user chooses method tsne, this function provides the number of dubious cells for all provided perplexity choices.
}





