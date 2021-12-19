\name{chooseK}
\alias{chooseK}
\title{
Choose the number of PC based on an elbow plot
}
\usage{
chooseK(input_data)
}

\arguments{
  \item{input_data}{a Seurat object or a count matrix with rows as the genes}
}
\value{
an elbow plot for users to choose suitable dimensionality, in terms of numbers of Principal Components (num_pc), as the input of umap_tsne_process.
}
\description{
Generate an elbow plot(a ranking of principle components based on the percentage of variance explained by each one) for users to choose num_pc as the input of umap_tsne_process.
}


