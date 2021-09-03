\name{pbmc.data}
\alias{pbmc.data}
\docType{data}
\title{
pbmc3k
}
\description{
The dataset of Peripheral Blood Mononuclear Cells (PBMC) including 2,700 single cells that were sequenced on the Illumina NextSeq 500
}
\usage{data("pbmc.data")}
\format{
  The format is:
Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  ..@ i       : int [1:2286884] 70 166 178 326 363 410 412 492 494 495 ...
  ..@ p       : int [1:2701] 0 781 2133 3264 4224 4746 5528 6311 7101 7634 ...
  ..@ Dim     : int [1:2] 32738 2700
  ..@ Dimnames:List of 2
  .. ..$ : chr [1:32738] "MIR1302-10" "FAM138A" "OR4F5" "RP11-34P13.7" ...
  .. ..$ : chr [1:2700] "AAACATACAACCAC-1" "AAACATTGAGCTAC-1" "AAACATTGATCAGC-1" "AAACCGTGCTTCCG-1" ...
  ..@ x       : num [1:2286884] 1 1 2 1 1 1 1 41 1 1 ...
  ..@ factors : list()
}

\source{
https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
}
\references{
10X Genomics
}
\examples{
data(pbmc.data)
## maybe str(pbmc.data) ; plot(pbmc.data) ...
}
\keyword{datasets}
