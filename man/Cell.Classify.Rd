\name{Cell.Classify}
\alias{Cell.Classify}

\title{
Classification of cells using similarity scores 
}

\usage{
Cell.Classify(rho_original, rho_permuted, dubious_cutoff = 0.05, trustworthy_cutoff = 0.95)
}

\arguments{
  \item{rho_original}{
A vector containing the similarity scores calculated on the original data
}
  \item{rho_permuted}{
A vector containing the similarity scores calculated on the permuted data
}
  \item{dubious_cutoff}{
The cutoff for dubious cells (default = 0.05). Cells with scores worse (lower) than the dubious_cutoff percentile of null scores will be considered dubious. A lower dubious_cutoff means that to be considered dubious, cells will have to have lower scores. A higher dubious_cutoff means that cells can score higher and still be considered dubious. It is similar to significance level in hypothesis testing. 
}
  \item{trustworthy_cutoff}{
The cutoff for trustworthy cells (default = 0.95). Cells with scores better (higher) than the trustworthy_cutoff percentile of null scores will be considered trustworthy. A lower trustworthy_cutoff means that to be considered trustworthy, cells will not have to score as high. A higher trustworthy_cutoff means that cells will need to score higher in order to be considered trustworthy. It is similar to significance level in hypothesis testing. 
}

   \item{similarity_percent}{
The percentage of cells to consider in the similarity score calculations (default = 0.5). scDEED uses the nearest floor(number of cells * similarity_percent) neighbors in the similarity percent calculations. Intuitively, a higher similarity score considers more cells as neighbors (emphasis on global preservation) while a lower similarity score considers less (emphasis on local preservation) 
}

}

\value{
Returns a list with three vectors, (1) 'dubious_cells': the cell embeddings with dubious classification, as determined by the \code{dubious_cutoff} on the null data (2) 'trustworthy_cells': the cell embeddings with trustworthy classification, as determined by the \code{trustworthy_cutoff} on the null data (3) 'intermediate_cells': the cells with intermediate classification, as determined by the \code{dubious_cutoff} and \code{trustworthy_cutoff} on the null data

}

\description{
This function classifies the cell embeddings as dubious, trustworthy, or intermediate. It uses the cell similarity scores calculated on permuted data (\code{rho_permuted}) to classify the cell similarity scores calculated on the original data (\code{rho_original}). The classifications are based on the dubious_cutoff and trustworthy cutoff. 

}





