\name{getConsensusClass}
\alias{getConsensusClass}
\title{getConsensusClass}
\description{Nearest-centroid single sample classifier according to the consensus molecular subtypes of muscle-invasive bladder cancer, based on log2-scaled gene expression profile.}
\usage{getConsensusClass(D, minCor=.2, minDelta=0, return.values=F)}
\arguments{
\item{D}{A dataframe with log2-scaled gene expression profile. Each column is a sample gene expression profile (at least one is required), with one gene value by row. Rownames must be Entrez gene IDs. }
\item{minCor}{Correlation threshold between sample gene expression profile and consensus centroids profiles. A sample showing no correlation above this threshold will remain unclassifed (NA). }
\item{minDelta}{Threshold for the correlation difference between a sample profile and its two nearest centroids. A sample whose two highest correlations with centroids differ by less than this value will remain unclassified (NA).}
\item{return.values}{If set to TRUE the function will return a data frame with detailed correlation values and metrics for each sample. Only resulting sample consensus labels are returned if set to FALSE (Default).}
}
\value{If return.values is TRUE a dataframe is returned with correlation values for each samples (one sample by row). Pearson correlation values are given between each sample gene expression profile and each centroid. The predicted consensus labels are summarized in the consensusClass column and correspond to the nearestCentroid column when the minCor and minDelta conditions are verified.  
}
\author{Aurelie Kamoun}
\keyword{methods}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
