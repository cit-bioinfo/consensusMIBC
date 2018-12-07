\name{getConsensusClass}
\alias{getConsensusClass}
\title{getConsensusClass}
\description{Nearest-centroid single sample classifier according to the consensus molecular subtypes of muscle-invasive bladder cancer, based on log2-scaled gene expression profile.}
\usage{getConsensusClass(D, minCor=.2)}
\arguments{
\item{D}{A dataframe with log2-scaled gene expression profile. Each column is a sample gene expression profile (at least one is required), with one gene value by row. Rownames must be Entrez gene IDs. }
\item{minCor}{Correlation threshold between sample gene expression profile and consensus centroids profiles. A sample showing no correlation above this threshold will remain unclassifed (NA). }
}
\value{Dataframe with classification results. Predicted consensus class labels are given in the first consensusClass column. The confidence level of the predictions is given in the second column and is computed for each sample based on the correlation difference between the two highest correlated centroids.
The Pearson's correlation values for each sample and each centroid are detailed in the additional columns.
consensusClass predictions are set to NA if the minCor condition is not verified.
}
\author{Aurelie Kamoun}
\keyword{methods}
\examples{
data(tcgadat)
getConsensusClass(tcga.dat)
}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
