\name{getConsensusClass}
\alias{getConsensusClass}
\title{getConsensusClass}
\description{Nearest-centroid single sample classifier according to the consensus molecular subtypes of muscle-invasive bladder cancer, based on log2-scaled gene expression profile.}
\usage{getConsensusClass(x, minCor=.2)}
\arguments{
\item{x}{A dataframe with log2-scaled gene expression profile. Each column is a sample gene expression profile (at least one is required), with one gene value by row. Rownames must be Entrez gene IDs. }
\item{minCor}{Correlation threshold between sample gene expression profile and consensus centroids profiles. A sample showing no correlation above this threshold will remain unclassifed (NA). }
}
\value{Dataframe with classification results. The consensusClass column returns the predicted consensus class label(s) of the sample(s). The adjusted_pval column returns the Bonferroni corrected (n=6) p-value(s) associated to the Pearson's correlation of the sample(s) with the nearest centroid. The separationLevel ranges from 0 to 1 and gives a measure of how a sample is representative of its consensus class, with 0 meaning the sample is too close to other consensus classes to be confidently assigned its consensus class label, and 1 meaning the sample is very representative of its consensus class and very different from the other consensus classes. This separationLevel is measured as follows : (correlation to nearest centroid - correlation to second nearest centroid) / median difference of sample-to-centroid correlation.
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
