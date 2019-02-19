\name{getConsensusClass}
\alias{getConsensusClass}
\title{getConsensusClass}
\description{Nearest-centroid single sample classifier according to the consensus molecular subtypes of muscle-invasive bladder cancer, based on log2-scaled gene expression profile.}
\usage{getConsensusClass(x, minCor = .2, gene_id = c("entrezgene", "ensembl_gene_id", "hgnc_symbol")[1])}
\arguments{
\item{x}{Either a single named vector of gene expression values or a dataframe formatted according to the example data sets provided (unique genes in row, samples in column). Gene names (vector names or dataframe rownames) may be supplied as Entrez IDs, Ensembl gene IDs, or HUGO gene symbols. RNA-seq data needs to be log-transformed, for example using log2(normalized counts + 1). }
\item{minCor}{Numeric value specifying a confidence minimal threshold for best Pearson's correlation between sample gene expression profile and consensus centroids profiles. A sample showing no correlation above this threshold will remain unclassifed and prediction results will be set to NA. Default minCor value is 0.2. }
}
\item{gene_id}{Character value specifying the type of gene identifiers used for the names/rownames of x : entrezgene for Entrez IDs, ensembl_gene_id for Ensembl gene IDs, or hgnc_symbol for HUGO gene symbols. Default value is entrezgene. }
}
\value{Dataframe with classification results. The consensusClass column returns the predicted consensus class label(s) of the sample(s). The cor_pval column returns the p-value(s) associated to the Pearson's correlation of the sample(s) with the nearest centroid. The separationLevel ranges from 0 to 1 and gives a measure of how a sample is representative of its consensus class, with 0 meaning the sample is too close to other consensus classes to be confidently assigned its consensus class label, and 1 meaning the sample is very representative of its consensus class and very different from the other consensus classes. This separationLevel is measured as follows : (correlation to nearest centroid - correlation to second nearest centroid) / median difference of sample-to-centroid correlation.
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
