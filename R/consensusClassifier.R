#' Bladder Cancer Consensus Class Inference
#' 
#' Nearest-centroid single sample classifier according to the consensus molecular subtypes 
#' of muscle-invasive bladder cancer, based on log2-scaled gene expression profile.
#'
#' @param x Either a single named vector of gene expression values or a dataframe 
#' formatted according to the example data sets provided (unique genes in row, samples in column). 
#' Gene names (vector names or dataframe rownames) may be supplied as Entrez IDs, Ensembl gene IDs, 
#' or HUGO gene symbols. RNA-seq data needs to be log-transformed, for example using 
#' log2(normalized counts + 1). 
#' @param minCor Numeric value specifying a confidence minimal threshold for best Pearson's correlation 
#' between sample gene expression profile and consensus centroids profiles. A sample showing no correlation 
#' above this threshold will remain unclassifed and prediction results will be set to NA. Default minCor 
#' value is 0.2.
#' @param gene_id Character value specifying the type of gene identifiers used for the names/rownames of 
#' x : entrezgene for Entrez IDs, ensembl_gene_id for Ensembl gene IDs, or hgnc_symbol for HUGO 
#' gene symbols. Default value is entrezgene.
#'
#' @return a Dataframe with classification results. The consensusClass column returns the predicted 
#' consensus class label(s) of the sample(s). The cor_pval column returns the p-value(s) associated 
#' to the Pearson's correlation of the sample(s) with the nearest centroid. The separationLevel ranges 
#' from 0 to 1 and gives a measure of how a sample is representative of its consensus class, with 0 
#' meaning the sample is too close to other consensus classes to be confidently assigned its consensus 
#' class label, and 1 meaning the sample is very representative of its consensus class and very different 
#' from the other consensus classes. This separationLevel is measured as follows : (correlation to 
#' nearest centroid - correlation to second nearest centroid) / median difference of sample-to-centroid 
#' correlation. The Pearson's correlation values for each sample and each centroid are detailed in 
#' the additional columns. consensusClass predictions are set to NA if the minCor condition is not verified.
#'
#' @aliases getConsensusClass
#' @author Aurelie Kamoun
#' 
#' @examples 
#' data(tcgadat)
#' getConsensusClass(tcga.dat)
#'
#' @note This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): 
#' \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}
#' 
#' @importFrom graphics layout mtext
#' @importFrom stats cor cor.test median setNames
#' @importFrom utils data
#' @export

getConsensusClass <- function(x, minCor = .2, gene_id = c("entrezgene", "ensembl_gene_id", "hgnc_symbol")[1]){
  
  centroids <- consensusMIBC::centroids
  lev.cs <- c("LumP", "LumNS", "LumU", "Stroma-rich", "Ba/Sq", "NE-like")
  
  if(is.vector(x)) {
    if(is.null(names(x))) stop("Input vector of gene expression is missing names.\n The names must be the type of gene identifiers specified by the gene_id argument.")
    x <- data.frame(ss = x, row.names = names(x))
  }
  
  gkeep <- intersect(centroids[, gene_id], rownames(x))
  if (length(gkeep) == 0) stop("Empty intersection between profiled genes and the genes used for consensus classification.\n Make sure that gene names correspond to the type of identifiers specified by the gene_id argument")
  if (length(gkeep) < 0.5 * nrow(centroids)) warning("Input gene expression profile(s) include(s) less than half of the genes used for consensus classification. Results may not be relevant") 
  cor.dat <- as.data.frame(cor(x[gkeep, ], centroids[match(gkeep, centroids[, gene_id]), lev.cs], use = "complete.obs"), row.names = colnames(x))

  # Best correlated centroid
  cor.dat$nearestCentroid <- apply(cor.dat, 1, function(y){lev.cs[which.max(y)]})
  cor.dat$corToNearest <- apply(cor.dat[, lev.cs], 1, max)
  cor.dat$cor_pval <- sapply(colnames(x), function(smp){
    cor.test(x[gkeep, smp], centroids[match(gkeep, centroids[, gene_id]), cor.dat[smp, "nearestCentroid"]])$p.value
    })

  # Separation level metrics
  cor.dat$deltaSecondNearest <- apply(cor.dat$corToNearest - cor.dat[, lev.cs], 1, function(x){sort(x)[2]})
  cor.dat$deltaMed <- apply(cor.dat$corToNearest - cor.dat[, lev.cs], 1, median)
  cor.dat$separationLevel <- cor.dat$deltaSecondNearest/cor.dat$deltaMed
  
  cor.dat$consensusClass <- cor.dat$nearestCentroid
  
  # Set to NA if best correlation < minCor
  try(cor.dat[which(cor.dat$corToNearest < minCor), "consensusClass"] <-  NA)
  try(cor.dat[which(cor.dat$corToNearest < minCor), "separationLevel"] <-  NA)
  
  cor.dat <- cor.dat[, c("consensusClass" , "cor_pval", "separationLevel", lev.cs)]
  return(cor.dat)
}

