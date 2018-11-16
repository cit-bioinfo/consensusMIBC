# Rd
# description >> Nearest-centroid single sample classifier according to the consensus molecular subtypes of muscle-invasive bladder cancer, based on log2-scaled gene expression profile.
# argument
# item >> D >> A dataframe with log2-scaled gene expression profile. Each column is a sample gene expression profile (at least one is required), with one gene value by row. Rownames must be Entrez gene IDs. 
# item >> minCor >> Correlation threshold between sample gene expression profile and consensus centroids profiles. A sample showing no correlation above this threshold will remain unclassifed (NA). 
# item >> minDelta >> Threshold for the correlation difference between a sample profile and its two nearest centroids. A sample whose two highest correlations with centroids differ by less than this value will remain unclassified (NA).
# item >> return.values >> If set to TRUE the function will return a data frame with detailed correlation values and metrics for each sample. Only resulting sample consensus labels are returned if set to FALSE (Default).
# value >> If return.values is TRUE a dataframe is returned with correlation values for each samples (one sample by row). Pearson correlation values are given between each sample gene expression profile and each centroid. The predicted consensus labels are summarized in the consensusClass column and correspond to the nearestCentroid column when the minCor and minDelta conditions are verified.  
# author >> Aurelie Kamoun
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end

getConsensusClass <- function(D, minCor = .2, minDelta = 0, return.values = F){
  
  data(centroids)
  gkeep <- intersect(rownames(centroids), rownames(D))
  if (length(gkeep) == 0) stop("empty intersection between profiled genes and the genes used for consensus classification. Make sure that rownames(D) are Entrez gene IDs")
  if (length(gkeep) < 0.5 * nrow(centroids)) warning("input gene expression profile(s) include less than half of the genes used for consensus classification. Results may not be relevant") 
  cor.dat <- as.data.frame(cor(D[gkeep, ], centroids[gkeep, ], use = "complete.obs"))
  
  cor.dat$nearestCentroid <- apply(cor.dat, 1, function(x){colnames(centroids)[which.max(x)]})
  cor.dat$corToNearest <- apply(cor.dat[, colnames(centroids)], 1, max)
  cor.dat$deltaSecondNearest <- sapply(1:nrow(cor.dat), function(i){
    w <- setdiff(colnames(centroids), cor.dat[i, "nearestCentroid"])
    cor.dat[i, "corToNearest"] - max(cor.dat[i, w])
  })
  
  cor.dat$consensusClass <- cor.dat$nearestCentroid
  try(cor.dat[which(cor.dat$corToNearest < minCor | cor.dat$deltaSecondNearest < minDelta), "predicted"] <-  NA)
  
  if(return.values) return(cor.dat) else return(setNames(cor.dat$consensusClass, rownames(cor.dat))) 
}

