getConsensusClass <- function(x, minCor = .2){
  
  data(centroids)
  
  if(is.vector(x)) x <- data.frame(ss = x, row.names = names(x))
  
  gkeep <- intersect(rownames(centroids), rownames(x))
  if (length(gkeep) == 0) stop("empty intersection between profiled genes and the genes used for consensus classification. Make sure that rownames(D) are Entrez gene IDs")
  if (length(gkeep) < 0.5 * nrow(centroids)) warning("input gene expression profile(s) include less than half of the genes used for consensus classification. Results may not be relevant") 
  cor.dat <- as.data.frame(cor(x[gkeep, ], centroids[gkeep, ], use = "complete.obs"), row.names = colnames(x))

  # Best correlated centroid
  cor.dat$nearestCentroid <- apply(cor.dat, 1, function(y){colnames(centroids)[which.max(y)]})
  cor.dat$corToNearest <- apply(cor.dat[, colnames(centroids)], 1, max)
  cor.dat$adjusted_pval <- ncol(centroids) * sapply(colnames(x), function(smp){
    cor.test(x[gkeep, smp], centroids[gkeep, cor.dat[smp, "nearestCentroid"]])$p.value
    })

  # Separation level metrics
  cor.dat$deltaSecondNearest <- apply(cor.dat$corToNearest - cor.dat[, colnames(centroids)], 1, function(x){sort(x)[2]})
  cor.dat$deltaMed <- apply(cor.dat$corToNearest - cor.dat[, colnames(centroids)], 1, median)
  cor.dat$separationLevel <- cor.dat$deltaSecondNearest/cor.dat$deltaMed
  
  cor.dat$consensusClass <- cor.dat$nearestCentroid
  
  # Set to NA if best correlation < minCor
  try(cor.dat[which(cor.dat$corToNearest < minCor), "consensusClass"] <-  NA)
  try(cor.dat[which(cor.dat$corToNearest < minCor), "separationLevel"] <-  NA)
  
  cor.dat <- cor.dat[, c("consensusClass" , "adjusted_pval", "separationLevel", colnames(centroids))]
  return(cor.dat)
}

