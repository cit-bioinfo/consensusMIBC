getConsensusClass <- function(D, minCor = .2){
  
  data(centroids)
  data(minDelta)
  
  if(is.vector(D)) D <- data.frame(ss = D, row.names = names(D))
  
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
  cor.dat$confidence <- c("High", "Medium")[match(cor.dat$deltaSecondNearest >= minDelta[cor.dat$nearestCentroid], c(TRUE, FALSE))]
  cor.dat$confidence[which(cor.dat$deltaSecondNearest < mean(minDelta))]
  
  try(cor.dat[which(cor.dat$corToNearest < minCor), "consensusClass"] <-  NA)
  try(cor.dat[which(cor.dat$corToNearest < minCor), "confidence"] <-  NA)
  
  cor.dat <- cor.dat[, c("consensusClass" , "confidence" , colnames(centroids))]
  return(cor.dat)
}

