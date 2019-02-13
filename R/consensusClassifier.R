getConsensusClass <- function(x, minCor = .2, gene_id = c("entrezgene", "ensembl_gene_id", "hgnc_symbol")[1]){
  
  data(centroids)
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

