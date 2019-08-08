#' Centroids for the BLCA Consensus Classes
#' 
#' A data.frame containing the centroid coordinates for the six consensus classes of
#' muscle-invasive bladder cancer
#'
#' @author Aurelie Kamoun
#' @references \url{https://www.biorxiv.org/content/10.1101/488460v2}
"centroids"

#' Gene expression data derived from RNA-seq of the TCGA-BLCA cohort
#' 
#' A matrix of gene expression values for 17987 ENTREZ gene IDs and 406
#' samples derived from the TCGA-BLCA cohort.
#'
#' @references \url{https://portal.gdc.cancer.gov/projects/TCGA-BLCA}
"tcga.dat"

#' Gene expression data derived from the Sjodahl et. al. (2017) cohort
#' 
#' A matrix of gene expression values for 22228 ENTREZ gene IDs and 243
#' samples derived from the Sjodahl 2017 cohort.
#'
#' @references \url{https://onlinelibrary.wiley.com/doi/full/10.1002/path.4886}
"sjo.dat"

#' consensusMIBC
#'
#' This package implements a nearest-centroid transcriptomic classifier, that assigns class 
#' labels according to the consensus molecular classification of 
#' Muscle-Invasive Bladder Cancer. For more information, please read the package vignette by
#' calling `vignette("consensusMIBC")`.
#' 
#' @seealso \link[consensusMIBC]{getConsensusClass}, \link[consensusMIBC]{plotCorrelations} 
#' @author Aurelie Kamoun
#' @note This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): 
#' \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}
"_PACKAGE"