#' Radarchart visualization of correlations between samples and consensus classes
#' 
#' Plot function using radarcharts to visualize global correlation profiles between 
#' samples and consensus classes
#' 
#' @param xres  The results returned by getConsensusClass, either the full dataset or a subset of rows. 
#' @param col.samp Color used to plot sample correlation profiles
#' @param col.med Color used to plot the median correlation profile by consensus class
#' @param axistype see \link[fmsb]{radarchart} documentation
#' @param seg see \link[fmsb]{radarchart} documentation
#' @param pty see \link[fmsb]{radarchart} documentation
#' @param plty see \link[fmsb]{radarchart} documentation
#' @param plwd see \link[fmsb]{radarchart} documentation
#' @param cglty see \link[fmsb]{radarchart} documentation
#' @param cglwd see \link[fmsb]{radarchart} documentation
#' @param cglcol see \link[fmsb]{radarchart} documentation
#' @param centerzero see \link[fmsb]{radarchart} documentation
#' @param ... further arguments passed to \link[fmsb]{radarchart}
#' 
#' @return Plot a single radar chart for a single sample result. 
#' If multiple sample results are used, sample profiles are grouped by consensus class prediction 
#' on 6 radar charts, also showing the median profile by consensus class.
#' 
#' @author Gordon Robertson and Aurelie Kamoun
#' 
#' @seealso \link[fmsb]{radarchart}
#' 
#' @examples 
#'   data(tcgadat)
#'   xres <- getConsensusClass(tcga.dat)
#'   
#'   #-- Example for a single sample plot, return a single radarchart
#'   plotCorrelations(xres[1, ]) 
#'   
#'   #-- Example for multiple samples, returns one radarchart per consensus class
#'   plotCorrelations(xres) 
#' @importFrom fmsb radarchart
#' @export
plotCorrelations <- function(xres, col.samp = "#1414BE", col.med = "red", 
                             axistype = 2, seg = 4, pty = 32, plty = 1, plwd = 2, 
                             cglty = 1, cglwd = 0.5, cglcol = "grey80",
                             centerzero = T, ...){
  #-- Set up classes
  cl <- c("LumP", "LumNS", "LumU", "Stroma-rich", "Ba/Sq", "NE-like")
  
  #-- Graphical parameters
  max.val <- setNames(rep(1, length(cl)), cl)
  min.val <- setNames(rep(0, length(cl)), cl)
  cex.lab <- setNames(rep(0.8, length(cl)), cl)
  
  #-- For single sample
  if(nrow(xres) == 1){
    
    df <- xres[, cl]
    df <- rbind(rep(1, length(cl)), rep(0, length(cl)), xres[, cl])
    fmsb::radarchart(df,
               axistype = axistype, seg = seg, pty = pty, plty = plty, plwd = plwd,
               paxislabels = round(xres[, cl], digits = 2), 
               pcol = col.samp, cglty = cglty, cglwd = cglwd, cglcol = cglcol, 
               vlcex = cex.lab + .4 * as.numeric(cl == xres$consensusClass), # external label sizes
               palcex = cex.lab + .4 * as.numeric(cl == xres$consensusClass), # external label sizes
               centerzero = centerzero, ...)
    mtext(rownames(xres), side = 3, line = 1)
    mtext(paste("Separation Level : ", round(xres$separationLevel, 2), sep = ""), side = 1, line = 1, cex = .8)
    
    #-- For multiple samples
  } else if (nrow(xres) > 1) {
    
    layout(matrix(1:length(cl), nrow = 2, byrow = T))
    for(cln in cl){
      
      med.val <- apply(xres[which(xres$consensusClass == cln), cl], 2, median, na.rm = T)
      med.sep <- median(xres$separationLevel[which(xres$consensusClass == cln)], na.rm = T)
      
      df <- rbind(max.val, min.val, xres[which(xres$consensusClass == cln), cl], med.val)
      if(all(is.na(med.val))) df[nrow(df), ] <- 0
      colnames(df) <- cl
      #alpha <- c(150, 50, 25)[match(findInterval(nrow(df), c(0, 10, 100)), 1:3)]
      alpha <- c("CC", "99", "33")[match(findInterval(nrow(df), c(0, 10, 100)), 1:3)]
      fmsb::radarchart(df, 
                 axistype = axistype, seg = seg, pty = pty, plty = plty, plwd = plwd,
                 paxislabels = round(med.val, digits = 2), axislabcol = col.med,
                 vlcex = cex.lab + .4 * as.numeric(cl == cln), 
                 palcex = cex.lab + .4 * as.numeric(cl == cln), 
                 pcol = c(rep(paste(col.samp, alpha, sep = ""), nrow(df) - 3), col.med),
                 #pcol = c(rep(rgb(20, 20, 190, alpha = alpha , maxColorValue = 255), nrow(df) - 3), "red"),
                 cglty = cglty, cglwd = cglwd, cglcol = cglcol, 
                 centerzero = centerzero, ...)
      mtext(paste(cln, " predictions (n=", nrow(df), ")", sep = ""), side = 3, line = 1)
      mtext(paste("Separation Level (Median): ", round(med.sep, 2) , sep = ""), side = 1, line = 1, cex = .8)
    }
  }
}
  



