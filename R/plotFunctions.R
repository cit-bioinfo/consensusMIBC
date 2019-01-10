plotCorrelations <- function(xres, col.samp = "#1414BE", col.med = "red", 
                             axistype = 2, seg = 4, pty = 32, plty = 1, plwd = 2, 
                             cglty = 1, cglwd = 0.5, cglcol = "grey80",
                             centerzero = T, ...){
  
  cl <- c("LumP", "LumNS", "LumU", "Stroma-rich", "Ba/Sq", "NE-like")
  max.val <- setNames(rep(1, length(cl)), cl)
  min.val <- setNames(rep(0, length(cl)), cl)
  cex.lab <- setNames(rep(0.8, length(cl)), cl)
  
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
      mtext(paste(cln, " predictions (n=", nrow(xres), ")", sep = ""), side = 3, line = 1)
      mtext(paste("Separation Level (Median): ", round(med.sep, 2) , sep = ""), side = 1, line = 1, cex = .8)
    }
  }
}
  



