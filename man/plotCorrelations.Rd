\name{plotCorrelations}
\alias{plotCorrelations}
\title{Radarchart visualization of correlations between samples and consensus classes}
\description{
Plot function using radarcharts to visualize global correlation profiles between samples and consensus classes
}
\usage{
plotCorrelations(xres, col.samp = "#1414BE", col.med = "red", axistype = 2, seg = 4, pty = 32, plty = 1, plwd = 2, cglty = 1, cglwd = 0.5, cglcol = "grey80", centerzero = T, ...)
}
\arguments{
  \item{xres}{
The results returned by getConsensusClass, either the full dataset or a subset of rows. 
}
  \item{col.samp}{
Color used to plot sample correlation profiles
}
  \item{col.med}{
Color used to plot the median correlation profile by consensus class
}
  \item{axistype}{
see \link[fmsb]{radarchart} documentation
}
  \item{seg}{
see \link[fmsb]{radarchart} documentation
}
  \item{pty}{
see \link[fmsb]{radarchart} documentation
}
  \item{plty}{
see \link[fmsb]{radarchart} documentation
}
  \item{plwd}{
see \link[fmsb]{radarchart} documentation
}
  \item{cglty}{
see \link[fmsb]{radarchart} documentation
}
  \item{cglwd}{
see \link[fmsb]{radarchart} documentation
}
  \item{cglcol}{
see radarchart function documentation
}
  \item{centerzero}{
see \link[fmsb]{radarchart} documentation
}
  \item{\dots}{
Further arguments passed to radarchart()}
}
\value{
Plot a single radar chart for a single sample result. If multiple sample results are used, sample profiles are grouped by consensus class prediction on 6 radar charts, also showing the median profile by consensus class.
}
\author{
Gordon Robertson and Aurelie Kamoun
}
\seealso{
\link[fmsb]{radarchart}
}
\examples{
data(tcgadat)
xres <- getConsensusClass(tcga.dat)

#Example for a single sample plot, return a single radarchart
plotCorrelations(xres[1, ]) 

#Example for multiple samples, returns one radarchart per consensus class
plotCorrelations(xres) 

}
\keyword{dplot}% use one of  RShowDoc("KEYWORDS")
