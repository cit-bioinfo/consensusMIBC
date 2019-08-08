## ----message=FALSE-------------------------------------------------------
#-- Get TCGA data
library(consensusMIBC)
data("tcgadat")

#-- Classify samples
sample_classes <- getConsensusClass(tcga.dat, minCor = .2, gene_id = "entrezgene")

## ----eval=FALSE----------------------------------------------------------
#  head(sample_classes)

## ----echo=FALSE, message=FALSE-------------------------------------------
dt <- head(sample_classes)
dt$cor_pval <- ifelse(dt$cor_pval < 0.001, "< 0.001", dt$cor_pval)
kableExtra::kable_styling(knitr::kable(dt, digits = 3, format = "html"))

## ----fig.width=4, fig.height=4-------------------------------------------
plotCorrelations(sample_classes[1,])

## ----fig.width=8, fig.height=8-------------------------------------------
plotCorrelations(sample_classes)

