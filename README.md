# consensusMIBC
This package implements a nearest-centroid transcriptomic classifier that assigns class labels according to the consensus molecular classification of Muscle-Invasive Bladder Cancer (Manuscript in preparation). The consensus classification identifies 6 molecular classes : Luminal Papillary (LumP), Luminal Non Specified (LumNS), Luminal Unstable (LumU), Stroma-rich, Basal/Squamous (Ba/Sq), Neuroendocrine-like (NE-like).
Two example data sets are provided to run the classifier.

## Citation
For now, please provide a link to this github repository:
<https://github.com/cit-bioinfo/consensusMIBC>

## Install
You may install this package with [devtools]:

[devtools]: https://github.com/hadley/devtools

```{r}
library(devtools)
devtools::install_github("cit-bioinfo/consensusMIBC")
library(consensusMIBC)
```
## Example
The classifier expect either a single named vector of gene expression values or a dataframe formatted according to the example data sets provided (genes in row, samples in column). Gene names (vector names or dataframe rownames) must be Entrez IDs.

```{r}
data(tcgadat)
res <- getConsensusClass(tcga.dat)

head(res)

#                 consensusClass confidence      LumP     LumNS      LumU Stroma-rich     Ba/Sq   NE-like
#TCGA-2F-A9KO-01A           LumP        Low 0.6176298 0.5735781 0.5641684   0.5890560 0.5754460 0.1820276
#TCGA-2F-A9KP-01A           LumP       High 0.7284300 0.6806556 0.6938315   0.5608757 0.4695459 0.2734213
#TCGA-2F-A9KQ-01A           LumP       High 0.7196396 0.6443850 0.6345979   0.5290542 0.4500750 0.2041000
#TCGA-2F-A9KR-01A          Ba/Sq        Low 0.5980346 0.4938276 0.4900237   0.5126923 0.6274487 0.2116246
#TCGA-2F-A9KT-01A          Ba/Sq       High 0.5772637 0.5144653 0.5335134   0.5393338 0.6615954 0.2865093
#TCGA-2F-A9KW-01A           LumU       High 0.5136967 0.5527899 0.5963426   0.5371054 0.4725626 0.3087334

table(res$consensusClass)

##      Ba/Sq       LumNS        LumP        LumU     NE-like Stroma-rich 
##        153          21         128          53           6          45 
```
