# consensusMIBC
This package implements a nearest-centroid transcriptomic classifier that assigns class labels according to the consensus molecular classification of Muscle-Invasive Bladder Cancer (Manuscript in preparation). The consensus classification identifies 6 molecular classes : Luminal Papillary (LumP), Luminal Non Specified (LumNS), Luminal Unstable (LumU), Stroma-rich, Basal, Neuroendocrine-like (NE-like).
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
cl <- getConsensusClass(tcga.dat)

head(cl)
## TCGA-2F-A9KO-01A TCGA-2F-A9KP-01A TCGA-2F-A9KQ-01A TCGA-2F-A9KR-01A TCGA-2F-A9KT-01A TCGA-2F-A9KW-01A 
##       "LumP"           "LumP"           "LumP"          "Basal"          "Basal"           "LumU"

table(cl)
## cl
##      Basal       LumNS        LumP        LumU     NE-like Stroma-rich 
##        153          21         128          53           6          45 

all.dat <- getConsensusClass(tcga.dat, return.values = T)
## Basal      LumP     LumNS      LumU   NE-like Stroma-rich nearestCentroid
## TCGA-2F-A9KO-01A 0.5754460 0.6176298 0.5735781 0.5641684 0.1820276   0.5890560            LumP
## TCGA-2F-A9KP-01A 0.4695459 0.7284300 0.6806556 0.6938315 0.2734213   0.5608757            LumP
## TCGA-2F-A9KQ-01A 0.4500750 0.7196396 0.6443850 0.6345979 0.2041000   0.5290542            LumP
## TCGA-2F-A9KR-01A 0.6274487 0.5980346 0.4938276 0.4900237 0.2116246   0.5126923           Basal
## TCGA-2F-A9KT-01A 0.6615954 0.5772637 0.5144653 0.5335134 0.2865093   0.5393338           Basal
## TCGA-2F-A9KW-01A 0.4725626 0.5136967 0.5527899 0.5963426 0.3087334   0.5371054            LumU
##                  corToNearest deltaSecondNearest consensusClass predicted
## TCGA-2F-A9KO-01A    0.6176298         0.02857384           LumP        NA
## TCGA-2F-A9KP-01A    0.7284300         0.03459846           LumP        NA
## TCGA-2F-A9KQ-01A    0.7196396         0.07525455           LumP        NA
## TCGA-2F-A9KR-01A    0.6274487         0.02941414          Basal        NA
## TCGA-2F-A9KT-01A    0.6615954         0.08433170          Basal        NA
## TCGA-2F-A9KW-01A    0.5963426         0.04355272           LumU        NA
```
