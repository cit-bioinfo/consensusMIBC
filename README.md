# consensusMIBC
This package implements a nearest-centroid transcriptomic classifier, that assigns class labels according to the consensus molecular classification of Muscle-Invasive Bladder Cancer (Manuscript submitted). The consensus classification identifies 6 molecular classes : Luminal Papillary (LumP), Luminal Non Specified (LumNS), Luminal Unstable (LumU), Stroma-rich, Basal/Squamous (Ba/Sq), Neuroendocrine-like (NE-like).  
Two example data sets are provided to run the classifier.

## Citation
For now, you can cite the following bioRxiv preprint: 
bioRxiv 488460; doi: https://doi.org/10.1101/488460

## Install
You may install this package with [devtools]:

[devtools]: https://github.com/hadley/devtools

```{r}
library(devtools)
devtools::install_github("cit-bioinfo/consensusMIBC")
library(consensusMIBC)
```

## Usage
```{r}
getConsensusClass(x, minCor = .2)
```
where `x` is either a single named vector of gene expression values or a dataframe formatted according to the example data sets provided (unique genes in row, samples in column). Gene names (vector names or dataframe rownames) must be Entrez IDs. RNA-seq data needs to be log-transformed.

and `minCor` is a confidence minimal threshold for best Pearson's correlation. Classifier predictions relying on a correlation lower than `minCor` are set to NA. Default is `0.2`.

## Example
```{r}
data(tcgadat)

# Single sample classification

getConsensusClass(tcga.dat[, 1])

#   consensusClass adjusted_pval separationLevel      LumP     LumNS      LumU Stroma-rich    Ba/Sq   NE-like
#ss           LumP  1.698954e-90       0.6626931 0.6176298 0.5735781 0.5641684    0.589056 0.575446 0.1820276

# Classification of a matrix of samples

res <- getConsensusClass(tcga.dat)
head(res)

#                 consensusClass adjusted_pval separationLevel      LumP     LumNS      LumU Stroma-rich     Ba/Sq   NE-like
#TCGA-2F-A9KO-01A           LumP  1.698954e-90       0.6626931 0.6176298 0.5735781 0.5641684   0.5890560 0.5754460 0.1820276
#TCGA-2F-A9KP-01A           LumP 8.546007e-142       0.3213549 0.7284300 0.6806556 0.6938315   0.5608757 0.4695459 0.2734213
#TCGA-2F-A9KQ-01A           LumP 8.031287e-137       0.5460606 0.7196396 0.6443850 0.6345979   0.5290542 0.4500750 0.2041000
#TCGA-2F-A9KR-01A          Ba/Sq  3.288756e-94       0.2368501 0.5980346 0.4938276 0.4900237   0.5126923 0.6274487 0.2116246
#TCGA-2F-A9KT-01A          Ba/Sq 3.171797e-108       0.6737278 0.5772637 0.5144653 0.5335134   0.5393338 0.6615954 0.2865093
#TCGA-2F-A9KW-01A           LumU  6.990305e-83       0.6139238 0.5136967 0.5527899 0.5963426   0.5371054 0.4725626 0.3087334

table(res$consensusClass)

##      Ba/Sq       LumNS        LumP        LumU     NE-like Stroma-rich 
##        153          21         128          53           6          45 
```
The classifier returns a dataframe with 9 columns :  

`consensusClass` returns the consensus calls for each sample. Calls are set to NA for low confidence predictions (maximal correlation is below the given `minCor` parameter).  

`adjusted_pval` returns the Bonferroni corrected (n=6) p-value(s) associated to the Pearson's correlation of the sample(s) with the nearest centroid.  

`separationLevel` gives a measure of how a sample is representative of its consensus class. It ranges from 0 to 1, with 0 meaning the sample is too close to other consensus classes to be confidently assigned one consensus class label, and 1 meaning the sample is highly representative of its consensus class and well separated from the other consensus classes. The separationLevel is measured as follows : (correlation to nearest centroid - correlation to second nearest centroid) / median difference of sample-to-centroid correlation.  

The 6 other columns return the Pearson's correlation between each sample and each consensus class.

