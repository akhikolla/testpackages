---
title: "Using the PAC package"
date: "2018-03-03"
output: rmarkdown::html_vignette
bibliography: references.bib
vignette: >
  %\VignetteIndexEntry{Using the PAC package}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

This vignette illustrates the basic usage of the PAC package for R. 

Biology Example
---------------------------

The PAC-MAN data analysis pipeline can be applied to mass-cytometry (CyTOF) data analysis. In this case, the user reads in the example data files (already saved as the Rdata format) subsetted from Bendall et al., 2011 and goes through the data analysis pipeline.



Load the required R packages


```r
library(PAC)
```

Construct the sampleIDs vector to analyze the data


```r
sampleIDs<-c("Basal", "BCR", "IL7")
```

Partition, cluster into desired number of subpopulations, and output subpopulation mutual information networks


```r
samplePass(sampleIDs, dim_subset=NULL, hyperrectangles=35, num_PACSupop=25, num_networkEdge=25, max.iter=50)
```

```
## Input Data: 2650 by 18
## Partition method: Discrepancy based partition
## Maximum level: 35
## partition completed
## [1] "Initial Clustering..."
## [1] "Merging..."
```

```
## Input Data: 3537 by 18
## Partition method: Discrepancy based partition
## Maximum level: 35
## partition completed
## [1] "Initial Clustering..."
## [1] "Merging..."
```

```
## Input Data: 3813 by 18
## Partition method: Discrepancy based partition
## Maximum level: 35
## partition completed
## [1] "Initial Clustering..."
## [1] "Merging..."
```

Multiple Alignments of Networks


```r
clades_network_only<-MAN(sampleIDs, num_PACSupop=25, smallSubpopCutoff=100, k_clades=5)
```

Refine the PAC labels with multiple alignments of networks representative labels for clades


```r
refineSubpopulationLabels(sampleIDs,dim_subset=NULL, clades_network_only, expressionGroupClamp=5)
```

Draw clade/representative mutual information networks


```r
getRepresentativeNetworks(sampleIDs, dim_subset=NULL, SubpopSizeFilter=200, num_networkEdge=25)
```

Obtain annotations of subpopulations


```r
aggregateMatrix_withAnnotation<-annotateClades(sampleIDs, topHubs=4)
head(aggregateMatrix_withAnnotation)
```

```
##                      Annotation ClusterID SampleID   pPLCgamma2    pSTAT5
## 1         pNFkB-pSTAT3-pH3-Ki67    clade1    Basal  0.868677090 1.4210433
## 2    IkBalpha-pP38-Ki67-pERK1.2    clade2    Basal  0.417650152 0.7512015
## 3    pNFkB-Ki67-pMAPKAPK2-pSHP2    clade3    Basal -0.009067766 0.1439336
## 4 pP38-pSrcFK-pSTAT3-pZAP70.Syk    clade1      BCR  1.414513278 1.4852767
## 5     pSrcFK-pCREB-pSHP2-pSTAT3    clade2      BCR  0.331270112 0.7172607
## 6    pSrcFK-pNFkB-pBtk.Itk-Ki67    clade4      BCR  0.912350993 1.5864291
##        Ki67     pSHP2   pERK1.2 pMAPKAPK2 pZAP70.Syk    pSTAT3       pSLP
## 1 1.6045171 1.0443365 1.4661011 1.5628676 1.43355321 2.1289258 0.92427696
## 2 0.7591787 0.6073356 0.9891308 0.9838982 0.68705268 1.4899395 0.46984991
## 3 0.5970280 0.0675446 0.1806816 0.2359121 0.05943879 0.1938544 0.06218395
## 4 1.6055084 1.3134441 1.5878684 1.6812649 1.74982011 2.3305502 1.16167123
## 5 0.2062699 0.5510241 0.7372884 1.0448433 0.73709113 1.4776139 0.33423706
## 6 3.5343868 0.9396143 1.1513465 1.7927573 1.14747844 2.1468960 0.69638449
##       pNFkB IkBalpha       pH3      pP38  pBtk.Itk       pS6    pSrcFK
## 1 2.5118414 1.930353 2.0931644 2.3369804 2.6054058 1.5696688 2.9274172
## 2 1.7826376 1.426169 0.9924652 1.5957836 1.7199279 0.6538710 1.8988397
## 3 0.4461451 0.206577 0.2719665 0.2130849 0.5528641 0.1526701 0.2939679
## 4 2.5682406 1.835245 2.6952983 2.8117867 2.5400931 1.8602256 3.1488323
## 5 1.7574310 1.692790 1.2809057 1.7151455 1.9010833 0.9869888 2.2340015
## 6 2.4373202 1.627335 2.1106601 2.3534429 5.0717325 2.3080344 2.1823834
##        pCREB      pCrkL count
## 1 1.80653721 1.03088589  1481
## 2 0.86030693 0.57582728   507
## 3 0.08732748 0.06956507   458
## 4 1.91621956 1.04601166   990
## 5 1.26005371 0.47350443  1341
## 6 1.48173361 0.66764625   437
```


Obtain heatmap input and plot heatmap


```r
cladeProportionMatrix<-heatmapInput(aggregateMatrix_withAnnotation)
heatmap(as.matrix(cladeProportionMatrix))
```

![plot of chunk Heatmap input and plot heatmap](figure/Heatmap input and plot heatmap-1.png)


Append subpopulation proportions for each sample in the annotation matrix


```r
annotationMatrix_prop<-annotationMatrix_withSubpopProp(aggregateMatrix_withAnnotation)
head(annotationMatrix_prop)
```

```
##                      Annotation ClusterID SampleID   pPLCgamma2    pSTAT5
## 1         pNFkB-pSTAT3-pH3-Ki67    clade1    Basal  0.868677090 1.4210433
## 2    IkBalpha-pP38-Ki67-pERK1.2    clade2    Basal  0.417650152 0.7512015
## 3    pNFkB-Ki67-pMAPKAPK2-pSHP2    clade3    Basal -0.009067766 0.1439336
## 4 pP38-pSrcFK-pSTAT3-pZAP70.Syk    clade1      BCR  1.414513278 1.4852767
## 5     pSrcFK-pCREB-pSHP2-pSTAT3    clade2      BCR  0.331270112 0.7172607
## 6    pSrcFK-pNFkB-pBtk.Itk-Ki67    clade4      BCR  0.912350993 1.5864291
##        Ki67     pSHP2   pERK1.2 pMAPKAPK2 pZAP70.Syk    pSTAT3       pSLP
## 1 1.6045171 1.0443365 1.4661011 1.5628676 1.43355321 2.1289258 0.92427696
## 2 0.7591787 0.6073356 0.9891308 0.9838982 0.68705268 1.4899395 0.46984991
## 3 0.5970280 0.0675446 0.1806816 0.2359121 0.05943879 0.1938544 0.06218395
## 4 1.6055084 1.3134441 1.5878684 1.6812649 1.74982011 2.3305502 1.16167123
## 5 0.2062699 0.5510241 0.7372884 1.0448433 0.73709113 1.4776139 0.33423706
## 6 3.5343868 0.9396143 1.1513465 1.7927573 1.14747844 2.1468960 0.69638449
##       pNFkB IkBalpha       pH3      pP38  pBtk.Itk       pS6    pSrcFK
## 1 2.5118414 1.930353 2.0931644 2.3369804 2.6054058 1.5696688 2.9274172
## 2 1.7826376 1.426169 0.9924652 1.5957836 1.7199279 0.6538710 1.8988397
## 3 0.4461451 0.206577 0.2719665 0.2130849 0.5528641 0.1526701 0.2939679
## 4 2.5682406 1.835245 2.6952983 2.8117867 2.5400931 1.8602256 3.1488323
## 5 1.7574310 1.692790 1.2809057 1.7151455 1.9010833 0.9869888 2.2340015
## 6 2.4373202 1.627335 2.1106601 2.3534429 5.0717325 2.3080344 2.1823834
##        pCREB      pCrkL count subpop_proportion
## 1 1.80653721 1.03088589  1481             60.55
## 2 0.86030693 0.57582728   507             20.73
## 3 0.08732748 0.06956507   458             18.72
## 4 1.91621956 1.04601166   990             29.43
## 5 1.26005371 0.47350443  1341             39.86
## 6 1.48173361 0.66764625   437             12.99
```




