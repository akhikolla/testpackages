## ----setup, include=FALSE------------------------------------------------
library(GenomicTools)
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(base.dir=getwd())
knitr::opts_chunk$set(echo = TRUE,
                      dev=c("png"))

## ----eval=FALSE----------------------------------------------------------
#  install.packages("GenomicTools")

## ----eval=FALSE----------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("snpStats")

## ----eval=FALSE----------------------------------------------------------
#  library("devtools")
#  install_github("fischuu/GenomicTools")

## ---- warning=FALSE, error=FALSE,message=FALSE---------------------------
library("GenomicTools")

## ------------------------------------------------------------------------
data("annotTrack")

## ---- comment=NA---------------------------------------------------------
annotTrack

## ----eval=FALSE----------------------------------------------------------
#  ensGTF <- importGTF(file="Homo_sapiens.GRCh38.85.gtf.gz")

## ------------------------------------------------------------------------
data("genotData")

## ---- comment=NA---------------------------------------------------------
genotData

## ----eval=FALSE----------------------------------------------------------
#  ownGenotypes <- importPED(file="myGenotypes.ped", snps="myGenotypes.map")

## ----eval=FALSE----------------------------------------------------------
#  ownGenotypes <- importVCF(file="myGenotypes.vcf")

## ------------------------------------------------------------------------
data("geneEXP")

## ---- comment=NA---------------------------------------------------------
geneEXP[1:5,1:4]

## ------------------------------------------------------------------------
# Make the example data available
  data("annotTrack")   # Standard gtf file, imported with importGTF
  data("geneEXP")      # Matrix with gene expression
  data("genotData")    # An imported Ped/Map filepair, using importPED
  # data("genotDataVCF") # An imported vcf file, using importVCF (too large for Cran)

## ------------------------------------------------------------------------
# Transform gtf to bed format (not necessarily required)
annot.bed <- gtfToBed(annotTrack)

## ----eval=FALSE----------------------------------------------------------
#  
#  # cis-eQTL
#  ###############################################
#  
#  # Most basic cis-eQTL runs:
#  EQTL1 <- eQTL(gex=geneEXP[,1:10], xAnnot = annotTrack, geno= genotData)
#  
#  # Same run, if gtf has been transformed to bed previously
#  EQTL1.1   <- eQTL(gex=geneEXP[,1:10], xAnnot = annot.bed, geno= genotData)
#  
#  # Same run, when the genotype data wasn't loaded and should be loaded
#  # here instead
#  EQTL1.2 <- eQTL(gex=geneEXP[,1:10], xAnnot = annotTrack,
#  +               geno= file.path("Datasets","genotypes.ped"))
#  
#  # Full set of genes, this time filtered with column names
#  EQTL2   <- eQTL(gex=geneEXP, xAnnot = annot.bed, geno= genotData,
#  +               which = colnames(geneEXP)[1:20])
#  
#  # Single vector of gene expression values, underlying gene is specified
#  # in the which option
#  EQTL3   <- eQTL(gex=as.vector(geneEXP[,1]), xAnnot = annot.bed,
#  +               geno= genotData, which="ENSG00000223972")
#  
#  # Same call, but instead of the name the row number in the gtf/bed
#  # file is provided
#  EQTL3.2 <- eQTL(gex=geneEXP[,1], xAnnot = annot.bed, geno= genotData,
#  +               which=1)
#  
#  # The same expression values are now assigned to three different genes
#  EQTL4   <- eQTL(gex=as.vector(geneEXP[,1]), xAnnot = annot.bed,
#  +               geno= genotData, which=1:3)
#  

## ---- eval=FALSE---------------------------------------------------------
#  EQTL1.vcf <- eQTL(gex=geneEXP[,1:10], xAnnot = annotTrack, geno= genotDataVCF)
#  
#  # Same run, if gtf has been transformed to bed previously
#  EQTL1.1.vcf   <- eQTL(gex=geneEXP[,1:10], xAnnot = annot.bed, geno= genotDataVCF)
#  
#  # Same run, when the genotype data wasn't loaded and should be loaded
#  # here instead
#  EQTL1.2.vcf <- eQTL(gex=geneEXP[,1:10], xAnnot = annotTrack,
#  +               geno= file.path("Datasets","genotypes.vcf"))
#  
#  # Full set of genes, this time filtered with column names
#  EQTL2.vcf   <- eQTL(gex=geneEXP, xAnnot = annot.bed, geno= genotDataVCF,
#  +               which = colnames(geneEXP)[1:20])
#  
#  # Single vector of gene expression values, underlying gene is specified
#  # in the which option
#  EQTL3.vcf   <- eQTL(gex=as.vector(geneEXP[,1]), xAnnot = annot.bed,
#  +               geno= genotData.vcf, which="ENSG00000223972")
#  
#  # Same call, but instead of the name the row number in the gtf/bed
#  # file is provided
#  EQTL3.2.vcf <- eQTL(gex=geneEXP[,1], xAnnot = annot.bed, geno= genotDataVCF,
#  +               which=1)
#  
#  # The same expression values are now assigned to three different genes
#  EQTL4.vcf   <- eQTL(gex=as.vector(geneEXP[,1]), xAnnot = annot.bed,
#  +                 geno= genotData.vcf, which=1:3)

## ----message=FALSE, eval=FALSE-------------------------------------------
#  # Same call, but this time is the corresponding column not casted
#  EQTL3.1 <- eQTL(gex=geneEXP[,1] , xAnnot = annot.bed, geno= genotData,
#                  which="ENSG00000223972")

## ----eval=FALSE----------------------------------------------------------
#  # Trans-eQTL
#  ######################################
#  
#  # Trans eQTL for the first and the last gene in our expression matrix
#  EQTL5   <- eQTL(gex=geneEXP[,c(1,1000)] , xAnnot = annot.bed,
#  +               geno= genotData, windowSize = NULL)
#  
#  # Same call, this time distributed to 8 cores (ony available on
#  # Linux computers)
#  EQTL5   <- eQTL(gex=geneEXP[,c(1,1000)] , xAnnot = annot.bed,
#  +               geno= genotData, windowSize = NULL, mc=8)

## ---- eval=FALSE---------------------------------------------------------
#  # Expression values from the first gene are used to test the 100st
#  # gene for trans-eQTL
#  EQTL6   <- eQTL(gex=as.vector(geneEXP[,1]) , xAnnot = annot.bed, geno= genotData, windowSize = NULL, which=100)

## ----fig.width=10, fig.height=10, fig.dev='png', eval=FALSE--------------
#  #png(file="cisEQTL.png", width=685, height=685)
#  plot(EQTL3.1)
#  #dev.off()

## ---- fig.retina = NULL, fig.cap="Example for a cis-eQTL", echo=FALSE----
knitr::include_graphics("./cisEQTL.png")

## ----fig.width=10, fig.height=10, fig.dev='png', eval=FALSE--------------
#  #png(file="transEQTL.png", width=685, height=685)
#  plot(EQTL6)
#  #dev.off()

## ---- fig.retina = NULL, fig.cap="Example for a trans-eQTL", echo=FALSE----
knitr::include_graphics("./transEQTL.png")

## ------------------------------------------------------------------------
# Make the example data available
  data("phenoData")
  data("genotData")

## ---- eval=FALSE---------------------------------------------------------
#  qtl1 <- QTL(pheno=phenoData[,2:3], geno=genotData)

## ----eval=FALSE----------------------------------------------------------
#  # The most basic approach
#    qtl1 <- QTL(pheno=phenoData, geno=genotData)
#  
#  # Use only a named subset of phenotypes
#    qtl2 <- QTL(pheno=phenoData, geno=genotData, which = c("Pheno1", "Pheno4"))
#  
#  # Use a numbers subset of genotypes, distributed to 3 cores
#    qtl2.1 <- QTL(pheno=phenoData, geno=genotData, which = 3:4, mc=3)
#  
#  # Use a single phenotype only
#    qtl2.2 <- QTL(pheno=phenoData, geno=genotData, which = 7)
#  
#  # Same thing, but filtering applied directly to the data
#    qtl3 <- QTL(pheno=phenoData[,5], geno=genotData)
#  
#  # Also a vector input isntead of a matrix is possible
#    qtl3.1 <- QTL(pheno=as.vector(phenoData[,5]), geno=genotData)
#  
#  # The genotype data can be loaded in runtime, without previous step
#    qtl4 <- QTL(pheno=phenoData[,5], geno=file.path("Datasets","genotypes.ped"))

## ----fig.width=10, fig.height=10, fig.dev='png', eval=FALSE--------------
#  # Visualize e.g. the 1st phenotype from previous runs
#  # png(file="QTL1.png", width=685, height=685)
#    plot(qtl1, which=1)
#  # dev.off()

## ---- fig.retina = NULL, fig.cap="Example 1 for a QTL", echo=FALSE-------
knitr::include_graphics("./QTL1.png")

## ----fig.width=10, fig.height=10, fig.dev='png', eval=FALSE--------------
#  # Visualize e.g. the 1st phenotype from previous runs
#  # png(file="QTL2.png", width=685, height=685)
#    plot(qtl1, which=1, genome = "Human68")
#  # dev.off()

## ---- fig.retina = NULL, fig.cap="Example 2 for a QTL", echo=FALSE-------
knitr::include_graphics("./QTL2.png")

## ------------------------------------------------------------------------
data(mdrExample)
mdrSNP <- mdrExample[,1:20]
fit.mdr <- mdr(mdrSNP, mdrExample$Class, fold=3, top=5)
fit.mdr
fit.mdr <- mdr(mdrSNP, mdrExample$Class)
fit.mdr

## ----eval=TRUE-----------------------------------------------------------
data(mdrExample)
mdrSNP.train <- mdrExample[1:350,1:20]
mdrSNP.test <- mdrExample[351:400,1:20]
fit.mdr <- mdr(mdrSNP.train, mdrExample$Class[1:350], fold=2, top=20)
ensResult <- mdrEnsemble(fit.mdr, data = mdrSNP.test)
table(ensResult, mdrExample[351:400,21])

## ----fig.width=10, fig.height=10, fig.dev='png', eval=FALSE--------------
#  #png(file="./MDR.png", width=685, height=685)
#  plot(fit.mdr)
#  #dev.off()

## ---- fig.retina = NULL, fig.cap="Example for a MDR plot", echo=FALSE----
knitr::include_graphics("./MDR.png")

