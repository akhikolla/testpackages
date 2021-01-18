## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  screenshot.force = FALSE, 
  echo = TRUE,
  rows.print = 5,
  message = FALSE,
  warning = FALSE)

## ----data_load----------------------------------------------------------------
library(PLNmodels)
data(trichoptera)

## ----trichoptera_structure----------------------------------------------------
str(trichoptera, max.level = 1)

## ----covariates_overview------------------------------------------------------
names(trichoptera$Covariate)

## ----first_try, eval = FALSE--------------------------------------------------
#  PLN(Abundance ~ Wind + Pressure, data = trichoptera)

## ----prepare_data_first_look--------------------------------------------------
trichoptera2 <- prepare_data(counts     = trichoptera$Abundance, 
                             covariates = trichoptera$Covariate)
str(trichoptera2)

## ----compute_offset-----------------------------------------------------------
## same as compute_offset(trichoptera$Abundance, offset = "TSS")
compute_offset(trichoptera$Abundance) 

## ----other_offsets, warning=TRUE, error = TRUE, results='hide'----------------
compute_offset(trichoptera$Abundance, "CSS")
compute_offset(trichoptera$Abundance, "RLE")
compute_offset(trichoptera$Abundance, "GMPR")

## ----pseudocounts-------------------------------------------------------------
compute_offset(trichoptera$Abundance, "RLE", pseudocounts = 1)

## ----prepare_data_other_offset------------------------------------------------
str(prepare_data(trichoptera$Abundance, 
             trichoptera$Covariate, 
             offset = "RLE", pseudocounts = 1))

## ----trim_down_samples--------------------------------------------------------
nrow(prepare_data(trichoptera$Abundance[-1, ], ## remove first sample
                  trichoptera$Covariate[-49,]  ## remove last sample
                  ))

## ----import_biom,eval = FALSE-------------------------------------------------
#  ## If biomformat is not installed, uncomment the following lines
#  # if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  #   install.packages("BiocManager")
#  # }
#  # BiocManager::install("biomformat")
#  library(biomformat)
#  biomfile <- system.file("extdata", "rich_dense_otu_table.biom", package = "biomformat")
#  biom <- biomformat::read_biom(biomfile)
#  ## extract counts
#  counts <- as(biomformat::biom_data(biom), "matrix")
#  ## extract covariates (or prepare your own)
#  covariates <- biomformat::sample_metadata(biom)
#  ## prepare data
#  my_data <- prepare_data(counts = counts, covariates = covariates)
#  str(my_data)

## ----import_phyloseq, eval = FALSE--------------------------------------------
#  ## If biomformat is not installed, uncomment the following lines
#  # if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  #   install.packages("BiocManager")
#  # }
#  # BiocManager::install("phyloseq")
#  library(phyloseq)
#  data("enterotype")
#  ## extract counts
#  counts <- as(phyloseq::otu_table(enterotype), "matrix")
#  ## extract covariates (or prepare your own)
#  covariates <- phyloseq::sample_data(enterotype)
#  ## prepare data
#  my_data <- prepare_data(counts = counts, covariates = covariates)
#  str(my_data)

