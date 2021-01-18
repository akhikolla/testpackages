## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BinaryDosage)

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
bd1afile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
bd1bfile <- system.file("extdata", "vcf1b.bdose", package = "BinaryDosage")
bd1file <- tempfile()

bdmerge(mergefiles = bd1file, bdfiles = c(bd1afile, bd1bfile))

bd1ainfo <- getbdinfo(bd1afile)
bd1binfo <- getbdinfo(bd1bfile)
bd1info <- getbdinfo(bd1file)

nrow(bd1ainfo$samples)
nrow(bd1binfo$samples)
nrow(bd1info$samples)

