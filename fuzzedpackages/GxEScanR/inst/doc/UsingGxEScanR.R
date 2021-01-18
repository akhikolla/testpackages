## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(GxEScanR)

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
covdatafile <- system.file("extdata", "covdata.rds", package = "GxEScanR")
covdata <- readRDS(covdatafile)

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(covdata[1:5,], caption = "First 5 Subjects")

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
bdinfofile <- system.file("extdata", "pdata_4_1.bdinfo", package = "GxEScanR")
bdinfo <- readRDS(bdinfofile)
bdinfo$filename <- system.file("extdata", "pdata_4_1.bdose", package = "GxEScanR")

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
modeldf <- readRDS(system.file("extdata", "models.rds", package = "GxEScanR"))
knitr::kable(modeldf, caption = "Models Fit")

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(modeldf[1,], caption = "Model Fit")

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
lingwas1 <- gwas(data = covdata,
                 bdinfo = bdinfo,
                 binary = FALSE)

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(lingwas1, caption = "Linear Regression GWAS")

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
outfile <- tempfile()
lingwas2 <- gwas(data = covdata,
                 bdinfo = bdinfo,
                 outfile = outfile,
                 binary = FALSE)
lingwas2
lingwas2 <- read.table(outfile, header = TRUE, sep ='\t')

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(lingwas2, caption = "Linear Regression GWAS")

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(modeldf[1:2,], caption = "Models Fit")

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
lingweis1 <- gweis(data = covdata,
                   bdinfo = bdinfo,
                   minmaf = 0.2,
                   binary = FALSE)

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(lingweis1, caption = "Linear Regression GWEIS")

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
skipfile = tempfile()
lingweis2 <- gweis(data = covdata,
                   bdinfo = bdinfo,
                   skipfile = skipfile,
                   minmaf = 0.2,
                   binary = FALSE)

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(lingweis2, caption = "Linear Regression GWEIS")
skipsnps <- read.table(skipfile, header = TRUE, sep = '\t')

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
knitr::kable(skipsnps, caption = "Skipped SNPs")

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
reasondf <- readRDS(system.file("extdata", "skipreason.rds", package = "GxEScanR"))
knitr::kable(reasondf, caption = "Skipped Reasons")

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(modeldf[1,], caption = "Model Fit")

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
loggwas1 <- gwas(data = covdata,
                 bdinfo = bdinfo,
                 blksize = 2,
                 binary = TRUE)

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(loggwas1, caption = "Logistic Regression GWAS", digits = 4)

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
defaultdf <- readRDS(system.file("extdata", "defaultblk.rds", package = "GxEScanR"))
knitr::kable(defaultdf, caption = "Default blksize")

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(modeldf, caption = "Models Fit")

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
loggweis1 <- gweis(data = covdata,
                   bdinfo = bdinfo,
                   snps = 1:2,
                   binary = TRUE)

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(loggweis1, caption = "Logistic Regression GWEIS", digits = 4)

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
covdata2 <- covdata
covdata2$e <- covdata2$e + 1
loggweis2 <- gweis(data = covdata2,
                   bdinfo = bdinfo,
                   snps = c("1:10001", "1:10002"))

## ---- eval = T, echo = F, message = F, warning = F, tidy = T------------------
knitr::kable(loggweis2, caption = "Logistic Regression GWEIS", digits = 4)

