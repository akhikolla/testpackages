## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BinaryDosage)

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
bd1afile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
bd1ainfo <- getbdinfo(bdfiles = bd1afile)

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
aaf <- unlist(bdapply(bdinfo = bd1ainfo, func = getaaf))

altallelefreq <- data.frame(SNP = bd1ainfo$snps$snpid, aafcalc = aaf)
knitr::kable(altallelefreq, caption = "Information vs Calculated aaf", digits = 3)

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
snp3 <- data.frame(getsnp(bdinfo = bd1ainfo, "1:12000:T:C", FALSE))

knitr::kable(snp3[1:20,], caption = "SNP 1:12000:T:C", digits = 3)

