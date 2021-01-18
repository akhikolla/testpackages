## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE--------------------------------------------------------
library(BinaryDosage)

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
bdfile1a_woinfo <- tempfile()
vcftobd(vcffiles = vcf1afile, bdfiles = bdfile1a_woinfo)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")
bdfile1a_winfo <- tempfile()
vcftobd(vcffiles = c(vcf1afile, vcf1ainfo), bdfiles = bdfile1a_winfo)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
bdinfo1a_woinfo <- getbdinfo(bdfiles = bdfile1a_woinfo)
bdinfo1a_woinfo$snpinfo

bdinfo1a_winfo <- getbdinfo(bdfiles = bdfile1a_winfo)
knitr::kable(data.frame(bdinfo1a_winfo$snpinfo), caption = "bdinfo1a_winfo$snpinfo")

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
vcf1afile_gz <- system.file("extdata", "set1a.vcf.gz", package = "BinaryDosage")
vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")

bdfile1a_woinfo_gz <- tempfile()
vcftobd(vcffiles = vcf1afile_gz, bdfiles = bdfile1a_woinfo_gz)

bdfile1a_winfo_gz <- tempfile()
vcftobd(vcffiles = c(vcf1afile_gz, vcf1ainfo), bdfiles = bdfile1a_winfo_gz)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------

bdinfo1a_woinfo_gz <- getbdinfo(bdfiles = bdfile1a_woinfo_gz)
bdinfo1a_winfo_gz <- getbdinfo(bdfiles = bdfile1a_winfo_gz)

aaf1a_woinfo <- unlist(bdapply(bdinfo = bdinfo1a_woinfo, getaaf))
aaf1a_winfo <- unlist(bdapply(bdinfo = bdinfo1a_winfo, getaaf))
aaf1a_woinfo_gz <- unlist(bdapply(bdinfo = bdinfo1a_woinfo_gz, getaaf))
aaf1a_winfo_gz <- unlist(bdapply(bdinfo = bdinfo1a_winfo_gz, getaaf))

aaf1a <- data.frame(SNPID = bdinfo1a_woinfo$snps$snpid,
                    aaf1a_woinfo = aaf1a_woinfo,
                    aaf1a_winfo = aaf1a_winfo,
                    aaf1a_woinfo_gz = aaf1a_woinfo_gz,
                    aaf1a_winfo_gz = aaf1a_winfo_gz)

knitr::kable(aaf1a, caption = "Alternate Allele Frequencies", digits = 4)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
vcf2afile <- system.file("extdata", "set2a.vcf", package = "BinaryDosage")
bdfile2a <- tempfile()

vcftobd(vcffiles = vcf2afile, bdfiles = bdfile2a)

bdinfo2a <- getbdinfo(bdfiles = bdfile2a)
snp1_2a <- data.frame(getsnp(bdinfo = bdinfo2a, snp = 1L, dosageonly = FALSE))

snp1 <- cbind(SubjectID = bdinfo2a$samples$sid, snp1_2a)

knitr::kable(snp1[1:10,], caption = "Dosage and Genotype Probabilities")

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
vcf1brsfile <- system.file("extdata", "set1b_rssnp.vcf", package = "BinaryDosage")
bdfile1b.snpid0 <- tempfile()
bdfile1b.snpid1 <- tempfile()
bdfile1b.snpid2 <- tempfile()
bdfile1b.snpidm1 <- tempfile()

vcftobd(vcffiles = vcf1brsfile, bdfiles = bdfile1b.snpid0)
vcftobd(vcffiles = vcf1brsfile, bdfiles = bdfile1b.snpid1, snpidformat = 1)
vcftobd(vcffiles = vcf1brsfile, bdfiles = bdfile1b.snpid2, snpidformat = 2)
vcftobd(vcffiles = vcf1brsfile, bdfiles = bdfile1b.snpidm1, snpidformat = -1)

bdinfo1b.snpid0 <- getbdinfo(bdfiles = bdfile1b.snpid0)
bdinfo1b.snpid1 <- getbdinfo(bdfiles = bdfile1b.snpid1)
bdinfo1b.snpid2 <- getbdinfo(bdfiles = bdfile1b.snpid2)
bdinfo1b.snpidm1 <- getbdinfo(bdfiles = bdfile1b.snpidm1)

snpnames <- data.frame(format0 = bdinfo1b.snpid0$snps$snpid,
                       format1 = bdinfo1b.snpid1$snps$snpid,
                       format2 = bdinfo1b.snpid2$snps$snpid,
                       formatm1 = bdinfo1b.snpidm1$snps$snpid)

knitr::kable(snpnames, caption = "SNP Names by Format")


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
bdfile1a_calcinfo <- tempfile()
vcftobd(vcffiles = vcf1afile, bdfiles = bdfile1a_calcinfo, bdoptions = c("aaf", "maf", "rsq"))

bdcalcinfo <- getbdinfo(bdfile1a_calcinfo)

snpinfo <- data.frame(aaf_info = bdinfo1a_winfo$snpinfo$aaf,
                      aaf_calc = bdcalcinfo$snpinfo$aaf,
                      maf_info = bdinfo1a_winfo$snpinfo$maf,
                      maf_calc = bdcalcinfo$snpinfo$maf,
                      rsq_info = bdinfo1a_winfo$snpinfo$rsq,
                      rsq_calc = bdcalcinfo$snpinfo$rsq)

knitr::kable(snpinfo, caption = "Information vs Calculated Information", digits = 3)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
vcfinfo <- getvcfinfo(vcf1afile, index = TRUE)
aaf <- unlist(vcfapply(vcfinfo = vcfinfo, getaaf))

altallelefreq <- data.frame(SNP = vcfinfo$snps$snpid, aafinfo = aaf1a_winfo, aafcalc = aaf)
knitr::kable(altallelefreq, caption = "Information vs Calculated aaf", digits = 3)


