## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE--------------------------------------------------------
library(BinaryDosage)

## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
gen3bchrfile <- system.file("extdata", "set3b.chr.imp", package = "BinaryDosage")
sample3bfile <- system.file("extdata", "set3b.sample", package = "BinaryDosage")
bdfile3b_chr <- tempfile()
gentobd(genfiles = c(gen3bchrfile, sample3bfile), bdfiles = bdfile3b_chr)

bdinfo3b_chr <- getbdinfo(bdfiles = bdfile3b_chr)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
gen1bfile <- system.file("extdata", "set1b.imp", package = "BinaryDosage")
bdfile1b <- tempfile()
gentobd(genfiles = gen1bfile,
        bdfiles = bdfile1b,
        snpcolumns = c(1L, 3L, 2L, 4L, 5L),
        header = TRUE)

bdinfo1b <- getbdinfo(bdfiles = bdfile1b)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
gen3bfile <- system.file("extdata", "set3b.imp", package = "BinaryDosage")
sample3bfile <- system.file("extdata", "set3b.sample", package = "BinaryDosage")
bdfile3b <- tempfile()
gentobd(genfiles = c(gen3bfile, sample3bfile),
        bdfiles = bdfile3b,
        snpcolumns = c(0L,2L:5L))

bdinfo3b <- getbdinfo(bdfiles = bdfile3b)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
gen4bfile <- system.file("extdata", "set4b.imp", package = "BinaryDosage")
sample4bfile <- system.file("extdata", "set4b.sample", package = "BinaryDosage")
bdfile4b <- tempfile()
gentobd(genfiles = c(gen4bfile, sample4bfile),
        bdfiles = bdfile4b,
        snpcolumns = c(1L,2L,4L,5L,6L),
        startcolumn = 7L,
        impformat = 2L)

bdinfo4b <- getbdinfo(bdfiles = bdfile4b)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
gen2bfile <- system.file("extdata", "set2b.imp", package = "BinaryDosage")
sample2bfile <- system.file("extdata", "set2b.sample", package = "BinaryDosage")
bdfile2b <- tempfile()
gentobd(genfiles = c(gen2bfile, sample2bfile),
        bdfiles = bdfile2b,
        snpcolumns = c(1L,3L,2L,4L,5L),
        impformat = 1L)

bdinfo2b <- getbdinfo(bdfiles = bdfile2b)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
gen3bfile <- system.file("extdata", "set3b.imp", package = "BinaryDosage")
sample3bfile <- system.file("extdata", "set3b.sample", package = "BinaryDosage")
bdfile3bm1 <- tempfile()
gentobd(genfiles = c(gen3bfile, sample3bfile),
        bdfiles = bdfile3bm1,
        snpcolumns = c(-1L,2L:5L),
        chromosome = "1")

bdinfo3bm1 <- getbdinfo(bdfiles = bdfile3bm1)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
gen3bfile <- system.file("extdata", "set3b.imp", package = "BinaryDosage")
sample3bnhfile <- system.file("extdata", "set3bnh.sample", package = "BinaryDosage")
bdfile3bnh <- tempfile()
gentobd(genfiles = c(gen3bfile, sample3bnhfile),
        bdfiles = bdfile3bnh,
        snpcolumns = c(0L,2L:5L),
        header = c(FALSE, FALSE))

bdinfo3bnh <- getbdinfo(bdfiles = bdfile3bnh)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
gen4bgzfile <- system.file("extdata", "set4b.imp.gz", package = "BinaryDosage")
sample4bfile <- system.file("extdata", "set4b.sample", package = "BinaryDosage")
bdfile4bgz <- tempfile()
gentobd(genfiles = c(gen4bgzfile, sample4bfile),
        bdfiles = bdfile4bgz,
        snpcolumns = c(1L,2L,4L,5L,6L),
        startcolumn = 7L,
        impformat = 2L,
        gz = TRUE)

bdinfo4bgz <- getbdinfo(bdfiles = bdfile4bgz)


## ---- eval = T, echo = T, message = F, warning = F, tidy = T------------------
gen3bchrfile <- system.file("extdata", "set3b.chr.imp", package = "BinaryDosage")
sample3bfile <- system.file("extdata", "set3b.sample", package = "BinaryDosage")
geninfo <- getgeninfo(genfiles = c(gen3bchrfile, sample3bfile), index = TRUE)

aaf <- unlist(genapply(geninfo = geninfo, getaaf))

altallelefreq <- data.frame(SNP = geninfo$snps$snpid, aafcalc = aaf)
knitr::kable(altallelefreq, caption = "Calculated aaf", digits = 3)


