## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
library(knitr)
# Set so that long lines in R will be wrapped:
opts_chunk$set(tidy.opts=list(width.cutoff=80),tidy=TRUE)

## ----mega2rtutDump, message=FALSE---------------------------------------------
library(Mega2R)
dump_mega2rtutorial_data()

## ----mega2rtutDumpls, message=FALSE-------------------------------------------
list.files(where_mega2rtutorial_data())

## ---- mega2rtutClean,eval=FALSE-----------------------------------------------
#  clean_mega2rtutorial_data()

## ---- install,eval=FALSE------------------------------------------------------
#  install.packages("Mega2R")

## ---- seqsimrdb,eval=FALSE----------------------------------------------------
#  library(Mega2R)
#  # Before issuing the next command, make sure you have issued
#  # this command `dump_mega2rtutorial_data()` first
#  # as instructed in the "Tutorial Data' section above.
#  db = file.path(where_mega2rtutorial_data(), "seqsimr.db")
#  ENV = read.Mega2DB(db, verbose = TRUE)

## ----seqsimrdb2---------------------------------------------------------------
# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
db = file.path(where_mega2rtutorial_data(), "seqsimr.db")
ENV = read.Mega2DB(db, verbose = TRUE)

## ---- ls----------------------------------------------------------------------
ls(ENV)

## ---- showMega2ENV------------------------------------------------------------
showMega2ENV() 

## ---- str---------------------------------------------------------------------
str(ENV$locus_table)

## ----setrange0----------------------------------------------------------------
    dim(ENV$refRanges)
    head(ENV$refRanges)

## ----setrange1----------------------------------------------------------------
    ranges = matrix(c(1, 2240000, 2245000,
                      1, 2245000, 2250000,
                      1, 3760000, 3761000,
                      1, 3761000, 3762000,
                      1, 3762000, 3763000,
                      1, 3763000, 3764000,
                      1, 3764000, 3765000,
                      1, 3765000, 3763760,
                      1, 3763760, 3767000,
                      1, 3767000, 3768000,
                      1, 3768000, 3769000,
                      1, 3769000, 3770000),
                      ncol = 3, nrow = 12, byrow = TRUE)

    setRanges(ranges, 1:3)

    dim(ENV$refRanges)
    head(ENV$refRanges)

## ----setrange2----------------------------------------------------------------
    ranges = matrix(c(1, 2240000, 2245000,
                      1, 2245000, 2250000,
                      1, 3760000, 3761000,
                      1, 3761000, 3762000,
                      1, 3762000, 3763000,
                      1, 3763000, 3764000,
                      1, 3764000, 3765000,
                      1, 3765000, 3763760,
                      1, 3763760, 3767000,
                      1, 3767000, 3768000,
                      1, 3768000, 3769000,
                      1, 3769000, 3770000),
                      ncol = 3, nrow = 12, byrow = TRUE)
    ranges = data.frame(ranges)
    ranges$name = LETTERS[1:12]
    names(ranges) = c("chr", "start", "end", "name")

    setRanges(ranges, 1:4)
    dim(ENV$refRanges)
    head(ENV$refRanges)

## ----callback-----------------------------------------------------------------
    show = function(m, r, e) {
       print("rrrrrrrrrr")
       print(r)
       print("mmmmmmmmmm")
       print(m)
       print("g6g6g6g6g6")
       print(head(getgenotypes(m, envir = e)))
    }

## ----applyfntoranges1,eval=TRUE-----------------------------------------------
    # apply function "show" to all the ranges
    # ranges
    ENV$verbose = FALSE
    applyFnToRanges(show)

## ----applyfntoranges2,eval=TRUE-----------------------------------------------
    # apply function "show" to all genotypes on chromosomes 1 
    # ranges
    applyFnToRanges(show, 
                    ranges_arg =
                    matrix(c(1, 4000000, 5000000, "range4m",
                             1, 5000000, 6000000, "range5m",
                             1, 6000000, 7000000, "range6m",
                             1, 7000000, 8000000, "range7m",
                             1, 8000000, 9000000, "range8m",
                             1, 9000000,10000000, "range9m"),
                            ncol = 4, nrow = 6, byrow = TRUE),
                    indices_arg = 1:4)

## ----setannotations0----------------------------------------------------------
    ENV$txdb
    ENV$entrezGene

## ----applyfntogenes0,eval=TRUE------------------------------------------------
    # apply function "show" to all transcripts on genes ELL2 and CARD15
    applyFnToGenes(show, genes_arg = c("CEP104"))

## ----applyfntogenes19,eval=TRUE-----------------------------------------------
    setAnnotations("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db")
    applyFnToGenes(show, genes_arg = c("CEP104"))

## ----applyfntogenes,eval=TRUE-------------------------------------------------
    # apply function "show" to all genotypes on chromosomes 1 for two base
    # pair ranges
    applyFnToGenes(show, ranges_arg = 
                            matrix(c(1,  5000000, 10000000,
                                     1, 10000000, 15000000),
                                    ncol = 3, nrow = 2, byrow = TRUE))

    # apply function "show" to all genotypes for first marker 
    # in each chromosome (We only have data for chromosome 1.)
    # NOTE: Since we are using an arbitrary collection of markers, the
    # range is not available.
    applyFnToGenes(show, markers_arg = ENV$markers[! duplicated(ENV$markers$chromosome), 3])

    # apply function "show" to all genotypes on chromosomes 24 and 26.
    # remember our example database is only chr 1
    applyFnToGenes(show, chrs_arg=c(24, 26))

## ----showstat,eval=TRUE-------------------------------------------------------
show.stat = function(m, r, e, fn) {
      print(r) 
# collect genotypes for the set of markers "m"
      mm = getgenotypes(m, envir=e) 
# apply xxx.test of trait vs marker (accumulating samples)
      pv = apply(mm, 2, fn)
      names(pv) = m$MarkerName
      print(pv)
}

show2 = function(m, r, e) {
    f = function(x) {
          tryCatch(fisher.test(table(e$fam$trait, x)),
                   error = function(e) {list(p.value = NA)}) $p.value }

    show.stat(m, r, e, f)
}

show3 = function(m, r, e) {
    f = function(x) {
          tryCatch(chisq.test(table(e$fam$trait, x)),
                   error = function(e) {list(p.value = NA)}) $p.value }

    show.stat(m, r, e, f)
}

## ----xshow2,eval=TRUE---------------------------------------------------------
    applyFnToGenes(show2, genes_arg = c("CEP104"))

## ----xshow3,eval=TRUE---------------------------------------------------------
    applyFnToGenes(show3, genes_arg = c("CEP104"))

## ----getgenotypes,eval=TRUE---------------------------------------------------
    genotype = getgenotypes(ENV$markers[1:10,])
    dim(genotype)
    head(genotype)

## ----getgenotypesraw,eval=TRUE------------------------------------------------
# two ints in upper/lower half integer representing allele 
    raw = getgenotypesraw(ENV$markers[1:10,])
    dim(raw)
    head(raw)

## ----init_pedgene,message=FALSE-----------------------------------------------
# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
db = file.path(where_mega2rtutorial_data(), "seqsimr.db")
ENV = init_pedgene(db)

## ----mega2pedgene-------------------------------------------------------------
ENV$verbose=TRUE
Mega2pedgene(gs=50:60)

## ---- FnToRanges,eval=FALSE---------------------------------------------------
#  # we will skip this line for the Rmd document production because it takes too long
#  applyFnToRanges(DOpedgene, ENV$refRanges, ENV$refIndices, envir = ENV)

## ----eval=FALSE---------------------------------------------------------------
#  applyFnToGenes(DOpedgene, genes_arg = c('DISP1', 'KIF26B', 'AK5', 'ST7L'), envir = ENV)

## ----FnToGenes,eval=TRUE------------------------------------------------------
applyFnToGenes(DOpedgene, genes_arg = c('DISP1', 'AK5'), envir = ENV)

## ---- FnToGenesall,eval=FALSE-------------------------------------------------
#  # we will skip this line for the Rmd document production because it takes too long
#  applyFnToGenes(DOpedgene, genes_arg = '*', envir = ENV)

## ----init_SKAT,message=FALSE--------------------------------------------------
# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
db = file.path(where_mega2rtutorial_data(), "seqsimr.db")
ENV = init_SKAT(db, verbose = F, allMarkers = F)

## ----init_SKAT2,eval=FALSE----------------------------------------------------
#  # Before issuing the next command, make sure you have issued
#  # this command `dump_mega2rtutorial_data()` first
#  # as instructed in the "Tutorial Data' section above.
#  db = file.path(where_mega2rtutorial_data(), "seqsimr.db")
#  ENV = init_SKAT(db, verbose = T, allMarkers = F)

## ---- FnToRangesSKAT,eval=TRUE------------------------------------------------
ENV$verbose = FALSE
ENV$SKAT_results = ENV$SKAT_results[0, ]
Mega2SKAT(ENV$phe[, 3] - 1 ~ 1, "D", kernel = "linear.weighted", weights.beta=c(0.5,0.5))

## ---- FnToRangesa,eval=TRUE---------------------------------------------------
print(ENV$SKAT_results)

## ---- FnToRange2s,eval=FALSE--------------------------------------------------
#  # we will skip this line for the Vignette document production because it takes too long
#  ENV$verbose = FALSE
#  ENV$SKAT_results = ENV$SKAT_results[0, ]
#  Mega2SKAT(ENV$phe[, 3] - 1 ~ 1, "D", kernel = "linear.weighted", weights.beta=c(0.5,0.5), gs=1:nrow(ENV$refRanges))
#  print(ENV$SKAT_results)

## ----FnToGenesSKAT,eval=TRUE--------------------------------------------------
ENV$SKAT_results = ENV$SKAT_results[0, ]
Mega2SKAT(ENV$phe[, 3] - 1 ~ 1, "D", kernel = "linear.weighted", weights.beta=c(0.5,0.5), genes = c('DISP1', 'AK5', 'KIF26B', 'ST7L'), envir = ENV)

## ---- FnToGenesa,eval=TRUE----------------------------------------------------
print(ENV$SKAT_results)

## ----FnToGene2s,eval=TRUE-----------------------------------------------------
ENV$verbose = FALSE
ENV$SKAT_results = ENV$SKAT_results[0, ]
Mega2SKAT(ENV$phe[, 3] - 1 ~ 1, "D", kernel = "linear.weighted", weights.beta=c(0.5,0.5), genes = '*', envir = ENV)

## ---- FnToGene2a,eval=TRUE----------------------------------------------------
print(ENV$SKAT_results)

## ----init_famSKAT,message=FALSE-----------------------------------------------
# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
db = file.path(where_mega2rtutorial_data(), "seqsimr.db")
ENV = init_famSKATRC(db, verbose = F)

## ----init_famSKATb------------------------------------------------------------
# The special hack below reduces the samples to 20% of the original, so the run will
# finish in reasonable time. There are 20 different pedigrees.
setfam(ENV$fam[(ENV$fam$pedigree_link %in% 0:3),], ENV)

ENV$phe = mkphenotype(ENV)
ENV$phe[ENV$phe == 0] = NA

## ----init_famSKAT2,eval=FALSE-------------------------------------------------
#  # Before issuing the next command, make sure you have issued
#  # this command `dump_mega2rtutorial_data()` first
#  # as instructed in the "Tutorial Data' section above.
#  db = file.path(where_mega2rtutorial_data(), "seqsimr.db")
#  ENV = init_famSKATRC(db, verbose = T)

## ---- FnToRangesfamSKAT,eval=TRUE---------------------------------------------
ENV$verbose = FALSE
Mega2famSKATRC(pheno=3,gs=1:60)

## ---- FnToRangesfama,eval=TRUE------------------------------------------------
print(ENV$famSKATRC_results)

## ---- FnToRange2fams,eval=FALSE-----------------------------------------------
#  # we will skip this line for the Vignette document production because it takes too long
#  ENV$verbose = FALSE
#  Mega2famSKATRc(pheno=3, gs=1:nrow(ENV$refRanges))
#  print(ENV$famSKATRC_results)

## ----FnToGenesfamSKAT,eval=TRUE-----------------------------------------------
ENV$verbose = TRUE
Mega2famSKATRC(pheno=3, genes = c('DISP1', 'AK5', 'KIF26B', 'ST7L'), envir = ENV)
print(ENV$famSKATRC_results)

## ----FnToGene2fams,eval=FALSE-------------------------------------------------
#  ENV$verbose = FALSE
#  Mega2famSKATRC(pheno=3, genes = '*', envir = ENV)
#  print(ENV$famSKATRC_results)

## ----mega2vcflib--------------------------------------------------------------
# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
vcfdir = file.path(where_mega2rtutorial_data(), "vcfr")
if (! dir.exists(vcfdir)) dir.create(vcfdir)

## ----mega2vcf-----------------------------------------------------------------
# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
vcffile = file.path(where_mega2rtutorial_data(), "vcfr", "vcf.01")
Mega2VCF(vcffile, ENV$markers[ENV$markers$chromosome==1 ,])

## ---- lsvcf-------------------------------------------------------------------
# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
vcfdir = file.path(where_mega2rtutorial_data(), "vcfr")
list.files(vcfdir)

## ----seqsimgwaa, eval=FALSE---------------------------------------------------
#  require("GenABEL")
#  # Before issuing the next command, make sure you have issued
#  # this command `dump_mega2rtutorial_data()` first
#  # as instructed in the "Tutorial Data' section above.
#  db = file.path(where_mega2rtutorial_data(), "seqsimr.db")
#  ENV = read.Mega2DB(db)
#  
#  # This line converts the database to a gwaa.data-class object. The intermediate
#  # files are in tempdir() and begin with "Mega2GenABEL"
#  seqsimgwaa = Mega2GenABEL()

## ----strseqsimgwaa, eval=FALSE------------------------------------------------
#  str(seqsimgwaa)

## ----mega2gds-----------------------------------------------------------------
# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
db = file.path(where_mega2rtutorial_data(), "seqsimr.db")
ENV = read.Mega2DB(db)

## ----mega2gdsfx1--------------------------------------------------------------
# NOTE: the gds file to be created must be closed with the function below,
# or by using an on.exit(closefn.gds(<name>))
showfile.gds(closeall=T, verbose=F)

# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
gdsfile = file.path(where_mega2rtutorial_data(), "foo.gds")
gdsn = Mega2gdsfmt(gdsfile, ENV$markers[ENV$markers$chromosome==1 ,], SeqArray=TRUE)

## ---- mega2gdspr1-------------------------------------------------------------
print(gdsn)
closefn.gds(gdsn)

## ----mega2gdsfx2,message=FALSE,warning=FALSE----------------------------------
# NOTE: the gds file to be created must be closed with the function below,
# or by using an on.exit(closefn.gds(<name>))
showfile.gds(closeall=T, verbose=F)

## ----mega2gdsfx3,message=FALSE,warning=FALSE----------------------------------
# Before issuing the next command, make sure you have issued
# this command `dump_mega2rtutorial_data()` first
# as instructed in the "Tutorial Data' section above.
gdsfile = file.path(where_mega2rtutorial_data(), "foo.gds")
gdsn = Mega2gdsfmt(gdsfile, ENV$markers[ENV$markers$chromosome==1 ,], SeqArray=FALSE)

## ---- mega2gdspr2-------------------------------------------------------------
print(gdsn)
closefn.gds(gdsn)

## ----srdta, eval=FALSE--------------------------------------------------------
#  GotGenABEL = require("GenABEL", quietly=FALSE)
#  if (GotGenABEL) data(srdta) else srdta=NULL

## ----dmpPed,results="hide",warning=FALSE,eval=FALSE---------------------------
#  GotGenABEL = require("GenABEL", quietly=FALSE)
#  
#  # Before issuing the next command, make sure you have issued
#  # this command `dump_mega2rtutorial_data()` first
#  # as instructed in the "Tutorial Data' section above.
#  
#  srdtafile = file.path(where_mega2rtutorial_data(), "srdta")
#  if (GotGenABEL) export.plink(srdta, transpose = FALSE, filebasename = srdtafile,
#               phenotypes = names(srdta@phdata)[-(1:2)])

## ---- srdta.db,eval=FALSE-----------------------------------------------------
#  GotGenABEL = require("GenABEL", quietly=FALSE)
#  # Before issuing the next command, make sure you have issued
#  # this command `dump_mega2rtutorial_data()` first
#  # as instructed in the "Tutorial Data' section above.
#  sdb = file.path(where_mega2rtutorial_data(), "srdta.db")
#  ENV = read.Mega2DB(sdb)
#  
#  mega = Mega2GenABEL()

## ----strs, eval=FALSE---------------------------------------------------------
#  GotGenABEL = require("GenABEL", quietly=FALSE)
#  str(mega)
#  if (GotGenABEL) str(srdta)

## ----mega2genabeltst, eval=FALSE----------------------------------------------
#  GotGenABEL = require("GenABEL", quietly=FALSE)
#  options(max.print = 30)
#  Mega2GenABELtst(mega_ = mega, gwaa_ = srdta)

## ---- fin,eval=TRUE,echo=FALSE------------------------------------------------
clean_mega2rtutorial_data()

