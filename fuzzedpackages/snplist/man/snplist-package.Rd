\name{snplist-package}
\alias{snplist-package}
\alias{snplist}
\docType{package}
\title{
Tools to Create Gene Sets
}
\description{
A set of functions to create SQL tables of gene and SNP information and compose them into a SNP Set, for example for use with the \code{RSNPset} package, or to export to a PLINK set.
}
\details{
\tabular{ll}{
Package: \tab snplist\cr
Type: \tab Package\cr
Version: \tab 0.18.1\cr
Date: \tab 2017-12-11\cr
License: \tab GPL-3\cr
}
Please see the example function calls below, or refer to the individual function documentation or the included vignette for more information.
}
\author{
Authors: Chanhee Yi, Alexander Sibley, and Kouros Owzar
Maintainer: Alexander Sibley <alexander.sibley@dm.duke.edu>
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%%~~ the R documentation directory ~~
\keyword{ package }
\seealso{
    \code{RSQLite}, \code{\link{Rcpp}}
}
\examples{
chromosome <- c(1,5,22,"X","Y","MT")

geneNum <- 5
snpNum <- 1200
annoDataNum <- 500

chrLength <- 1000
geneLength <- 100

gene <- paste("gene",1:geneNum,sep="")
chr <- sample(chromosome,geneNum,replace=TRUE)
start <- sample(chrLength,geneNum,replace=TRUE)
d <- sample(geneLength,geneNum,replace=TRUE)
end <- start+d
geneInfo <- data.frame(gene,chr,start,end)

rsid <- paste("rs",1:snpNum,sep="")
chr <- sample(chromosome,snpNum,replace=TRUE)
pos <- sample(chrLength+geneLength,snpNum,replace=TRUE)
snpInfo <- data.frame(rsid,chr,pos)

annoInfo <- data.frame("rsid"=sample(rsid,annoDataNum))

dim(geneInfo)
dim(snpInfo)
dim(annoInfo)

setGeneTable(geneInfo)
setSNPTable(snpInfo)
geneset <- makeGeneSet(annoInfo)
exportPLINKSet(geneset,"geneSet.set")
file.show("geneSet.set")
}
