\name{collUtils-package}
\alias{collUtils-package}
\alias{collUtils}
\docType{package}
\title{
A auxiliary package for CollapsABEL.
}
\description{
This package includes some low level functions for processing PLINK input and output files written in Java or C++.
Normally you shouldn't need to directly use functions from this package. 
}
\details{
\tabular{ll}{
Package: \tab collUtils\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-06-12\cr
License: \tab GPL-3\cr
}
}
\author{
Kaiyin Zhong
Maintainer: Kaiyin Zhong <kindlychung@gmail.com>
}
\references{
To be updated.
}
\keyword{ 
GWAS
Compound heterozygosity
Genomics
}
\examples{
\dontrun{
require(collUtils)
rbed_obj = rBed("test.bed")
geno = rbed_obj$readBed()
geno = getJArray(geno)
print(geno)
fn = tempfile()
f = file(fn, "wb")
writeBin("a", f)
writeBin("b", f)
writeBin("c", f)
close(f)
file.info(fn)$size == 6
truncateEndOfFile(fn, 1)
}
}