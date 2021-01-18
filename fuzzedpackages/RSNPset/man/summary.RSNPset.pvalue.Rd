\name{summary.RSNPset.pvalue}
\alias{summary.RSNPset.pvalue}
\title{RSNPset P-value Summary Function}
\description{
    Summary function to sort and display p-values resulting from \code{rsnpset.pvalue()}.
}
\usage{
\method{summary}{RSNPset.pvalue}(object, sort="p", decreasing=FALSE, 
        nrows=10, dropcols=c(""), verbose=FALSE, ...) 
}

\arguments{
	\item{object}{Result from \code{rsnpset.pvalue()}, an "RSNPset.pvalue" S3 class object. Required.}
	\item{sort}{Character string indicating column by which to sort results. If not one of \code{c("W",} \code{"rank",} \code{"m",} \code{"p",} \code{"pB",} \code{"PB",} \code{"q",} \code{"qB")} results will be sorted by the row names (i.e. the names of the SNP sets) instead. Default is \code{"p"}.}
	\item{decreasing}{Boolean indicating if the sort column should be arranged in decreasing order. Default is \code{FALSE}.}
	\item{nrows}{Integer indicating number of rows to display. Default is 10.}
	\item{dropcols}{Character vector corresponding names of columns of to be suppressed from the summary. Default is none.}
	\item{verbose}{Boolean indicating if additional information about the p-value calculations should be reported. Default is \code{FALSE}.}
	\item{...}{Additional arguments affecting the summary produced.}
}

\details{
As a typical GWAS study may span thousands of SNPs and SNP sets, this function allows for the succinct reporting of p-values for the most significant results. For more information about the different columns reported, see the documentation for \code{rsnpset.pvalue()}. If \code{verbose=TRUE}, a note will be printed with the total number of SNP sets and replications used in the calculations, as well as the value of the \code{pval.transform} argument from \code{rsnpset.pvalue()}.
}

\value{
    A \code{data.frame} object subset from \code{object}, the result of \code{rsnpset.pvalue()}. Rows are selected based on the \code{sort}, \code{decreasing}, and \code{nrows} arguments, and columns are selected based on the \code{dropcols} argument.
}

\seealso{
    The function \code{\link{rsnpset.pvalue}} provides a description of the different p-values computed, as well as the other columns in the results.
}

\examples{
n <- 200    # Number of patients
m <- 1000   # Number of SNPs

set.seed(123)
G <- matrix(rnorm(n*m), n, m)   # Normalized SNP expression levels
rsids <- paste0("rs", 1:m)      # SNP rsIDs 
colnames(G) <- rsids
 
K <- 15                         # Number of SNP sets
genes <- paste0("XYZ", 1:K)     # Gene names 
gsets <- lapply(sample(3:50, size=K, replace=TRUE), sample, x=rsids)
names(gsets) <- genes

# Survival outcome
time <- rexp(n, 1/10)           # Survival time
event <- rbinom(n, 1, 0.9)      # Event indicator

\dontrun{
# Optional parallel backend
library(doParallel)
registerDoParallel(cores=8) }

# B >= 1000 is typically recommended
res <- rsnpset(Y=time, delta=event, G=G, snp.sets=gsets, score="cox", 
               
               B=50, r.method="permutation", ret.rank=TRUE)
pvals <- rsnpset.pvalue(res, pval.transform=TRUE)

summary(pvals)

summary(pvals, sort="W", decreasing=TRUE, nrows=5, dropcols=c("p","rank"), verbose=TRUE)
}
