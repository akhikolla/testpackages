\name{rsnpset.pvalue}
\alias{rsnpset.pvalue}
\title{RSNPset P-value Function}
\description{
    Calculate observed, resampling, FWER-adjusted, and FDR-adjusted p-values for statistics from the function \code{rsnpset()}.
}
\usage{
    rsnpset.pvalue(result, pval.transform=FALSE, qfun=function(x)\{qvalue(x)$qvalue\}) 
}

\arguments{
	\item{result}{Result from \code{rsnpset()}, an "RSNPset" S3 class object. Required.}
	\item{pval.transform}{Boolean indicating if the resampling p-values should be computed by comparing the observed p-value to the resampling p-values (\code{TRUE}). If not (\code{FALSE}), they are computed by comparing the observed statistics to the resampling statistics (may not be appropriate for \code{r.method="permutation"}). Note that \code{rsnpset()} must be run with \code{B > 0} in order to use \code{pval.transform=TRUE}. Default is \code{FALSE}.}
	\item{qfun}{Function used to calculate false discovery rate adjusted p-values. See below. Default is \code{function(x)\{qvalue(x)$qvalue\}}.}
}

\details{
See below.
}

\value{
An S3 class \code{RSNPset.pvalue} object that extends \code{data.frame}, with one row for each of the \code{K} SNP sets in \code{result}, columns \code{W, rank, m,} and two or more additional columns of p-values. Two columns, \code{p,} and \code{q} are always returned. If \code{rsnpset()} was run with \code{B > 0}, the columns \code{pB} and \code{qB} are included as well. If \code{pval.transform=TRUE}, the returned p-value columns will be \code{p, pB, PB, q,} and \code{qB}. 
    
\tabular{ll}{
\bold{Column} \tab \bold{Definition} \cr
\code{W} \tab Observed statistic \cr
\code{rank} \tab Rank of the variance matrix for the observed data \cr
\code{m} \tab Number of SNPs analyzed in the SNP set \cr
}

\tabular{lll}{
\bold{Column} \tab \bold{P-value} \tab \bold{Definition} \cr
\code{p} \tab Asymptotic* \tab \code{pchisq(W,rank,lower.tail=FALSE)} \cr
\code{pB} \tab Resampling** \tab See below. \cr
\code{PB} \tab Family-wise error adjusted*** \tab See below. \cr
\code{q} \tab False discovery rate adjusted \tab \code{qfun(p)} \cr
\code{qB} \tab Resampling FDR adjusted \tab \code{qfun(pB)} \cr 
}

\bold{*} For \code{W} and \code{rank} from \code{rsnpset()}.

\bold{**} By default, the unadjusted resampling p-values are computed by comparing the observed statistics to the replication statistics. Note that a large number of replications may be required in order to account for multiple testing. For each SNP set, the value for \code{pB} is \code{sum(W <= Wb)/B}, where \code{W} is the observed statistic for the SNP set, \code{Wb} is a vector of resampling statistics, and \code{B} is the number of replications. If \code{pval.transform=TRUE}, then for each SNP set, the value for \code{pB} is \code{sum(p > pb)/B} where \code{p} is the observed p-value, and \code{pb} is a vector of the p-values of the \code{B} resampling statistics. It is possible that \code{pB} may be 0 for some SNP sets. To prevent this, \code{pmax(pB,1/B)} is returned instead.

\bold{***} The column \code{PB} is only returned if \code{pval.transform=TRUE}. For each SNP set, the value for \code{PB} is \code{sum(p > Zb)/B}, where \code{Zb} a vector of length \code{B}. Each element of \code{Zb} is the smallest resampling p-value across all \code{K} SNP sets for the \code{b}th replication. It is possible that \code{PB} may be 0 for some SNP sets. To prevent this, \code{pmax(PB,1/B)} is returned instead.

}

\note{
The \code{qvalue()} function, used by default in \code{qfun}, can fail for small numbers of replications/SNP sets. To overcome this, the \code{qfun} argument can be used to define a new q-value function, or to assign arguments for the \code{qvalue()} function. For example: 

\code{qfun=function(x)\{qvalue(x, robust=TRUE)$qvalue\}}.
}

\seealso{
This function computes p-values for the statistics from the function \code{\link{rsnpset}}. 

For sorting and reviewing the p-values, see \code{\link{summary.RSNPset.pvalue}}.

More information on \code{\link{qvalue}}.
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
rsnpset.pvalue(res, pval.transform=TRUE)
}

