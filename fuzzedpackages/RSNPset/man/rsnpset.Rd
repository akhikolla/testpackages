\name{rsnpset}
\alias{rsnpset}
\title{RSNPset Analysis Function}
\description{
    Compute the (binomial, Gaussian, or Cox) efficient score for genome-wide SNP set analysis with binary, uncensored quantitative, or right-censored time-to-event traits.
}
\usage{
    rsnpset(Y, delta=NULL, G, X=NULL, snp.sets, 
            score=c("cox", "binomial", "gaussian"), B=0, 
            r.method="monte carlo", v.method="empirical",
            v.permute=FALSE, ret.rank=FALSE, pinv.check=FALSE, 
            pinv.method="specdecomp", pinv.tol=7.8e-8) 
}

\arguments{
	\item{Y}{Vector of outcomes. Can be binary, continuous, or non-negative (time-to-event), dictating the type of score to request in argument \code{score}. Required.}
	\item{delta}{Vector of event indicators for survival outcomes. Must be binary: 1 indicates event, 0 indicates censored. Required when \code{score="cox"}, unused otherwise.}
	\item{G}{Matrix of SNP allele counts or expression levels. One row per patient and one column per SNP. Must be type \code{double}. Required.}
	\item{X}{Optional matrix of covariates. One row per patient and one column per cofactor. Must be numeric with non-colinear columns. If used, an intercept column will be appended automatically. Currently implemented only for \code{score="gaussian"}.}
	\item{snp.sets}{List of sets of SNPs to be analyzed together. Each element should be a vector of values (usually rsIDs) corresponding to names of columns of \code{G}. Required.}
	\item{score}{Score statistic to be calculated. (Must be appropriate for the type of outcome given in \code{Y}). Required to be one of \code{c("cox", "binomial", "gaussian")}.}
	\item{B}{Number of replication sets used to calculate the empirical distribution of the score statistics and empirical p-values. Default is 0 (no replications will be generated).}
	\item{r.method}{Resampling method. Default is \code{"monte carlo"}. \code{r.method="permutation"} is allowed when the model does not include covariates (\code{X=NULL}).}
	\item{v.method}{Method for calculating the variance of the score statistic. Default is \code{"empirical"} for all values of \code{score}, and is the only method currently implemented for \code{score="gaussian"} or \code{score="cox"}. \code{score="binomial"} also allows \code{"asymptotic"}.}
	\item{v.permute}{Boolean indicating whether or not to recalculate the variance of the score statistics for each permutation replicate. If \code{FALSE}, variance matrices are computed only for the observed data and then reused for the permutation replicates. \code{FALSE} when \code{r.method="monte carlo"} (does not apply). \code{r.method="permutation"} also allows \code{TRUE}.}
	\item{ret.rank}{Boolean indicating whether or not to return the rank of the variance matrices for the permutation replicates. \code{FALSE} when \code{r.method="monte carlo"} (does not apply). \code{r.method="permutation"} also allows \code{TRUE}.}
	\item{pinv.check}{Boolean. If \code{TRUE}, the function returns several diagnostic measures (see below) of the Penrose-Moore inverse of the variance matrices. Default is \code{FALSE}.}
	\item{pinv.method}{Method for calculating the inverse of the variance matrices. Currently, only the default of \code{"specdecomp"} is implemented.}
	\item{pinv.tol}{Number indicating the tolerance for determining the rank of the variance matrix (see below). Default is \code{7.8e-8}.}
}

\details{
For each SNP set, the function computes the statistic \emph{W} as
\deqn{W=\bm{U}_{\bullet}^T \Sigma^{+} \bm{U}_{\bullet}}{W = t(Uvec)\%*\%\Sigma\%*\%Uvec} 
where \eqn{\Sigma^{+}}{\Sigma} is the Penrose-Moore inverse of the variance matrix \eqn{\Sigma = \bm{U}^T \bm{U}}{V = t(U)\%*\%U}, and \eqn{\bm{U}_{\bullet} = \bm{U}^T \bm{1}}{Uvec = t(colSums(U))}. Here, \eqn{\bm{U}}{U} is an \emph{n} by \emph{J} matrix where entry \emph{i,j} corresponds to the \emph{i}th patient's contribution to the score statistic for SNP \emph{j}. Statistical performance is improved by centering the values of \code{G} for each SNP prior to calculating \eqn{\bm{U}}{U}.

Under suitable regularity conditions, the distribution of \emph{W} can be approximated by a chi-squared distribution with degrees of freedom equal to the rank of \eqn{\Sigma}{V}. The rank is determined as the number of eigenvalues greater than \code{pinv.tol}. For more information on SNP sets and the efficient score, see the package vignette.

If \code{B > 0} and \code{r.method="monte carlo"}, \code{B} resampling replicates of \emph{W} are obtained by replacing \eqn{\bm{U}_{\bullet} = \bm{U}^T \bm{1}}{Uvec = t(colSums(U))} with \eqn{\bm{U}_{\bullet} = \bm{U}^T \bm{Z}}{Uvec = t(colSums(U*Z))}, where \eqn{\bm{Z}}{Z} is a vector of \emph{n} normal random values. Replications are executed in parallel, if a backend is available.

If \code{B > 0} and \code{r.method="permutation"}, \code{B} permutation replicates of \emph{W} are obtained by permuting the values of \code{Y}, or, in the case of \code{score="cox"}, by permuting \code{(Y,delta)} pairs. Permutation replicates are executed in parallel, if a backend is available. Note that \code{r.method="permutation"} is not appropriate when the model includes covariates.
    
If \code{pinv.check=TRUE}, the following diagnostic measures of the Penrose-Moore inverse are calculated.
\tabular{ll}{
\bold{Column} \tab \bold{Absolute largest element of:} \cr
\code{d0} \tab \eqn{\Sigma - \bm{QDQ}}{\Sigma - QDQ} \cr
\code{d1} \tab \eqn{\Sigma \Sigma^{+} \Sigma - \Sigma}{V\%*\%\Sigma\%*\%V-V} \cr
\code{d2} \tab \eqn{\Sigma^{+} \Sigma \Sigma^{+} - \Sigma^{+}}{\Sigma\%*\%V\%*\%\Sigma-\Sigma} \cr
\code{d2} \tab \eqn{(\Sigma \Sigma^{+})^T - \Sigma \Sigma^{+}}{t(V\%*\%\Sigma)-V\%*\%\Sigma} \cr
\code{d4} \tab \eqn{(\Sigma^{+} \Sigma)^T - \Sigma^{+} \Sigma }{t(\Sigma\%*\%V)-\Sigma\%*\%V} \cr
}
where \eqn{\bm{QDQ}}{QDQ} is the spectral decomposition of \eqn{\Sigma}{V}. Departure of these values from zero indicates poor performance of the Penrose-Moore inverse.
}

\value{
An S3 class \code{RSNPset} object consisting of a list of objects of class \code{data.frame}. The first data frame in the list describes the observed statistics, with \code{B} additional data frames corresponding to the requested resampling replicates. The rows of each data frame correspond to the elements of the \code{snp.sets} argument. The first column is \code{W}, the calculated efficient score for that set, and the second is \code{rank}, the rank of the variance matrix of the computed statistic. If \code{ret.rank=FALSE}, the ranks are not returned for the replicates. (They are assumed to be identical). The first data frame in the list also includes a third column, \code{m}, giving the number of SNPs analyzed in that SNP set.

If \code{pinv.check=TRUE}, a list of \code{B+1} data frames, each with one row per SNP set and five columns of diagnostic measures of the calculated Penrose-Moore inverse (see above), is returned as an attribute. This and other attributes of the function's execution can be accessed using the class's \code{summary} function.
}

\note{
    I. This function does not require that all entries of an element of \code{snp.sets} be present in the matrix \code{G}. If an element contains column names that are not present in \code{G}, the function will execute without objection and return a value based on the subset of columns that \emph{are} present.
    
    II. No statistics are returned for SNP sets which include SNPs with missing values (i.e. \code{NA}s in the selected columns of the matrix \code{G}).

    III. As the Cox score statistic (\code{method="cox"}) is not a sum of independent patient level scores, some level of pruning of SNPs is recommended. The permutation resampling facilities of the package can be utilized to assess the performance of the asymptotic inference. The development of robust methods for calculating the score statistics and approximating the covariance matrix is under way.
}

\seealso{
    The function \code{\link{rsnpset.pvalue}} provides options for computing p-values for the returned statistics.
    
    The function \code{\link{summary.RSNPset}} provides diagnostics and information about the function's execution.
}

\examples{
n <- 200    # Number of patients
m <- 1000   # Number of SNPs

set.seed(123)

G <- matrix(rnorm(n*m), n, m)   # Normalized SNP expression levels
rsids <- paste0("rs", 1:m)      # SNP rsIDs 
colnames(G) <- rsids
 
K <- 10                         # Number of SNP sets
genes <- paste0("XYZ", 1:K)     # Gene names 
gsets <- lapply(sample(3:50, size=K, replace=TRUE), sample, x=rsids)
names(gsets) <- genes

# Binary outcome
Yb <- rbinom(n, 1, 0.5) 
head(rsnpset(Y=Yb, G=G, snp.sets=gsets, score="binomial", v.method="empirical"))
head(rsnpset(Y=Yb, G=G, snp.sets=gsets, score="binomial", v.method="asymptotic"))

# Quantitative outcome
Yq <- rbinom(n, 1, 0.5)
head(rsnpset(Y=Yq, G=G, snp.sets=gsets, score="gaussian"))

# Survival outcome
time <- rexp(n, 1/10)           # Survival time
event <- rbinom(n, 1, 0.9)      # Event indicator
head(rsnpset(Y=time, delta=event, G=G, snp.sets=gsets, score="cox"))

\dontrun{
# Optional parallel backend
library(doParallel)
registerDoParallel(cores=8)

res <- rsnpset(Y=Yb, G=G, snp.sets=gsets, score="binomial", B=1000)
length(res) # = 1001 }
}
