\name{TestCor-package}
\alias{TestCor-package}
\alias{TestCor}
\docType{package}
\title{
  FWER and FDR controlling procedures for multiple correlation tests
}
\description{
   The package compiles some multiple testing procedures which theoretically control asymptotically the FWER in the framework of correlation testing.
    Four tests statistics can be considered: the empirical correlation, the Student statistics, the Fisher's z-transform and the usual Gaussian statistics considering random variables \eqn{(X_i-mean(X_i))(X_j-mean(X_j))}.
    Four methods are implemented: Bonferroni (1935)'s, Šidák (1967)'s, Romano & Wolf (2005)'s bootstrap and (Drton & Perlman (2007)'s procedure based on the asymptotic distributions of the test statistics, called MaxTinfty.
    The package also includes some multiple testing procedures which are related to the control of the FDR : Cai & Liu (2016)'s procedures called LCT-N and LCT-B -which have been proven to control the FDR for correlation tests- and Benjamini & Hochberg (1995)'s -which has no theoretical results in correlation testing.
}
\details{
  Consider \eqn{\lbrace \bold{X}_\ell = (X_{1\ell},\dots X_{p\ell}),\; \ell=1,...,n\rbrace} a set of \eqn{n} independent and identically distributed \eqn{R^p}-valued random variables. Denote \code{data} the array containing \eqn{\lbrace\mathbf{X_\ell},\; \ell=1,\dots,n\rbrace}, with observation indexes \eqn{l} in row. The aim is to test simultaneously \deqn{(H_{0ij})~ Cor(X_i,X_j)=0 {~~against~~} (H_{1ij})~ Cor(X_i,X_j)\neq 0,~~ i,j=1,..., p,~ i<j.}
  Four tests statistics are implemented: the empirical correlation, the Student statistics, the Fisher's z-transform and the usual test statistics on expectancy considering the product of random variables. They are available in function \code{eval_stat}.
  Next, two main types of procedures are available:
\describe{
  \item{Asymptotically FWER controlling procedures:}{Bonferroni (1935)'s method, Šidák (1967)'s procedure, Romano & Wolf (2005)'s bootstrap procedure and Drton & Perlman (2007)'s procedure. A description of these methods can be found in Chapter 5 of Roux (2018). To apply these procedures, function \code{ApplyFwerCor} can be used as follows:

\code{ApplyFwerCor(data,alpha,stat_test,method)}, with \code{alpha} the desired level of control for FDR and \code{stat_test}, \code{method} respectively the kind of test statistic and the FDR controlling method. The function returns the list of indexes \eqn{\lbrace (i,j),  i < j \rbrace} for which null hypothesis \eqn{(H_{0ij})} is rejected.}
   \item{Asymptotically FDR controlling procedures:}{Cai & Liu(2016)'s two procedures and Benjamini & Hochberg (1995)'s procedure (with no theoretical proof for the latest). To apply these procedures, use function \code{ApplyFdrCor} as follows:
\code{ApplyFdrCor(data,alpha,stat_test,method)} with \code{alpha} the desired level of control for FWER and \code{stat_test}, \code{method} respectively the kind of test statistic and the FDR controlling method. The function returns the list of indexes \eqn{\lbrace (i,j),  i < j \rbrace} for which null hypothesis \eqn{(H_{0ij})} is rejected. }
 }
Functions \code{SimuFwer} and \code{SimuFdr} provide simulations of Gaussian random variables for a given correlation matrix and return estimated FWER, FDR, Power and true discovery rate obtained applying one of the procedure above. Some example of results obtained can be found in Chapter 6 of Roux (2018). 
}
\author{Irene Gannaz

Maintainer: Irene Gannaz <irene.gannaz@insa-lyon.fr>
}
\references{
Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300, \url{https://doi.org/10.1111/j.2517-6161.1995.tb02031.x}.

Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in onore del professore salvatore ortu carboni, 13-60. 

Cai, T. T., & Liu, W. (2016). Large-scale multiple testing of correlations. Journal of the American Statistical Association, 111(513), 229-240, \url{https://doi.org/10.1080/01621459.2014.999157}.

Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449, \url{https://doi.org/10.1214/088342307000000113}.

Romano, J. P., & Wolf, M. (2005). Exact and approximate stepdown methods for multiple hypothesis testing. Journal of the American Statistical Association, 100(469), 94-108, \url{https://doi.org/10.1198/016214504000000539}.

Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, \url{https://tel.archives-ouvertes.fr/tel-01971574v1}.

Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. Journal of the American Statistical Association, 62(318), 626-633.
}
\concept{multiple testing} \concept{FWER} \concept{FDR} \concept{correlation testing} \concept{graph inference}

\examples{
\donttest{
# Parameters for simulations
Nsimu  <- 100                # number of Monte-Carlo simulations
seqn   <- seq(100,400,100)   # sample sizes
p      <- 10                 # number of random variables considered
rho    <- 0.3                # value of non-zero correlations
seed   <- 156724
 
corr_theo <- diag(1,p)       # the correlation matrix
corr_theo[1,2:p] <- rho
corr_theo[2:p,1] <- rho               

# Parameters for multiple testing procedure
stat_test <- 'empirical'     # test statistics for correlation tests
method <- 'BootRW'           # FWER controlling procedure
SD <- FALSE                  # logical determining if stepdown is applied
alpha  <- 0.05               # FWER threshold 
Nboot  <- 100                # number of bootstrap or simulated samples

# Simulations and application of the chosen procedure
res <- matrix(0,nrow=length(seqn),ncol=5)
for(i in 1:length(seqn)){
    temp <- SimuFwer(corr_theo,n=seqn[i],Nsimu=Nsimu,alpha=alpha,stat_test=stat_test,
           method='BootRW',Nboot=Nboot,stepdown=SD,seed=seed)
    res[i,] <- temp
}
rownames(res) <- seqn
colnames(res) <- names(temp)

# Display results
par(mfrow=c(1,2))
plot(seqn,res[,'fwer'],type='b',ylim=c(0,max(alpha*1.1,max(res[,'fwer']))),
    main='FWER',ylab='fwer',xlab='number of observations')
plot(seqn,res[,'sensitivity'],type='b',ylim=c(0,1.1),
    main='Power',ylab='sensitivity',xlab='number of observations')
}
}

