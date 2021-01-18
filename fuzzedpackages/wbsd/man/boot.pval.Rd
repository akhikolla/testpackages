\name{boot.pval}

\alias{boot.pval}

\title{Bootstrap p-value}
\encoding{latin1}
\description{
This function computes a bootstrap p-value as discussed in the article \enc{Pötscher}{Poetscher} and Preinerstorfer (2020).
}

\usage{
boot.pval(y, X, R, r, hcmethod, restr.cov, wilddist, wildmult, 
wildmult.restr, boot.res.restr, boot.center.restr,    
tol = 1e-07, comp.meth = "exact", Boot.supp = NULL, 
checks = TRUE, cores = 1)
}

\arguments{
\item{y}{Matrix (n rows)}
\item{X}{Matrix (n times k, rank k).}
\item{R}{Matrix (q times k, rank q).}
\item{r}{Vector (q dimensional).}
\item{hcmethod}{Integer in [-1, 4]. Determines the method used in the construction of the covariance estimator
used in the test statistic. The value -1 corresponds to uncorrected F statistic without df adjustment; the value 0 
corresponds to the HC0 estimator; ...; the value 4 corresponds to the HC4 estimator. Cf. \enc{Pötscher}{Poetscher} 
and Preinerstorfer (2020) and the references there for details.}
\item{restr.cov}{TRUE or FALSE. Covariance matrix estimator based on null-restricted (TRUE) or unrestricted (FALSE) residuals.}
\item{wilddist}{Auxiliary distribution used to generate the bootstrap sample. Either ``rademacher'' or ``mammen''. This is only used in case comp.meth = exact.}
\item{wildmult}{Integer in [0, 4]. Determines the wild-bootstrap multiplier weights used. The value 0 corresponds to constant multipliers of 1;
the value 1 corresponds to multipliers that are obtained by taking the square of the HC1 weights; ...; the value
4 corresonds to multipliers that are obtained by taking the square of the HC4 weights.}
\item{wildmult.restr}{TRUE or FALSE. Bootstrap multiplier weights computed on null-restricted (TRUE) or unrestricted (FALSE) projection matrix.}
\item{boot.res.restr}{TRUE or FALSE. Bootstrap residuals based on null-restricted sample residuals (TRUE) or unrestricted sample residuals (FALSE)}
\item{boot.center.restr}{TRUE or FALSE. Bootstrap sample centered at null-restricted predictor (TRUE) or at the unrestricted predictor (FALSE).}
\item{tol}{Real number. Tolerance parameter used in checking invertibility of the covariance matrix in the test statistic. Default is 1e-07. A negativ input value will be converted to the machine epsilon. To check the stability of the numerical computations, it
is recommended to try different values of tol and compare the obtained value of theta. Larger values of tol will lead to larger theta.}
\item{comp.meth}{``exact'' or ``approximation''. Determines whether probabilities over the bootstrap distribution are computed exactly or approximately (the latter requires input Boot.supp). Exact computations are only feasible for small sample sizes.}
\item{Boot.supp}{Bootstrap approximation distribution (e.g., obtained by drawing vectors each with n independent entries distributed according to a Rademacher or Mammen distribution); (matrix, n rows, columns = number of bootstrap samples). This is only used if comp.meth = ``approximation''; note that these vectors are further multiplied by wildmult weights.}
\item{checks}{TRUE (default) or FALSE. Determines whether input checks are conducted (TRUE) or not (FALSE).}
\item{cores}{Positive integer. Maximal number of cores used in the computations.}
}

\details{For details see the relevant sections in \enc{Pötscher}{Poetscher} 
and Preinerstorfer (2020).
}

\value{
\item{p}{Bootstrap p-value.}
}

\examples{
y <- matrix(rnorm(10), ncol = 1)
X <- cbind(rep(1, 10), rnorm(10))                         
R <- matrix(c(0, 1), nrow = 1, ncol = 2)
r <- 0

boot.pval(y, X, R, r, 3, FALSE, "rademacher", 3, FALSE, TRUE, TRUE)
}


\references{
\enc{Pötscher}{Poetscher}, B. M. and Preinerstorfer, D. (2020). How Reliable are Bootstrap-based Heteroskedasticity Robust Tests? <arXiv:2005.04089>
}             