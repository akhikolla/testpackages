\name{theta}

\alias{theta}

\title{theta}
\encoding{latin1}
\description{
This function computes theta as discussed in the article \enc{Pötscher}{Poetscher} and Preinerstorfer (2020).
}

\usage{
theta(
X, R, r, hcmethod, restr.cov, wilddist, wildmult, 
wildmult.restr, boot.res.restr,  boot.center.restr,    
tol = 1e-07, as.tol = 1e-07, in.tol = 1e-05, comp.meth = "exact", 
Boot.supp = NULL, checks = TRUE, cores = 1)
}

\arguments{
\item{X}{Matrix (n times k, rank k).}
\item{R}{Matrix (q times k, rank q).}
\item{r}{Vector (q dimensional).}
\item{hcmethod}{Integer in [-1, 4]. Determines the method used in the construction of the covariance estimator
used in the test statistic. The value -1 corresponds to uncorrected F statistic without df adjustment; the value 0 
corresponds to the HC0 estimator; ...; the value 4 corresponds to the HC4 estimator. Cf. \enc{Pötscher}{Poetscher} 
and Preinerstorfer (2020) and the references there for details.}
\item{restr.cov}{TRUE or FALSE. Covariance matrix estimator based on null-restricted (TRUE) or unrestricted (FALSE) residuals.}
\item{wilddist}{Distribution used to generate the bootstrap sample. Either ``rademacher'' or ``mammen''. This is only used in case comp.meth = exact.}
\item{wildmult}{Integer in [0, 4]. Determines the wild-bootstrap multiplier weights used. The value 0 corresponds to constant multiplieres of 1;
the value 1 corresponds to multipliers that are obtained by taking the square of the HC1 weights; ...; the value
4 corresonds to multipliers that are obtained by taking the square of the HC4 weights. }
\item{wildmult.restr}{TRUE or FALSE. Bootstrap multiplier weights computed on null-restricted (TRUE) or unrestricted (FALSE) projection matrix.}
\item{boot.res.restr}{TRUE or FALSE. Bootstrap residuals based on null-restricted sample residuals (TRUE) or unrestricted sample residuals (FALSE)}
\item{boot.center.restr}{TRUE or FALSE. Bootstrap sample centered at null-restricted predictor (TRUE; or FALSE)}
\item{tol}{Real number. Tolerance parameter used in checking invertibility of the covariance matrix in the test statistic. This tolerance parameter is also used in checking \eqn{R\hat{\beta}(e_i) \neq 0} in the computation of \eqn{\vartheta_{2, Het}} and \eqn{\vartheta_{2, uc}} (it is checked whether \eqn{\|R\hat{\beta}(e_i)\|_{\infty} > tol}. Default is 1e-07. A negativ input value will be converted to the machine epsilon. To check the stability of the numerical computations, it
is recommended to try different values of tol and compare the obtained value of theta. Larger values of tol will lead to larger theta.}
\item{as.tol}{Real number. Tolerance parameter used in checking Assumptions 1 or 2. Default is 1e-07. 
A negativ input value will be converted to the machine epsilon. To check the stability of the numerical computations, it
is recommended to try different values of as.tol and compare the obtained value of theta. Larger values of as.tol will lead to larger theta.}
\item{in.tol}{Real number. Tolerance parameter used in checking a strict inequality in the theta computation. Default is 1e-05. 
A negativ input value will be converted to the machine epsilon. To check the stability of the numerical computations, it
is recommended to try different values of in.tol and compare the obtained value of theta. Larger values of in.tol will lead to larger theta.}
\item{comp.meth}{``exact'' or ``approximation''. Determines whether probabilities over the bootstrap distribution are computed exactly
or approximately. Exact computations are only feasible for small sample sizes.}
\item{Boot.supp}{Bootstrap approximation distribution (e.g., obtained by drawing vectors each with n independent entries distributed according to a Rademacher distribution, or a Mammen distribution); (matrix, n rows, columns = number of bootstrap samples). This is only used if comp.meth = ``approximation''; note that these vectors are further multiplied by wildmult weights.}
\item{checks}{TRUE (default) or FALSE. Determines whether input checks are conducted (TRUE) or not (FALSE).}
\item{cores}{Positive integer. Maximal number of cores used in the computations.}
}

\details{For details see the relevant sections in \enc{Pötscher}{Poetscher} 
and Preinerstorfer (2020).
}

\value{
\item{theta}{The value of theta; in case Assumptions are not satisfied NA is returned.}
\item{Aspt.sat}{TRUE if checks passed, FALSE if checks did not pass.}
\item{Max.ind}{The index that led to the minimum in the computation of theta. If Astp.sat = FALSE, then Max.ind = NA. Furthermore, Max.ind = NA in case theta was determined by the sufficient condition (k = q, boot.res.restr =       
boot.center.restr = TRUE) for theta = 1 discussed in the paper}
}

\examples{
v1 <- rnorm(10)
X <- cbind(rep(1, 10), rnorm(10), rnorm(10), sign(v1)*exp(v1))                         
R <- rbind(c(0, 0, 1, 0), c(0, 0, 0, 1))
r <- c(0, 0)

theta(X, R, r, 3, FALSE, "rademacher", 3, FALSE, TRUE, TRUE)

}

\references{
\enc{Pötscher}{Poetscher}, B. M. and Preinerstorfer, D. (2020). How Reliable are Bootstrap-based Heteroskedastiity Robust Tests? <arXiv:2005.04089>
}             