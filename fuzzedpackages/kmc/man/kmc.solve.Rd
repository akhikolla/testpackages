\name{kmc.solve}
\alias{kmc.solve}
\title{
Calculate NPMLE with constriants for right censored data
    %%  ~~function to do ... ~~
}
\description{
This function calculate the Kaplan-Meier estimator with mean constraints recursively. 
\deqn{El(F)=\prod_{i=1}^{n}(\Delta F(T_i))^{\delta_i}(1-F(T_i))^{1-\delta_i}}
with constraints
\deqn{\sum_i g(T_i)\Delta F(T_i)=0,\quad,i=1,2,...}
It uses Lagrange multiplier directly.
}
\usage{
kmc.solve(x, d, g, em.boost = T, using.num = T, using.Fortran =
                 T, using.C = F, tmp.tag = T, rtol = 1e-09, control =
                 list(nr.it = 20, nr.c = 1, em.it = 3),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     Non-negative real vector. The observed time.
}
  \item{d}{
      0/1 vector. Censoring status indictator, 0: right censored; 1 uncensored
}
  \item{g}{
list of contraint functions. It should be a list of functions list(f1,f2,...) 
}
  \item{em.boost}{
  	 A logical value. It determines whether the EM algorithm is used to get the initial value, default=TRUE. See 'Details' for EM control.
}
\item{using.num}{
  	 A logical value. It determines whether the numeric derivatives is used in iterations, default=TRUE.
}
\item{using.Fortran}{
  	 A logical value. It determines whether Fortran is used in root solving, default=F.
}
\item{using.C}{
	A logical value. It determines whether to use Rcpp in iteraruib, default=T. This option will promote the computational efficiency of the KMC algorithm. Development version works on one constraint only, otherwise it will generate a Error information. It won't work on using.num=F.
}
\item{tmp.tag}{
	Development version needs it, keep it as TRUE.
}
\item{rtol}{
	Tolerance used in rootSolve(multiroot) package, see 'rootSolve::multiroot'.
}
\item{control}{
  A list. The entry nr.it controls max iterations allowed in N-R algorithm default=20; nr.c is the scaler used in N-R algorithm default=1; em.it is max iteration if use EM algorithm (em.boost) to get the initial value of lambda, default=3.
}
\item{...}{
Unspecified yet.
}
}
\details{
The function check_G_mat checks whether the solution space is null or not under the constraint. But due to the computational complexity, it will detect at most two conditions.
}
\value{
 A list with the following components:
 \item{loglik.ha}{
 	The log empirical  likelihood without constraints
 }
 \item{loglik.h0}{
 	The log empirical likelihood with constraints
 }
 \item{"-2LLR"}{
 	The -2 Log empirical likelihood ratio
 }
 \item{phat}{
 	\deqn{\Delta F(T_i)}
 }
 \item{pvalue}{
  The p-value of the test
 }
 \item{df}{
  Degree(s) of freedom. It equals the number of constraints.
 }
 \item{lambda}{
  The lambda is the Lagrangian multiplier described in reference.
 }
}
\references{
Zhou, M. and Yang, Y. (2015).
A recursive formula for the Kaplan-Meier estimator with mean constraints and its application to empirical likelihood
\emph{Computational Statistics}
\bold{Online ISSN} 1613-9658.
}

\author{
    Mai Zhou(mai@ms.uky.edu), Yifan Yang(yifan.yang@uky.edu)
}


\seealso{
	  	\code{\link{plotkmc2D}}.
}
\examples{
x <- c( 1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3) # positive time
d <- c( 1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)               # status censored/uncensored

#################
# dim =1
#################

f<-function(x) { x-3.7}                                       # \sum f(ti) wi ~ 0 
g=list( f=f) ;                                                #define constraint as a list

kmc.solve( x,d,g) ;                                           #using default
kmc.solve( x,d,g,using.C=TRUE) ;                              #using Rcpp

#################
# dim =2
#################

myfun5 <- function( x)  { 
 x^2-16.5
} 

g = list( f1=f,f2=myfun5) ;                                     #define constraint as a list

re0 <- kmc.solve( x,d,g);

#################
# Print Estimation and other information 
# with option: digits = 5
#################
f_print <- function(x, digits = 5){
        cat("\n---------------------------------------------------------------------------------\n")
        cat("A Recursive Formula for the Kaplan-Meier Estimator with Constraint\n")
        cat("Information:\n")
        cat("Number of Constraints:\t", length(x$g), "\n")
        cat("lamda(s):\t",x$lambda,'\n');
        cat("\n---------------------------------------------------------------------------------\n")
        names <- c("Log-likelihood(Ha)", "Log-likelihood(H0)",
        "-2LLR", paste("p-Value(df=", length(x$g), ")",sep = ""))
        re <- matrix(c(x[[1]], x[[2]], x[[3]], 1 - pchisq(x[[3]],
        length(x$g))), nrow = 1)
        colnames(re) <- names
        rownames(re) <- "Est"
        print.default(format(re, digits = digits), print.gap = 2,
        quote = FALSE, df = length(x$g))
        cat("---------------------------------------------------------------------------------\n")
}
f_print(re0)

}