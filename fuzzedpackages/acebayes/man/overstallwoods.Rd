\name{overstallwoods}
\alias{utilcomp18bad}
\alias{optdescomp18bad}
\alias{utilcomp15bad}
\alias{optdescomp15bad}
\alias{optdescomp15sig}
\alias{utilcomp15sig}
\alias{utilcomp15sigDRS}
\alias{optdescomp15sigDRS}

\alias{utillrbad}
\alias{optdeslrbad}
\alias{utillrsig}
\alias{inideslrsig}
\alias{optdeslrsig}
\alias{utillrbaa}
\alias{optdeslrbaa}
\alias{utillrnsel}
\alias{inideslrnsel}
\alias{optdeslrnsel}
\alias{utilhlrbad}
\alias{optdeshlrbad}
\alias{utilhlrsig}
\alias{inideshlrsig}
\alias{optdeshlrsig}
\alias{utilhlrbaa}
\alias{optdeshlrbaa}
\alias{utilhlrnsel}
\alias{inideshlrnsel}
\alias{optdeshlrnsel}

\alias{utilbeetle}
\alias{optdesbeetle}

\alias{utillinmod}
\alias{optdeslinmod}

\title{
Functions implementing the examples of Overstall & Woods (2017).
}
\encoding{UTF-8}
\description{
This suite of functions implement the examples in Overstall & Woods (2017).
}
\usage{

######## Compartmental model #################################

utilcomp18bad(d, B)
optdescomp18bad(type = "ACE")
utilcomp15bad(d, B)
optdescomp15bad()
utilcomp15sig(d, B)
optdescomp15sig()
utilcomp15sigDRS(d, B)
optdescomp15sigDRS()

######## Logistic regression model ###########################

utillrbad(d, B)
optdeslrbad(n, type = "ACE")
utillrsig(d, B)
inideslrsig(n, rep)
optdeslrsig(n)

utilhlrbad(d, B)
optdeshlrbad(n)
utilhlrsig(d, B)
inideshlrsig(n, rep)
optdeshlrsig(n)

utillrbaa(d, B)
optdeslrbaa(n)
utillrnsel(d, B)
inideslrnsel(n, rep)
optdeslrnsel(n)

optdeshlrbaa(n)
utilhlrbaa(d, B)
utilhlrnsel(d, B)
inideshlrnsel(n, rep)
optdeshlrnsel(n)

######## Beetle mortality experiment #########################

utilbeetle(d, B)
optdesbeetle(n)

######## Linear model ########################################

utillinmod(d, B)
optdeslinmod(n, type = "ACE")

##############################################################
}

\arguments{
  \item{d}{
An \code{n} by \code{k} matrix specifying the design matrix, where \code{n} and \code{k} denote the number of runs and factors, respectively. See \bold{Details} and \bold{Value} for the values that \code{n} and \code{k} can take for each of the examples. Each element of \code{d} is scaled to the interval [-1,1].
}
  \item{n}{
The number of runs ins the experiment.
}
  \item{rep}{
A scalar integer in the set \eqn{\left\{1,\dots,20\right\}}{{1,...,20}} specifying the initial design.
}
  \item{B}{
A scalar integer specifying the Monte Carlo sample size.
}
  \item{type}{
An optional character argument specifying which design to return. 

For \code{optdeslrbad}, possible values are \code{c("ACE","Gotwalt","Atkinson")}. If \code{"ACE"} (the default) then the design found by the ACE algorithm will be returned. If \code{"Gotwalt"} then the design published in Gotwalt et al (2009) is returned. If \code{"Atkinson"} then the design found by Atkinson et al (1993) is returned.

For \code{optdeslinmod}, possible values are \code{c("ACE","BoxDraper")}. If \code{"ACE"} (the default) then the design found by the ACE algorithm will be returned. If \code{"BoxDraper"} then the true optimal design, as found by Box & Draper (1971), will be returned.
}
}
\details{
This section provides details on the examples considered and the functions implemented in \code{acebayes}.

\bold{Compartmental model}

Compartmental models are used in Pharmokinetics to study how materials 
flow through an organism. A drug is administered to an individual or animal and then the 
amount present at a certain body location is measured at a set of \code{n} pre-determined sampling 
times (in hours). There is one design variable: sampling time. Therefore the design matrix \code{d} is 
an \code{n} by 1 matrix with elements controlling the \code{n} sampling times, i.e. the number of factors is
\code{k=1}.

In Overstall & Woods (2017), two different compartmental model examples are considered. The first (in the Supplementary Material) 
comes from Atkinson et al (1993) and Gotwalt et al (2009) where there are \code{n = 18} sampling times and interest 
lies in finding a Bayesian D-optimal design. The functions whose name includes \code{"comp18"} refer to this example.

The second example (in Section 3.2) comes from Ryan et al (2014), where there are \code{n = 15} sampling times and 
the ultimate interest lies in finding an optimal design under the Shannon information gain utility. Also considered is the
Bayesian D-optimal design. The functions whose name includes \code{"comp15"} refer to this example. Note that Ryan et al (2014) used a dimension reduction scheme (DRS) to find optimal designs. The function whose name is suffixed by \code{"DRS"} refer to this situation.

\bold{Logistic regression model}

In Section 3.3 of Overstall & Woods (2017), a first-order logistic regression model in \code{k = 4} factors and \code{n} runs is considered. Woods et al (2006) and Gotwalt et al (2009) considered generating Bayesian D-optimal designs for \code{n = 16} and \code{n = 48}. Overstall & Woods (2017) extended this example by considering Bayesian A-optimal, Shannon information gain (SIG) and negative squared error loss (NSEL) utility functions, a range of number of runs from 6 to 48, and "random effects" to form a hierarchical logistic regression model.

\bold{Beetle mortality experiment}

Overstall & Woods (2017, Section 3.4) considers generating an optimal design for a follow-up experiment. The original design and data (Bliss, 1935) involves administering different doses of poison to \eqn{N = 8} groups of beetles. The number of beetles that die in each group are recorded. Six different models are considered formed from the combination of three link functions and two linear predictors (following the analysis of O'Hagan & Forster, 2004). Interest lies in the quantity known as lethal dose 50 (LD50) which is the dose required to kill 50\% of the beetles and is a function of the model parameters for a given model. Consider finding an optimal design for estimating LD50 under the negative squared error loss (NSEL) function for \code{n} new doses of poison (i.e. \code{k = 1} factor). The prior distribution is equivalent to the posterior distribution arising from the original data and includes model uncertainty.

\bold{Linear model}

In the supplementary material, Overstall & Woods (2017) considers finding D-optimal designs for a second-order (i.e. \code{k = 2}) response surface in \code{n=6,7,8,9} runs. Note that the D-optimal design is equivalent to the optimal design under Shannon information gain and a non-informative prior distribution.

The expected utility function in this case is available in closed form, i.e. it does not require approximation. Box & Draper (1971) found optimal designs analytically for the number of runs considered here. Overstall & Woods (2017) use this example to demonstrate the efficacy of the ACE algorithm.

}

\value{

\bold{Compartmental model}

For the example in the Supplementary Material;
\itemize{
\item{
The function \code{utilcomp18bad} will return a vector of length \code{B} where each element is the value of the Bayesian D-optimal utility function (i.e. the log determinant of the Fisher information matrix) evaluated at a sample of size \code{B} generated from the prior distribution of the model parameters.}
\item{
The function \code{optdescomp18bad} will return an 18 by 1 matrix giving the optimal design (specified by the argument \code{type}). The elements will be scaled to be in the interval [-1, 1], i.e. a -1 corresponds to a sampling times of 0 hours, and 1 corresponds to a time of 24 hours.}}
For the example in Section 3.2;
\itemize{
\item{
The function \code{utilcomp15bad} will return a vector of length \code{B} where each element is the value of the Bayesian D-optimal utility function (i.e. the log determinant of the Fisher information matrix) evaluated at a sample of size \code{B} generated from the prior distribution of the model parameters.}
\item{
The function \code{optdescomp15bad} will return an 15 by 1 matrix giving the optimal design (found using ACE) under Bayesian D-optimality. The elements will be scaled to be in the interval [-1, 1], i.e. a -1 corresponds to a sampling times of 0 hours, and 1 corresponds to a time of 24 hours.}
\item{
The function \code{utilcomp15sig} will return a vector of length \code{B} where each element is the value of the SIG utility function evaluated at a sample of size \code{B} generated from the joint distribution of model parameters and unobserved responses.}
\item{
The function \code{optdescomp15sig} will return an 18 by 1 matrix giving the optimal design (found using ACE) under the SIG utility. The elements will be scaled to be in the interval [-1, 1], i.e. a -1 corresponds to a sampling times of 0 hours, and 1 corresponds to a time of 24 hours.}
\item{
The function \code{utilcomp15sigDRS} will return a vector of length \code{B} where each element is the value of the SIG utility function (where a DRS has been used) evaluated at a sample of size \code{B} generated from the joint distribution of model parameters and unobserved responses. Here the Beta DRS (see Overstall & Woods, 2017) has been used so \code{d} should be a 2 by 1 matrix containing the positive beta parameters.}
\item{
The function \code{optdescomp15sigDRS} will return a 2 by 1 matrix giving the optimal design (found using ACE) under the SIG utility, where a DRS has been used. The elements correspond to the parameters of a beta distribution.}
}

\bold{Logistic regression model}

A function whose name includes \code{"lr"} refers to standard logistic regression, whereas \code{"hlr"} refers to hierarchical logistic regression. Under standard logistic regression the possible values for the argument \code{n} can be any even integer between 6 and 48. For hierarchical logistic regression, \code{n} can be any integer divisible by 6 between 12 and 48. The function name also indicates the utility function:
\itemize{
\item{\code{"bad"}}{
Bayesian D-optimal}
\item{\code{"baa"}}{
Bayesian A-optimal}
\item{\code{"sig"}}{
Shannon information gain}
\item{\code{"nsel"}}{
Negative squared error loss}}
The functions prefixed by \code{"util"} will return a vector of length \code{B} where each element is the utility function evaluated at a sample generated from the prior distribution of model parameters (for Bayesian D- and A-optimality) or the joint distribution of model parameters and unobserved responses (for SIG and NSEL).

The functions prefixed by \code{"optdes"} will return an \code{n} by \code{k = 4} matrix giving the optimal design found by ACE. The designs given by this function are those reported on in Overstall & Woods (2017). The function \code{optdeslrbad} will also return designs (for \code{n = 16, 48}) found by Woods et al (2006) and Gotwalt et al (2009) by specifying the \code{type} argument appropriately.

The functions prefixed by \code{"inides"} will return an \code{n} by \code{k = 4} matrix giving an initial design for ACE to find the optimal designs under the SIG and NSEL utility functions. These are 20 designs found using ACE under approximations to the Bayesian A- and D-optimal utility functions, respectively. The argument \code{rep} specifies which of these 20 designs to use.

\bold{Beetle mortality experiment}

The function \code{utilbeetle} will return a vector of length \code{B} where each element is the value of the utility function for a sample generated from the joint distribution of the model parameters, model and unobserved responses.

The function \code{optdesbeetle} will return an \code{n} by 1 matrix giving the optimal design under the NSEL utility function (found using ACE) for estimating the LD50. The elements will be scaled to be in the interval [-1, 1], where -1 corresponds to a dose of 1.6907, 0 to a dose of 1.7873 and 1 to a dose of 1.8839. The designs given by this function are those reported in Overstall & Woods (2017) for \code{n} = 1, \ldots, 10.

\bold{Linear model}

The function \code{utillinmod} will return a vector of length \code{B} where each element is a realisation of a stochastic approximation to the expected utility.

The function \code{optdeslinmod} will return an \code{n} by 2 matrix giving the D-optimal design (specified by the argument \code{type}). If \code{type = "ACE"}, the designs returned by this function are those found using the ACE algorithm and are reported in the Supplementary Material of Overstall & Woods (2017), and if \code{type = "BoxDraper"}, the designs returned are the exact D-optimal designs.

}

\references{
Atkinson, A., Chaloner, K., Herzberg, A., & Juritz, J. (1993). Experimental Designs for Properties
of a Compartmental Model. \emph{Biometrics}, \bold{49}, 325-337.

Bliss, C. (1935). The calculation of the dosage-mortality curve. \emph{Annals of Applied Biology},
\bold{22}, 134-167.

Box, M. & Draper, N. (1971). Factorial designs, the \eqn{|F^T F|} criterion and some related matters. 
\emph{Techometrics}, \bold{13}, 731-742.

Gotwalt, C., Jones, B. & Steinberg, D. (2009). Fast Computation of Designs Robust to Parameter
Uncertainty for Nonlinear Settings. \emph{Technometrics}, \bold{51}, 88-95.

O'Hagan, A, & Forster, J.J. (2004). \emph{Kendall's Advanced Theory of Statistics, Volume 2B: Bayesian
Inference. 2nd edition}. John Wiley & Sons.

Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, \bold{59}, 458-470.

Ryan, E., Drovandi, C., Thompson, M., Pettitt, A. (2014). Towards Bayesian experimental design
for nonlinear models that require a large number of sampling times. \emph{Computational
Statistics and Data Analysis}, \bold{70}, 45-60.

Woods, D.C., Lewis, S., Eccleston, J., Russell, K. (2006). Designs for Generalized Linear Models
With Several Variables and Model Uncertainty. \emph{Technometrics}, \bold{48}, 284-292.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}

\seealso{
\code{\link{ace}}, \code{\link{pace}}.
}
\examples{

######## Compartmental model #################################

set.seed(1)
## Set seed for reproducibility

d<-optimumLHS(n = 18, k = 1) * 2 - 1
## Generate an 18-run design.

u<-utilcomp18bad(d = d, B = 20000)
## Calculate the D-optimal utility function for a 
## sample of size 20000. 

u[1:5]
## Look at first 5 elements.
#[1] 14.283473 10.525390  4.126233  7.061498 12.793245

d0<-optdescomp18bad( )
u0<-utilcomp18bad(d = d0, B = 20000)
## Optimal design found by ACE and calculate the D-optimal 
## utility function for a sample of size 20000.

u0[1:5]
## Look at first 5 elements.
#[1] 15.04721 15.37141 16.84287 14.06750 14.01523

mean(u)
mean(u0)
## Calculate expected Bayesian D-optimal utility.

d<-matrix(runif(2), ncol = 1)
## Generate two beta parameters.

u<-utilcomp15sigDRS(d = d, B = 5)
u
## Calculate the SIG utility function for a 
## sample of size 5.
#[1] 17.652044  4.878998 19.919559 22.017760  5.600473

######## Logistic regression model ###########################

set.seed(1)
## Set seed for reproducibility

d<-optimumLHS(n = 16,k = 4) * 2 - 1
## Generate an 16-run design.

u<-utillrbad(d = d, B = 20000)
## Calculate the D-optimal utility function for a 
## sample of size 20000. 

u[1:5]
## Look at first 5 elements.
#[1] -11.630683  -5.748912  -9.554413 -10.150132  -7.940938

d0<-optdeslrbad(16)
u0<-utillrbad(d = d0, B = 20000)
## Optimal design found by ACE and calculate the D-optimal 
## utility function for a sample of size 20000.

u0[1:5]
## Look at first 5 elements.
#[1] -4.644116 -2.411431 -4.999891 -2.906558 -2.282687

mean(u)
mean(u0)
## Calculate expected Bayesian D-optimal utility.
#[1] -9.38253
#[1] -2.992012

######## Beetle mortality experiment #########################

set.seed(1)
## Set seed for reproducibility

d<-optimumLHS(n = 10, k = 1)*2-1
## Generate a design of 10 doses with elements in [-1, 1]

utilbeetle(d = d, B = 5)
## Calculate the utility function for a sample of size 5.

#-4.720491e-06 -1.198955e-05 -1.681380e-05 -3.123498e-06 -1.412722e-05

d0<-optdesbeetle(10)
d0
## Print out optimal design from Overstall & Woods (2017) for 10 doses 

0.5*( d0 + 1)*( 1.8839 - 1.6907 ) + 1.6907
## On original scale.
#          [,1]
# [1,] 1.769957
# [2,] 1.769520
# [3,] 1.768641
# [4,] 1.777851
# [5,] 1.768641
# [6,] 1.769520
# [7,] 1.777851
# [8,] 1.765997
# [9,] 1.768641
#[10,] 1.768641

######## Linear model ########################################

set.seed(1)
## Set seed for reproducibility

d<-cbind(rep(c( -1, 0, 1), each = 3), rep(c( -1, 0, 1), 3))
## Create a 9-run design which is the true D-optimal design

utillinmod(d = d, B = 5)
## Calculate the approximation to the true expected D-optimal utility 
## function for a sample of size 5. 

#[1]  7.926878  8.736976  7.717704 10.148613  8.882840

d0<-optdeslinmod(9)
## Optimal D-optimal design found using ACE

X<-cbind(1, d, d^2, d[,1] * d[,2])
X0<-cbind(1, d0, d0^2, d0[,1] * d0[,2])
## Calculate model matrices under both designs

detX<-determinant( t(X) \%*\% X)$modulus[1]
detX0<-determinant( t(X0) \%*\% X0)$modulus[1]
## Calculate true expected D-optimal utility function for both designs

100 * exp( 0.2 * ( detX0 - detX ))
## Calculate D-efficiency of ACE design.

# 99.93107
}
