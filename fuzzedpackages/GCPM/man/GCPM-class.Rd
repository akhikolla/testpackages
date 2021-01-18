\name{GCPM-class}
\Rdversion{1.1}
\docType{class}
\alias{GCPM-class}
\title{Class \code{"GCPM"}}
\description{
The class represents a generalized credit portfolio framework. Users which are
not familiar with credit portfolio models in general and the CreditRisk+ model
as well as the CreditMetrics model in particular should refer to the references
given below. Models can be simulative or analytical (in case of a CreditRisk+
type model). The link function can be chosen to be either of the CreditRisk+ or
the CreditMetrics type. Counterparties' default distribution can be specified
to be either Bernoulli or Poisson, which is the default
distribution in the basic CreditRisk+ framework. 
}
\section{Objects from the Class}{
Objects can be created via the \code{init} function (see \code{\link{init}})
}
\section{Slots}{
  \describe{
    \item{\code{model.type}:}{Character value, specifying the model type.
    One can choose between \dQuote{simulative} and \dQuote{CRP} which
    corresponds to the analytical version of the CreditRisk+ model
    (see First Boston Financial Products, 1997)}
    \item{\code{default}:}{Character vector specifying the counterparties'
    default distribution (either \dQuote{Bernoulli} or \dQuote{Poisson})}
    \item{\code{link.function}:}{character value, specifying the type of the
    link function. One can choose between \dQuote{CRP}, which corresponds to
    \eqn{\overline{PD}=PD\cdot (w^Tx)} and \dQuote{CM} which corresponds to
    \eqn{\overline{PD}=\Phi\left(\frac{\Phi^{-1}PD-w^Tx}{\sqrt{1-w^T\Sigma w}}\right)},
    where PD is the original PD from portfolio data, x is the vector of sector
    drawings, \eqn{\Phi} is the CDF of the standard normal distribution, w is
    the vector of sector weights given in the portfolio data and \eqn{\Sigma}
    is the correlation matrix of the sector variables estimated from
    \code{random.numbers}. \dQuote{CRP} will be used automatically if
    \code{model.type} == "CRP".
}
    \item{\code{loss.unit}:}{numeric value used to discretize potential losses.}
    \item{\code{NS}:}{number of sectors}
    \item{\code{NC}:}{number of counterparties}
    \item{\code{name}:}{counterparties' names defined in the portfolio}
    \item{\code{NR}:}{counterparties' identification numbers defined in the
    portfolio}
    \item{\code{EAD}:}{counterparties' exposure at default defined in the
    portfolio}
    \item{\code{LGD}:}{counterparties' loss given default defined in the
    portfolio}
    \item{\code{PL}:}{counterparties' potential loss (\eqn{EAD*LGD})}
    \item{\code{PD}:}{counterparties' probability of default defined in the
    portfolio}
    \item{\code{business}:}{counterparties' business line defined in the
    portfolio}
    \item{\code{country}:}{counterparties' country defined in the portfolio}
    \item{\code{EL.analyt}:}{Expected loss calculated from portfolio data
    (without discretization)}
    \item{\code{EL}:}{Expected loss derived from loss distribution}
    \item{\code{nu}:}{multiples of \code{loss unit} representing discretized
    potential losses within an analytical CreditRisk+ type model}
    \item{\code{PL.disc}:}{counterparties' potential loss (\eqn{EAD*LGD})
    after discretization}
    \item{\code{PD.disc}:}{counterparties' probability of default defined in
    the portfolio after discretization}
    \item{\code{sec.var}:}{sector variances within an analytical CreditRisk+
    type model}
    \item{\code{sector.names}:}{sector names}
    \item{\code{SD.div}:}{diversifiable part of portfolio risk (measured by
    standard deviation) in case of a CreditRisk+ type model}
    \item{\code{SD.syst}:}{Non-diversifiable part of portfolio risk (measured
    by standard deviation) in case of a CreditRisk+ type model}
    \item{\code{SD.analyt}:}{portfolio standard deviation derived from
    portfolio data in case of a CreditRisk+ type model}
    \item{\code{SD}:}{portfolio standard deviation derived from loss
    distribution}
    \item{\code{W}:}{counterparties' sector weights}
    \item{\code{idiosyncr}:}{counterparties idiosyncratic weight in case of a
    CreditRisk+ type model}
    \item{\code{alpha.max}:}{maximum level of CDF of the loss distribution
    within an analytical CreditRisk+ type model}
    \item{\code{a}:}{internal parameter used to calculate risk contributions
    in case of an analytical CreditRisk+ type model}
    \item{\code{PDF}:}{probability density function of portfolio losses}
    \item{\code{CDF}:}{cumulative distribution function of portfolio losses}
    \item{\code{B}:}{internal parameter used to calculate risk contributions
    in case of an analytical CreditRisk+ type model}
    \item{\code{loss}:}{portfolio losses corresponding to \code{PDF} and
    \code{CDF}}
    \item{\code{random.numbers}:}{sector drawing in case of a simulative model}
    \item{\code{LHR}:}{likelihood ration of sector drawing in case of a
    simulative model}
    \item{max.entries}{
    numeric value defining the maximum number of loss scenarios stored to
    calculate risk contributions.}
    \item{\code{N}:}{number of simulations in case of a simulative model}
    \item{\code{scenarios}:}{scenarios (rows) of \code{random.numbers} used
    within the simulation of portfolio losses}
    \item{\code{seed}:}{parameter used to initialize the random number
    generator. If \code{seed} is not provided a value based on current system 
    time will be used.}
    \item{\code{loss.thr}:}{specifies a lower bound for portfolio losses to be
    stored in order to derive risk contributions on counterparty level. Using
    a lower value needs a lot of memory but will be necessary in order to
    calculate risk contributions on lower CDF levels. This parameter is used
    only if \code{model.type} == "simulative".}
    \item{\code{sim.losses}:}{simulated portfolio losses in case of a
    simulative model}
    \item{\code{CP.sim.losses}:}{simulated losses on counterparty level when
    the overall portfolio loss is greater or equal to \code{loss.thr}}
  }
}
\references{
Jakob, K. & Fischer, M. "GCPM: A flexible package to explore credit portfolio risk" Austrian Journal of Statistics 45.1 (2016): 25:44\cr
Morgan, J. P. "CreditMetrics-technical document." JP Morgan, New York, 1997\cr 
First Boston Financial Products, "CreditRisk+", 1997 \cr
Gundlach & Lehrbass, "CreditRisk+ in the Banking Industry", Springer, 2003\cr
}
\author{
Kevin Jakob
}
\seealso{
\code{\link{GCPM-package}}, \code{\link{init}}, \code{\link{analyze}}
}
\keyword{classes}
