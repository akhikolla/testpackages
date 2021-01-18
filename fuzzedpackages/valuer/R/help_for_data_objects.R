#'BBM2010 financial processes
#'@description
#'List of parameters to initialize a va_sde_engine object to
#'simulate the interest rate, volatility and log price processes
#'according to the stochastic differential equations specified
#'in BBM2010 -  See \bold{References}.
#'@format A list with elements:
#'\describe{
#' \item{[[1]]}{List of parameters for \code{\link[yuima]{simulate}}}
#' \item{[[2]]}{List of parameters for \code{\link{setModel}}}
#' \item{[[3]]}{Vector with indices indicating the interest rate and log price
#' in solve.variable \code{\link{setModel}}}
#'}
#'@references
#'\enumerate{
#' \item{[BBM2010]}{
#'  \cite{Bacinello A.R., Biffis E. e Millossovich P.
#'        "Regression-based algorithms
#'         for life insurance contracts with surrender guarantees".
#'         In: Quantitative Finance 10.9 (2010), pp. 1077-1090.}
#'        }
#'  }
"financials_BBM2010"



#'BBM2010 demographic processes
#'@description
#'List of parameters to initialize a va_sde_engine object to
#'simulate the intensity of mortality process
#'according to the stochastic differential equation specified
#'in BBM2010 -  See \bold{References}.
#'@format A list with elements:
#'\describe{
#' \item{[[1]]}{List of parameters for \code{\link[yuima]{simulate}}}
#' \item{[[2]]}{List of parameters for \code{\link{setModel}}}
#' \item{[[3]]}{Vector with indices indicating the intensity of mortality
#' in solve.variable \code{\link{setModel}}}
#'}
#'@references
#'\enumerate{
#' \item{[BBM2010]}{
#'  \cite{Bacinello A.R., Biffis E. e Millossovich P.
#'        "Regression-based algorithms
#'         for life insurance contracts with surrender guarantees".
#'         In: Quantitative Finance 10.9 (2010), pp. 1077-1090.}
#'        }
#'  }
"mortality_BBM2010"


#'BMOP2011 financial processes
#'@description
#'List of parameters to initialize a va_sde_engine object to
#'simulate the interest rate, volatility and log price processes
#'according to the stochastic differential equations specified
#'in BMOP2011 -  See \bold{References}.
#'@format A list with elements:
#'\describe{
#' \item{[[1]]}{List of parameters for \code{\link[yuima]{simulate}}}
#' \item{[[2]]}{List of parameters for \code{\link{setModel}}}
#' \item{[[3]]}{Vector with indices indicating the interest rate and log price
#' in solve.variable \code{\link{setModel}}}
#'}
#'@references
#'\enumerate{
#' \item{[BMOP2011]}{
#'  \cite{Bacinello A.R., Millossovich P., Olivieri A. e Pitacco E.
#'        "Variable annuities: a unifying valuation approach."
#'         In: Insurance: Mathematics and Economics 49 (2011), pp. 285-297.}
#'       }
#'  }
"financials_BMOP2011"



#'BMOP2011 demographic processes
#'@description
#'List of parameters to initialize a va_sde_engine object to
#'simulate the intensity of mortality process
#'according to the stochastic differential equation specified
#'in BMOP2011 - See \bold{References}.
#'@format A list with elements:
#'\describe{
#' \item{[[1]]}{List of parameters for \code{\link[yuima]{simulate}}}
#' \item{[[2]]}{List of parameters for \code{\link{setModel}}}
#' \item{[[3]]}{Vector with indices indicating the intensity of mortality
#' in solve.variable \code{\link{setModel}}}
#'}
#'@references
#'\enumerate{
#' \item{[BMOP2011]}{
#'  \cite{Bacinello A.R., Millossovich P., Olivieri A. e Pitacco E.
#'        "Variable annuities: a unifying valuation approach."
#'         In: Insurance: Mathematics and Economics 49 (2011), pp. 285-297.}
#'       }
#'  }
"mortality_BMOP2011"


#'BZ2016 financial processes
#'@description
#'List of parameters to initialize a \code{\link{va_sde_engine2}} object to
#'simulate the interest rate and log price processes being the
#'volatility constant. The interest rate and fund processes
#'follow the stochastic differential equations specified
#'in BMOP2011 -  See \bold{References}. The volatility is constant with
#'default value 0.2
#'@format A list with elements:
#'\describe{
#' \item{[[1]]}{List of parameters for \code{\link[yuima]{simulate}}}
#' \item{[[2]]}{List of parameters for \code{\link{setModel}}}
#' \item{[[3]]}{Vector with indices indicating the interest rate and log price
#' in solve.variable \code{\link{setModel}}}
#'}
#'@references
#'\enumerate{
#' \item{[BMOP2011]}{
#'  \cite{Bacinello A.R., Millossovich P., Olivieri A. e Pitacco E.
#'        "Variable annuities: a unifying valuation approach."
#'         In: Insurance: Mathematics and Economics 49 (2011), pp. 285-297.}
#'       }
#'  }
#'@examples
#'  #Sets the constant volatility to 0.3
#'  financials_BZ2016[[1]]$K <- 0.3 ^ 2
#'
"financials_BZ2016"

#'BZ2016bis financial processes
#'@description
#'List of parameters to initialize a \code{\link{va_sde_engine3}} object to
#'simulate the  log price and volatility processes which
#'follow the stochastic differential equations specified
#'in BMOP2011 -  See \bold{References}. The interest rate is constant
#'with default value 0.03.
#'@format A list with elements:
#'\describe{
#' \item{[[1]]}{List of parameters for \code{\link[yuima]{simulate}}}
#' \item{[[2]]}{List of parameters for \code{\link{setModel}}}
#' \item{[[3]]}{Vector with indices indicating the log price
#' in solve.variable \code{\link{setModel}}}
#'}
#'@references
#'\enumerate{
#' \item{[BMOP2011]}{
#'  \cite{Bacinello A.R., Millossovich P., Olivieri A. e Pitacco E.
#'        "Variable annuities: a unifying valuation approach."
#'         In: Insurance: Mathematics and Economics 49 (2011), pp. 285-297.}
#'       }
#'  }
#'@examples
#' #Sets the interest rate to 2%
#' financials_BZ2016bis[[1]]$r <- 0.02
#'
"financials_BZ2016bis"

