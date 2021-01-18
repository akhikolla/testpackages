\name{WeeklyReturnsData}
\alias{WeeklyReturnsData}
\docType{data}
\title{
Logarithmic Returns of Exchange Rates of 7 Currencies Against US Dollar
}

\description{
This data set has logarithmic returns of exchange rates of 7 currencies against US dollar extracted from the International Monetary Fund's (IMF) database. These currencies are Australian Dollar (AUD), Canadian Dollar (CAD), Norwegian Kroner (NOK), Singapore Dollar (SGD), Swedish Kroner (SEK), Swiss Franc (CHF) and British Pound (GBP).
}

\usage{
data("WeeklyReturnsData")
}

\format{
  An object of class ts with 605 observations on the following 7 variables.
  \describe{
    \item{\code{AUD}}{The weekly logarithmic returns \eqn{\mathbf r_{AUD, t}}{r_(AUD,t)} of the exchange rates of AUD against US Dollar.}
    \item{\code{CAD}}{The weekly logarithmic returns \eqn{\mathbf r_{CAD, t}}{r_(CAD,t)} of the exchange rates of CAD against US Dollar.}
    \item{\code{NOK}}{The weekly logarithmic returns \eqn{\mathbf r_{NOK, t}}{r_(NOK,t)} of the exchange rates of NOK against US Dollar.}
    \item{\code{SGD}}{The weekly logarithmic returns \eqn{\mathbf r_{SGD, t}}{r_(SGD,t)} of the exchange rates of SGD against US Dollar.}
    \item{\code{SEK}}{The weekly logarithmic returns \eqn{\mathbf r_{SEK, t}}{r_(SEK,t)} of the exchange rates of SEK against US Dollar.}
    \item{\code{CHF}}{The weekly logarithmic returns \eqn{\mathbf r_{CHF, t}}{r_(CHF,t)} of the exchange rates of CHF against US Dollar.}
    \item{\code{GBP}}{The weekly logarithmic returns \eqn{\mathbf r_{GBP, t}}{r_(GBP,t)} of the exchange rates of GBP against US Dollar.}
    }
}

\details{
The daily exhange rates of the currencies against US Dollar from March 22, 2000 to October
26, 2011 are extracted from the International Monetary Fund's (IMF) Exchange Rates database from \url{https://www.imf.org/external/np/fin/ert/GUI/Pages/CountryDataBase.aspx}. These rates are representative rates (currency units per US Dollar), which are reported daily to the IMF by the issuing central bank.

The weekly averages of these exchange rates are then calculated. The logarithmic returns of the average weekly exchange rates are calculated for the currency \eqn{j} as follows.

Let \eqn{\mathbf x_{j, t}}{x_(j,t)} be the exchange rates of \eqn{j} against US Dollar. Then
\deqn{\mathbf r_{j, t} = \textrm{ log }(\mathbf x_{j, t}) - \textrm{ log } (\mathbf x_{j, t-1}),}{r_(j,t) = log x_(j, t)) - log(x_(j,t-1)),}
where \eqn{t = 1, \ldots, 605} and \eqn{j = AUD, CAD, NOK, SGD, SEK, CHF, GBP}.
The six missing values in \eqn{\mathbf r_{SEK, t}}{r_(SEK,t)} are changed to 0. The assumption used here is that there has not been any change from the previous week.

The weekly returns data is then changed to a multivariate time series object of class \code{\link{ts}}. The resulting \code{\link{ts}}
object is then dataset \code{WeeklyReturnsData}.

An example analysis of the data is given in Miettinen et al. (2018). Same data has also been used in Hu and Tsay (2014).
}

\source{
International Monetary Fund (2017), \emph{IMF Exchange Rates}, \url{https://www.imf.org/external/np/fin/ert/GUI/Pages/CountryDataBase.aspx}

For IMF Copyrights and Usage, and special terms and conditions pertaining to the use of IMF data, see \url{https://www.imf.org/external/terms.htm}
}

\references{
Miettinen, M., Matilainen, M., Nordhausen, K. and Taskinen, S. (2020), \emph{Extracting Conditionally Heteroskedastic Components Using Independent Component Analysis}, Journal of Time Series Analysis,41, 293--311.

Hu and Tsay (2014), \emph{Principal Volatility Component Analysis}, Journal of Business & Economic Statistics, 32(2), 153--164.
}

\examples{
plot(WeeklyReturnsData)

res <- gSOBI(WeeklyReturnsData)
res

coef(res)
plot(res)
head(bss.components(res))
}

\keyword{datasets}