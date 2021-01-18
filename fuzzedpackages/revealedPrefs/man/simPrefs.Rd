\name{simPrefs}
\alias{simPrefs}

\title{
Generate random data according to a given matrix of direct preferences.
}
\description{
Function for generating random data (prices and quantities) consistent with a given direct preferences matrix.
}
\usage{
simPrefs(pref.mat, ngoods, afriat.par= 1, 
         qmin= 0, qmax= 1, pmin= 0, pmax= 1,
         maxit= 100, verbose= FALSE)
}
\arguments{
  \item{pref.mat}{the desired matrix of direct preferences (must be a square matrix, see details),}
  \item{ngoods}{the number of goods in the dataset (number of columns in the quantities and prices matrices),}
  \item{afriat.par}{the Afriat parameter, a real number in [0,1], which allows a certain level of error in the optimization of choices; default is 1, ie. no optimization error allowed,}
  \item{qmin}{minimum quantities for each good,}
  \item{qmax}{maximum quantities for each good,}
  \item{pmin}{minimum prices for each good,}
  \item{pmax}{maximum prices for each good,}
  \item{maxit}{maximum number of iterations of the optimization algorithm (default to 100),}
  \item{verbose}{logical value: if \code{TRUE} a trace of the optimization algorithm is printed.}
}
\value{
  \item{\code{x}}{numeric matrix of generated quantities, or \code{NULL} if optimization fails,}
  \item{\code{p}}{numeric matrix of generated prices, or \code{NULL} if optimizaiton fails.}
}
\details{
The input \code{pref.mat} must be a square matrix, with number of rows equal to the number of desired observations. All entries must be either 0, 1, or 2. The interpretation of the matrix is the same as in function \code{directPrefs}:

\code{pref.mat[i, j] = 0} if bundle \code{i} is not revealed prefered to bundle \code{j}

\code{pref.mat[i, j] = 1} if bundle \code{i} is revealed prefered to bundle \code{j}

\code{pref.mat[i, j] = 2} if bundle \code{i} is strictly revealed prefered to bundle \code{j}.
  
All diagonal entries of \code{pref.mat} must be 1 (each bundle is revealed prefered to itself), except when \code{afriat.par} is strictly less than 1.

The simulated data (quantities and prices) are obtained by particle swarm optimization (of package \code{pso}). Fitness must reach 0 for the data to be consistent with the preference matrix. If optimization fails, a warning is issued.
}
\references{
Varian, H. R. (1982) The Nonparametric Approach to Demand Analysis, 
\emph{Econometrica}, 50(4):945-973.

Varian, H. R. (1984) \emph{Microeconomic Analysis}. New York/London:
Norton, 2nd edition, pp 141-143.
}
\author{
Julien Boelaert \email{jubo.stats@gmail.com}
}
\seealso{
See \link{directPrefs} for computation of preferences.
}
\examples{
# Generate 3 observations of 5 goods, according to a given preference matrix
pref.mat <- matrix(c(1,0,0,2,1,0,2,2,1), 3)
simdata <- simPrefs(pref.mat = pref.mat, ngoods = 5)
}
