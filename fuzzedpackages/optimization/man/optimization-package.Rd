\name{optimization-package}
\alias{optimization-package}
\docType{package}
\title{
\packageTitle{optimization}
}
\description{
\packageDescription{optimization}
}
\details{
Package: optimization
Type:    Package
Version: 1.0-6
Date:    2017-09-23
License: GPL-2
}
\author{
\packageAuthor{optimization}

Maintainer: \packageMaintainer{optimization}
}
\references{
Corana, A., Marchesi, M., Martini, C. and Ridella, S. (1987), Minimizing Multimodal Functions of Continuous Variables with the 'Simulated Annealing' Algorithm. ACM Transactions on Mathematical Software, 13(3):262-280.

Gao, F. and Han, L. (2012). Implementing the nelder-mead simplex algorithm with adaptive parameters. Computational Optimization and Applications, 51(1):259 277.

Geiger, C. and Kanzow, C. (1999). Das nelder-mead-verfahren. Numerische Verfahren zur Loesung unregestrier er Optimierungsaufgaben.

Kirkpatrick, S., Gelatt, C. D. and Vecchi, M. P. (1983). Optimization by Simulated Annealing. Science, 220(4598): 671-680.

Nelder, J. and Mead, R. (1965). A simplex method for function minimization. Computer Journal, 7(4).

Pronzato, L., Walter, E., Venot, A. and Lebruchec, J.-F. (1984). A general-purpose global optimizer: Implimentation and applications. Mathematics and Computers in Simulation, 26(5):412-422.
}

\seealso{
\code{\link{optim_nm}}, \code{\link{optim_sa}}, \code{\link{optim}},  \code{\link{plot}}
}


\examples{
hi <- function(x){(x[1]**2 + x[2] - 11)**2 + (x[1] + x[2]**2 -7)**2}
optim_nm(fun = hi, k = 2)
optim_sa(fun = hi, start = c(runif(2, min = -1, max = 1)),
  trace = FALSE,
  lower = c(-4, -4),
  upper = c(4, 4),
  control = list(dyn_rf = FALSE,
    rf = 1.2,
    t0 = 10,
    nlimit = 100,
    r = 0.6,
    t_min = 0.1
  )
)
}
