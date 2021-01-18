\name{rplot}
\alias{rplot}

\title{
Recombination Plot}
\description{
This function creates a plot which shows the sum of all recombination events across a half-sib family.}
\usage{
rplot(x, distance, start = 1, end = ncol(x), maximum = 100,
overwrite = FALSE, method = "constant")
}

\arguments{
  \item{x}{\code{matrix} of half-sib genotypes  (one half-sib per row, with SNP ordered by mapping position in the columns. Data should be numeric. Use 0, 1 and 2 for respectively AA, AB and BB. Use 9 for missing data).
}
\item{distance}{\code{integer} of physical distances between markers
}


  \item{start}{\code{integer}
 first marker selected for the plot}
 
  \item{end}{\code{integer}
 last marker selected for the plot}
 
  \item{maximum}{\code{integer}
maximum number of recombinations to show (higher recombination rates will be omitted from the plot)}

  \item{overwrite}{\code{logical}
draw a diagram over the current diagram  (default \code{FALSE})}

 \item{method}{\code{character} please refer to the \code{\link{pm}} document}
 
}
\examples{
genotype <- matrix(c(
0,2,0,1,0,
2,0,1,2,2,
2,2,1,0,2,
2,2,1,1,1,
0,0,2,1,0), ncol = 5, byrow = TRUE) # each row contains the SNP of individuals

rplot(genotype, c(1,2,3,4,8)) 
}


\keyword{snp}
\keyword{halfsib}
\keyword{phase}
\keyword{recombination}
