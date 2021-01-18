\name{add.noise}
\alias{add.noise}
\alias{add.house.noise}
\alias{add.candle.noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Apply Noise to Discrete-Valued Tables
}
\description{
    The function can apply two types of noise to contingency tables of discrete values. A house noise model is designed for ordinal variables; a candle noise model is for categorical variables. Noise is applied independently for each data point in a table.
}
\usage{
add.noise(tables, u, noise.model, margin=0)
add.house.noise(tables, u, margin=0)
add.candle.noise(tables, u, margin=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tables}{
    a list of tables or one table. A table can be either a matrix or a data frame of integer values.
  }
  \item{u}{
    a numeric value between 0 and 1 to specify the noise level to be applied to the input tables. See Details.
  }
  \item{noise.model}{
    a character string indicating the noise model of either \code{"house"} for ordinal variables or \code{"candle"} for categorical variables. See Details.
  }
  \item{margin}{
    a value of either 0, 1, or 2. Default is 0.

    0: noise is applied along both rows and columns in a table. The sum of values in the table is the same before and after noise application.

    1: noise is applied along each row. The sum of each row is the same before and after noise application.

    2: noise is applied along each column. The sum of each column is the same before and after noise application.
  }
}

\details{
  Each noise model defines a conditional probability function of a noisy version given an original discrete value and a noise level. In the house noise model for ordinal variables, defined in \insertCite{zhang2015chinet}{FunChisq}, the probability decreases as the noisy version deviates from the original ordinal value. The shape of the function is like a pitched house roof. In the candle noise model for categorical variables, the probability of the noisy version for any value other than the original categorical value is the same given the noise level. The function shape is like a candle.

  At a minimum level of 0, no noise is applied on the input table(s). A maximum level of 1 indicates that the original sample will be changed to some other values with a probability of 1. For a discrete random variable of two possible values, a noise level of 1 will flip the values and create a non-random pattern; a noise level of 0.5 creates the most random pattern.
}

\value{
    If \code{tables} is a list, the function returns a list of tables with noised applied. If \code{tables} is a numeric matrix or a data frame, the function returns one table with noise applied.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{
%% ~put references to the literature/web site here ~
%Zhang, Y., Liu, Z. L. and Song, M. (2015) ChiNet uncovers rewired transcription subnetworks in tolerant yeast for advanced biofuels conversion. \emph{Nucleic Acids Research} \bold{43}(9), 4393--4407. Retrieved from \url{https://nar.oxfordjournals.org/content/43/9/4393.long}

\insertAllCited{}

}

\author{
    Hua Zhong, Yang Zhang, and Joe Song.
}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{simulate_tables}}.
}
\examples{
  # Example 1. Add house noise to a single table

  # Create a 4x4 table
  t <- matrix(c(3,0,0,0,
                0,2,2,0,
                0,0,0,4,
                3,3,2,0),
              nrow=4, ncol=4, byrow=TRUE)
  # Two ways to apply house noise at level 0.1 along both rows
  # and columns of the table:
  add.noise(t, 0.1, "house", 0)
  add.house.noise(t, 0.1, 0)


  # Example 2. Add candle noise to a list of tables

  # Create a list of tables
  t.list <- list(t+5, t*10, t*2)
  # Two ways to apply candle noise at level 0.2 along the rows
  # of the table:
  add.noise(t.list, 0.2, "candle", 1)
  add.candle.noise(t.list, 0.2, 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{ datagen }
\concept{ noise model }
