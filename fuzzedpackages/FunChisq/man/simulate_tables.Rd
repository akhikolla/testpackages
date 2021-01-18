\name{simulate_tables}
\alias{simulate_tables}
\title{
Simulate Noisy Contingency Tables to Represent Diverse Discrete Patterns
}

\description{
Generate random contingency tables representing various functional, non-functional, dependent, or independent patterns, without specifying a parametric model for the patterns.
}

\usage{
simulate_tables(
  n = 100, nrow = 3, ncol = 3,
  type = c("functional", "many.to.one",
           "discontinuous", "independent",
           "dependent.non.functional"),
  n.tables = 1,
  row.marginal = NULL,
  col.marginal = NULL,
  noise = 0.0, noise.model = c("house", "candle"),
  margin = 0
)
}

\arguments{

   \item{n}{
   a positive integer specifying the sample size to be distributed in each table. For \code{"functional"},  \code{"many.to.one"}, and \code{"discontinuous"} tables, \code{n} must be no less than \code{nrow}. For \code{"dependent.non.functional"} tables, \code{n} must be no less than \code{nrow}*\code{ncol}. For \code{"independent"} tables, \code{n} must be a positive integer.
  }

  \item{nrow}{
  a positive integer specifying the number of rows in each table. The value must be no less than 2. For \code{"many.to.one"} tables, \code{nrow} must be no less than 3.
  }

  \item{ncol}{
  a positive integer specifying the number of columns in output table. \code{ncol} must be no less than 2.
  }

  \item{type}{
  a character string to specify the type of pattern underlying the table. The options are \code{"functional"} (default),  \code{"many.to.one"}, \code{"discontinuous"}, \code{"independent"}, and \code{"dependent.non.functional"}. See Details.
%The options are \code{"functional"} (default, where y=f(x) but x may or may not be function of y),  \code{"many.to.one"} (where y=f(x) but x != f(y)), \code{"independent"} (x and y are statistically independent),  \code{"dependent.non.functional"} (x and y are statistically dependent but y!=f(x) and x!=f(y)).
  }
  \item{n.tables}{
  a positive integer value specifying the number of tables to be generated.
  }
  \item{row.marginal}{
  a non-negative numeric vector of length \code{nrow} specifying row marginal probabilities. The vector is linearly scaled so that the sum is 1. The default is a uniform distribution.
  %For \code{"functional table"}, \code{"dependent.non.functional"} and \code{"independent"} tables, \code{row.marginal} must be greater than or equal to 2.
%For \code{"many.to.one"} tables, the length of \code{row.marginal} vector must be no less than 3.
  }
  \item{col.marginal}{
  a non-negative numeric vector of length \code{ncol} specifying column marginal probabilities. The vector is linearly scaled so that the sum is 1.
This argument is ignored by \code{"dependent.non.functional"} tables.
  % The vector is used only in generating \code{"Functional"} and \code{"independent"} tables. The default is a uniform distribution.
  }

   \item{noise}{
  a numeric value between 0 and 1 specifying the noise level to be added to a table using function \code{\link{add.noise}}. The noise can be applied along row, column, or both, which can be specified by the \code{margin} argument. See \code{\link{add.noise}} for details.
  }

   \item{noise.model}{
  a character string indicating the noise model of either \code{"house"} for ordinal variables \insertCite{zhang2015chinet}{FunChisq} or \code{"candle"} for categorical variables. See \code{\link{add.noise}} for details.
  }

  \item{margin}{
  a numeric value of either 0, 1 or 2. Default is 0.
  0: noise is applied along both rows and columns.
  1: noise is applied along each row.
  2: noise is applied along each column.
  See \code{\link{add.noise}} for details.
  }
}

\details{
This function generates five types of table representing different interaction patterns between row and column discrete random variables \eqn{X}{X} and \eqn{Y}{Y}. Three of the five types are non-constant functional patterns (\eqn{Y}{Y} is a non-constant function of \eqn{X}{X}):

\code{type="functional"}: \eqn{Y}{Y} is a function of \eqn{X}{X} but \eqn{X}{X} may or may not be a function of \eqn{Y}{Y}. The samples are distributed using either the given row or column marginal probabilities. If both marginals are provided, randomly select one for each table. If neither is provided, randomly select one of them with uniform marginal for each table.

\code{type="many.to.one"}: \eqn{Y}{Y} is a many-to-one function of \eqn{X}{X} but \eqn{X}{X} is not a function of \eqn{Y}{Y}. The samples are distributed using either the given row or column marginal probabilities. If both marginals are provided, randomly select one for each table. If neither is provided, randomly select one of them with uniform marginal for each table.

\code{type="discontinuous"}: \eqn{Y}{Y} is a function of \eqn{X}{X}, where the function value of X must differ from its neighbors. \eqn{X}{X} may or may not be a function of \eqn{Y}{Y}. The samples are distributed using either the given row or column marginal probabilities. If both marginals are provided, randomly select one for each table. If neither is provided, randomly select one of them with uniform marginal for each table. A discontinuous function forms a contrast with those that are close to constant functions.

The fourth \code{type} \code{"dependent.non.functional"} is non-functional patterns where \eqn{X}{X} and \eqn{Y}{Y} are statistically dependent but not function of each other. The samples are distributed according to \code{row.marginal} probabilities.

The fifth \code{type} \code{"independent"} represents patterns where \eqn{X}{X} and \eqn{Y}{Y} are statistically independent whose joint probability mass function is the product of their marginal probability mass functions.

Random noise can be optionally applied to the tables using either the house or the candle noise model. See \code{\link{add.noise}} for details.
%(Zhang et al., 2015) by specifying a non-zero value (\eqn{\le 1}{<= 1}) for argument \code{noise}.

\insertCite{sharma2017simulating;textual}{FunChisq} provide full mathematical and statistical details of the simulation strategies for the above table types except the \code{"discontinuous"} type which was introduced after the publication.
}

\value{
A list containing the following components:

\item{pattern.list}{a list of tables containing binary patterns in 0's and 1's. Each table is created by setting all non-zero entries in the corresponding sampled contingency table from \code{sample.list} to 1. Each table strictly satisfies the  mathematical relationship required for a given pattern \code{type} requested, but it does not meet the statistical requirements. As each table represents the truth regarding the mathematical relationship between the row and column variables, they can be used as the ground truth or gold standard for benchmarking.}

\item{sample.list}{a list of tables satisfying both the mathematical and statistical requirements. These tables are noise free.}

\item{noise.list}{a list of tables after applying noise to the corresponding tables in \code{sample.list}. Each table is the noisy version of the corresponding sampled contingency table. Due to the added noise, each table may no longer strictly satisfy the required mathematical or statistical relationships. These tables are the main output to be used for the evaluation of a discrete pattern discovery algorithm.}

\item{pvalue.list}{a list of p-values reporting the statistical significance of the generated tables for the required type. When the pattern type specifies a functional relationship, the p-values are computed by the functional chi-square test \insertCite{zhang2013deciphering}{FunChisq}; otherwise, the Pearson's chi-square test of independence is used to calculate the p-value.}
}

\references{
\insertAllCited{}
}

\author{
Ruby Sharma, Sajal Kumar, Hua Zhong, and Joe Song
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{add.noise}} for details of the noise model.
}

\examples{
\donttest{
# In all examples, x is the row variable and y is the column
#    variable of a table.

# Example 1. Simulating a noisy function where y=f(x),
#            x may or may not be g(y) with given row.marginal.

tbls <- simulate_tables(n=100, nrow=4, ncol=5, type="functional",
                noise=0.2, n.tables = 1,
                row.marginal = c(0.3,0.2,0.3,0.2))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 1. Functional pattern")
plot_table(tbls$sample.list[[1]], main="Ex 1. Sampled pattern (noise free)")
plot_table(tbls$noise.list[[1]], main="Ex 1. Sampled pattern with 0.2 noise")
plot.new()

# Example 2. Simulating a noisy functional pattern where
#            y=f(x), x may or may not be g(y) with given row.marginal.

tbls <- simulate_tables(n=100, nrow=4, ncol=5, type="functional",
                noise=0.5, n.tables = 1,
                row.marginal = c(0.3,0.2,0.3,0.2))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 2. Functioal pattern", col="seagreen2")
plot_table(tbls$sample.list[[1]], main="Ex 2. Sampled pattern (noise free)", col="seagreen2")
plot_table(tbls$noise.list[[1]], main="Ex 2. Sampled pattern with 0.5 noise", col="seagreen2")
plot.new()


# Example 3. Simulating a noisy many.to.one function where
#            y=f(x), x!=f(y) with given row.marginal.

tbls <- simulate_tables(n=100, nrow=4, ncol=5, type="many.to.one",
                noise=0.2, n.tables = 1,
                row.marginal = c(0.4,0.3,0.1,0.2))
par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 3. Many-to-one pattern", col="limegreen")
plot_table(tbls$sample.list[[1]], main="Ex 3. Sampled pattern (noise free)", col="limegreen")
plot_table(tbls$noise.list[[1]], main="Ex 3. Sampled pattern with 0.2 noise", col="limegreen")
plot.new()

# Example 4. Simulating noisy discontinuous
#   pattern where y=f(x), x may or may not be g(y) with given row.marginal.

tbls <- simulate_tables(n=100, nrow=4, ncol=5,
                type="discontinuous", noise=0.2,
                n.tables = 1, row.marginal = c(0.2,0.4,0.2,0.2))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 4. Discontinuous pattern", col="springgreen3")
plot_table(tbls$sample.list[[1]], main="Ex 4. Sampled pattern (noise free)", col="springgreen3")
plot_table(tbls$noise.list[[1]], main="Ex 4. Sampled pattern with 0.2 noise", col="springgreen3")
plot.new()


# Example 5. Simulating noisy dependent.non.functional
#            pattern where y!=f(x) and x and y are statistically
#            dependent.

tbls <- simulate_tables(n=100, nrow=4, ncol=5,
                type="dependent.non.functional", noise=0.3,
                n.tables = 1, row.marginal = c(0.2,0.4,0.2,0.2))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 5. Dependent.non.functional pattern",
col="sienna2", highlight="none")
plot_table(tbls$sample.list[[1]], main="Ex 5. Sampled pattern (noise free)",
col="sienna2", highlight="none")
plot_table(tbls$noise.list[[1]], main="Ex 5. Sampled pattern with 0.3 noise",
col="sienna2", highlight="none")
plot.new()

# Example 6. Simulating a pattern where x and y are
#            statistically independent.

tbls <- simulate_tables(n=100, nrow=4, ncol=5, type="independent",
                noise=0.3, n.tables = 1,
                row.marginal = c(0.4,0.3,0.1,0.2),
                col.marginal = c(0.1,0.2,0.4,0.2,0.1))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 6. Independent pattern",
col="cornflowerblue", highlight="none")
plot_table(tbls$sample.list[[1]], main="Ex 6. Sampled pattern (noise free)",
col="cornflowerblue", highlight="none")
plot_table(tbls$noise.list[[1]], main="Ex 6. Sampled pattern with 0.3 noise",
col="cornflowerblue", highlight="none")
plot.new()


# Example 7. Simulating a noisy function where y=f(x),
#            x may or may not be g(y), with given column marginal


tbls <- simulate_tables(n=100, nrow=4, ncol=5, type="functional",
                noise=0.2, n.tables = 1,
                col.marginal = c(0.2,0.1,0.4,0.2,0.1))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 7. Functional pattern")
plot_table(tbls$sample.list[[1]], main="Ex 7. Sampled pattern (noise free)")
plot_table(tbls$noise.list[[1]], main="Ex 7. Sampled pattern with 0.2 noise")
plot.new()


# Example 8. Simulating a noisy many.to.one function where
#            y=f(x), x!=f(y) with given column marginal.

tbls <- simulate_tables(n=100, nrow=4, ncol=4, type="many.to.one",
                noise=0.2, n.tables = 1,
                col.marginal = c(0.4,0.3,0.1,0.2))
par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 8. Many-to-one pattern", col="limegreen")
plot_table(tbls$sample.list[[1]], main="Ex 8. Sampled pattern (noise free)", col="limegreen")
plot_table(tbls$noise.list[[1]], main="Ex 8. Sampled pattern with 0.2 noise", col="limegreen")
plot.new()


# Example 9. Simulating noisy discontinuous
#   pattern where y=f(x), x may or may not be g(y) with given column marginal

tbls <- simulate_tables(n=100, nrow=4, ncol=4,
                type="discontinuous", noise=0.2,
                n.tables = 1, col.marginal = c(0.1,0.4,0.2,0.3))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 9. Discontinuous pattern", col="springgreen3")
plot_table(tbls$sample.list[[1]], main="Ex 9. Sampled pattern (noise free)", col="springgreen3")
plot_table(tbls$noise.list[[1]], main="Ex 9. Sampled pattern with 0.2 noise", col="springgreen3")
plot.new()




}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{ datagen }
