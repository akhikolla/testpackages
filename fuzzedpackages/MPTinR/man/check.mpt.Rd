\name{check.mpt}
\alias{check.mpt}

\title{
Check construction of MPT models.
}
\description{
A helper function which can aid in the process of constructing a MPT model file for MPTinR. It will check if the probabilities in each trees sum to 1 (if so, a tree is well constructed). If probabilities do not sum to 1, \code{check.mpt} will return for which trees. Furthermore, it will return the number of parameters and their names (helpful in spotting typos), the number of categories and the number of dfs the model provides. Finally, you can also pass restrictions as an argument and will receive the number and names of the parameters after restrictions are applied.
}
\usage{
check.mpt(model.filename, restrictions.filename = NULL, model.type = c("easy", "eqn"))
}
\arguments{
  \item{model.filename}{
A character \code{vector} specifying the location and name of the model file. 
}
  \item{restrictions.filename}{
\code{NULL} or a character \code{vector} specifying the location and name of the restrictions file. Default is \code{NULL} which corresponds to no restrictions.
}
  \item{model.type}{
Character vector specifying whether the model file is formatted in the \emph{easy} format (\code{"easy"}; i.e., each line represents all branches corresponding to a response category) or the traditional EQN syntax (\code{"eqn"}; see e.g., Stahl & Klauer, 2007). If the model filename ends with ".eqn" or ".EQN" the model is automatically treated as an EQN file.
}
}
\details{
As default, \code{check.mpt} expects a model file in the easy format, but if the filename ends with .eqn or .EQN \code{check.mpt} will expect the EQN format.\cr

In case of inequality restrictions, the original parameters which are inequality restricted are replaced  with dummy parameters starting with \code{hankX}. When using \code{\link{fit.mpt}} you will not notice this, as the output only shows the original parameters. In contrast, \code{check.mpt} removes the original parameters and shows the dummy parameters called \code{hankX}. Note that this does not change the number of parameters in the model.
}
\value{
A list with 
\item{probabilites.eq.1}{A logical value indicating whether or not the probabilities in each tree sum to 1. If \code{FALSE}, a warning is shown indicating in which trees the probabilities do not sum to 1.}
\item{n.trees}{Number of trees in the model.}
\item{n.model.categories}{Total number of categories expected in a dataset for that model.}
\item{n.independent.categories}{Number of independent response categories (i.e., independent data points) the model provides (i.e., \code{n.model.categories - n.trees}). The number of parameters can not be higher than this value for a model to be identifiable.}
\item{n.params}{Number of parameters in the model.}
\item{parameters}{Names of parameters in the model.}

If restrictions are present, the \code{n.params} and \code{parameters} are displayed for the unrestricted model (\code{orig.model}) as well as for the restricted model (\code{restr.model}).

}
\seealso{
\code{\link{fit.mpt}}
}
\examples{
# model of example 1 from example(fit.mpt)
model1 <- system.file("extdata", "rb.fig1.model", package = "MPTinR")
check.mpt(model1)

#model 1 in eqn format
model1.eqn <- system.file("extdata", "rb.fig1.model.eqn", package = "MPTinR")
check.mpt(model1.eqn)



#models of example 2 from example(fit.mpt)
model2 <- system.file("extdata", "rb.fig2.model", package = "MPTinR")
check.mpt(model2)

model2r.r.eq <- system.file("extdata", "rb.fig2.r.equal", package = "MPTinR")
check.mpt(model2, model2r.r.eq)

model2r.c.eq <- system.file("extdata", "rb.fig2.c.equal", package = "MPTinR")
check.mpt(model2, model2r.c.eq)

}


\keyword{models}
\keyword{tree}
