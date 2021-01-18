\name{make.mpt.cf}
\alias{make.mpt.cf}
\alias{lbmpt.to.mpt}
\title{
Functions to transform MPT models.
}
\description{
(Helper) functions that takes an MPT model file and transforms it into a representation in the context-free language of MPT models L-BMPT (Purdy & Batchelder, 2009) or takes the representation in LBMPT and returns the model equations.
}

\usage{
make.mpt.cf(model.filename, restrictions.filename = NULL, 
            model.type = c("easy", "eqn"), treewise = FALSE)

lbmpt.to.mpt(model.list, outfile = NULL, category.names= TRUE)
}

\arguments{
  \item{model.filename}{
A character \code{vector} specifying the location and name of the model file. 
}
  \item{restrictions.filename}{
\code{NULL} or a character \code{vector} or a \code{list} of characters. The default is \code{NULL} which corresponds to no restrictions. A character \code{vector} specifies the location or name of the restrictions file. A \code{list} of characters contains the restrictions directly. Inequality/order restrictions are silently ignored.
}
  \item{model.type}{
Character vector specifying whether the model file is formatted in the \code{easy} format (\code{"easy"}; i.e., each line represents all branches corresponding to a response categories) or the traditional EQN syntax (\code{"eqn"} or \code{"eqn2"}). See Details in \code{\link{fit.mpt}}.
}
  \item{treewise}{
  logical. Should the model be concatenated to one tree before transforming to LBMPT? Default is \code{FALSE}.
}
  \item{model.list}{
  A list of character \code{vector}s representing a model in LBMPT. Each element of the vector corresponds to either a parameter or category. Each list element corresponds to one tree. Can be returned from \code{make.mpt.cf}.
}
  \item{outfile}{
  Name of the file the model equation should be saved in (in easy format). If \code{NULL} (the default) prints it to the console instead (\code{stdout()}).
}
  \item{category.names}{
  logical. Should category names (e.g., "category 1") be printed at the end of each line?
}
}
\details{
Purdy and Batchelder (2009) provide a new way of how binary multinomial processing tree (MPT) models can be represented, a context free language called L-BMPT. This function takes a model file that consists of the equations defining a model and returns a character vector representing this model in L-BMPT.

There are three important things to know about this function:

1. L-BMPT distinguishes between observable categories (C) and parameters (theta). As MPTinR allows parameters to have any name that is legal for a variable in R (with the only restriction that parameters should NOT start with \code{hank}) the L-BMPT representation of the model's parameters can also consist of any name that is a legal variable name in R. To distinguish parameters from categories, categories are represented as integers (i.e., numbers) (a number is not a legal variable name in R). Furthermore, as legal variable names may end with a number, concatenating parameters and categories into one string could lead to an ambiguous representation of the model. Therefore, the returned representation in L-BMPT is a character vector with each element representing either a parameter (any legal variable name in R) or a category (an integer).\cr
Theta: Names that are legal variable names in R.\cr
C: Integers.

2. If a model consists of more than n > 1 trees, this function per default concatenates the trees into a single binary MPT model by adding n-1 parameters (named \code{hank.join.x} with \code{x} be a integer starting at 1; see Wu, Myung & Batchelder, 2010). This can be turned off, by setting \code{treewise} to TRUE.

3. It is absolutely necessary that the representation of the model via equations in the model file exactly maps on the structure of the binary tree. In other words, equations in the model file can NOT be simplified in any way. The equations in the model file may only consist of the parameters and their negations (e.g. \code{d} and \code{(1-d)}). Simplifications and aggregations in the model file (e.g., from \code{u * (1-u) + u * (1-u)} to \code{2*u*(1-u}) will lead to erroneous results! Similarly, reparameterizations for inequality constraints (which can be done by \code{fit.mpt}) can not be represented in L-BMPT.

}
\value{
\code{make.mpt.cf}: A character vector with each element representing either a parameter or a category (categories are represented by integer numbers). In case of multiple trees and \code{treewise = TRUE}, a list of such vectors.

\code{lbmpt.to.mpt}: Either prints the model to the screen or returns nothing and saves the model equations in the specified file.
}
\references{
Purdy, B. P., & Batchelder, W. H. (2009). A context-free language for binary multinomial processing tree models. \emph{Journal of Mathematical Psychology}, 53, 547-561.

Riefer, D. M., & Batchelder, W. H. (1988). Multinomial modeling and the measurement of cognitive processes. \emph{Psychological Review}, 95, 318-339.

Wu, H., Myung, J., I., & Batchelder, William, H. (2010). Minimum description length model selection of multinomial processing tree models. \emph{Psychonomic Bulletin & Review}, 17, 275-286.
}
\author{
Henrik Singmann (make.mpt.cf)

Quentin Gronau and Franz Dietrich (lbmpt.to.mpt, using a function from Akhil S Bhel, LinearizeNestedList)
}
\note{
It is absolutely necessary that the model file exactly maps on the structure of the binary tree. See Details.
}

\seealso{
\code{\link{get.mpt.fia}} and \code{link{bmpt.fia}} for functions calling \code{make.mpt.cf} to obtain the FIA of a MPT model.

See \code{\link{fit.mpt}} for the main function of MPTinR which also calls \code{make.mpt.cf} for obtaining the FIA.

\code{\link{prepare.mpt.fia}} will provide the code needed for obtaining the Minimum Description Length of a MPT model using Matlab (Wu, Myung & Batchelder, 2010) and calls \code{make.mpt.cf} for obtaining the L-BMPT representation.\cr
}
\examples{
model2 <- system.file("extdata", "rb.fig2.model", package = "MPTinR")

make.mpt.cf(model2)

make.mpt.cf(model2, treewise = TRUE)
  
lbmpt.to.mpt(make.mpt.cf(model2, treewise = TRUE))
}

