\name{prepare.mpt.fia}
\alias{prepare.mpt.fia}
\title{
Provides MATLAB command to get FIA
}
\description{
This function needs data and a model files and outputs the exact command needed to obtain the minimum description length measure for MPT models using the procedure by Wu, Myung, and Batchelder (2010) for MATLAB. It can be considered an extended wrapper for \code{\link{make.mpt.cf}}.
}
\usage{
prepare.mpt.fia(data, model.filename, restrictions.filename = NULL,
	outfile = "clipboard", Sample = 2e+05, model.type = c("easy", "eqn", "eqn2"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Either a \emph{numeric} \code{vector} for individual fit or a \emph{numeric} \code{matrix} or \code{data.frame} for multi-individual fit. The data on each position (column for multi individual fit) must correpsond to the relevant line in the model file.
}
  \item{model.filename}{
A character \code{vector} specifying the location and name of the model file.
}
  \item{restrictions.filename}{
\code{NULL} or a character \code{vector} specifying the location and name of the restrictions file. Default is \code{NULL} which corresponds to no restrictions.
}
  \item{outfile}{
A character vector specifying the name of the file where the MATLAB code is saved. Default is \code{"clipboard"} which will copy the output to the clipboard and will not write it to a file (Windows only). Actually, this parameter is directly passed to \code{\link{writeLines}} which interprets character vectors as filenames, so any other legal \code{\link{connection}} can be used.
}
  \item{Sample}{
Number of Monte Carlo samples to be used by the procedure of Wu, Myung, and Batchelder (2010). Default is 200.000.
}
  \item{model.type}{
Character vector specifying whether the model file is formatted in the \code{easy} format (\code{"easy"}; i.e., each line represents all branches corresponding to a response categories) or the traditional EQN syntax (\code{"eqn"} or \code{"eqn2"}). See Details in \code{\link{fit.mpt}}.
}
}
\details{
This function uses \code{\link{make.mpt.cf}} to create the representation in the L-BMPT. Therefore, it is necessary that the representation of the model via equations in the model file exactly maps on the structure of the binary tree (see \code{\link{make.mpt.cf}} for more details).

Whereas \code{\link{fit.mpt}} can reparameterize MPT models for fitting inequality constraints, Wu, Myung, and Batchelder (2010) have used another method to deal with these issues that is also adopted here. Our function does not report a reparameterized version of the MPT model that satisfies the inequality contraints, but modifies the appropriate argument in the call to the function by Wu et al (2010).

Note that MATLAB needs the statistics toolbox to run the script by Wu, Myung, and Batchelder (2010).
}
\value{
The most important value is the output to a file or clipboard (Windows only) of the MATLAB code to get the minimum description length. For multiple individuals multiple outputs are generated which only differ if the ns of the data differ. Furthermore, each argument is returned in a list:
\item{s}{The string representation of the model.}
\item{parameters}{A list of the numbers representing the parameters.}
\item{param.codes}{A vector describing which number corresponds to which parameter in the parameters vector.}
\item{category}{The numbers representing the categories.}
\item{ineq}{The matrix representing the inequality constraints.}
\item{n}{The n of the data.}
\item{internal}{The L-BMPT representation as returned by \code{make.mpt.cf}.}
}
\references{
Wu, H., Myung, J., I., & Batchelder, William, H. (2010). Minimum description length model selection of multinomial processing tree models. \emph{Psychonomic Bulletin & Review}, 17, 275-286.
}
\author{
Henrik Singmann
}

\seealso{
Since we ported the original BMPTFIA function by Wu, Myung, & Batchelder (2010) to R (\code{\link{bmpt.fia}}), this function is a little bit outdated. However, getting the FIA in Matlab is (still) faster than getting it in R.

See also \code{\link{get.mpt.fia}} which takes the same arguments but will then compute the FIA using the function provided by Wu et al. (2010) ported to R.

\code{\link{make.mpt.cf}}
}

\examples{
\dontrun{
# This example produces the code for the first example of how to use the 
# function by Wu, Myung & Batchelder (2010, pp. 280):
# Value should be around 12.61 and 12.62

model.1htm <- system.file("extdata", "wmb.ex1.model", package = "MPTinR")
model.1htm.restr <- system.file("extdata", "wmb.ex1.restr", package = "MPTinR")

prepare.mpt.fia(c(250,0,0,250,0,0,500,0,0), model.1htm, model.1htm.restr)
}
}

