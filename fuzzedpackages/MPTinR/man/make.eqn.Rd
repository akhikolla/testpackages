\name{make.eqn}
\alias{make.eqn}
\alias{make.mdt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates an EQN model file oir MDT data file
}
\description{
\code{make.eqn} takes a model file in the "easy" format and creates a model file in the EQN format. \code{make.mdt} takes a data vector and produces an .mdt data file.
}
\usage{
make.eqn(model.filename, eqn.filename)
make.mdt(data, mdt.filename, index, prefix = "dataset")
}
\arguments{
  \item{model.filename}{
A character \code{vector} specifying the location and name of the model file in the \code{easy} format.
}
  \item{eqn.filename}{
A character vector specifying the location and name of the target .eqn file.
}
  \item{data}{A vector, matrix or data.frame containing an individual data set to write to a .mdt file}
  \item{mdt.filename}{character vector specifying name and location of mdt file to be written.}
  \item{index}{index or second word written to the header of the mdt file. Ignored if \code{data} is a matrix or data.frame. }
  \item{prefix}{first word or prefix written to the mdt file. Default is \code{dataset}}
  
}
\details{
eqn and mdt files are the usual files used for programs to fit MPTs. You can use these functions to compare the results of MPTinR with other prgrams such as HMMTree or multiTree.

Note that these function do not add the endings .eqn or .mdt to the filename.

Since the MPTinR 0.9.4 \code{make.mdt} writes a single mdt file from a matrix or data.frame separating the participants via \code{===}.
}
\value{
Nothing
}
\references{
More information on the .eqn format in e.g.:\cr
Stahl, C., & Klauer, K. C. (2007). HMMTree: A computer program for latent-class hierarchical multinomial processing tree models. \emph{Behavior Research Methods}, 39, 267-273.
}
\author{
Henrik Singmann
}
\seealso{
\code{\link{fit.mpt}} for the main function of MPTinR
}

