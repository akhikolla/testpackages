\name{SFS_read}
\alias{read}
\alias{SFS_read}
\title{Read similarity or dissimilarity input data}
\description{
  Read the similarity (or dissimilarity) information between the objects that one wants to order and build a 3-columns \code{data frame}, where each row \eqn{(i, j, A_{ij})} represents the (dis)similarity \eqn{A_{ij}} between objects \eqn{i} and \eqn{j}. In case of symmetric data (i.e., \eqn{A_{ij} = A_{ji}}), only the entries for pairs \eqn{(i,j)} with \eqn{i<j} are listed.
}
\usage{
read(data, zero_epsilon = 1e-200, symmetric = TRUE, identical_val = FALSE)
}

\arguments{
  \item{data}{a representation of the similarity (or dissimilarity) between pairs of objects.}
  \item{zero_epsilon}{a numerical value which determines that values in \code{data} below this threshold are considered to be \code{0}.}
  \item{symmetric}{a boolean value equal to \code{TRUE} if the input data is a symmetric matrix (i.e., \eqn{A_{ij} = A_{ji}} for all \eqn{i} and \eqn{j}).}
  \item{identical_val}{a boolean value equal to \code{TRUE} if the data is given as a 3-columns \code{data frame} and entries at both positions \eqn{(i,j)} and \eqn{(j,i)} are included.}
}
\value{
  Returns a 3-columns \code{data frame} representation of the original data listing all the pairwise (dis)similarities \eqn{(i, j, A_{ij})} between objects and selecting only the entries \eqn{A_{ij}} with \eqn{i<j} when the data is a symmetric matrix \eqn{A}.
}
\details{
  The input data can be a weighted adjacency matrix (represented by the objects: \code{matrix}, \code{dist} or \code{data frame}), or a list of all the weighted edges of a weighted graph (represented by a 3-col \code{data frame}) where each row \eqn{(i, j, A_{ij})} represents the (dis)similarity \eqn{A_{ij}} between objects \eqn{i} and \eqn{j} with \eqn{i<j}). If not specified, the data is assumed to be symmetric (i.e., same entry at positions \eqn{(i,j)} and \eqn{(j,i)}). Since by default the data is assumed to be symmetric, if it is represented by a 3-columns \code{data frame}, then it is assumed that symmetric pairs are not listed, and thus by default \code{identical_val = FALSE}. The reason for this choice is that for large symmetric data, it is more efficient to list the symmetric entries only once. 
 However, note that if \code{symmetric = FALSE} then \code{identical_val = TRUE} automatically.
 }  
 
\author{Matteo Seminaroti (SFS) and Utz-Uwe Haus (R wrapping)}

