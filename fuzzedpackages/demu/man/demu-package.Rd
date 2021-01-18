\name{demu-package}
\alias{demu-package}
\docType{package}
\title{
  \code{demu} is an open-source R package implementing a Gaussian process optimal design emulator based on Determinantal point processes.
}
\description{
  \code{demu} implements a determinantal point process emulator for probabilistically sampling optimal designs for Gaussian process (GP) regression models. Currently, \code{demu} is a proof of concept implementation that implements basic DPP sampling, conditional DPP sampling for drawing designs of fixed size \code{n}, sequential DPP sampling to build designs iteratively and a faster C++ implementation of the conditional DPP sampler using sparse matrices.  The package supports popular stationary correlation functions commonly used in GP regression models, including the Gaussian and Wendland correlation functions.
}
\details{
  The main model fitting functions in the package include \code{sim.dpp.modal()} for dense correlation matrices and \code{sim.dpp.modal.fast()} for sparse correlation matrices. These functions use a grid-based approximation to sample from the relevant DPP model.
}
\author{
Matthew T. Pratola <mpratola@stat.osu.edu> [aut, cre, cph]
}
\references{
  Pratola, Matthew T., Lin, C. Devon, and Craigmile, Peter. (2018)
  Optimal Design Emulators: A Point Process Approach. 
  \emph{arXiv:1804.02089}.

}
\keyword{ package }
\seealso{
  \code{\link{sim.dpp.modal},\link{sim.dpp.modal.fast},\link{sim.dpp.modal.seq},\link{sim.dpp.modal.fast.seq}}
}
