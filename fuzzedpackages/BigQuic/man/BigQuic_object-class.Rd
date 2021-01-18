\name{BigQuic_object-class}
\Rdversion{1.1}
\docType{class}
\alias{BigQuic_object-class}

\title{Class \code{"BigQuic_object"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
Reference Class that holds all the relevant results of the BigQuic computation.
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\examples{
showClass("BigQuic_object")
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{precision_matrices}:}{Object of class \code{list} ~~ }
    \item{\code{X}:}{Object of class \code{matrix} ~~ }
    \item{\code{inputFileName}:}{Object of class \code{character} ~~ }
    \item{\code{lambda}:}{Object of class \code{numeric} ~~ }
    \item{\code{numthreads}:}{Object of class \code{numeric} ~~ }
    \item{\code{maxit}:}{Object of class \code{numeric} ~~ }
    \item{\code{epsilon}:}{Object of class \code{numeric} ~~ }
    \item{\code{k}:}{Object of class \code{numeric} ~~ }
    \item{\code{memory_size}:}{Object of class \code{numeric} ~~ }
    \item{\code{verbose}:}{Object of class \code{numeric} ~~ }
    \item{\code{isnormalized}:}{Object of class \code{numeric} ~~ }
    \item{\code{seed}:}{Object of class \code{numeric} ~~ }
    \item{\code{use_ram}:}{Object of class \code{logical} ~~ }
    \item{\code{clean}:}{Object of class \code{logical} ~~ }
    \item{\code{output_file_names}:}{Object of class \code{character} ~~ }
    \item{\code{opt.lambda}:}{Object of class \code{numeric} ~~ }
    \item{\code{inFlag}:}{Object of class \code{logical} ~~ }
    \item{\code{outFlag}:}{Object of class \code{logical} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{\code{cleanFiles(verbose)}:}{ ~~ }
    \item{\code{setOptLambda(optLambda)}:}{ ~~ }
    \item{\code{setX(inputX)}:}{ ~~ }
    \item{\code{setSeed(inputSeed)}:}{ ~~ }
  }
}
