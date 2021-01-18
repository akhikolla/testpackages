\name{shapeSlot}
\alias{shapeSlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extracts the Slots of a Structure
}
\description{

  Extract the slot of a structure.

}
\usage{
shapeSlot(object, slotName = "par", type = "all", as = "vector")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{

    An object to extract from, typically a covariance kernel.

  }
  \item{slotName}{

    Name of the slot to be extracted.

  }
  \item{type}{

    Type of slot to be extracted. Can be either a type of parameter,
    \code{"var"} or \code{"all"}.

  }
  \item{as}{

    Type of result wanted. Can be \code{"vector"}, \code{"list"}
    or \code{"matrix"}.

  }
}
%\details{}
\value{

  A vector, list or matrix containing the extraction.

}
%% \references{}
%%\author{Yves}
\note{

  This function is for internal use only.

}
%% \seealso{}
%% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
