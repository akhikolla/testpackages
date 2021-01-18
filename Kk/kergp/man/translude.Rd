\name{translude}
\alias{translude}
\title{Make Translucent colors}
\usage{

  translude(colors, alpha = 0.6)

}
\arguments{
  \item{colors}{

    A vector of colors in a format that can be understood by
    \code{\link{col2rgb}}.

  }
  \item{alpha}{

    Level of opacity ("0" means fully transparent and "max" means
    opaque).  After recycling to reach the required length, this value
    or vector is used as \code{alpha} in \code{\link{rgb}}.

  }
}
\value{
  
  A vector of translucent (or semi-transparent) colors.

}
\description{

  Make translucent colors.

}
