\name{plot.simulate.gp}
\alias{plot.simulate.gp}
\title{

  Plot Simulations from a \code{gp} Object

}
\usage{
\method{plot}{simulate.gp}(x, y,
        col = list(sim = "SpringGreen3", trend = "orangered"),
        show = c(sim = TRUE, trend = TRUE, y = TRUE),
        ...)
}
\arguments{
  
  \item{x}{

    An object containing simulations, produced by 'simulate' with
    \code{output = "list"}.

  }
  \item{y}{

    Not used yet.

  }
  \item{col}{

    Named list of colors to be used, with elements \code{"sim"} and
    \code{"trend"}.

  }
  \item{show}{

    A logical vector telling which elements must be shown.

  }
  \item{...}{

    Further argument passed to \code{plot}.

  }
  
}
\value{

  Nothing.

}
\description{

  Function to plot simulations from a \code{gp} object.

}
\note{

  For now, this function can be used only when the number of inputs is
one.

}
\seealso{
  
\code{\link{simulate.gp}}.

}
