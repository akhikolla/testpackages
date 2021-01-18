\name{plot.functional}
\alias{plot.functional}
\alias{lines.functional}
\alias{points.functional}
\title{
Plot functions for the Functional Data
}
\description{
Plots the functional data given in the form which is described in the topic \code{\link{dataf.*}}.
}
\usage{
\method{plot}{functional}(x, 
      main = "Functional data", xlab = "args", ylab = "vals", 
      colors = c("red", "blue", "green", "black", "orange", "pink"), ...)
      
\method{lines}{functional}(x, 
      colors = c("red", "blue", "green", "black", "orange", "pink"), ...)
      
\method{points}{functional}(x, 
      colors = c("red", "blue", "green", "black", "orange", "pink"), ...)
}

\arguments{
  \item{x}{
The functional data as in the topic \code{\link{dataf.*}}. 
Note, that the in order to use s3 methods the data must be of class "functional". 
}
  \item{main}{
an overall title for the plot: see \code{\link{title}}
}
  \item{xlab}{
a title for the x axis: see \code{\link{title}}
}
  \item{ylab}{
a title for the y axis: see \code{\link{title}}
}
  \item{colors}{
the colors for the classes of the data. The colors are applied to the classes sorted in alphabetical order. Use the same set of classes to ensure that the same colours are selected in \code{lines} and \code{points} as in \code{plot} (do not remove entire classes). 
}
  \item{\dots}{
  additional parameters
}
}

\seealso{
\code{\link{dataf.*}} for functional data description
}
\examples{

\dontrun{
  ## load the Growth dataset
  dataf = dataf.growth()
  
  labels = unlist(dataf$labels)
  plot(dataf, 
       main = paste("Growth: girls red (", sum(labels == "girl"), "),", 
                    " boys blue (", sum(labels == "boy"), ")", sep=""),
       xlab="Year", ylab="Height, cm",
       colors = c("blue", "red")   # in alphabetical order of class labels   
  )
  
  # plot an observation as a line
  observation = structure(list(dataf = list(dataf$dataf[[1]])), class = "functional")
  lines(observation, colors = "green", lwd = 3)
  
  # plot hight at the age of 14 
  indexAge14 = which(observation$dataf[[1]]$args == 14)
  hightAge14 = observation$dataf[[1]]$vals[indexAge14]
  atAge14 = structure(list(
                      dataf = list(dataf = list(args = 14, vals = hightAge14))
                      ), class = "functional")
  points(atAge14, colors = "yellow", pch = 18)
}

}
\keyword{ visualization }
\keyword{ functional }
