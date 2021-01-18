\name{dataf.medflies}
\alias{dataf.medflies}
\alias{medflies}
\docType{data}
\title{
Relationship of Age Patterns of Fecundity to Mortality for Female Medflies.
}
\description{
The data set consists of number of eggs laid daily for each of 1000 medflies (Mediterranean fruit flies, Ceratitis capitata) until time of death. Data were obtained in Dr. Carey's laboratory. The main questions are to explore the relationship of age patterns of fecundity to mortality, longevity and lifetime reproduction.
 
A basic finding was that individual mortality is associated with the time-dynamics of the egg-laying trajectory. An approximate parametric model of the egg laying process was developed and used in Muller et al. (2001). Non-parametric approaches which extend principal component analysis for curve data to the situation when covariates are present have been developed and discussed in  Chiou, Muller and Wang (2003) and Chiou et al. (2003).
}
\usage{
dataf.medflies()
}
\format{
  The functional data as a data structure.
  \describe{
    \item{\code{dataf}}{
      The functional data as a list of objects. Each object is characterized by two coordinates.
      \describe{
        \item{\code{args}}{\bold{day} - a numeric vector of the days numbers}
        \item{\code{vals}}{\bold{#eggs} - a numeric vector of numbers of eggs laid daily}
       }
    }
    \item{\code{labels}}{The classes of the objects: long-lived, short-lived}
  }
}
\source{
Carey, J.R., Liedo, P., Muller, H.G., Wang, J.L., Chiou, J.M. (1998). Relationship of age patterns of fecundity to mortality, longevity, and lifetime reproduction in a large cohort of Mediterranean fruit fly females. J. of Gerontology --Biological Sciences 53, 245-251. 

Muller, H.G., Carey, J.R., Wu, D., Liedo, P., Vaupel, J.W. (2001). Reproductive potential predicts longevity of female Mediterranean fruit flies. Proceedings of the Royal Society B 268, 445-450. 

Chiou, J.M., Muller, H.G., Wang, J.L. (2003). Functional quasi-likelihood regression models with smooth random effects. J. Royal Statist. Soc. B65, 405-423. 

Chiou, J.M., Muller, H.G., Wang, J.L., Carey, J.R. (2003). A functional multiplicative effects model for longitudinal data, with application to reproductive histories of female medflies. Statistica Sinica 13, 1119-1133. 
 
Chiou, J.M., Muller, H.G., Wang, J.L. (2004).Functional response models. Statistica Sinica 14,675-693.
 
}
\seealso{
\code{\link{dataf.*}} for other functional data sets

\code{\link{plot.functional}} for building plots of functional data
}
\examples{
## load the dataset
dataf = dataf.medflies()

## view the classes
unique(dataf$labels)

## access the 5th point of the 2nd object
dataf$dataf[[2]]$args[5]
dataf$dataf[[2]]$vals[5]

## plot the data
\dontrun{
labels = unlist(dataf$labels)
plot(dataf,
  xlab="Day", ylab="# eggs", 
  main=paste("Medflies (training time):\n short-lived red (", sum(labels == "short-lived"), "),", 
                    " long-lived blue (", sum(labels == "long-lived"), ")", sep=""),
  colors = c("blue", "red")   # in alphabetical order of class labels
  )

}
}
\keyword{datasets}
\keyword{functional}
