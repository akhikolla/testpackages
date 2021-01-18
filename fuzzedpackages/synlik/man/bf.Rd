\name{bf}
\alias{bf}
\alias{bf1}
\alias{bf2}
\alias{bf3}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Nicholson's 1954 blowfly data}
\description{
Data from figures 3 and 4 of Nicholson, 1954. 
} %- end description

\usage{
data(bf1)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{bf1}{the dataset name}

}

\details{ \code{bf1} is Nisbet and Gurney's run 1, and Nicholson's (1954) figure 3 (adult food limitation). The 
data are actually from the global population dynamics database at Silwood. They are daily: Nicholson's figure 3 
plots data every other day, but the text says that measurements were taken daily. However elsewhere they are 
reported every other day. Probably best to assume that they have been interpolated to daily.

\code{bf2} and \code{bf3} are digitized from Nicholson's (1954) figure 4. \code{bf2} is the upper series: 
larval food limitation, with 50g per day of larval food provided. \code{bf3} is the lower series: same set up, half as much food. 
These are Nisbet and Gurney's runs 2 and 3, respectively.
}

\value{matrix of replicate data series
}
\references{
Alexander J Nicholson. An outline of the dynamics of animal populations. Australian
Journal of Zoology, 2(1):9--65, 1954.
}
\author{ Simon N. Wood, maintainer Matteo Fasiolo <matteo.fasiolo@@gmail.com>}


\seealso{  
blowfly
}

\examples{
  library(synlik)
  data(bf1)
  data(bf2)
  data(bf3)
  par(mfrow=c(3,1),mar=c(4,4,1,1))
  with(bf1,plot(day,pop,type="l"))
  with(bf1,points(day,pop,pch=20,cex=.8))
  abline(mean(bf1$pop),0,col=2); abline(median(bf1$pop),0,col=3);
  with(bf2,plot(day,pop,type="l"))
  with(bf2,points(day,pop,pch=20,cex=.8))
  abline(mean(bf2$pop),0,col=2); abline(median(bf2$pop),0,col=3);
  with(bf3,plot(day,pop,type="l"))
  with(bf3,points(day,pop,pch=20,cex=.8))
  abline(mean(bf3$pop),0,col=2); abline(median(bf3$pop),0,col=3);
}







