
\name{wordlayout}
\alias{wordlayout}
\title{
  Word Layout
}
\description{
  finds text plot layout coordinates such that no text overlaps
}
\usage{
wordlayout(x, y, words, cex=1, rotate90 = FALSE,
		xlim=c(-Inf,Inf), ylim=c(-Inf,Inf), tstep=.1, rstep=.1, ...)
}
\arguments{
  \item{x}{x coordinates}
  \item{y}{y coordinates}
  \item{words}{the text to plot}
  \item{cex}{font size}
  \item{rotate90}{a value or vector indicating whether words should be rotated 90 degrees}
  \item{xlim}{x axis bounds for text}
  \item{ylim}{y axis bounds for text}
  \item{tstep}{the angle (theta) step size as the algorithm spirals out}
  \item{rstep}{the radius step size (in standard deviations) as the algorithm spirals out}
  \item{...}{Additional parameters to be passed to strwidth and strheight.}
}
\value{
A matrix with columns representing x, y width and height.
 }
\examples{

#calculate standardized MDS coordinates
dat <- sweep(USArrests,2,colMeans(USArrests))
dat <- sweep(dat,2,sqrt(diag(var(dat))),"/")
loc <- cmdscale(dist(dat))
x <- loc[,1]
y <- loc[,2]
w <- rownames(loc)

#plot with no overlap and all words visible
plot(x,y,type="n",xlim=c(-3,3),ylim=c(-3,2))
lay <- wordlayout(x,y,w,xlim=c(-3,3),ylim=c(-3,2))
text(lay[,1]+.5*lay[,3],lay[,2]+.5*lay[,4],w)

#notice north dakota is only partially visible
textplot(x,y,w)

}
