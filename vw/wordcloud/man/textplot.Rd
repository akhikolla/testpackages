\name{textplot}
\alias{textplot}
\title{
  Text Plot
}
\description{
  An x y plot of non-overlapping text
}
\usage{
textplot(x, y, words, cex=1,new=TRUE, show.lines=TRUE, ...)
}
\arguments{
  \item{x}{x coordinates}
  \item{y}{y coordinates}
  \item{words}{the text to plot}
  \item{cex}{font size}
  \item{new}{should a new plot be created}
  \item{show.lines}{if true, then lines are plotted between x,y and the word, for those words not covering their x,y coordinates}
  \item{...}{Additional parameters to be passed to wordlayout and text.}
}
\value{
nothing
 }
\examples{

#calculate standardized MDS coordinates
dat <- sweep(USArrests,2,colMeans(USArrests))
dat <- sweep(dat,2,sqrt(diag(var(dat))),"/")
loc <- cmdscale(dist(dat))

#plot with no overlap
textplot(loc[,1],loc[,2],rownames(loc))

#scale by urban population size
textplot(loc[,1],loc[,2],rownames(loc),cex=USArrests$UrbanPop/max(USArrests$UrbanPop))

#x limits sets x bounds of plot, and forces all words to be in bounds
textplot(loc[,1],loc[,2],rownames(loc),xlim=c(-3.5,3.5)) 

#compare to text (many states unreadable)
plot(loc[,1],loc[,2],type="n")
text(loc[,1],loc[,2],rownames(loc))

}