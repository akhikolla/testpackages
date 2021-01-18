\name{plot.cumhaz.l1mstate}
\alias{plot.cumhaz.l1mstate}
\title{
 Plot the estimated cumulative hazard rates of the multi-state model. 
}
\description{
 Plot the estimated cumulative hazard rates of the multi-state model.
}
\usage{
\method{plot}{cumhaz.l1mstate}(x,type=c("single","separate"),cols,
     xlab="Years since transplant",ylab="Cumulative hazard",
     ylim,lwd=3,lty,legend,legend.pos,bty="o",...)
}
\arguments{
  \item{x}{fitted cumhaz.l1mstate object}
  \item{type}{One of "single" or "separate"; in case of "single", all estimated cumulative hazards are drawn in a single plot, in case of "separate", separate plots are shown for the estimated cumulative hazards.}
  \item{cols}{A vector specifying colors for the different transitions}
  \item{xlab}{A title for the x-axis; default is "Years since transplant"}
  \item{ylab}{A title for the y-axis; default is "Cumulative hazard"}
  \item{ylim}{The y limits of the plot}
  \item{lwd}{The line width; default is 3}
  \item{lty}{The line type}
  \item{legend}{The transition numbers; if missing, these will be taken from the transition matrix contained in cumhaz.l1mstate object.}
  \item{legend.pos}{The position of the legend; default is "topleft"}
  \item{bty}{The box type of the legend}
  \item{\dots}{Further arguments to plot}
}
\value{No return value}
\details{
A plot is produced, and nothing is returned.
}

\author{
Xuan Dang \email{xuandang11289@gmail.com}
}

\examples{
#### EBMT data
data(ebmt4)

x1 <- ifelse(ebmt4$match=="gender mismatch", 1, 0)
x2 <- ifelse(ebmt4$proph=="yes", 1, 0)
x3 <- ifelse(ebmt4$year=="1990-1994", 1, 0)
x4 <- ifelse(ebmt4$year=="1995-1998", 1, 0)
x5 <- ifelse(ebmt4$agecl=="20-40", 1, 0)
x6 <- ifelse(ebmt4$agecl==">40", 1, 0)
dt <- data.frame(ebmt4[,1:11],x1,x2,x3,x4,x5,x6)
col_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
tmat <- transMat(x = list(c(2,3,5,6),c(4,5,6),c(4,5,6),c(5,6),c(),c()),  
                 names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
dat <- msprep(data = dt, trans = tmat, time = c(NA, "rec", "ae", 
                "recae","rel", "srv"),status = c(NA, "rec.s", "ae.s", 
                "recae.s","rel.s", "srv.s"), keep = col_names)
# select the first 450 lines for a sake of illustration
dat <- dat[1:450,]
out <- cv.l1mstateR(dat, nlambda = 100, nfolds = 10)

### Compute the cumulative hazard and the transition probabilities
## patient A: include the highest number of patients (287)
ptsA <- dat[which(dat$x1==0 & dat$x2==0 & dat$x3==1 & dat$x4==0 & dat$x5==1 & dat$x6==0),]
## observed transitions (ground truth)
# predicted time = 0
events(ptsA)
# use models to predict the transition probabilities
# L1MSTATE
ptA <- ptsA[which(ptsA$id == unique(ptsA$id)[1]),]
ptA <- ptA[which(ptA$trans \%in\% c(1:4)),c(4,9:14)]

cumhazA <- cumhaz.l1mstate(object = out, longdt = dat, newdata = ptA, cv.method = "pcvl")
plot.cumhaz.l1mstate(cumhazA, main = "L1MSTATE", type = "single", lwd=2, las=1, 
ylim = c(0,1))
legend("topleft",inset=.01,legend=c("1","2","3","4"),col=c("black","red","green",
       "blue"),lty=1,cex=0.8,title="Transitions",text.font=2,bg='white')
}

\keyword{ multi-state models}% use one of  RShowDoc("KEYWORDS")
\keyword{ regularization }% __ONLY ONE__ keyword per line
