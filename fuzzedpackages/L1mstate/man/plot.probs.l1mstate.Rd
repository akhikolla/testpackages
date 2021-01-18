\name{plot.probs.l1mstate}
\alias{plot.probs.l1mstate}
\title{
Plot the transition probabilities
}
\description{
Plot the transition probabilities produced by probs.l1mstate
}
\usage{
\method{plot}{probs.l1mstate}(x,from,type=c("stacked","filled","single","separate"),
     ord,cols,xlab="Years since transplant",ylab="Probability",
     xlim,ylim,lwd,lty,cex,legend,legend.pos,bty="o",...)
}
\arguments{
  \item{x}{object produced from probs.l1mstate}
  \item{from}{the starting state from which the probabilities are used to plot}
  \item{type}{with "stacked" type, the transition probabilities are stacked and the distance between two adjacent curves indicates the probability; with "filled" type, it is the same but the space between adjacent curves are filled; with "single" type, the probabilities are shown as different curves in a single plot; with "separate" type, separate plots are shown for the estimated transition probabilities}
  \item{ord}{A vector indicates the order of plotting in case "stacked" or "filled"}
  \item{cols}{A vector specifying colors for the different transitions}
  \item{xlab}{A title for the x-axis; default is "Years since transplant"}
  \item{ylab}{A title for the y-axis; default is "Probability"}
  \item{xlim}{The x limits of the plot; default is range of time}
  \item{ylim}{The y limits of the plot}
  \item{lwd}{The line width; default is 3}
  \item{lty}{The line type}
  \item{cex}{Character size}
  \item{legend}{The transition numbers; if missing, these will be taken from the transition matrix contained in cumhaz.l1mstate object.}
  \item{legend.pos}{The position of the legend; default is "topleft"}
  \item{bty}{The box type of the legend}
  \item{\dots}{Further arguments to plot}
}
\value{No return value.}
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
probA_0 <- probs.l1mstate(cumhazA, longdt = longdt, tmat, predt = 0, direction = "forward")
statecols <- heat_hcl(6, c = c(90, 10), l = c(20, 80), power = c(1/5, 2))[c(6, 5, 3, 4, 2, 1)]
ord <- c(1, 2, 3, 4, 5, 6)
plot.probs.l1mstate(probA_0, main = "L1MSTATE", from = 1, ord = ord,
                    las=1, cex=0, type = "filled", col = statecols[ord])
}
\keyword{ multi-state models}% use one of  RShowDoc("KEYWORDS")
\keyword{ regularization }% __ONLY ONE__ keyword per line
