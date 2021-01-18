\name{probs.l1mstate}
\alias{probs.l1mstate}
\title{
Compute subject-specific or overall transition probabilities
}
\description{
Compute subject-specific or overall transition probabilities.
}
\usage{
probs.l1mstate(object, longdt, tmat, predt, direction=c("forward","fixedhorizon"))
}
\arguments{
  \item{object}{fitted cumhaz.l1mstate object}
  \item{longdt}{long-format data input}
  \item{tmat}{a transition matrix for multi-state model}
  \item{predt}{a prediction time}
  \item{direction}{indicates whether prediction is forward or for a fixed horizon}
}

\value{
   A list with each element [[s]] containing a data frame with the estimated transition probabilities from state s. It also includes transition and tmat information for plotting purpose
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
}

\keyword{ multi-state models}% use one of  RShowDoc("KEYWORDS")
\keyword{ regularization }% __ONLY ONE__ keyword per line
