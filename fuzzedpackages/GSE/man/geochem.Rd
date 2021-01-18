\name{geochem}
\alias{geochem}
\title{Geochemical Data}
\usage{ data(geochem) }
\description{Geochemical data  analyzed by Smith et al (1984). The variables in the data
measures the contents (in parts per million) for 20 chemical elements (e.g., Copper and Zinc)  in 53 samples of rocks in Western Australia.
}
\format{The data contains 53 observations on 20 variables corresponding to the 20 chemical elements.}
\references{  
  Smith, R.E., Campbell, N.A., Licheld, A. (1984).
  Multivariate statistical techniques applied to pisolitic laterite geochemistry at Golden Grove, 
  Western Australia. 
  \emph{Journal of Geochemical Exploration}, \bold{22},
  193--216.

  Agostinelli, C., Leung, A. , Yohai, V.J., and Zamar, R.H. (2014) Robust estimation of 
  multivariate location and scatter in the presence of cellwise and casewise contamination. 
  arXiv:1406.6031[math.ST]
}  
\keyword{datasets}
\examples{
\dontrun{
library(ICSNP)
library(rrcov)

data(geochem)
n <- nrow(geochem)
p <- ncol(geochem)

# MLE
res.ML <- list(mu=colMeans(geochem), S=cov(geochem))

# Tyler's M
geochem.med <- apply(geochem,2,median,na.rm=TRUE)
res.Tyler <- tyler.shape(geochem, location=geochem.med)
res.Tyler <- res.Tyler*(median(mahalanobis( geochem, geochem.med, res.Tyler))/qchisq(0.5, df=p) )
res.Tyler <- list(mu=geochem.med, S=res.Tyler)

# Rocke's Covariace 
res.Rock <- CovSest(geochem, method="rocke")
res.Rock <- list(mu=res.Rock@center, S=res.Rock@cov)

# Fast-MCD
res.FMCD <- CovMcd( geochem)
res.FMCD <- list(mu=res.FMCD@center, S=res.FMCD@cov)

# MVE
res.MVE <- CovMve( geochem)
res.MVE <- list(mu=res.MVE@center, S=res.MVE@cov)

# S-estimator with bisquare rho function
res.S <- CovSest(geochem, method="bisquare")
res.S <- list(mu=res.S@center, S=res.S@cov)

# Fast-S
res.FS <- CovSest(geochem)
res.FS <- list(mu=res.FS@center, S=res.FS@cov)

# 2SGS
res.2SGS <- TSGS( geochem, seed=999 )
res.2SGS <- list(mu=res.2SGS@mu, S=res.2SGS@S)

# Combine all the results
geochem.res <- list(ML=res.ML, Tyler=res.Tyler, Rocke=res.Rock, MCD=res.FMCD, 
	MVE=res.MVE, FS=res.FS, MVES=res.S, TSGS=res.2SGS)

## Compare LRT distances between different estimators
res.tab <- data.frame( 	LRT.to.2SGS=c(slrt( res.ML$S, res.2SGS$S),
			slrt( res.Tyler$S, res.2SGS$S), 
			slrt( res.Rock$S, res.2SGS$S),
			slrt( res.FMCD$S, res.2SGS$S),
			slrt( res.MVE$S, res.2SGS$S),
			slrt( res.FS$S, res.2SGS$S),
			slrt( res.S$S, res.2SGS$S),
			slrt( res.2SGS$S, res.2SGS$S) ))
row.names(res.tab) <- c("ML","Tyler","Rocke","MCD","MVE","FS","MVES","TSGS")

# Calculate proportion of outliers cellwise
pairwise.mahalanobis <- function(x, mu, S){
	# function that computes pairwise mahalanobis distances
	p <- ncol(x)
	pairs.md <- c()
	for(i in 1:(p-1)) for(j in (i+1):p)
		pairs.md <- c(pairs.md, mahalanobis( x[,c(i,j)], mu[c(i,j)], S[c(i,j),c(i,j)]))
	pairs.md
}
res.tab$Full <- res.tab$Pairs <- res.tab$Cell <- NA
for(i in names(geochem.res) ){
	## Identify cellwise outliers
	uni.dist <- sweep(sweep(geochem, 2, geochem.res[[i]]$mu, "-"), 2, 
		sqrt(diag(geochem.res[[i]]$S)), "/")^2
	uni.dist.stat <- mean(uni.dist > qchisq(.99^(1/(n*p)), 1))
	res.tab$Cell[ which( row.names(res.tab) == i)] <- round(uni.dist.stat,3)

	## Identify pairwise outliers
	pair.dist <- pairwise.mahalanobis( geochem, geochem.res[[i]]$mu, geochem.res[[i]]$S)
	pair.dist.stat <- mean(pair.dist > qchisq(0.99^(1/(n*choose(p,2))), 2))
	res.tab$Pairs[ which( row.names(res.tab) == i)] <- round(pair.dist.stat,3)

	## Identify any large global MD
	full.dist <- mahalanobis( geochem, geochem.res[[i]]$mu, geochem.res[[i]]$S)
	full.dist.stat <- mean(full.dist > qchisq(0.99^(1/n), p))
	res.tab$Full[ which( row.names(res.tab) == i)] <- round(full.dist.stat,3)
}
res.tab
}
}
