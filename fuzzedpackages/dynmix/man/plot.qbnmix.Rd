\name{plot.qbnmix}
\alias{plot.qbnmix}

\title{Plots Selected Outcomes from \code{qbnmix} Object.}

\description{The function plots selected outcomes from \code{qbnmix} object.}

\usage{\method{plot}{qbnmix}(x, ...)}

\arguments{
\item{x}{an object of \code{qbnmix} class}
\item{...}{not used}
}

\details{
The function plots a few outcomes from \code{qbnmix} object. First, it plots means for each cluster. Then, it plots posterior probabilities for each cluster. Finally, estimates of mixing weights for each cluster.
}

\examples{
R <- list(matrix(c(1,0.3,0,
                   0.3,0.3,0,
                   0,0,0.15),3,3),
          matrix(c(1,0,0,
                   0,0.5,0,
                   0,0,0.2),3,3))
data <- rbind(MASS::mvrnorm(n=180,c(5,2,3),R[[1]]),
              MASS::mvrnorm(n=20,c(1,2,3),R[[2]]))
data <- data[sample(nrow(data)),]
                  
mu0 <- list(matrix(c(4.8689,1.9417,3.0175),nrow=1,ncol=3),
            matrix(c(1.0182,1.9903,2.8847),nrow=1,ncol=3))     
est <- qbnmix(y=data,mu0=mu0)
plot(est)
}
