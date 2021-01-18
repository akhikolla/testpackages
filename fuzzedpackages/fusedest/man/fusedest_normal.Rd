\name{fusedest_normal}
\alias{fusedest_normal}
\title{
The block splitting algorithm for linear regression estimation with the fused group lasso penalty function
}
\description{
A function for computing linear regression estimation with the fused group lasso penalty function
}
\value{
Return a list of output, e.g. the solution, runtime and iteration error, for the block splitting algorithm. For more details, please see the example below.

}
\examples{

library(fusedest)
library(igraph)

####### Functions for data generation #########


generating_normal_data <- function(beta.true, N, m, sigma2.y){

  p <- dim(beta.true)[2]
  M <- dim(beta.true)[1]

  label.list <- sample(c(1:M), m, replace = TRUE)
  n.list <- rpois(m, N)
  X <- matrix(rnorm(sum(n.list)*p, 0, 1), nrow = sum(n.list), ncol = p)

  ind.strt <- c(1, cumsum(n.list[1:(m-1)])+1)
  ind.end <- cumsum(n.list)
  label.dc <- rep(c(1:m), n.list)

  y <- as.numeric(unlist(sapply(c(1:m),
        function(i){
          X[ind.strt[i]:ind.end[i],]\%*\%as.numeric(beta.true[label.list[i],]) +
          rnorm(n.list[i], 0, sqrt(sigma2.y))
          })))

  label.true <- rep(label.list, n.list)

  results <- list(X, y, n.list, label.dc, label.true)
  names(results) <- c("X", "y", "n.list", "label.dc", "label.true")
  return(results)
}

generatingEdgelistID03 <- function(m, deg){

  c1 <- NULL
  c2 <- NULL

  if(deg < m-deg){

    c1 <- rep(0, (m-deg)*deg)
    c2 <- rep(0, (m-deg)*deg)

    for(i in 1:(m-deg)){

      ind.i <- c(((i-1)*deg+1):(i*deg))

      c1[ind.i] <- rep(i, deg)
      c2[ind.i] <- c((i+1):(i+deg))
    }

    if(deg > 1){
      c3 <- rep(0, deg*(deg-1)/2)
      c4 <- rep(0, deg*(deg-1)/2)
      l <- 0
      for(i in (m-deg+1):(m-1)){

        c3[c((l+1):(l+m-i))] <- rep(i, m-i)
        c4[c((l+1):(l+m-i))] <- c((i+1):m)
        l <- l + (m-i)
      }

    }
  }

  return(cbind(c(c1,c3),c(c2,c4)))
}



RcppInvGram <- function(X, w, lambda) {
    .Call('_fusedest_RcppInvGram', PACKAGE = 'fusedest', X, w, lambda)
}

RcppXtwy <- function(X, y, w) {
    .Call('_fusedest_RcppXtwy', PACKAGE = 'fusedest', X, y, w)
}

RcppWolsSolver03 <- function(invXtwX, Xtwy, b) {
    .Call('_fusedest_RcppWolsSolver03', PACKAGE = 'fusedest', invXtwX, Xtwy, b)
}



############ Setting true parameters ##########

p.star <- 10

beta.true <- t(matrix(
  c(rep(c(-2,2), p.star),
    rep(c(2,-2), p.star),
    c(rep(2, p.star),rep(-2,5)),
    c(rep(-2,p.star),rep(2,5)),
    rep(c(-1,3), p.star)), nrow = p.star, ncol = 5
))
N <- 100
m <- 10
p <- dim(beta.true)[2]

########## Generating data ###################

strt <- Sys.time()

mydata <- generating_normal_data(beta.true, N, m, sigma2.y = 1)

end <- Sys.time()
difftime(end, strt, units="sec")

y <- mydata$y
X <- mydata$X
label_dc <- mydata$label.dc
label.true <- mydata$label.true
n.list <- mydata$n.list

sum(n.list)
length(n.list)
length(y)
dim(X)
min(n.list)
max(n.list)
sum(n.list)

###### Run simulation #########################################

no.cores <- 1
m.total <- 10
m.list <- 10
ind.strt <- c(1, cumsum(n.list[1:(m.total-1)])+1)
ind.end <- cumsum(n.list)

no_lambda <- 1
lambda_list <- 0.01

u <- 1
H <- generatingEdgelistID03(m = m.list[u], deg = 2)
q_H <- sum(degree(graph_from_edgelist(H, directed = FALSE)))/2

max_iter <- 10
tol_err <- 10^(-100)
rho <- 1

set.seed(2, kind = NULL, normal.kind = NULL)

##### Computing initial values ####################################

beta_ini <- t(parallel::mcmapply(function(i){
  W <- rep(1, n.list[i])
  inv_XTX_i <- RcppInvGram(X[ind.strt[i]:ind.end[i],], W, 0)
  XTy_i <- RcppXtwy(X[ind.strt[i]:ind.end[i], ],y[ind.strt[i]:ind.end[i]],W)
  RcppWolsSolver03(inv_XTX_i, XTy_i, rep(0, p))
}, c(1:m.total), mc.cores = no.cores))

beta_ini_norm <- sqrt(apply(beta_ini^2, 1, sum))

####### Running the proposed method ##########################

result.uv <- fusedest_normal(X = X[ind.strt[1]:ind.end[m.list[u]],],
                             y = y[ind.strt[1]:ind.end[m.list[u]]],
                             label_dc = label_dc[ind.strt[1]:ind.end[m.list[u]]], H = H,
                             rho = rho, no_lambda = no_lambda, lambda_list = lambda_list,
                             beta_ini = beta_ini[1:m.list[u],], max_iter = max_iter,
                             tol_err = tol_err, no.cores = no.cores)


result.BS <- result.uv$alg.matrix

}

