## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
require(lhs)

## ----q1-----------------------------------------------------------------------
a <- (1:10) 
b <- (20:30) 
dataGrid <- expand.grid(a, b)

## ----a1-----------------------------------------------------------------------
X <- randomLHS(22, 2) 
X[,1] <- 1 + 9*X[,1] 
X[,2] <- 20 + 10*X[,2] 

# OR 

Y <- randomLHS(22, 2) 
Y[,1] <- qunif(Y[,1], 1, 10) 
Y[,2] <- qunif(Y[,2], 20, 30) 

head(X)
head(Y)

## ----a12----------------------------------------------------------------------
a <- c(sample(1:3,1), sample(4:6, 1), sample(7:9, 1)) 
b <- c(sample(21:24,1), sample(25:28, 1), sample(29:32,1))

## ----a13----------------------------------------------------------------------
integerLHS <- function(n, intGroups) 
{ 
  stopifnot(all(lapply(intGroups, function(X) length(X) %% n) == 0)) 
  stopifnot(require(lhs)) 
  stopifnot(is.list(intGroups)) 
  ranges <- lapply(intGroups, function(X) max(X) - min(X)) 
  A <- matrix(nrow = n, ncol = length(intGroups)) 
  for (j in 1:length(ranges)) 
  { 
    sequ <- order(runif(n)) 
    if (length(intGroups[[1]]) > 1) 
    { 
      spacing <- intGroups[[j]][2] - intGroups[[j]][1] 
    } else stop("must have more than 1 intGroup") 
    for (k in 1:n) 
    { 
      i <- sequ[k] 
      a <- min(intGroups[[j]]) + (i - 1)*(ranges[[j]] + spacing)/n 
      b <- min(intGroups[[j]]) + i*(ranges[[j]] + spacing)/n - 1 
      if (a < b) 
      { 
        A[k,j] <- sample(seq(a,b,spacing), 1) 
      } else if (a == b) 
      { 
        A[k,j] <- a 
      } else stop("error") 
    } 
  } 
  return(A) 
} 

integerLHS(10, list(1:10, 31:40)) 
integerLHS(5, list(1:10, 31:40)) 
integerLHS(2, list(1:10, 31:40)) 
integerLHS(5, list(1:20, 31:60, 101:115)) 
integerLHS(5, list(seq(2,20,2), 31:60, 101:115)) 

## ----a21----------------------------------------------------------------------
x <- seq(0.05, 0.95, length = 10) 
y <- 1 - x 
all.equal(x + y, rep(1, length(x))) 
hist(x, main = "") 
hist(y, main = "") 

## ----a22----------------------------------------------------------------------
x <- seq(0.05, 0.95, length = 10) 
y <- runif(length(x), 0, 1 - x) 
z <- 1 - x - y 
hist(x, main = "") 
hist(y, main = "") 
hist(z, main = "") 

## ----a3, fig.width=5, fig.height=5--------------------------------------------
N <- 1000 
x <- randomLHS(N, 3) 
y <- x 
y[,1] <- qnorm(x[,1], 1, 2) 
y[,2] <- qnorm(x[,2], 3, 4) 
y[,3] <- qunif(x[,3], 5, 10) 

par(mfrow = c(2,2)) 
dummy <- apply(x, 2, hist, main = "") 

par(mfrow = c(2,2)) 
dummy <- apply(y, 2, hist, main = "") 

## ----a24, fig.width=5, fig.height=5-------------------------------------------
x <- randomLHS(N, 5) 
y <- x 
y[,1] <- x[,1]/rowSums(x[,1:3]) 
y[,2] <- x[,2]/rowSums(x[,1:3]) 
y[,3] <- x[,3]/rowSums(x[,1:3]) 
y[,4] <- x[,4] 
y[,5] <- x[,5] 

par(mfrow = c(2,3)) 
dummy <- apply(x, 2, hist, main = "") 

par(mfrow = c(2,3)) 
dummy <- apply(y, 2, hist, main = "") 

all.equal(rowSums(y[,1:3]), rep(1, nrow(y))) 

## ----a25----------------------------------------------------------------------
par(mfrow = c(1,1)) 
pairs(x) 
pairs(y, col = "red") 

## ----a26----------------------------------------------------------------------
N <- 10 
x <- randomLHS(N, 5) 
y <- x 
y[,1] <- x[,1]/rowSums(x[,1:3]) 
y[,2] <- x[,2]/rowSums(x[,1:3]) 
y[,3] <- x[,3]/rowSums(x[,1:3]) 
y[,4] <- x[,4] 
y[,5] <- x[,5] 

pairs(x) 
pairs(y, col = "red") 

## ----qdirichlet---------------------------------------------------------------
qdirichlet <- function(X, alpha) 
{ 
  # qdirichlet is not an exact quantile function since the quantile of a 
  #  multivariate distribtion is not unique 
  # qdirichlet is also not the quantiles of the marginal distributions since 
  #  those quantiles do not sum to one 
  # qdirichlet is the quantile of the underlying gamma functions, normalized 
  # This has been tested to show that qdirichlet approximates the dirichlet 
  #  distribution well and creates the correct marginal means and variances 
  #  when using a latin hypercube sample 
  lena <- length(alpha) 
  stopifnot(is.matrix(X)) 
  sims <- dim(X)[1] 
  stopifnot(dim(X)[2] == lena) 
  if (any(is.na(alpha)) || any(is.na(X))) 
    stop("NA values not allowed in qdirichlet") 

  Y <- matrix(0, nrow = sims, ncol = lena) 
  ind <- which(alpha != 0) 
  for (i in ind) 
  { 
    Y[,i] <- qgamma(X[,i], alpha[i], 1) 
  } 
  Y <- Y / rowSums(Y) 
  return(Y) 
} 

X <- randomLHS(1000, 7) 
Y <- qdirichlet(X, rep(1,7)) 
stopifnot(all(abs(rowSums(Y) - 1) < 1E-12)) 
range(Y) 

ws <- randomLHS(1000, 7) 
wsSums <- rowSums(ws) 
wss <- ws / wsSums 
stopifnot(all(abs(rowSums(wss) - 1) < 1E-12)) 
range(wss)

## ----custom, fig.width=5, fig.height=5----------------------------------------
require(lhs) 

# functions you described 
T1 <- function(t) t*t 
WL1 <- function(T1, t) T1*t 
BE1 <- function(WL1, T1, t) WL1*T1*t 

# t is distributed according to some pdf (e.g. normal) 
# draw a lhs with 512 rows and 3 columns (one for each function) 
y <- randomLHS(512, 3) 
# transform the three columns to a normal distribution (these could be any 
# distribution) 
t <- apply(y, 2, function(columny) qnorm(columny, 2, 1)) 
# transform t using the functions provided 
result <- cbind( 
  T1(t[,1]), 
  WL1(T1(t[,2]), t[,2]), 
  BE1(WL1(T1(t[,3]), t[,3]), T1(t[,3]), t[,3]) 
) 
# check the results 
# these should be approximately uniform 
par(mfrow = c(2,2)) 
dummy <- apply(y, 2, hist, breaks = 50, main = "") 
# these should be approximately normal 
par(mfrow = c(2,2)) 
dummy <- apply(t, 2, hist, breaks = 50, main = "") 
# these should be the results of the functions 
par(mfrow = c(2,2)) 
dummy <- apply(result, 2, hist, breaks = 50, main = "") 

## ----q6, fig.height=5, fig.width=5--------------------------------------------
N <- 1000 
set.seed(1919) 

x <- randomLHS(N, 4) 
y <- x 
# uniform on 1-10 
y[,1] <- ceiling(qunif(x[,1], 0, 10)) 
# three colors 1,2,3 
y[,2] <- ceiling(qunif(x[,2], 0, 3)) 
# other distributions 
y[,3] <- qunif(x[,3], 5, 10) 
y[,4] <- qnorm(x[,4], 0, 2) 

par(mfrow=c(2,2)) 
dummy <- apply(x, 2, hist, main="") 

par(mfrow=c(2,2)) 
plot(1:10, c(table(y[,1])), type="h", col="blue", lwd=2, ylim=c(0,120), 
ylab="Frequency", xlab="y[,1]") 
plot(1:3, c(table(y[,2])), type="h", col="blue", lwd=2, ylim=c(0,400), 
ylab="Frequency", xlab="y[,2]") 
hist(y[,3], main="") 
hist(y[,4], main="") 

# change to color names 
z <- as.data.frame(y) 
z[,2] <- factor(y[,2], labels=c("R","G","B")) 
z[1:10,] 

