# Copyright Lilia Leticia Ramirez Ramirez, llramirezramirez@gmail.com


###################### NETWORK CONSTRUCTION FUNCTIONS ##########################


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

li <- function(delta, b, lim = 10000) {
  # b can be a vector
  res <- rep(0, length(b))
  if (length(delta) != 1)
    stop("delta debe ser escalar") else res[b == 1] <- VGAM_zeta(delta)
    # else if(length(delta)==length(b))res[b==1]<-VGAM_zeta(delta[b==1])
    a <- which(b != 1)
    if (length(a))
      for (j in a) res[j] <- sum(b[j]^(1:lim)/(1:lim)^delta)
    res
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

resample <- function(x, size, rep = FALSE, prob = NULL) {
  # Suggested function to improve the sample function in R
  if (length(x) <= 1) {
    if (!missing(size) && size == 0)
      x[FALSE] else x
  } else sample(x, size, replace = rep, prob = prob)
  # browser()
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

NDIM <- function(a, as.row = FALSE) {
  # I created NDIM,to use with either matrices or vectors takes a vector as
  # column matrix if as.row=FALSE and row matrix if as.row=TRUE
  if ((is.vector(a) | is.list(a)) & as.row)
    n.dim <- c(1, length(a)) else n.dim <- c(NROW(a), NCOL(a))
    n.dim
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#


elements.matrix <- function(matrix, esc, cases = 0, rows.mat = 0,
                            as.row = FALSE, cols = FALSE) {
  # "matrix" can be a matrix or a vector representation of a matrix with rows=rows.mat.
  # If cols=TRUE This function returns:
  # a) The (row,column) position of the elements in "matrix" with value "esc", or
  # b)transform the element "case"-th in the vector "matrix" into (row, column)
  # for "matrix" as a matrix if cols=FALSE only returns the rows
  # When "matrix" is a vector, it can be considered column (as.row=F)
  # or row (as.row=TRUE) cases is the position in the matrix of the elements
  # equal to esc (by columns)
  if (cases == 0)
    cases <- which(matrix == esc)
  if (rows.mat == 0)
    rows.mat <- NDIM(matrix, as.row)[1]
  rows <- (cases - 1)%%rows.mat + 1
  result <- rows
  if (cols) {
    columns <- (cases - 1)%/%rows.mat + 1
    result <- list(rows = rows, cols = columns)
  }
  result
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

order.edges <- function(edges, ord.col = TRUE) {
  # this function is to order the network matrix (or vector). Output is always a
  #  matrix order the network.edges so the rows are increasing ord.col=T
  #  if we want the program to order columns as well, so in the first column
  #  appears the smaller ID's
  if (is.matrix(edges)) {
    if (ord.col) {
      keep <- edges[, 1] < edges[, 2]
      edges[keep == 0, ] <- cbind(edges[keep == 0, 2], edges[keep == 0, 1])
    }
    edges <- edges[order(edges[, 1], edges[, 2]), ]
  } else edges <- t(as.matrix(sort(edges)))
  dimnames(edges) <- NULL
  edges
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

compare.vectors <- function(y, x) {
  # I define this function us use it with "sapply"
  a <- sum(y == x)
  a
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

concatenate <- function(edges.hist, vec.mat, constant = NULL) {
  # Function to paste vec.mat (a matrix or vector) as new rows (or row)
  # to a matrix edges.hist constant is a vector of
  # constants (usually: time,arc.width,color) edges.hist is data.frame
  # vec.mat is not empty (is.null is not good)
  if (sum(NDIM(vec.mat) > 0) == 2) {
    if (is.vector(vec.mat))
      vec.mat <- t(as.matrix(vec.mat))  # vector (scalars are vectors) is now
                                        # a matrix with one row
    vec.mat <- cbind(vec.mat, matrix(rep(constant, nrow(vec.mat)), nrow(vec.mat),
                                     length(constant), byrow = TRUE))
    vec.mat <- as.data.frame(vec.mat)
    names(vec.mat) <- names(edges.hist)  #to be able to use rbind next
    edges.hist <- rbind(edges.hist, vec.mat)
  }
  edges.hist
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

rztpois <- function(n, lambda) {
  # Function so simulate random numbers from: Zero-truncated Poisson distribution
  # P(N=k)=(lambda^k*exp^{-lambda})/(k!*(1-e^{-lambda})), for k=1,2,3,..
  maxval <- max(50, lambda + 10 * sqrt(lambda))
  cdf <- (cumsum(stats::dpois(1:maxval, lambda)))/(1 - exp(-lambda))
  cut(stats::runif(n), unique(c(0, cdf, 1)), labels = FALSE, right = FALSE)
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

rlog <- function(n, param) {
  # Function to generate random numbers from
  # the Logarithmic distribution (poly-logarithmic when alpha=1)
  # where param is a positive real number.
  if (length(param) != 1 | param < 0)
    stop("the Logarithmic parameter must be a positive real number") else {
    pdf.s <- (1:500)^(-1) * exp(-(1:500)/param)/(-log(1 - exp(-1/param)))
    sim <- cut(stats::runif(n), unique(c(0, cumsum(pdf.s))), labels = FALSE, right = FALSE)
  }
  sim
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

rpower.law <- function(n, param) {
  # Function to generate random numbers from the power law k^{-param}
  # where param a real number >1
  if (length(param) != 1 | param <= 1)
    stop("the power law parameter must be a real number grater than one") else {
    pdf.s <- (1:1000)^(-param) * (VGAM_zeta(param))^(-1)
    sim <- cut(stats::runif(n), unique(c(0, cumsum(pdf.s))), labels = FALSE, right = FALSE)
  }
  sim
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

dpower.law <- function(x, param) {
  # Function to generate the distribution for the power law k^{-param},
  # where param is real number
  if (length(param) != 1)
    stop("the power law parameter dimension is wrong (it must be of length 1)
         or the parameter value is incorrect") else {
               pdf.s <- x^(-param) * (VGAM_zeta(param))^(-1)
         }
  pdf.s
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

rpoly.log <- function(n, param) {
  # Function to generate random numbers with
  # the distribution Gutenberg-Richter law ck^{-delta}e^{-k/lambda},
  # with: c=(li_lambda(e^{-1/lambda})) param=c(delta,lambda)
  if (length(param) != 2)
    stop("the Gutenber-Richter law parameter is wrong") else {
          sim <- cut(stats::runif(n), unique(c(0, cumsum(dpoly.log(1:1000, param)))),
                     labels = FALSE, right = FALSE)
    }
  sim
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

dpoly.log <- function(x, param, lim = 10000) {
  # Function to generate the distribution for the Gutenberg-Richter law
  # ck^{-delta}e^{-k/lambda}
  # c=(li_lambda(e^{-1/lambda})) param=c(delta,lambda)
  if (length(param) != 2)
    stop("the Gutenber-Richter law parameter is wrong") else {
          pdf.s <- x^(-param[1]) * exp(-x/param[2])/li(param[1], exp(-1/param[2]), lim = lim)
    }
  pdf.s
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

sdegree <- function(n, distrib, param = 0) {
  # rutine to generate the degrees based on the specified distribution n is
  # the number of indivudual to connect distrib is
  # the name of the distribution param is the unlist distribution parameter
  degree <- 3
  while (sum(degree)%%2 == 1) {
    # first step to realizable graph (Newman, Stogatz and Watts, 2001)
    if (distrib == "fixed") {
      if (length(param) != n | sum(param)%%2 == 1)
        stop("the length of the degrees must be equal to the number of nodes
             and the sum of degree must to be even")
      degree <- param
    } else if (distrib == "pois")
      degree <- stats::rpois(n, param) else if (distrib == "ztpois")
      degree <- rztpois(n, param) else if (distrib == "geom")
      degree <- stats::rgeom(n, 1/(param + 1))  # sometimes called exponential graph.
                                         # param is then the expected value
 else if (distrib == "ztgeom")
      degree <- stats::rgeom(n, 1/(param + 1)) + 1 else if (distrib == "nbinom" && length(param) == 2)
      degree <- stats::rnbinom(size = param[1], prob = param[2]) else if (distrib == "poly.log")
      degree <- rpoly.log(n, param) else if (distrib == "logarithmic")
      degree <- rlog(n, param) else if (distrib == "power.law")
      degree <- rpower.law(n, param) else if (distrib == "full")
      degree <- rep(n - 1, n) else if (distrib == "none")
      degree <- rep(0, n) else stop("incorrect distribution specification")
  }
  degree
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

edge.check <- function(edges, new.edge) {
  new.edge <- sort(new.edge)
  a <- which(edges == new.edge[1])
  b <- elements.matrix(edges, new.edge[2])
  if (length(a) == 0 | length(b) == 0)
    comp <- 0 else comp <- sapply(a, compare.vectors, b)
  comp
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

table.mult.degree <- function(vector, place, m, freq = TRUE) {
  if (freq) {
    if (place < m){
      res <- c(vector[place] * (vector[place] - 1)/2,
               vector[place] * vector[(place + 1):m])
      } else {
            res <- vector[place] * (vector[place] - 1)/2
      }
  } else {
    res <- vector[place] * vector[place:m]
  }
  res
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
