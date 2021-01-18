
k1FunExp <- function(x1, x2, par) {
  res <- .Call(k1ExpC, x1, x2, par)
}

k1FunGauss <- function(x1, x2, par) {
  res <- .Call(k1GaussC, x1, x2, par)
}

k1FunPowExp <- function(x1, x2, par) {
  res <- .Call(k1PowExpC, x1, x2, par)
}

k1FunMatern3_2<- function(x1, x2, par) {
  res <- .Call(k1Matern3_2C, x1, x2, par)
}

k1FunMatern5_2 <- function(x1, x2, par) {
  res <- .Call(k1Matern5_2C, x1, x2, par)
}

k1Exp <- new("covMan",            
             kernel = k1FunExp,
             hasGrad = TRUE, 
             acceptMatrix = FALSE,   
             label = "exponential",
             d = 1L,
             inputNames = "x",
             parLower = c("range" = 1e-8, "var" = 1e-8),
             parUpper = c("range" = Inf, "var" = Inf),
             par = c("range" = 1.0, "var" = 1.0),
             parN = 2L,
             kernParNames = c("range", "var"))

k1Gauss <- new("covMan",            
               kernel = k1FunGauss,
               hasGrad = TRUE, 
               acceptMatrix = FALSE,   
               label = "gaussian",
               d = 1L,
               inputNames = "x",
               parLower = c("range" = 1e-8, "var" = 1e-8),
               parUpper = c("range" = Inf, "var" = Inf),
               par = c("range" = 1.0, "var" = 1.0),
               parN = 2L,
               kernParNames = c("range", "var"))

k1Matern3_2 <- new("covMan",            
                   kernel = k1FunMatern3_2,
                   hasGrad = TRUE, 
                   acceptMatrix = FALSE,   
                   label = "Matern nu = 3/2",
                   d = 1L,
                   inputNames = "x",
                   parLower = c("range" = 1e-8, "var" = 1e-8),
                   parUpper = c("range" = Inf, "var" = Inf),
                   par = c("range" = 1.0, "var" = 1.0),
                   parN = 2L,
                   kernParNames = c("range", "var"))

k1Matern5_2 <- new("covMan",            
                   kernel = k1FunMatern5_2,
                   hasGrad = TRUE, 
                   acceptMatrix = FALSE,   
                   label = "Matern nu = 5/2",
                   d = 1L,
                   inputNames = "x",
                   parLower = c("range" = 1e-8, "var" = 1e-8),
                   parUpper = c("range" = Inf, "var" = Inf),
                   par = c("range" = 1.0, "var" = 1.0),
                   parN = 2L,
                   kernParNames = c("range", "var"))

k1PowExp <- new("covMan",            
                kernel = k1FunPowExp,
                hasGrad = TRUE, 
                acceptMatrix = FALSE,   
                label = "power exponential",
                d = 1L,
                inputNames = "x",
                parLower = c("range" = 1e-8, "shape" = 0.0, "var" = 1e-8),
                parUpper = c("range" = Inf, "shape" = 2.0, "var" = Inf),
                par = c("range" = 1.0, "shape" = 1.5, "var" = 1.0),
                parN = 3L,
                kernParNames = c("range", "shape", "var"))

setAs("covMan", "function", function(from) from@kernel)

GPkernNames <- c("k1Exp", "k1Matern3_2", "k1Matern5_2", "k1PowExp", "k1Gauss")
