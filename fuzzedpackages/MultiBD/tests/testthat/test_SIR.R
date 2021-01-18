library(testthat)

test_that("Transition probabilities of SIR model", {
  
  data(Eyam)
  data <- Eyam
  
  alpha = 3.204
  beta = 0.019
  k = 1
  N <- data$S[1] + data$I[1] + data$R[1]
  
  p1 <- bbd_prob(
    t  = data$time[k + 1] - data$time[k], # Time increment
    a0 = data$R[k], b0 = data$I[k],       # From: R(t_k), I(t_k)                                      
    brates1 <- function(a,b){ 0 }, 
    brates2 <- function(a,b){ beta  * max(N - a - b, 0)  * b }, 
    drates2 <- function(a,b){ 0 }, 
    trans21 <- function(a,b){ alpha * b },
    A = data$R[k + 1], B = data$R[k + 1] + data$I[k] - data$R[k],
    computeMode = 4, nblocks = 80         # Compute using 4 threads
    )[data$R[k + 1] - data$R[k] + 1, data$I[k + 1] + 1]
  
  p2 <- dbd_prob(  # Compute the transition probability matrix
    t  = data$time[k + 1] - data$time[k], # Time increment
    a0 = data$S[k], b0 = data$I[k],       # From: S(t_k), I(t_k)                                      
    drates1 <- function(a, b) { 0 }, 
    brates2 <- function(a, b) { 0 }, 
    drates2 <- function(a, b) { alpha * b }, 
    trans12 <- function(a, b) { beta  * a * b },
    a = data$S[k + 1], B = data$S[k] + data$I[k] - data$S[k + 1],
    computeMode = 4, nblocks = 80         # Compute using 4 threads
    )[1, data$I[k + 1] + 1]
  
  p3 <- SIR_prob( 
    t  = data$time[k + 1] - data$time[k], # Time increment
    alpha = alpha, beta = beta,
    S0 = data$S[k], I0 = data$I[k],       # From: R(t_k), I(t_k)                                      
    nSI = data$S[k] - data$S[k + 1], nIR = data$R[k + 1] - data$R[k],
    direction = "Forward",
    computeMode = 4, nblocks = 80         # Compute using 4 threads
    )[data$S[k] - data$S[k + 1] + 1, data$R[k + 1] - data$R[k] + 1]
  
  expect_equal(0.0, abs(p1-p2), 1E-5)
  expect_equal(0.0, abs(p2-p3), 1E-8)
   
})