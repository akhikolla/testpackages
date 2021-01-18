########################################################################################################################
####################### The Nelder-Mead Algorithm for numerical optimization ###########################################
########################################################################################################################

# Necessary imputs:
#    fun := Function to optimize
#      k := number of parameters (or give alternitavely the starting vector)
# optional inputs:
#  start := start vector
#maximum := if maximization or minimization
#  trace := tracking the function values and parameters in every iteration step for later use e.g. in graphics
# Parameters for approximation (fine tuning)
#  alpha := 1, reflection
#   beta := 2, expansion
#  gamma := 1/2, contraction
#  delta := 1/2, shrinking
#    tol := 0.00001, tolerance of algorithm
#   exit := 500, max number of iterations
#   edge := 1, edge length of initial simplex

#------------------------#
# Starting the Algorithm #
#------------------------#

optim_nm <- function(fun, k = 0, start, maximum = FALSE, trace = FALSE, alpha = 1, beta = 2,
                       gamma = 1/2, delta = 1/2, tol = 0.00001, exit = 500, edge = 1){

  # checking for inputs ----------------------------------------------------------------------

  if (missing(fun)) {
    stop("You have to specify a function to optimize")
  }
  if (k == 0 & missing(start)) {
    stop("You have to specify at least starting values or the number of parameters")
  }
  if (missing(start)) {
    start <- c(runif(k, min = -1, max = 1))
  }
  if (is.vector(start) == FALSE) {
    stop("start has to be a vector of numbers")
  }
  if (k == 0) {
    k <- as.numeric(length(start))
  }
  if (alpha <= 0 | beta <= 0 | gamma <= 0 | delta <= 0 | edge <= 0 | tol <= 0 | exit <= 0) {
    stop("parameter have to be > 0")
  }

  # creating start simplex ---------------------------------------------------------------------

  verteces <- function(x, k){
    Simplex <- matrix(NA, nrow = k + 1, ncol = k) # creating (k+1 x k) matrix
    Simplex[1, ] <- x # first row is determined by starting values
    p1 <- function(x){
      edge * ((sqrt(k + 1) + k - 1)/sqrt(2) * k)
    }
    p2 <- function(x){
      edge * ((sqrt(k + 1) - 1)/sqrt(2) * k)
    }
    p_1 <- p1(x)
    p_2 <- p2(x)
    p <- as.numeric(t(lapply(x, p2)))
    identityv <- diag(k) # creating identity matrix (k x k)
    for (i in 1:k) { # loop to create the several vertices
      Simplex[i + 1, ] <- x + p + (p_1 - p_2) * identityv[, i]
    }
    return(Simplex)
  }

  initial_simplex <- verteces(start, k) # apply verteces function to starting vector
  function_value <- apply(initial_simplex, 1, fun) # creating function values of initial simplex
  simplex <- cbind(function_value, initial_simplex) # combines function and x values
  fun_length = length(start)

  # Maximizing Function ------------------------------------------------------------------------------

  if (maximum == TRUE) {

    simplex <-  simplex[order(simplex[, 1], decreasing = TRUE), ] # ordering simplex according to function values
    param <- simplex[, -1] # extracting the parameters from simplex
    kparam <- param[-(k + 1), ] # extracting the k best vertices
    M <- colMeans(kparam)# calculating the mean of the k best verticess
    End <- sqrt(sum((simplex[, 1] - fun(M)) ^ 2) / (k + 1)) # determine the exit out of the loop if a certain accuracy is reache
    counter <- 0

    # creating trace matrix
    if (trace == TRUE) {
      trace_array <- matrix(ncol = 1 + length(simplex[1, ]), nrow = 1,
                            c(counter, simplex[1, ]), dimnames = list(character(0),
                                                                      c("iteration", "function_value", paste("x", c(1:fun_length), sep = "_"))))
    } else {
      trace_array <- NULL # if trace = FALSE, trace matrix is empty
    }

    while (tol < End) { # loop is executed until a certain auccuaricy or maximum number of iterations are reached.

      simplex <-  simplex[order(simplex[, 1], decreasing = TRUE), ] # ordering simplex according to function values
      param <- simplex[, -1] # extracting the parameters from simplex
      kparam <- param[-(k + 1), ] # extracting the k best vertices
      M <- colMeans(kparam) # calculating the mean of the k best vertices
      counter <- counter + 1 # counts number of iterations

      if (trace == TRUE & counter >= 2) {
        trace_array <- rbind(trace_array, c(counter + 1, simplex[1], simplex[1, -1])) # writing trace for visulaization
      }

      R <- M + alpha * (M - simplex[k + 1, c(2:(k + 1))])  # Reflection
      Ry <- fun(R)

      if (Ry > simplex[k, 1]) {  # case 1: Reflection or Expension (correct direction)

        if (simplex[1, 1] >= Ry) { # using reflection if Ry is not better than the best
          simplex[k + 1, ] <- c(Ry, R)
        } else {  # Expansion if Ry is better than the best
          E <- M + beta * (R - M)
          Ey <- fun(E)
          if (Ey > Ry) { # checking if expansion leads to better values or not
            simplex[k + 1, ] <- c(Ey, E)
          } else {
            simplex[k + 1, ] <- c(Ry, R)
          }
        }

      } else { # case 2: Contraction or shrinking

        if (simplex[k + 1, 1] < Ry & Ry <= simplex[k, 1]) {  # Outside Contraction: Between worst and 2nd worst
          OC <- M + gamma * (R - M)
          OCy <- fun(OC)
          if (OCy >= simplex[k + 1, 1]) {  # using OCy if it's better than worst
            simplex[k + 1, ] <- c(OCy, OC)
          } else if (OCy < simplex[k + 1, 1]) { # shrinking if OCy didn't lead to better points
            for (i in 1:k) {
              S <- delta * (simplex[1, c(2:(k + 1))] + simplex[i + 1, c(2:(k + 1))])
              simplex[i + 1, ] <- c(fun(S), S)
            }
          }
        } else if (Ry <= simplex[k + 1, 1]) {    # Inside Contraction; worse than the worst
          IC <- M - gamma * (R - M)
          ICy <- fun(IC)
          if (ICy > simplex[k + 1, 1]) {  # using ICy if it's better than worst
            simplex[k + 1, ] <- c(ICy, IC)
          } else if (ICy <= simplex[k + 1, 1]) { # shrinking if ICy didn't lead to better points
            for (i in 1:k) {
              S <- delta * (simplex[1, c(2:(k + 1))] + simplex[i + 1,c(2:(k + 1))])
              simplex[i + 1, ] <- c(fun(S), S)
            }
          }
        }
      }
      End <- sqrt(sum((simplex[, 1] - fun(M)) ^ 2) / (k + 1)) # determine the exit out of the loop if a certain accuracy is reache
      if (counter == exit) { # emergency exit to avoid infinite loops
        warning("Maximum number of Iterations reached \n check your input parameters or change the 'exit' value")
        break
      }
    }

    # Minimizing Function -----------------------------------------------------------------------------------

  } else {

    simplex <-  simplex[order(simplex[, 1]), ] # ordering simplex according to function values
    param <- simplex[, -1] # extracting the parameters from simplex
    if(k > 1){
      kparam <- param[-(k + 1), ] # extracting the k best vertices from simplex
      M <- colMeans(kparam)  # calculating the mean of the k best points
    }else{
      kparam <- param[-(k + 1)] # extracting the k best vertices from simplex
      M <- kparam  # calculating the mean of the k best points
    }
    
    End <- sqrt(sum((simplex[, 1] - fun(M)) ^ 2) / (k + 1)) # determine the exit out of the loop if a certain accuracy is reache
    counter <- 0

    # creating trace matrix
    if (trace == TRUE) {
      trace_array <- matrix(ncol = 1 + length(simplex[1, ]), nrow = 1,
                            c(counter + 1, simplex[1, ]), dimnames = list(character(0),
                                                                          c("iteration", "function_value", paste("x", c(1:fun_length), sep = "_"))))
    } else {
      trace_array <- NULL # if trace = FALSE, trace matrix is empty
    }

    while (tol < End) {  # loop is executed until a certain auccuaricy or maximum number of iterations are reached.

      simplex <- simplex[order(simplex[, 1]), ] # ordering simplex according to function values
      param <- simplex[, -1] # extracting the parameters from simplex
      if(k > 1){
        kparam <- param[-(k + 1), ] # extracting the k best vertices from simplex
        M <- colMeans(kparam)  # calculating the mean of the k best points
      }else{
        kparam <- param[-(k + 1)] # extracting the k best vertices from simplex
        M <- kparam  # calculating the mean of the k best points
      }
      counter <- counter + 1 # counts number of iterations

      if (trace == TRUE & counter >= 2) {
        trace_array <- rbind(trace_array, c(counter, simplex[1], simplex[1, -1])) # writing trace for visulaization
      }

      R <- M + alpha * (M - simplex[k + 1, c(2:(k + 1))])  # Reflection
      Ry <- fun(R)

      if (Ry < simplex[k, 1]) {  # case 1: Reflection or Expension (correct direction)

        if (simplex[1, 1] <= Ry) { # using reflection if Ry is not better than the best
          simplex[k + 1, ] <- c(Ry, R)
        } else {  # Expansion if Ry is better than the best
          E <- M + beta * (R - M)
          Ey <- fun(E)
          if (Ey < Ry) { # checking if expansion leads to better values or not
            simplex[k + 1, ] <- c(Ey, E)
          } else {
            simplex[k + 1, ] <- c(Ry, R)
          }
        }

      } else { # case 2: Contraction or shrinking

        if (simplex[k + 1, 1] > Ry & Ry >= simplex[k, 1]) {  # Outside Contraction: Between worst and 2nd worst
          OC <- M + gamma * (R - M)
          OCy <- fun(OC)
          if (OCy <= simplex[k + 1, 1]) {  # using OCy if it's better than worst
            simplex[k + 1, ] <- c(OCy, OC)
          } else if (OCy > simplex[k + 1, 1]){ # shrinking if OCy didn't lead to better points
            for (i in 1:k) {
              S <- delta * (simplex[1, c(2:(k + 1))] + simplex[i + 1, c(2:(k + 1))])
              simplex[i + 1, ] <- c(fun(S), S)
            }
          }
        } else if (Ry >= simplex[k + 1, 1]){    # Inside Contraction: worse than the worst
          IC <- M - gamma * (R - M)
          ICy <- fun(IC)
          if (ICy < simplex[k + 1, 1]) {  # using ICy if it's better than worst
            simplex[k + 1, ] <- c(ICy, IC)
          } else if (ICy >= simplex[k + 1, 1]) { # shrinking if ICy didn't lead to better points
            for (i in 1:k) {
              S <- delta * (simplex[1, c(2:(k + 1))] + simplex[i + 1, c(2:(k + 1))])
              simplex[i + 1, ] <- c(fun(S), S)
            }
          }
        }
      }
      End <- sqrt(sum((simplex[, 1] - fun(M)) ^ 2) / (k + 1)) # determine the exit out of the loop if a certain accuracy is reache
      if (counter == exit) { # emergency exit to avoid infinite loops
        warning("Maximum number of Iterations reached \n check your input parameters or change the 'exit' value")
        break
      }
    }
  }

  # Creating the Output ----------------------------------------------------------------------------------------

  output <- list(par            = simplex[1, -1],
                 function_value = simplex[1],
                 trace          = trace_array,
                 fun            = fun,
                 start          = start,
                 upper          = NA,
                 lower          = NA,
                 control        = list(
                   k          = k,
                   iterations = counter)
  )

  #class(output) <- append(class(output), "optim_nmsa")
  class(output) <- "optim_nmsa"
  return(output)
}
