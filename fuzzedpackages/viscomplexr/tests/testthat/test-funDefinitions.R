# -----------------------------------------------------------------------------
# Test that phasePortrait works with all different ways of defining
# Input functions and their arguments
# -----------------------------------------------------------------------------

# The function makeFunctionFromInput handles the user input
# if it returns an object which is a function (instead of NULL),
# this function has also be tested for producing a useful
# result with a random complex number (and additional arguments,
# if such are provided by the user).
test_that("All allowed kind of inputs are transformed into
          working functions", {
  # Single string input without more arguments
  out <- makeFunctionFromInput("(2-z)^2*(-1i+z)^3*(4-3i-z)/((2+2i+z)^4)")
  expect_equal(is.function(out), TRUE)
  # Single string input with more arguments
  out <- makeFunctionFromInput("sin(r/z)^(z-b)",
                               moreArgs = list(b = 1.2, r = 1))
  expect_equal(is.function(out), TRUE)
  # Single string input with vapply
  out <- makeFunctionFromInput(
    "vapply(z, function(z, a){
      return(prod(abs(a)/a * (a-z)/(1-Conj(a)*z)))
     },
     a = c(0.12152611+0.06171533i,  0.53730315+0.32797530i,
           0.35269601-0.53259644i, -0.57862039+0.33328986i,
          -0.94623221+0.06869166i, -0.02392968-0.21993132i),
    FUN.VALUE = complex(1))")
  expect_equal(is.function(out), TRUE)
  # Match.fun compatible definition with string
  out <- makeFunctionFromInput("cos")
  expect_equal(is.function(out), TRUE)
  # Match.fun compatible definition with function object
  out <- makeFunctionFromInput(tan)
  expect_equal(is.function(out), TRUE)
  # Match.fun compatible definition with user function object
  whatEver <- function(z) {
    (1.2 - z)^cos(z)
  }
  out <- makeFunctionFromInput(whatEver)
  expect_equal(is.function(out), TRUE)
  # Match.fun compatible definition with user function name
  strangeFun <- function(z) {
    (1.2 - z)/cos(z)
  }
  out <- makeFunctionFromInput("strangeFun")
  expect_equal(is.function(out), TRUE)
  # Match.fun compatible definition with user function object
  # and moreArgs
  jacobiTheta <- function(z, tau = 1i, nIter = 30) {
    k <- c(1:nIter)
    q <- exp(pi*1i*tau)
    g <- exp(2*pi*1i*z)
    return(1 + sum(q^(k^2)*g^k + q^(k^2)*(1/g)^k))
  }
  out <- makeFunctionFromInput(jacobiTheta, moreArgs = list(tau = 1i/2 - 1/4))
  expect_equal(is.function(out), TRUE)
  # Match.fun compatible definition with user function name
  # and moreArgs
  out <- makeFunctionFromInput("jacobiTheta", moreArgs = list(tau = 1i/2 - 1/4))
  expect_equal(is.function(out), TRUE)
  # Match.fun compatible definition with user function object
  # but without moreArgs
  out <- makeFunctionFromInput(jacobiTheta)
  expect_equal(is.function(out), TRUE)
  # Match.fun compatible definition with user function name
  # but without moreArgs
  out <- makeFunctionFromInput("jacobiTheta")
  expect_equal(is.function(out), TRUE)
})

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

