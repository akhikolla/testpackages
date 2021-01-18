Sys.setenv("R_TESTS" = "")
require(Rdtq)
require(Rcpp)
# Example 1:
# Define the drift function f(x) = -x and diffusion function g(x) = 1
# using C++ code:
sourceCpp(code = '#include <Rcpp.h>
  using namespace Rcpp;
  double drift(double& x)
  {
    return(-x);
  }
  double diff(double& x)
  {
    return(1.0);
  }
  typedef double (*funcPtr)(double& x);
  // [[Rcpp::export]]
  XPtr<funcPtr> driftXPtr()
  {
    return(XPtr<funcPtr>(new funcPtr(&drift)));
  }
  // [[Rcpp::export]]
  XPtr<funcPtr> diffXPtr()
  {
    return(XPtr<funcPtr>(new funcPtr(&diff)));
  }')

# Solve for the PDF (at final time fT=1) of the SDE with drift f,
# diffusion g, and deterministic initial condition X(0) = 0.
# First we solve using the grid specified by k and bigm.
# Then we solve using the grid specified by a, b, and bigm.
# We then check that we get the same PDF either way.
k = 0.01
M = 250
test1 = rdtq(h=0.1,k,bigm=M,init=0,fT=1,
             drift=driftXPtr(),diffusion=diffXPtr(),method="cpp")
test2 = rdtq(h=0.1,a=-2.5,b=2.5,bigm=501,init=0,fT=1,
             drift=driftXPtr(),diffusion=diffXPtr(),method="cpp")
print(k*sum(abs(test1$pdf-test2$pdf)))

