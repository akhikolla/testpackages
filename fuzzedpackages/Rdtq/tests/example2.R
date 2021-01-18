Sys.setenv("R_TESTS" = "")
require(Rdtq)
# Example 2:
# We again use the drift function f(x) = -x and diffusion function g(x) = 1.
# This time, we use the method="sparse" version of DTQ.
# This requires us to define the drift and diffusion functions in R:
mydrift = function(x) { -x }
mydiff = function(x) { rep(1,length(x)) }
test = rdtq(h=0.1,k=0.01,bigm=250,init=0,fT=1,
            drift=mydrift,diffusion=mydiff,method="sparse")
plot(test$xvec,test$pdf,type='l')

