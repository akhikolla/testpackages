# This file calculates the confidence interval for Laura's strength in
# series 6 of MasterChef, using profile likelihood.  It takes a very
# long time to run.


library("hyper2")
data("masterchef")

`profile_likelihood` <- function(lauras_strength,upper=TRUE,give=FALSE){

  n <- 13
  UI <- rbind(diag(n-1),-1) # p_i >= 0
  CI <- c(rep(0,n-1),-1)     # p_1+...+p_{n-1} <= 1


  if(upper){val <- 1.01} else {val <- 0.9}
  startp <-  # hot start
    c(0.084048569574384, 0.0558987436910824, 0.104876293529251,
      0.0193074053977768, 0.0866801340969936, 1.87723524092112e-09,
      0.0767382663925213, 0.0133022689718716, lauras_strength * val,
      0.0286590852735936, 0.0171995606825605, 1.4764744281375e-07)
  
  if(upper){
    UI <- rbind(UI,c(0,0,0,0,0,0,0,0,1,0,0,0))  
    CI <- c(CI,lauras_strength)
  } else {
    UI <- rbind(UI,c(0,0,0,0,0,0,0,0,-1,0,0,0)) 
    CI <- c(CI,-lauras_strength)
  }

  ans2 <-
    constrOptim(
        theta = startp,
        f = function(p){-like_series(p,masterchef_series6)},
        grad=NULL,
        ui = UI, ci=CI,
        control=list(trace=100,maxit=300)
    )
  if(give){
    return(ans2)
  } else {
    return(ans2$value)
  }
}

probs <-
  c(0.1, 0.14, 0.145, 0.146, 0.15, 0.2, 0.27, 0.3, 0.4, 0.46,
    0.465, 0.47, 0.5)

laura_likelihood <- probs + NA

f <- function(i){
  profile_likelihood(lauras_strength=probs[i],upper=probs[i] > 0.27,give=FALSE)
}

for(i in seq_along(probs)){
  laura_likelihood[i] <- f(i)
}

plot(probs,laura_likelihood-min(laura_likelihood))
abline(h=2)
