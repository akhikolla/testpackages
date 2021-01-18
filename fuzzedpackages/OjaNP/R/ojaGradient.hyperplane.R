`ojaGradient.hyperplane` <-
  function(d, x){
    return(rep(sign(round(d%*%c(x,1),digits=8)),length(x))*d[seq_along(x)])
  }

