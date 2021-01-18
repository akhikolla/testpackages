get.best <-
function(fit.grid){
  ## input:
  ## fit.grid is output of regmed.grid
    
  if(class(fit.grid) != "regmed.grid") {
      stop("input not regmed.grid class")
  }
    
  lambda <- fit.grid$grid.dat$lambda
  bic <- fit.grid$grid.data$bic
  index.best <- (1:length(lambda))[bic == min(bic)]
  fit.best <- fit.grid$fit.list[index.best][[1]]
  grid.best <- fit.grid$grid.data[index.best,]
  return(list(fit=fit.best, grid=grid.best))
}
