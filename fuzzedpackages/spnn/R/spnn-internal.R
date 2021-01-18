
.spnn.create <- function(){ # create empty list for spnn
  
  nn <- list(model = "Scale Invariant Probabilistic Neural Network", 
             set = NULL)
  
  return(nn)
}

.cspnn.create <- function(){ # create empty list for cspnn
  
  nn <- list(model = "Condensed Scale Invariant Probabilistic Neural Network",
             set = NULL)
  
  return(nn)
}
