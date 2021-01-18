checkinputs <- function(y,
                        K,
                        M,
                        smartinit,
                        nbinit,
                        tol,
                        nbKeep,
                        iterSmall,
                        nbcores
                      ){
  if (!is.list(y)) 
    stop("y must be a list, each element of the list corresponds to one subject")
  if (any(sapply(y, function(u) !is.numeric(u))))
    stop("Each element of y must be a numeric vector")
  if (any(sapply(y, function(u) any(na.omit(u)<0))))
    stop("Each element of y must be positive or zero")
  if (any(sapply(y, function(u) is.na(u[1]))))
    stop("The first element of an observed sequence cannot be a missing value")
  if (any(sapply(y, function(u) is.na(u[length(u)]))))
    stop("The last element of an observed sequence cannot be a missing value")
  
  if (!is.numeric(K))
    stop("K must be a numeric vector of length one")
  if (length(K)!=1)
    stop("K must be a numeric vector of length one")
  if ( (K!=ceiling(K)) || (K<1) )
    stop("K must be a positive integer")
  
  if (!is.numeric(M))
    stop("M must be a numeric vector of length one")
  if (length(M)!=1)
    stop("M must be a numeric vector of length one")
  if ( (M!=ceiling(M)) || (M<2) )
    stop("M must be a positive integer greater or equal to two")  
  if (all(sapply(y, function(u) length(u)<(2*M))))
    stop("M is to large compare to the length of the sequences")
  
  if (!is.logical(smartinit))
    stop("smartinit must be a logical")
  if (length(smartinit)!=1)
    stop("smartinit must be a logical of length one")  
  
  if (!is.numeric(nbinit))
    stop("nbinit must be a numeric vector")
  if (length(nbinit)!=1)
    stop("nbinit must be a numeric vector of length one")
  if ( (nbinit!=ceiling(nbinit)) || (nbinit<1) )
    stop("nbinit must be a positive integer greater")  
  
  if (!is.numeric(tol))
    stop("tol must be a numeric vector of length one")
  if (length(tol)!=1)
    stop("tol must be a numeric vector of length one")
  if ( (tol<=0) )
    stop("tol must be a positive number")   
  
  if (!is.numeric(nbKeep))
    stop("nbinit must be a numeric vector")
  if (length(nbKeep)!=1)
    stop("nbinit must be a numeric vector of length one")
  if ( (nbKeep!=ceiling(nbKeep)) || (nbKeep<1)  || (nbKeep>nbinit) )
    stop("nbinit must be a positive integer smaller than nbinit")  
  
  if (!is.numeric(iterSmall))
    stop("iterSmall must be a numeric vector of length one")
  if (length(iterSmall)!=1)
    stop("iterSmall must be a numeric vector of length one")
  if ( (iterSmall!=ceiling(iterSmall)) || (iterSmall<1) )
    stop("iterSmall must be a positive integer")
  
  if (!is.numeric(nbcores))
    stop("nbcores must be a numeric vector of length one")
  if (length(nbcores)!=1)
    stop("nbcores must be a numeric vector of length one")
  if ( (nbcores!=ceiling(nbcores)) || (nbcores<1) )
    stop("nbcores must be a positive integer")  
  


}