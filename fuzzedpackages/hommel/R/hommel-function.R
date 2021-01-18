hommel <- function(p, simes = TRUE)
{
  if(any(is.na(p))) 
    stop("missing values in input p-values")
  
  if(length(p)==0)
    stop("the pvalues vector is empty")    
  
  #save names in right order
  names <- names(p)
  
  perm <- order(p)
  sortedp <- p[perm]
  
  m <- length(p)
  
  simesfactor <- findsimesfactor(simes, m)
  
  jumpalpha <- findalpha(sortedp, m, simesfactor, simes)
  
  adjusted <- adjustedElementary(sortedp, jumpalpha, m, simesfactor)
  adjusted[perm] <- adjusted
  names(adjusted) <- names

  out <- new("hommel",
             p = p,
             jumpalpha = jumpalpha,
             sorter = perm,
             adjusted = adjusted,
             simesfactor = simesfactor,
             simes = simes)
  
  return(out)
  
}

