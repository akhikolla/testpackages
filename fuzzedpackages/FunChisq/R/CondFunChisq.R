cond.fun.chisq.test = function(x, y, z=NULL, data=NULL, log.p = FALSE, method = c("fchisq", "nfchisq"))
{
  # parameter check
  if(!is.null(data)){
    
    if(class(data)!="data.frame"){
      stop("data can only be a data.frame!")
    }
    
    x = data[,x]
    
    y = data[,y]
    
    if(!is.null(z))
      z = data[,z]
  }
  
  if(!is.vector(x) | any(x!=round(x))){
    stop("X can only be a discrete vector!")
  }
  if(!is.null(z) && (!is.vector(z) | any(z!=round(z)))){
    stop("Z can only be a discrete vector or NULL!")
  }
  if(!is.vector(y) | any(y!=round(y))){
    stop("Y can only be a discrete vector!")
  }
  
  method = match.arg(method)
  
  # if no conditional variable has been provided then return the standard FunChisq test.
  if(is.null(z))
  {
    return(fun.chisq.test(table(x,y), method = method, log.p = log.p))
  }
  
  # if a conditional variable is provided :
  # compute combined effect, that is, X_f (x, z -> y)
  combined_effect = fun.chisq.test(table(paste(x,z),y))
  
  # compute marginal effect, that is X_f (z -> y)
  marginal_effect = fun.chisq.test(table(z,y))
  
  # compute conditional statistic, 
  # X_f (x->y | z) = X_f (x,z -> y) (combined effect) - X_f (z->y) (marginal effect)
  statistic = combined_effect$statistic - marginal_effect$statistic
  parameter = (length(unique(x)))*(length(unique(y))-1)*(length(unique(z))-1)
  
  # compute conditional functional index
  # Given as : X_f (x->y | z) / (n*(s-1) - marginal effect)
  max.chisq = ((nrow(as.matrix(x)) * (length(unique(y))-1)) - marginal_effect$statistic)
  estimate = sqrt(statistic/max.chisq)
  
  names(statistic) = "statistic"
  names(parameter) = "parameter"
  names(estimate) = "conditional functional index"
  
  if(method == "fchisq") {
    
    p.value = pchisq(statistic, df=parameter, lower.tail = FALSE, log.p = log.p)
    return(structure(list(statistic = statistic, parameter = parameter, p.value = p.value, data.name = "Y = f(X) | Z",
                          method = "Conditional functional chi-square test", estimate = estimate),
                     class = "htest"))
  } else {
    
    statistic = (statistic - parameter)/sqrt(2*parameter)
    p.value = pnorm(statistic, lower.tail = FALSE, log.p = log.p)
    return(structure(list(statistic = statistic, parameter = parameter, p.value = p.value, data.name = "Y = f(X) | Z",
                          method = "Normalized conditional functional chi-square test", estimate = estimate),
                     class = "htest"))
  }
}