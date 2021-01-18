
print.mcmc_hsar <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("Type:", ' hsar ', "\n")
  
  cat("\n Coefficients:\n")
  print( x$Mbetas )
  
  cat("\n Spatial Coefficients:\n")
  print( cbind( rho= x$Mrho, lambda=x$Mlambda) )
  
  cat("\n Log likelihood:", x$Log_Likelihood, "\n")
  invisible(x)
}

print.mcmc_sar <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("Type:", ' sar ', "\n")
  
  cat("\n Coefficients:\n")
  print( x$Mbetas )
  
  rho<-x$Mrho
  names(rho)<-'rho'
  cat("\n Spatial Coefficients:\n")
  print( rho )
  
  cat("\n Log likelihood:", x$Log_Likelihood, "\n")
  invisible(x)
}

print.mcmc_hsar_rho_0 <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("Type:", ' hsar with rho = 0 ', "\n")
  
  cat("\n Coefficients:\n")
  print( x$Mbetas )
  
  lambda<-x$Mlambda
  names(lambda)<-'lambda'
  cat("\n Spatial Coefficients:\n")
  print( lambda )
  
  cat("\n Log likelihood:", x$Log_Likelihood, "\n")
  invisible(x)
}

print.mcmc_hsar_lambda_0 <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("Type:", ' hsar with lambda = 0 ', "\n")
  
  cat("\n Coefficients:\n")
  print( x$Mbetas )
  
  rho<-x$Mrho
  names(rho)<-'rho'
  cat("\n Spatial Coefficients:\n")
  print( rho )
  
  cat("\n Log likelihood:", x$Log_Likelihood, "\n")
  invisible(x)
}

