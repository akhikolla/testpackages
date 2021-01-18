
calculate_reduced_vara <- function(Zmat=NULL, X=NULL, varE=NULL, varG=NULL, invMMt=NULL, MMtsqrt=NULL, quiet=TRUE)
{
 ## internal function to AM
 ## Using var(\hat(a)) = simgaG - Cjj  where Cjj is the component from C^-1 (Henderson's 
 ##   mixed model equations coefficient matrix.   See Verbyla et al. TAG 2007.

 ##  Mixed model equations for the linear mixed model
 ##
 ##  X^T %*% R^-1  X                  X^T %*% R^-1 %*% Ze %*% Zmat
 ##
 ##
 ##  Zmat^t %*% Ze^t %*% R^-1 %*% X   Zmat %*% Ze^t %*% R^-1 %*% Ze %*% Zmat^t  +  G^-1
 ##
 ##  Ze = MMt^0.5
 ##  R  = (varE * I)^-1
 ##  G  = (varG * I)^-1
 ## 


  if(is.null(Zmat)){
     ## first principals
     Ze <- MMtsqrt
     #R1  <- solve( varE * diag(nrow(invMMt)))
     #G1  <- solve( varG * diag(nrow(invMMt)))
     A <- (1/varE) * crossprod(X)      ##old way:-  A <- t(X) %*% R1 %*% X
     B <- (1/varE) * t(X) %*% Ze       ##old way:-  B <- t(X) %*% R1 %*% Ze
     C <- (1/varE) * t(Ze) %*% X       ##old way :- C <-  t(Ze) %*% R1 %*% X
     D <- (1/varE) *  t(Ze) %*% Ze + (1/varG) * diag(nrow(invMMt))     ##   D <- t(Ze) %*% R1 %*% Ze + G1
  } else {
     Ze <- MMtsqrt
     A <- (1/varE) * crossprod(X)      
     B <- (1/varE) * t(X) %*% Zmat %*% Ze 
     C <- (1/varE) * t(Ze) %*% t(Zmat) %*%  X       
     D <- (1/varE) * t(Ze) %*% t(Zmat) %*% Zmat %*% Ze + (1/varG) * diag(nrow(invMMt))     
  }
  D1 <- solve(D)
  vars <- varG * diag(nrow(D1))  - ( D1 + D1 %*% C %*% solve(A - B %*% D1 %*% C) %*% B %*% D1 )
  return(vars )

}



