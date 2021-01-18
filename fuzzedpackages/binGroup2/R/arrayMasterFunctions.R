
# Start  MasterPool.Array.Measures() functions
###################################################################
#    Brianna Hitt - 05-01-17
#    Purpose: calculates descriptive measures for array testing with master pooling
#      calls: f, beta0, beta1, gamma0, gamma1, beta0.y, beta1.y, sum.beta.gamma,
#             mu.y, and nu.y - functions given in the equations for array testing
#             with master pooling in "Comparison of Group Testing Algorithms for Case
#             Identification in the Presence of Test Error" by Kim et al. (2007)
#      inputs: results = an object containing results from Array.Measures()
#              n = size of a row/column in the square array
#              pmat = matrix of individual probabilities
#              Se = sensitivity of the diagnostic test
#              Sp = specificity of the diagnostic test
#      outputs: list of the expected number of tests (ET), and measures of testing accuracy,
#               including PSe, PSp, PPPV, and PNPV

# Brianna Hitt - 03.18.2020
# Edited to allow Se/Sp to vary across stages of testing

MasterPool.Array.Measures <- function(results, n, pmat, Se, Sp){
  
  # extract the measures from the results for Array.Measures()
  ET.A2 <- results$ET
  PSe.A2 <- results$PSe
  PSp.A2 <- results$PSp
  PPV.A2 <- results$PPV
  NPV.A2 <- results$NPV
  
  p <- pmat[1,1]
  
  # calculate the measures for array testing with master pooling
  
  # ET.A2M <- 1/(n^2) + (1 - Se - Sp)*((1 - p)^(n^2))*(2/n + (1 - Sp)^2 + 2*(1 - Sp)*Sp^n) + Se*ET.A2
  # the code below allows Se/Sp to vary between row and column tests
  #   (master pool, row test, column test, individual test)
  # ET.A2M <- 1/(n^2) + (1 - Se[1] - Sp[1])*((1 - p)^(n^2))*(2/n + prod(1 - Sp[2:3]) + (1 - Sp[2])*Sp[3]^n + (1 - Sp[3])*Sp[2]^n) + Se[1]*ET.A2
  # the code below assumes the row/column test have the same Se/Sp
  #   (master pool, row/col test, individual test)
  ET.A2M <- 1/(n^2) + (1 - Se[1] - Sp[1])*((1 - p)^(n^2))*(2/n + (1 - Sp[2])^2 + 2*(1 - Sp[2])*Sp[2]^n) + Se[1]*ET.A2
  
  # PSe.A2M <- Se^4 + 2*(Se^3)*(1 - Se)*(1 - f(n=n, p=pmat, Se=Se, Sp=Sp))^(n-1)
  # the code below allows Se/Sp to vary between row and column tests
  #   (master pool, row test, column test, individual test)
  # PSe.A2M <- prod(Se[1:4]) + 
  #   prod(Se[c(1,4)])*(Se[2]*(1 - Se[3])*(1 - f(n=n, p=pmat, Se=Se[3], Sp=Sp[3]))^(n-1) + 
  #                       Se[3]*(1 - Se[2])*(1 - f(n=n, p=pmat, Se=Se[2], Sp=Sp[2]))^(n-1))
  # the code below assumes the row/column test have the same Se/Sp
  #   (master pool, row/col test, individual test)
  PSe.A2M <- prod(Se[1:3], Se[2]) +
    2*prod(Se[1:3])*(1 - Se[2])*(1 - f(n=n, p=pmat, Se=Se[2], Sp=Sp[2]))^(n-1)

  # PSp.A2M <- 1 - ((1 - Sp)*mu.y(n=n, p=pmat, Se=Se, Sp=Sp) + 2*(1 - Sp)*nu.y(n=n, p=pmat, Se=Se, Sp=Sp))
  # the code below allows Se/Sp to vary between row and column tests
  #   (master pool, row test, column test, individual test)
  # NOTE: mu.y() and nu.y() will need to be revised further if Se/Sp 
  #   differs for row and column tests
  # PSp.A2M <- 1 - ((1 - Sp[4])*mu.y(n=n, p=pmat, Se=Se, Sp=Sp) + 2*(1 - Sp[4])*nu.y(n=n, p=pmat, Se=Se, Sp=Sp))
  # the code below assumes the row/column test have the same Se/Sp
  #   (master pool, row/col test, individual test)
  PSp.A2M <- 1 - ((1 - Sp[3])*mu.y(n=n, p=pmat, Se=Se, Sp=Sp) + 2*(1 - Sp[3])*nu.y(n=n, p=pmat, Se=Se, Sp=Sp))
  
  PPV.A2M <- (pmat*PSe.A2M)/((1-pmat)*(1-PSp.A2M) + pmat*PSe.A2M)
  NPV.A2M <- ((1-pmat)*PSp.A2M)/(pmat*(1-PSe.A2M) + (1-pmat)*PSp.A2M)
  
  list("ET"=ET.A2M, "PSe"=PSe.A2M, "PSp"=PSp.A2M, "PPV"=PPV.A2M, "NPV"=NPV.A2M)
}

##################################################################
# Support Functions needed for MasterPool.Array.Measures()       #
##################################################################
f <- function(n, p, Se, Sp){
  (1 - Sp)*(1 - p)^n + Se*(1 - (1 - p)^n)
}

beta0 <- function(n, c, p){
  choose(n=n, k=c)*((1 - p)^(n^2 - n*c + c))*(1 - (1 - p)^(n-1))^c
}

beta1 <- function(n, c, p){
  choose(n=n, k=c)*((1 - p)^(n^2 - n*c))*(1 - (1 - p)^n)^c - beta0(n=n, c=c, p=p)
}

gamma0 <- function(n, c, Se, Sp){
  (1 - Sp)*((1 - Se)^c)*(Sp^(n-c))
}

gamma1 <- function(n, c, Se, Sp){
  Se*(Sp^(n-c))*(1 - Se)^c
}

beta0.y <- function(n, c, p){
  beta0(n=n, c=c, p=p)/(1 - p)
}

beta1.y <- function(n, c, p){
  choose(n=(n-1), k=c)*((1 - (1 - p)^n)^c)*((1 - p)^(n^2 - n*c - 1)) + choose(n=(n-1), k=(c-1))*((1 - (1 - p)^n)^(c-1))*((1 - p)^(n^2 - n*c))*(1 - (1 - p)^(n-1)) - beta0.y(n=n, c=c, p=p)
}

sum.beta.gamma <- function(n,p,Se,Sp){
  result <- 0
  for(c in 1:n){
    result <- result + beta0.y(n=n, c=c, p=p)*gamma0(n=n, c=c, Se=Se, Sp=Sp) + beta1.y(n=n, c=c, p=p)*gamma1(n=n, c=c, Se=Se, Sp=Sp)
  }
  result
}

# original code 
# mu.y <- function(n, p, Se, Sp){
#   f(n=(n^2 - 2*n + 1), p=p, Se=Se, Sp=Sp)*((1 - Sp)*(1 - p)^(n-1))^2 + (Se^2)*(1 - (1 - p)^(n-1))*((1 - Sp)*(1 - p)^(n-1) + f(n=(n - 1), p=p, Se=Se, Sp=Sp))
# }
# Brianna Hitt - 03.18.2020
# Revised to allow Se/Sp to vary across stages of testing
#   (master pool, row/col test, individual testing)
# would need to revise further if allowing Se/Sp to vary for row and column testing
mu.y <- function(n, p, Se, Sp){
  f(n=(n^2 - 2*n + 1), p=p, Se=Se[1], Sp=Sp[2])*((1 - Sp[2])*(1 - p)^(n-1))^2 + 
    prod(Se[1:2])*(1 - (1 - p)^(n-1))*((1 - Sp[2])*(1 - p)^(n-1) + f(n=(n - 1), p=p, Se=Se[2], Sp=Sp[2]))
}

# original code 
# nu.y <- function(n, p, Se, Sp){
#   (1 - Sp)*beta0.y(n=n, c=0, p=p)*gamma0(n=n, c=0, Se=Se, Sp=Sp) + Se*sum.beta.gamma(n=n, p=p, Se=Se, Sp=Sp)
# }
# Brianna Hitt - 03.18.2020
# Revised to allow Se/Sp to vary across stages of testing
#   (master pool, row/col test, individual testing)
# would need to revise further if allowing Se/Sp to vary for row and column testing
#   and would need to revise gamma0() and gamma1() appropriately
nu.y <- function(n, p, Se, Sp){
  # original code
  # (1 - Sp)*beta0.y(n=n, c=0, p=p)*gamma0(n=n, c=0, Se=Se, Sp=Sp) + Se*sum.beta.gamma(n=n, p=p, Se=Se, Sp=Sp)
  (1 - Sp[1])*beta0.y(n=n, c=0, p=p)*gamma0(n=n, c=0, Se=Se[2], Sp=Sp[2]) + Se[1]*sum.beta.gamma(n=n, p=p, Se=Se[2], Sp=Sp[2])
}

