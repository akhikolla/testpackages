"%^%" <- function(M,pow)
{
  if ((dim(M)[1]==1) & (dim(M)[2]==1)) return(as.matrix(M^pow))
  eM <- eigen(M); U <- eM$vectors
  D <- (eM$values + abs(eM$values))/2
  return(U %*% diag(D^pow) %*% t(U))
}

#----------------------------------------------------------#
# This function is of use in GOptim to form a pxp matrix A.#
# It is called repeatedly in a repeat loop                 #   
# We may consider writing its C++ code.                    #
#----------------------------------------------------------#
getA <- function(alpha, p, d) 
{
  A <- matrix(0, p, p)
  for (i in 1:d) 
  {
    for (j in (d + 1):p) 
    {
      Eij <- matrix(0, p, p)
      Eij[i, j] <- 1
      Eij[j, i] <- -1
      A <- A + alpha[i, j] * Eij
    }
  }
  return(round(A, digits = 4))
}

#----------------------------------------------------------------#
# This function is of use in Goptim. It computes the gradient of #
# the objective function numerically if the analytical expression#
# is not provided. Should be converted into C++                  #
#----------------------------------------------------------------#
getGradient <- function(objfun, W, fvalue, eps_grad) 
{
# SRM TESTING
  alpha <- objfun(W)$gradient
#   alpha <- NULL
 
  if (is.null(alpha)) 
  {
    Qt <- W$Qt
    p <- nrow(Qt)
    d <- W$dim[1]
    alpha <- matrix(0, nrow = d, ncol = (p - d))
    
    for (i in 1:d) 
    {
      for (j in (d + 1):p) 
      {
        Q_tilde <- Qt
        Q_tilde[, i] <- cos(eps_grad)*Qt[, i]-sin(eps_grad)*Qt[, j]
        W$Qt <- Q_tilde
        f_tilde <- round(objfun(W)$value, digits = 5)
        alpha[i, j - d] <- (f_tilde - fvalue)/eps_grad
      }
    }
  }
  return(alpha)
}

#---------------------------------#
# This function is used in GOptim.# 
# No need to be written in C++    # 
#---------------------------------#
func_delta <- function(delta, cW, matA)
{
  Expm <- matrix(attributes(expm(Matrix(-delta * matA)))$x, nrow = nrow(matA), ncol = ncol(matA))
  
  cW$Qt <- orthonorm(cW$Qt %*% t(Expm))

  objfun(cW)$value
}

#------------------------------------------------#
# The following functions fun_t and fun_YtMN     #
# are in use in SOptim. No need to turn into C++ #
#------------------------------------------------#
fun_t <- function(t, objfun, Yk, A, Q, R, dd)
{
  Yt <- fun_YtMN(t, Yk, A, Q, R, dd)$Yt
  return(objfun(Yt))
}

fun_YtMN <- function(t,Yk,A,Q,R,dd)
{
  zeros <- matrix(0,dd,dd)
  mat <- cbind(rbind(A,R), rbind(-t(R),zeros))
  
  MN <- Matrix::expm(t * mat)[,1:dd]
  
  if (dd > 1){ M <- MN[1:dd,];  N <- MN[(dd+1):(2*dd),] } else
  { M <- MN[1:dd]; N <- MN[(dd+1):(2*dd)] }
  
  Yt <- orthonorm(Yk %*% M + Q %*% N)
  
  return(list(Yt = Yt, M = as.matrix(M), N = as.matrix(N)))
}

