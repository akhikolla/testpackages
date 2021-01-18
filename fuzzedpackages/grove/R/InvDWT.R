#' @title Inverse discrete wavelet transform  
#'
#' @description This function performs the inverse discrete wavelet transform.
#'
#' @param grove.obj An object of class \code{grove}.
#' @param x A vector of the values of a predictor. 
#' @param include.C If \code{TRUE}, C is used for reconstructing
#' the function.
#' @param sample.C If \code{TRUE}, draws from C are used for recontructing
#' the function.

#' @return A matrix with each row representing a draw from the reconstructed signal.
#' @export
#' @examples
#' data <- wavethresh::DJ.EX(n = 512, noisy = TRUE, rsnr = 5)$doppler
#' W <- DWT(data)
#' ans <- Denoise(W)
#' denoised.data <- InvDWT(ans)
#' plot(data, type = "l")
#' lines(denoised.data[1, ], col = "red")

InvDWT <- function(grove.obj, 
                   x = NULL, # TODO: change for anova 
                   include.C = TRUE, 
                   sample.C = FALSE) {
  
  if (class(grove.obj) != "grove") {
    stop("Input should be a grove class object")
  }
  
  D <- grove.obj$samples$mean
  n_samp <- dim(D)[3]
  
  if (include.C) {
    C <- grove.obj$C_hat
  } else {
    C <- rep(0, length(grove.obj$C_hat))
  }
  
  m <- dim(D)[2] + 1 
  
  output <- matrix(NA, ncol = m, nrow = n_samp)
  temp <- wd(rep(1, m))
  
  y <- grove.obj$data$W$C
  X <- model.matrix(grove.obj$data$formula, grove.obj$data$X)
  
  s20 <- summary(lm(y ~ X))$sigma
  nu0 <- 10
  
  for (i in 1:n_samp) {
    if ((nrow(grove.obj$data$X) == 1) && (ncol(grove.obj$data$X) == 1)) {
      x <- 1
      temp$D <- rev(D[, , i])
    } else {
      if (is.null(x)) {
        x <- c(1, rep(0, dim(D)[1] - 1))
      }
      temp$D <- rev(t(D[, , i]) %*% x)
    }
    
    if (include.C && sample.C) { 
      # draw C from its posterior
      
      # This code segment is from Peter Hoff's textbook 
      # "A first course in Bayesian Statistics"
      g <- length(y)
      S <- 1
      n <- dim(X)[1] 
      p <- dim(X)[2]
      Hg <- ( g /( g+1)) * X%*% solve( t (X)%*%X)%*%t (X)
      SSRg <- t(y)%*%(diag(1,nrow=n)-Hg)%*%y
      s2 <- 1/rgamma(S,(nu0+n )/2 , ( nu0* s20+SSRg)/2 )
      Vb <- g* solve(t(X)%*%X)/(g+1)
      Eb <- Vb%*%t (X)%*%y
      E <- matrix ( rnorm(S*p , 0 , sqrt(s2)),S,p)
      
      C <- t(t(E%*%chol (Vb) ) +c(Eb) )
      #
    }
    
    temp$C[length(temp$C)] <- sum(C * x)
    output[i, ] <- wr(temp)    
  }
  return(output)
}
