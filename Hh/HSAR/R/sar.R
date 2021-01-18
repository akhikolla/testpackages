#
#

sar <- function(formula, data= NULL, W, burnin=5000, Nsim=10000, thinning=1,
                parameters.start = NULL ) {
  
    ## check input data and formula
    frame <- check_formula(formula, data)
    X <- get_X_from_frame(frame)
    y <- get_y_from_frame(frame)
    
    if (any(is.na(y))) stop("NAs in dependent variable", call. = FALSE)
    if (any(is.na(X))) stop("NAs in independent variable", call. = FALSE)
    
    n <- nrow(X)
    
    check_matrix_dimensions(W,n,'Wrong dimensions for matrix W' )
    
    detval <- lndet_imrw(W)
    
    #start parameters
    if (! is.null(parameters.start)){
      if(is_there_parameter(parameters.start, "rho")) rho <- parameters.start$rho else rho<-0.5
      if(is_there_parameter(parameters.start, "sigma2e")) sigma2e<- parameters.start$sigma2e else sigma2e <-2.0
      if(is_there_parameter(parameters.start, "betas")) {
        betas <- parameters.start$betas
        if (dim(X)[2]!= length(betas) ) stop("Starting values for Betas have got wrong dimension", call. = FALSE)
      } 
      else betas <- coef(lm(formula,data))
    }
    else{
      rho<-0.5
      sigma2e <-2.0
      betas <- coef(lm(formula,data))
    }
    
    result<- sar_cpp_arma(X, y, W, detval, burnin, Nsim, thinning, rho, sigma2e, betas ) 
      #.Call("HSAR_sar_cpp_arma", PACKAGE = 'HSAR', X, y, W, detval, 
       #            burnin, Nsim, thinning, rho, sigma2e, betas )
    
    class(result) <- "mcmc_sar"
    
    result$cbetas<-put_labels_to_coefficients(result$cbetas, colnames(X))
    result$Mbetas<-put_labels_to_coefficients(result$Mbetas, colnames(X))
    result$SDbetas<-put_labels_to_coefficients(result$SDbetas, colnames(X))
    
    result$labels <- colnames(X)
    result$call <-match.call()
    result$formula <- formula
    
    result
}
