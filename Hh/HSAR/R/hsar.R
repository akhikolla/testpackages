#
#

hsar <- function(formula, data = NULL, W=NULL, M=NULL, Delta, burnin=5000, Nsim=10000
                 , thinning=1, parameters.start = NULL) {
    
    ## check input data and formula
    frame <- check_formula(formula, data)
    X <- get_X_from_frame(frame)
    y <- get_y_from_frame(frame)
    
    if (any(is.na(y))) stop("NAs in dependent variable", call. = FALSE)
    if (any(is.na(X))) stop("NAs in independent variable", call. = FALSE)
    
    if( is.null(W) & is.null(M) ) stop('Both weight matrices can not be NULL', call. = FALSE)
    
    n <- nrow(Delta)
    p <- ncol(Delta)
    
    if( !is.null(W) ) check_matrix_dimensions(W,n,'Wrong dimensions for matrix W' )
    if( !is.null(M) ) check_matrix_dimensions(M,p,'Wrong dimensions for matrix M' )
    
    Unum <- apply(Delta,2,sum)
    
    #start parameters
    if (! is.null(parameters.start)){
      if(is_there_parameter(parameters.start, "rho")) rho <- parameters.start$rho else rho<-0.5
      if(is_there_parameter(parameters.start, "sigma2e")) sigma2e<- parameters.start$sigma2e else sigma2e <-2.0
      if(is_there_parameter(parameters.start, "rho")) lambda <- parameters.start$lambda else lambda<-0.5
      if(is_there_parameter(parameters.start, "sigma2e")) sigma2u<- parameters.start$sigma2u else sigma2u <-2.0
      if(is_there_parameter(parameters.start, "betas")) {
        betas <- parameters.start$betas
        if (dim(X)[2]!= length(betas) ) stop("Starting values for Betas have got wrong dimension", call. = FALSE)
      } 
      else betas <- coef(lm(formula,data))
    }
    else{
      rho <- 0.5
      lambda <- 0.5
      sigma2e <- 2.0
      sigma2u <- 2.0
      betas <- coef(lm(formula,data))
    }
    
    ## Call various models
    # Special case where rho =0 ; dependent regional effect 
    if (is.null(W)){
      detval <- lndet_imrw(M)
      result <- hsar_cpp_arma_rho_0(X, y, M, Delta, detval, Unum, burnin, Nsim, thinning, lambda, sigma2e, sigma2u, betas)
        #.Call("HSAR_hsar_cpp_arma_rho_0", PACKAGE = 'HSAR', X, y, M, Delta, detval, Unum, 
         #             burnin, Nsim, thinning, lambda, sigma2e, sigma2u, betas)
      class(result) <- "mcmc_hsar_rho_0"
    }
    # Special case where lamda =0 ; independent regional effect
    if ( is.null(M)){
      detval <- lndet_imrw(W)
      result <- hsar_cpp_arma_lambda_0(X, y, W, Delta, detval, Unum, burnin, Nsim, thinning, rho, sigma2e, sigma2u, betas) 
        #.Call("HSAR_hsar_cpp_arma_lambda_0", PACKAGE = 'HSAR', X, y, W, Delta, detval, Unum, 
         #             burnin, Nsim, thinning, rho, sigma2e, sigma2u, betas)
      class(result) <- "mcmc_hsar_lambda_0"
    }
    # Full HSAR model
    if ( (!is.null(M)) & (!is.null(W))){
      detval <- lndet_imrw(W)
      detvalM <- lndet_imrw(M)
      result <- hsar_cpp_arma(X, y, W, M, Delta, detval, detvalM, Unum, burnin, Nsim, thinning, rho, lambda, sigma2e, sigma2u, betas) 
        #.Call("HSAR_hsar_cpp_arma", PACKAGE = 'HSAR', X, y, W, M, Delta, detval, detvalM, Unum, 
         #             burnin, Nsim, thinning, rho, lambda, sigma2e, sigma2u, betas)
      class(result) <- "mcmc_hsar"
    }
    
    result$cbetas<-put_labels_to_coefficients(result$cbetas, colnames(X))
    result$Mbetas<-put_labels_to_coefficients(result$Mbetas, colnames(X))
    result$SDbetas<-put_labels_to_coefficients(result$SDbetas, colnames(X))
    
    result$labels <- colnames(X)
    result$call <- match.call()
    result$formula <- formula
    
    return(result)
}
