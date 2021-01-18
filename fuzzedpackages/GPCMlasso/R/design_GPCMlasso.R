## create responses for acat from ordinal values
createResponse <- function(Y){
  Y <- as.factor(Y)
  model.matrix(~0+Y)[,-length(levels(Y))]
}

design_GPCMlasso <- function(formula = formula, Y=Y, data = data, RSM = RSM, 
                             GPCM = GPCM, DSF = DSF, all.dummies = TRUE, main.effects = main.effects){

  ## extract covariates
  if(all.dummies){
    term.labels <- attr(terms(formula), "term.labels")
    X <- matrix(rep(1,nrow(data)))
    for(ij in 1:length(term.labels)){
      x.now <- data[,term.labels[ij], drop = FALSE]
      if(is.factor(x.now[,1])){
        if(nlevels(x.now[,1])==2){
          X <- cbind(X,model.matrix(~.,data=x.now)[,-1,drop=FALSE])
        }else{
          X <- cbind(X,model.matrix(~0+.,data=x.now))
        }
      }else{
        X <- cbind(X,model.matrix(~.,data=x.now)[,-1,drop=FALSE])
      }
    }
     X <- X[,-1,drop = FALSE]
  }else{
    X <- model.matrix(formula, data = data)
    if(ncol(X)>=1){
      if(colnames(X)[1]=="(Intercept)"){
        X <- X[,-1,drop = FALSE]
      }
    }
  }


  x.names <- colnames(X)
  
  ## factorize response vector 
 
  
  ## initialize basic parameters
  k <- apply(Y, 2, function(x){length(levels(as.factor(x)))})
  q <- k-1
  n <- nrow(Y)
  I <- ncol(Y)
  m <- ncol(X)
  n_sigma <- 1
  if(GPCM){
    n_sigma <- I
  }

  ## get final response vector
  if(all(q==1)){
    response <- as.numeric(as.factor(c(t(Y))))-1
  }else{
    response <- lapply(as.data.frame(Y),createResponse)
    response <- c(t(do.call("cbind",response)))
  }
  
  ## total number of parameters to be optimized
  px <- sum(q)+n_sigma
  if(RSM){
    px <- q[1]+I+n_sigma-1
  }
  
  ## scale X and save standard deviations for re-scaling
  # X <- scale(X, center = FALSE)
  X <- scale(X, center = FALSE, scale = apply(X, 2, sd, na.rm = TRUE))
  sd.vec <- attributes(X)$scale
  
  sd.vec <- create.sd.vec(sd.vec, DSF, px, n_sigma, I, q, main.effects) 
  
  ## create design matrix for basic item parameters
  if(!RSM){
    design <- -diag(px-n_sigma)
  }else{
    design <- -cbind(matrix(rep(diag(I),each=q[1]),ncol=I),
               matrix(rep(diag(q[1]),I),ncol=q[1],byrow = TRUE))[,-(I+1)]
  }
  
  ## create design matrix for covariate part
  if(m>=1){
    designX <- -get_designX(X, DSF, m, I, q, n)
    
    ## create penalization matrix
    acoefs <- get_acoefs(RSM, DSF, m, I, q, n_sigma)

        ## update number of parameters to be optimized
    px <- px + ncol(designX)
  }else{
    ## dummy matrix in case no covariates are specified
    designX <- matrix(0,0,0)
    acoefs <- matrix(0,nrow=px,ncol=1)
  }
  
  if(main.effects){
    design.main <- matrix(rep(X, each = sum(q)),ncol = ncol(X))
    designX <- cbind(design.main, designX)
    acoefs <- rbind(matrix(0, ncol = ncol(acoefs), nrow = ncol(design.main)), acoefs)
    px <- px + ncol(X)
  }
  
 
  ret.list <- list(q = q, I = I, m = m, px = px, n = n, response = response,
                   design = design, designX = designX, sd.vec = sd.vec,
                   acoefs = acoefs, n_sigma = n_sigma, x.names = x.names,
                   RSM = RSM, GPCM = GPCM, Y = Y, DSF = DSF)
  
  return(ret.list)
}

create.sd.vec <- function(sd.vec, DSF, px, n_sigma, I, q, main.effects){

  if(DSF){
    new_vec <- rep(sd.vec, sum(q))
  }else{
    new_vec <- rep(sd.vec, I)
  }
  if(main.effects){
    new_vec <- c(sd.vec, new_vec)
  }
  new_vec <- c(rep(1,px-n_sigma),new_vec,rep(1,n_sigma))
  new_vec
}
