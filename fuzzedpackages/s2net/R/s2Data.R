# Data interface

zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  # https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

s2Data <- function(xL, yL, xU = NULL, preprocess = T){
  
  # Check data
  n_l = nrow(xL)
  p = ncol(xL)
    
  #NA's
  if(any(is.na(xL)) | any(is.na(xU)) | any(is.na(yL))){
    stop("NA's not supported :( yet :)")
  }  
  if(!is.null(xU)){
    if(ncol(xU)!=p){
      stop("xL and xU have different number of columns")
    }
  }
  
  if(is.data.frame(xL)){
    if (!isFALSE(preprocess)) {
      X = rbind(xL, xU)
      X = stats::model.matrix(~., X)[,-1]
      xL = X[1:n_l, ]
      if(!is.null(xU)) xU = X[-(1:n_l), ]
    }
  }
  
  type = "regression"
  base = 0
  
  if(is.factor(yL)){
    if(length(levels(yL))!=2){
      stop("If yL is a factor, it should have exactly two levels")
    }
    y = as.numeric(yL) - 1
    type = "classification"
    base = levels(yL)[1]
    yL = y
  }
  
  if(!is.numeric(yL)) stop("yL is not numeric or factor")
  # Explicit convertion to matrix type
  xL = as.matrix(xL)
  yL = as.matrix(yL)
  if(!is.null(xU)) xU = as.matrix(xU)
  
  if(ncol(yL)!=1){
    stop("yL has more than one column. Multivariate response not supported :( yet :)")
  }
  # size
  if(nrow(yL)!=n_l){
    stop("xL and yL have different number of rows")
  }
  
  if(isTRUE(preprocess)){
    
    #remove constant columns 
    rm_cols = apply(xL, 2, zero_range)
    xL = xL[, !rm_cols]
    if(!is.null(xU)) xU = xU[, !rm_cols]
    
    # Scale, with respect to xL...
    xL = scale(xL)
    s_scale = attr(xL, "scaled:scale")
    s_center = attr(xL, "scaled:center")
    if(!is.null(xU)) xU = scale(xU, center = s_center, scale = s_scale)
    
    # center yL too
    if(type == "regression"){
      yL = scale(yL, center = T, scale = F)
      y_center = attr(yL, "scaled:center")
    }else{
      y_center = 0
    }
    
  }else{
    if(class(preprocess) == "s2Data"){
      rm_cols = attr(preprocess, "pr:rm_cols")
      xL = xL[, !rm_cols]
      if(!is.null(xU)) xU = xU[, !rm_cols]
      
      s_scale = attr(preprocess,  "pr:scale")
      s_center = attr(preprocess,  "pr:center")
      
      xL = scale(xL, center = s_center, scale = s_scale)
      if(!is.null(xU)) xU = scale(xU, center = s_center, scale = s_scale)
      
      y_center = attr(preprocess, "pr:ycenter")
      #y_scale =  attr(preprocess, "pr:yscale")
      
      if(type == "regression"){yL = scale(yL, center = y_center, scale = F)}
    }else{
      rm_cols = rep(FALSE, ncol(xL))
      s_scale = rep(1, ncol(xL))
      s_center = rep(0, ncol(xL))
      y_center = 0
    }
  }
  # s2P = list(
  #   rm_cols = rm_cols,
  #   center = s_center,
  #   scale = s_scale
  # )
  # class(s2P) = "s2Pr"
  
  s2D = list(
    xL = xL,
    yL = yL,
    xU = xU,
    type = type,
    base = base
  )
  class(s2D) = "s2Data"
  attr(s2D, "pr:rm_cols") = unname(rm_cols)
  attr(s2D, "pr:center") = unname(s_center)
  attr(s2D, "pr:scale") = unname(s_scale)
  attr(s2D, "pr:ycenter") = unname(y_center)
  
  return(s2D)
}

print.s2Data <- function(x, ...){
  plog("s2Data frame:")
  plog("Labeled data: ", nrow(x$xL), " ", ncol(x$xL))
  plog("Unlabeled data: ", nrow(x$xU), " ", ncol(x$xU))
  plog("Task ", x$type)
}

plog <- function(text,...){
  cat(paste0(text,..., "\n"))
}
