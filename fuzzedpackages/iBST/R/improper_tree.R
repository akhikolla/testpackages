improper_tree <- function(xdata, Y.names, P.names, T.names, method = "R2", args.rpart){
  methodSurvR2 = list(eval = .surveval, split = .survsplitR2, init = .survinit)
  methodSurvLR = list(eval = .surveval, split = .survsplitLR, init = .survinit)
  if(method == "R2") methodsurv = methodSurvR2 ;
  if(method == "LR") methodsurv = methodSurvLR ;
  YSTRATE <- exp(xdata[, P.names]) ;
  tree_surv <- rpart(as.formula(paste(paste('Surv(', paste(Y.names, collapse = ','),')', sep = ''), '~', paste(T.names, collapse = '+'), sep = '')), 
                     data = xdata, method = methodsurv,  control = args.rpart, weights = YSTRATE, y = FALSE ) ;
  return(tree_surv)
}