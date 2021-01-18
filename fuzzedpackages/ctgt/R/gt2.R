
## the global test
gt2 <- function (y, X, hyps,alpha = 0.05){
  sqrW = sqrt(mean(y)*(1-mean(y)) ) 
  WIHZ = sqrW *(sweep(X,2,colMeans(X))) ## W^{1/2}*(I-H)Z
  IHZ = WIHZ/sqrW ##(I-H)Z
  
  test = sum(colSums(y*IHZ[,hyps,drop=FALSE])^2) 
  lamt = round(eigen(tcrossprod(WIHZ[,hyps,drop=FALSE]),symmetric = TRUE,only.values = TRUE)$values,8)
  crt = criticalvalue(lamt,alpha = alpha)
  pva = pv(test,lamt)
  res = as.numeric(c(format(pva,scientific = TRUE,digits = 3), format(test,scientific = TRUE,digits = 3), 
                     format(crt,scientific = TRUE,digits = 3), format(length(hyps),digits = 0)) )
  names(res) = c("p-value", "Statistic", "Expected", "#Cov")
  return(res)
}


