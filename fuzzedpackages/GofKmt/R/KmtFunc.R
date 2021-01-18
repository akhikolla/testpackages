#' Obtain Critical Value Table 
#'
#' Obtain critical values of the Khmaladze martingale transformation test for various significane levels and sample sizes
#'@param strDist  the name of the null distribution for the hypothesis test: Normal, Cauchy, or Logistic. Other distributions such as Gumbel, Weibull and Frechet will be available in later versions. 
#'@param Modified a logical value which specifies whether or not to use the modeifed version of the test: False calls the original version while True calls the modified version.
#'@return A 10-by-6 table of critical values for various significance levles (0.1, 0.075, 0.05, 0.025, 0.01) and sample sizes (10,20,...,100). 
#'
#'@examples
#'
#'
#'
#'### Critical values of the original test for a normal distribution
#'
#'
#'strDist = "Normal"
#'Modified=FALSE
#'CritValue = GetCV(strDist, Modified)
#'
#'
#'
#'## Critical values of the modified test for the logistic distribution
#'strDist = "Logistic"
#'Modified=TRUE
#'CritValue = GetCV(strDist, Modified)
#'
#'
#'
#'##Critical values of the modified test for the Cauchy distribution
#'strDist = "Cauchy"
#'Modified=TRUE
#'CritValue = GetCV(strDist, Modified)
#'
#'@export
#'@seealso KhmaladzeTrans()
#'@useDynLib GofKmt


GetCV = function(strDist, Modified){
  
  
  if(Modified==FALSE){
    strName = "Original"
  }else{
    strName = "Modified"
  }
  
  if(strDist=="Normal"){
    strDist="Normal"
  }else if(strDist=="Logistic"){
    strDist="Logistic"
  }else if(strDist=="Cauchy"){
    strDist="Cauchy"
  }else{
    strDist="Normal"
  }
  
  targetstr = paste("Critical.Value.for.", strName, "Test.", strDist, sep="")
  
  print(targetstr)
  CritMat = Tables[[targetstr]]
  return(CritMat)
}




#' Implementing Khmaladze Martingale Transformation.
#'
#' Performs goodness-of-fit test through Khmaladze matringale transformation
#'@param X  a random sample of n observations
#'@param Modified a logical value which specifies whether or not to use the modeifed version of the test: False calls the original version while True calls the modified version.  
#'@param strDist  the name of the null distribution for the hypothesis test: Normal, Cauchy, or Logistic. Other distributions such as Gumbel, Weibull and Frechet will be available in later versions. 
#'@param bGraph  a logical value which specifies whether or not to get the graph of the objective function of the martingale transformation.
#'@param nNum  the number of ticks on each segmented interval when drawing the graph of the objective function. The default is 10. Bigger value will result in a smoother graph.
#'@return A list of the following values: 
#'\describe{
#'\item{opt.x}{ 
#'\itemize{
#'  \item When Modified is False, opt.x is the value of x where the optimum of the objective function - which is also the test statistic - occurs.
#'  \item When Modified is True, opt.x is the vector of the value of x's where the infimum and supremum of \eqn{U_{n}} occur. 
#'}
#'} 
#'
#'\item{test.stat}{ 
#'\itemize{
#'  \item When Modified is False, test.stat is the test statistic obtained through Khmaladze martingale transformation.
#'  \item When Modified is True, test.stat is the vector of the supremum of \eqn{U_{n}}, the infimum of \eqn{U_{n}}, and the difference of them. 
#'}
#'}
#'
#'\item{graph.data}{ a data frame which includes the information of the objective function.}
#'
#'\item{graph}{ a ggplot object which includes the graph of the objective function.}
#'
#'\item{intervals}{ a list of segmented intervals over which the graph of the objective function is defined.}
#'
#'\item{mu}{ the point estimate for the location parameter mu}
#'
#'\item{sigma}{ the point estimate for the scale parameter sigma}
#'}
#'
#'@examples
#'####################
#'n = 10
#'X = rnorm(n, 1,3)    # Generate a random sample of n observations from N(1,3)
#'strDist = "Normal"
#'Modified=FALSE
#'lResult = KhmaladzeTrans(X, Modified, strDist, bGraph=TRUE, nNum=10)
#'KMT_OptimalX = lResult$opt.x
#'KMT_TestStat = lResult$test.stat
#'KMT_DM = lResult$graph.data
#'KMT_Graph = lResult$graph
#'
#'#### Draw the graph of the objective function
#'KMT_Graph
#'
#'KMT_Intervals = lResult$intervals
#'KMT_Muhat = lResult$mu
#'KMT_Sigmahat = lResult$sigma
#'
#'
#'#####################
#'
#'#####################
#'n = 10
#'X = rlogis(n, 1,2)  # Generate a random sample of n observations from the logistic distribution
#'strDist = "Logistic"
#'Modified=TRUE
#'lResult = KhmaladzeTrans(X, Modified, strDist, bGraph=TRUE, nNum=10)
#'KMT_Optimal_Positive_X = lResult$opt.x[1]
#'KMT_Optimal_Negative_X = lResult$opt.x[2]
#'KMT_Postive_TestStat = lResult$test.stat[1]
#'KMT_Negative_TestStat = lResult$test.stat[2]
#'KMT_TestStat = lResult$test.stat[3]
#'KMT_DM = lResult$graph.data
#'KMT_Graph = lResult$graph
#'
#'#### Draw the graph of the objective function
#'KMT_Graph
#'
#'KMT_Intervals = lResult$intervals
#'KMT_Muhat = lResult$mu
#'KMT_Sigmahat = lResult$sigma
#'#####################
#'
#'
#'
#'
#'
#'#####################
#'n = 10
#'X = rcauchy(n, 0,1)  # Generate a random sample of n observations from Cauchy distribution
#'strDist = "Cauchy"
#'Modified=FALSE
#'lResult = KhmaladzeTrans(X, Modified, strDist, bGraph=TRUE, nNum=10)
#'KMT_OptimalX = lResult$opt.x
#'KMT_TestStat = lResult$test.stat
#'KMT_DM = lResult$graph.data
#'KMT_Graph = lResult$graph
#'
#'#### Draw the graph of the objective function
#'KMT_Graph
#'
#'
#'KMT_Intervals = lResult$intervals
#'KMT_Muhat = lResult$mu
#'KMT_Sigmahat = lResult$sigma
#'#####################





#'@references
#'[1] Khmaladze, E.V., Koul, H.L. (2004). Martingale transforms goodness-of-fit tests in regression models. Ann. Statist., 32. 995-1034
#'@references
#'[2] E.V. Khmaladze, H.L. Koul (2009). Goodness-of-fit problem for errors in nonparametric regression: distribution free approach. Ann. Statist., 37(6A) 3165-3185.
#'@references
#'[3] Kim, Jiwoong (2020). Implementation of a goodness-of-fit test through Khmaladze martingale transformation. Comp. Stat., 35(4): 1993-2017    
#'@export
#'@importFrom Rcpp evalCpp
#'@useDynLib GofKmt


KhmaladzeTrans = function(X, Modified=FALSE, strDist, bGraph=FALSE, nNum=10){
  

  fast.estimate = FALSE
  n = length(X)
  
  if(strDist == "Normal"){
    muhat = mean(X)
    sigmahat = sd(X)
    
  }else if(strDist == "Logistic"){
    
    if(fast.estimate == TRUE){
      muhat = median(X)
      sigmahat = sqrt(3)/pi*sd(X)
      
    }else{
      x01 = median(X)
      x02 = sqrt(3)/pi*sd(X)
      
      startVal = c(x01, x02)
      mle_method = optim(startVal, mleLogis(X))
      mle_par = mle_method$par
      
      muhat = mle_par[1]
      sigmahat = mle_par[2]
    }
    
    
  }else if(strDist == "Cauchy"){

    if(fast.estimate == TRUE){
      
      muhat = median(X)
      sighat = GetCauchyScale(X, muhat)
      
    }else{
      x01 = median(X)
      
      quan3 = quantile(X, 0.75)
      quan1 = quantile(X, 0.25)
      
      x02 = (quan3[[1]]-quan1[[1]])/2
      
      ####################
      x0 = c(x01, x02)
      LBV = c(-Inf, 0)
      UBV = c(Inf, Inf)
      
      MRe = solnp( pars=x0, fun=MLECauchyObjFunction(X), eqfun=NULL, eqB=NULL, ineqfun=NULL, ineqLB = NULL,
                   ineqUB = NULL, LB = LBV, UB = UBV)
      
      CauchyMle = MRe$pars
      
      muhat = CauchyMle[1]
      sigmahat = CauchyMle[2]
      
    }
      
    
  }else{
    message("Name of null distribution is not valid.")
    stop()
  }
  
  StandX = rep(0, times=n)
  
  for(i in 1:n){
    StandX[i] = (X[i]-muhat)/sigmahat
    
  }
  
  StandX = sort(StandX)
  
  NormalMat = Tables$Integration.Table.Normal
  LogisMat = Tables$Integration.Table.Logistic1
  ReMat = Tables$Integration.Table.Logistic2
  CauchyMat = Tables$Integration.Table.Cauchy
  
  

  if(Modified==FALSE){
    bModified=0
  }else{
    bModified=1
  }
  
  lst = KmtMain(StandX, bModified, NormalMat, LogisMat, ReMat, CauchyMat, strDist, bGraph, nNum)
  
  lineVec = lst[[3]]
  gVec = lst[[4]]
  nLen = length(lineVec)
  
  DM = matrix(0, nLen, 3)
  DM = data.frame(DM)
  
  lineVec = round(lineVec, digits=3)
  gVec=round(gVec, digits=3)
  
  colnames(DM, do.NULL = TRUE)
  colnames(DM) = c("x", "y", "Interval")
  
  DM[, "x"] = lineVec
  DM[, "y"] = gVec
  
  nGroup = n+1
  
  xIntVec = rep(0, times=n)
  
  for(i in 1:nGroup){
    
    Sn = (i-1)*nNum+1
    En = i*nNum
    
    if(i==1){
      Sval = "-Inf"
    }else{
      Sval = DM[Sn, "x"]
    }
    
    if(i==nGroup){
      Eval = "Inf"
    }else{
      Eval = DM[En+1, "x"]
      xIntVec[i] = Eval
    }
    
    Str = paste("", i, ". [", Sval, ", ", Eval, ")", sep="")
    print(Str)
    DM[Sn:En, "Interval"] = Str
  }
  
  DM[, "Interval"] = as.factor(DM[, "Interval"])
  
  Intervals = levels(DM[, "Interval"])
  
  gObj = ggplot(data=DM, mapping=aes(x=DM[,1], y=DM[,2], group = DM[,3]))+
    geom_line(color="red")+
    labs(title = "Graph of the Objective Function", x="x", y="y")
  
  
  if(Modified==FALSE){
    lResult = list(opt.x = lst[[1]], test.stat=lst[[2]], 
                   graph.data = DM, graph = gObj, intervals=Intervals, mu = muhat, sigma = sigmahat)
  }else{
    lResult = list(opt.x = c(lst[[5]], lst[[7]]), test.stat=c(lst[[6]], lst[[8]], (lst[[6]]-lst[[8]])), 
                   graph.data = DM, graph = gObj, intervals=Intervals, mu = muhat, sigma = sigmahat)
  }
  

  return(lResult)
}


GetCauchyScale = function(XVec, mu){
  
  al = 0.5
  quant = quantile(XVec, al )
  Qval = quant[[1]]
  
  sig = (al-mu)/( tan( pi*( pcauchy(Qval, 0,1) - 0.5 ) ) )
  
  return( abs(sig)+0.0001)
}


MLECauchyObjFunction = function(XVec){
  nLength = length(XVec)
  
  Dual = function(x){
    tempsum = -nLength*log(pi*x[2])
    for(i in 1:nLength){
      tempsum = tempsum - log(1+( (XVec[i]-x[1] ) / x[2]  )^2)
    }
    return(-tempsum)
  }
  return(Dual)
}



mleLogis = function(X){
  
  nLeng = length(X)
  Dual = function(param){
    tempsum=0
    for (i in 1:nLeng){
      normx = (X[i]-param[1])/param[2]
      tempsum = tempsum + normx - log(param[2]) -2*log(1+exp(normx))
    }
    return(-tempsum)
  }
  return(Dual)
  
}












