##################################################################
# bgtWidth() function                                            #
##################################################################

"bgtWidth" <-
  function(n, s, p, conf.level=0.95, alternative="two.sided", 
           method="CP"){
    
    if( any(n<=3) )
    {stop("the number of groups n allowed in calculations must be integers greater than 1")}
    
    if( any(s<1) ){stop("group size s must be specified as integers > 0")}
    
    if( length(conf.level)!=1 || conf.level<0 || conf.level>1)
    {stop("conf.level must be a positive number between 0 and 1")}
    
    if( length(p)!=1 || p>1 || p<0)
    {stop("true proportion p must be specified as a single number between 0 and 1")}
    
    method<-match.arg(method, choices=c("CP","Blaker","AC","score","Wald","soc"))
    
    alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))
    
    # calculations:
    
    matnsp <- cbind(n,s,p)
    matnsp <- cbind("ns"=matnsp[,1]*matnsp[,2], matnsp)
    power <- numeric(length=nrow(matnsp))
    bias <- numeric(length=nrow(matnsp))
    
    expCIwidth<-numeric(length=nrow(matnsp))
    
    for (i in 1:length(expCIwidth))
    {
      expCIwidth[i]<-bgtWidthI(n=matnsp[[i,2]], s=matnsp[[i,3]], p=matnsp[[i,4]], conf.level=conf.level, alternative=alternative, method=method)$expCIWidth
    }
    
    return(as.matrix(cbind(matnsp,expCIwidth)))
  }



# Brianna Hitt - 02-13-2020
# Changed class from "binWidth" to "gtWidth"

"bgtWidthI" <-
  function(n, s, p, conf.level=0.95, alternative="two.sided", method="CP")
  {
    
    # indicator function for the CI length at a special event
    # in one sided case: length is defined as absolute difference between estimator and confidence bound
    
    L.Ind<-function(y, n, s, p, conf.level, alternative, method)
      
    {
      
      if(method=="Wald"){int=bgtWald(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
      if(method=="Wilson"){int=bgtWilson(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
      if(method=="AC"){int=bgtAC(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
      if(method=="soc"){int=bgtSOC(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
      if(method=="CP"){int=bgtCP(y=y, n=n, s=s, conf.level=conf.level, alternative=alternative)}
      if(method=="Blaker"){int=bgtBlaker(y=y, n=n, s=s, conf.level=conf.level)}
      
      if(alternative=="less")
      {CIlength <- int[[2]]-p}
      
      if(alternative=="greater")
      {CIlength <- p-int[[1]]}
      
      if(alternative=="two.sided")
      {CIlength <- int[[2]]-int[[1]]}
      CIlength
    }
    
    # Probability of a single event, the group testing density:
    
    bgt.prob<-function(n,y,s,p.tr)
    {
      theta<-1-(1-p.tr)^s
      dbinom(x=y,size=n, prob=theta)
    }
    
    
    #  calculate this for all possible events: 
    
    yvec<-0:n
    
    Lvec<-numeric(length=length(yvec))   
    probvec<-numeric(length=length(yvec))
    
    for(i in 1:length(yvec))
    {Lvec[i]<-L.Ind(y=yvec[i], n=n, s=s, p=p, conf.level=conf.level, alternative=alternative, method=method)
    probvec[i]<-bgt.prob(y=yvec[i], n=n, s=s, p.tr=p)
    }
    expCILength=sum(Lvec * probvec)
    
    # E(X)= sum(Xi * prob(Xi))
    
    out<-list(expCIWidth=expCILength, alternative=alternative, n=n,s=s, p=p)
    
    class(out)<-"gtWidth"
    return(out)
  }




##################################################################
# gtWidth() function                                             #
##################################################################

#' @title Expected width of confidence intervals in group testing
#' 
#' @description Calculation of the expected value of the width of 
#' confidence intervals for one proportion in group testing. Calculations 
#' are available for the confidence interval methods in \code{\link{propCI}}.
#' 
#' @param n integer specifying the number of groups. A vector of 
#' integers is also allowed.
#' @param s integer specifying the common size of groups. A vector 
#' of integers is also allowed.
#' @param p the assumed true proportion of individuals showing 
#' the trait to be estimated. A vector is also allowed.
#' @param conf.level the required confidence level of the interval.
#' @param alternative character string specifying the alternative 
#' hypothesis, either \kbd{"two.sided"}, \kbd{"less"}, or \kbd{"greater"}.
#' @param method character string specifying the confidence 
#' interval method. Available options include those in 
#' \code{\link{propCI}}.
#' 
#' @details The two-sided (\kbd{method="two.sided"}) option calculates the 
#' expected width between the lower and upper bound of a two-sided 
#' \eqn{conf.level*100} percent confidence interval. See Tebbs & Bilder (2004) 
#' for expression. The one-sided (\kbd{method="less"} or \kbd{method="greater"}) 
#' options calculate the expected distance between the one-sided limit and the 
#' assumed true proportion \kbd{p} for a one-sided \eqn{conf.level*100} percent 
#' confidence interval.
#' 
#' @return A matrix containing the columns:
#' \item{ns}{the resulting total number of units, \eqn{n*s}.}
#' \item{n}{the number of groups.}
#' \item{s}{the group size.}
#' \item{p}{the assumed true proportion.}
#' \item{expCIWidth}{the expected value of the confidence 
#' interval width as defined under the argument \kbd{alternative}.}
#' 
#' @author This function was originally written as \code{bgtWidth} by Frank 
#' Schaarschmidt for the \code{binGroup} package. Minor modifications have 
#' been made for inclusion of the function in the \code{binGroup2} package.
#' 
#' @references 
#' \insertRef{Tebbs2004}{binGroup2}
#' 
#' @seealso \code{\link{propCI}} for confidence intervals in 
#' group testing.
#' 
#' @family estimation functions
#' 
#' @examples 
#' # Examine different group sizes to determine 
#' #   the shortest expected width.
#' gtWidth(n=20, s=seq(from=1, to=200, by=10),
#'         p=0.01, alternative="less", method="CP")
#' 
#' # Calculate the expected width of the confidence  
#' #   interval with a group size of 1 (individual testing).
#' gtWidth(n=20, s=1, p=0.005, alternative="less", method="CP")

gtWidth <- bgtWidth

#
