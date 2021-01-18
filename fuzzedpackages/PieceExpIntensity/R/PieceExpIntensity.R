#'  Runs the PieceExpIntensity sampler and returns posterior results.
#'
#'  Returns a list of posterior samples along with summaries for the most visited number of split points.
#' @importFrom Rcpp evalCpp
#' @importFrom stats quantile
#' @param X   Vector containing observed event times.
#' @param Y Vector containing poisson count intensities.
#' @param B Number of iterations to run the MCMC with half burned in.
#' @param Poi Prior mean number of split points.
#' @return A list of all posterior quantities and a summary of the most commonly visited model.
#' @references
#' Chapple (2017). Modeling ISIL terror attacks and their intensities via flexible Bayesian piecewise models. Currently Under Review.
#' @examples
#' B=1000
#' n=100
#' X=rexp(n,1)
#' Y=X
#' Y[X<.5]=rpois(sum(X<.5),20)
#' Y[X>.5]=rpois(sum(X>.5),3)
#' Poi=10
#' PieceExpIntensity(X,Y,B,Poi)
#' @export
PieceExpIntensity=function(X,Y,B,Poi){

  B1=B/2

  if(sum(Y%%1==0)==length(Y)){

    if(min(X)<0){
 cat("Some event times are less than 0, please fix

     ")

      return(NULL)
    }else{

      cat("Ok Lets Go")



G1=PieceExpIntensity2(X,Y,B,Poi)

##Now Get the mode and return some stuff

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

###Now find samples for the two most visited number of split points
NumSplit=Mode(G1[[4]])
J=NumSplit+2
J1=J-1

cat("Posterior Mode of Split Points:",NumSplit,"

    ")

s1=matrix(rep(NA,J*B1),ncol=J)
lam1=matrix(rep(NA,J1*B1),ncol=J1)
rate1=matrix(rep(NA,J1*B1),ncol=J1)


for(b in 1:B1){
  if(G1[[4]][b]==NumSplit){
    s1[b,]=G1[[1]][b,1:J]
    lam1[b,]=G1[[2]][b,1:J1]
    rate1[b,]=G1[[3]][b,1:J1]


  }


}



cat("Posterior Mean Location of Split Points (and Credible Interval) for J =",NumSplit,"

    ")


print(colMeans(s1,na.rm=TRUE))


print(apply(s1,2,function(z){quantile(z,probs=c(.025,.975),na.rm=TRUE)}))



cat("Posterior Mean Log-Hazards (and Credible Interval) for J =",NumSplit,"

    ")

##How did rates change
print(colMeans(lam1,na.rm=TRUE))

print(apply(lam1,2,function(z){quantile(z,probs=c(.025,.975),na.rm=TRUE)}))


##How much did they happen

cat("Posterior Mean Intensity Rates (and Credible Interval) for J =",NumSplit,"

    ")

print(colMeans(exp(rate1),na.rm=TRUE))

print(apply(exp(rate1),2,function(z){quantile(z,probs=c(.025,.975),na.rm=TRUE)}))



return(G1)





}


  }else{

    cat("Y does not contain count data

        ")
    return(NULL)
}


}
