### R package depcoeff

#' @import copula

#### function psi
psi<- function(x,id,par=-1){
  h<- abs(x)
  if (par<=0) {
    if (id==3) {par<- 1.5}
    if (id==4) {par<- 0.5}}
  if (id==1){return(x^2)}   # Spearman-coefficient
  if (id==2){return(h)}  # Spearman's footrule
  if (id==3){return(h^par)}  # power coefficient
  if (id==4){if (h<par) {return(0.5*h*h)} else {return(par*(h-0.5*par))}}  # Huberfunktion
}
#### normalisation factor psi_bar
psif<- function(id,par=-1){
  if (par<=0) {
    if (id==3) {par<- 1.5}
    if (id==4) {par<- 0.5}}
  if ((id==4)&(par>=1.0)) {par<- 0.5}
  if (id==1){return(0.166666666667)}
  if (id==2){return(0.333333333333)}
  if (id==3){return(2.0/(par*(par+3.0)+2.0))}
  if (id==4){return(par*(2.0-par)*(par*(par-2.0)+2.0)/12.0)}
}
### minimal value of zeta for the specific coefficients
zetamin<- function(id,par=-1){
  if (par<=0) {
    if (id==3) {par<- 1.5}
    if (id==4) {par<- 0.5}}
  if ((id==4)&(par>=1.0)) {par<- 0.5}
  if (id==1){return(-1.0)}
  if (id==2){return(-0.5)}
  if (id==3){return(-par/2.0)}
  if (id==4){return((par^3-2.0*par^2+2.0)/(par^3-4.0*par^2+6.0*par-4.0))}
}

#' Zeta dependence coefficient
#'
#' zetac is a function to evaluate the zeta dependence coefficient (one interval)
#' of two random variables x and y which is based on the copula. Four specific coefficients
#' are available: the Spearman coefficient, Spearman's footrule, the power coefficient
#' and the Huber function coefficient.
#'
#' Let \eqn{X_{1},\ldots ,X_{n}} be the sample of the \eqn{X} variable. Formulas
#' for the estimators of values \eqn{F(X_{i})} of the distribution function:
#'   methodF = 1 \eqn{\rightarrow \hat{F}(X_{i})=\frac{1}{n}\textrm{rank}(X_{i})}
#'   methodF = 2 \eqn{\rightarrow \hat{F}^{1}(X_{i})=\frac{1}{n+1}\textrm{rank}(X_{i})}
#'   methodF = 3 \eqn{\rightarrow \hat{F}^{2}(X_{i})=\frac{1}{\sqrt{n^{2}-1}}\textrm{rank}(X_{i})}
#' The values of the distribution function of \eqn{Y} are treated analogously.
#' @usage zetac(x,y,method="Spearman",methodF=1,parH=0.5,parp=1.5)
#' @param x,y data vectors of the two variables whose dependence is analysed.
#' @param method list of names of the coefficients: "Spearman" stands for the
#' Spearman coefficient, "footrule" means Spearman's footrule, "power" stands
#' for the power function coefficient, "Huber" means the Huber function
#' coefficient. If "all" is assigned to method then all methods are used.
#' @param methodF value 1,2 or 3 refers to several methods for computation of
#' the distribution function values, 1 is the default value.
#' @param parH parameter of the Huber function (default 0.5). Valid values for
#' parH are between 0 and 1.
#' @param parp parameter of the power function (default 1.5). The parameter has
#' to be positive.
#' @return zeta dependence coefficient of two random variables. This coefficient
#' is bounded by 1. The higher the value the stronger is the dependence.
#' @references Eckhard Liebscher (2014). Copula-based dependence measures. Dependence Modeling 2 (2014), 49-64
#' @examples
#' library(MASS)
#' data<- gilgais
#' zetac(data[,1],data[,2])
#' @export

zetac<- function(x,y,method="Spearman",methodF=1,parH=0.5,parp=1.5){
  n<- length(x)
  if (n!=length(y)){ stop("number of sample items in x and y are different")} #check the parameters
  if (n<4) {stop("not enough sample items")}
  if ((!is.vector(x))|(!is.vector(y))){ stop("data format wrong")}
  if (!methodF %in% 1:3) methodF<- 1

    ### compute normed difference of ranks
  if (methodF==1) {
    uu<- (rank(x)-rank(y))/n
  }else{
    if (methodF==2){
    uu<- (rank(x)-rank(y))/(n+1.0)} else {
    uu<- (rank(x)-rank(y))/sqrt(n^2-1.0)
    }}

  if (method=="all") {method<- c("Spearman","footrule","power","Huber")}
  ### Spearman coefficient
  if ("Spearman" %in% method) {
    s<- 1.0-sum(psi(uu,1))/(n*psif(1))
    outl<- list(Spearman=s)
  } else {outl<- list()}
  ### Spearman's footrule
  if ("footrule" %in% method) {
    s<- 1.0-sum(psi(uu,2))/(n*psif(2))
    outl<- append(outl,list(footrule=s))}
  ### power function coefficient
  if ("power" %in% method) {
    s<- 1.0-sum(psi(uu,3,parp))/(n*psif(3,parp))
    outl<- append(outl,list(power=s))}
  ### Huber function coefficient
  if ("Huber" %in% method) {
    s<- 0.0
    for (i in 1:(length(uu))){
      s<- s+ psi(uu[i],4,parH)
    }
    s<- 1.0-s/(n*psif(4,parH))
    outl<- append(outl,list(Huber=s))}

  return(outl)}  ## output: list of computed coefficients,

#' Zeta dependence coefficient of piecewise monotonicity
#'
#' zetapm is a function to evaluate the zeta dependence coefficients of piecewise
#' monotonicity of two random variables x and y which is based on the copula. The regressor
#' domain (domain of x) is split into two parts. The function searches
#' for the optimal splitting point to obtain maximum
#' depedence. The main part of the function is coded as C++ procedure
#'
#' Let \eqn{X_{1},\ldots ,X_{n}} be the sample of the \eqn{X} variable. Formulas
#' for the estimators of values \eqn{F(X_{i})} of the distribution function:
#'   methodF = 1 \eqn{\rightarrow \hat{F}(X_{i})=\frac{1}{n}\textrm{rank}(X_{i})}
#'   methodF = 2 \eqn{\rightarrow \hat{F}^{1}(X_{i})=\frac{1}{n+1}\textrm{rank}(X_{i})}
#'   methodF = 3 \eqn{\rightarrow \hat{F}^{2}(X_{i})=\frac{1}{\sqrt{n^{2}-1}}\textrm{rank}(X_{i})}
#' The values of the distribution function of \eqn{Y} are treated analogously.
#' @usage zetapm(x,y,amin=0.25,method="all",methodF=1,parp=1.5,parH=0.5)
#' @param x,y data vectors of the two variables whose dependence is analysed.
#' @param amin minimum fraction of sample items to be used for one split region
#' @param method vector of chosen special coefficients:
#'   Spearman...Spearman coefficient
#'   footrule...Spearman's footrule
#'   power...power coefficient
#'   Huber...Huber function coefficient,
#'   "all" refers to all coefficients
#' @param methodF value 1,2 or 3 refers to several methods for computation of the
#' distribution function values, 1 is the default value.
#' @param parH parameter of the Huber function (default 0.5). Valid values for parH
#' are between 0 and 1.
#' @param parp parameter of the power function (default 1.5). The parameter has
#' to be positive.
#' @return list of zeta dependence coefficients (plusminus coefficient and minusplus one)
#' of piecewise monotonicity of two random variables containing the
#' following elements or a subset of it in this order:
#' Spearman coefficient, footrule, power coefficient, Huber function coefficient.
#' position1 and position2 indicate the number of the sample items where the optimized
#' split point is located
#' @references Eckhard Liebscher (2017). Copula-based dependence measures for piecewise
#' monotonicity. Dependence Modeling 5 (2017), 198-220
#' @examples
#' library(MASS)
#' data<- gilgais
#' zetapm(data[,1],data[,2])
#' @export

zetapm<- function(x,y,amin=0.25,method="all",methodF=1,parp=1.5,parH=0.5){
  n<- length(x)
  if (n!=length(y)){ stop("number of sample items in x and y are different")} #check the parameters
  if ((parp<=0.0)|(parH<=0.0)|(parH>=1.0)){ stop("parameter(s) not valid")}
  if (!methodF %in% 1:3) methodF<- 1
  if (n<8) {  stop("not enough sample items")}
  if ((!is.vector(x))|(!is.vector(y))){ stop("data format wrong")}
  if ((amin>=0.5)|(amin<=0.0)) {amin<- 0.25}  ## amin--> default 0.25
  method0<- c("Spearman","footrule","power","Huber")    #vector of all methods
  if (length(method)==1) {if (method=="all") {method<- method0}}
  mv<- which(method == method0)  # vector with numbers of chosen methods

  u<- rank(x)
  u1<-u2<-v1<-v2<- vector()
  na<- max(ceiling(amin*n+0.5),4)  #number of items in region of minimum length
  ne<- n-na  #length of u1,u2
  # data for the left-side part
  u1<- u[u<=ne]  # x ranks of the left part
  v1<- rank(y[u<=ne])  # y ranks of the left part
  v3<- ne+1.0-v1   # y ranks for the minus-plus coefficient

  # data for the right-side part
  u2<- u[u>na]-na  # x ranks of the right part
  v2<- rank(y[u>na])  # y ranks of the leftright part
  v4<- ne+1.0-v2  # y ranks for the minus-plus coefficient

  ## plus-minus-coefficient
  h1<- coeffpml(u1,v1,u2,v2,amin,parp,parH,n,na,methodF)

   # data for the left-side part
  u1<- u[u<=ne]  # x ranks of the left part

  # data for the right-side part
  u2<- u[u>na]-na  # x ranks of the right part

  ### minus-plus coefficient
  h2<- coeffpml(u1,v3,u2,v4,amin,parp,parH,n,na,methodF)

  ###output of coefficients: Spearman, Spearman's footrule, power function, Huber function
  outl<- list(h1[mv],h1[(mv+4)],h2[mv],h2[(mv+4)])
  names(outl)<- c("plusminuscoeff","position1","minuspluscoeff","position2")
  return(outl)
}

#' Zeta coefficient of piecewise monotonicity with split domain
#'
#' The function zetaci evaluates the coefficient of piecewise monotonicity of variables x and y where the x-domain
#' is split into a fixed number of intervals.
#'
#' Let \eqn{X_{1},\ldots ,X_{n}} be the sample of the \eqn{X} variable. Formulas
#' for the estimators of values \eqn{F(X_{i})} of the distribution function:
#'   methodF = 1 \eqn{\rightarrow \hat{F}(X_{i})=\frac{1}{n}\textrm{rank}(X_{i})}
#'   methodF = 2 \eqn{\rightarrow \hat{F}^{1}(X_{i})=\frac{1}{n+1}\textrm{rank}(X_{i})}
#'   methodF = 3 \eqn{\rightarrow \hat{F}^{2}(X_{i})=\frac{1}{\sqrt{n^{2}-1}}\textrm{rank}(X_{i})}
#' The values of the distribution function of \eqn{Y} are treated analogously.
#' @usage zetaci(x,y,a,method="Spearman",methodF=1,parH=0.5,parp=1.5)
#' @param x,y data vectors of the two variables whose dependence is analysed.
#' @param a vector of fractions \eqn{a_{i},0<a_{i}<a_{i+1}<1} for the splitting. A fraction
#' of \eqn{a_{1},a_{2}-a_{1},a_{3}-a{2}}... of data points are in the corresponding split region.
#' The number of split regions is equal to the length of \eqn{a} plus 1.
#' @param method value (default "Spearman")
#' @param methodF value 1,2 or 3 refers to several methods for computation of
#' the distribution function values, 1 is the default value.
#' @param parH parameter of the Huber function (default 0.5). Valid values for
#' parH are between 0 and 1.
#' @param parp parameter of the power function (default 1.5). The parameter has
#' to be positive.
#' @return list of zeta dependence coefficients of piecewise monotonicity of two
#' random variables containing the following elements:
#'   Spearman...Spearman coefficient
#'   footrule...Spearman's footrule
#'   power...power coefficient
#'   Huber...Huber function coefficient
#' @references Eckhard Liebscher (2017). Copula-based dependence measures for piecewise
#' monotonicity. Dependence Modeling 5 (2017), 198-220
#' @examples
#' library(MASS)
#' data<- gilgais
#' zetaci(data[, 1], data[, 2], a=c(0.25, 0.5, 0.75))
#' @export


zetaci<- function(x,y,a,method="Spearman",methodF=1,parH=0.5,parp=1.5){
  if ((parp<=0.0)|(parH<=0.0)|(parH>=1.0)|(min(a)<=0.0)|(max(a)>=1.0)){
    stop("parameter(s) not valid")}  #check the parameters
  if (length(method)==1) { if (method=="all") {method<- c("Spearman","footrule","power","Huber")}}
  n<- length(x)
  if (n!=length(y)){ stop("number of sample items in x and y are different")}
  if ((!is.vector(x))|(!is.vector(y))){ stop("data format wrong")}
  if (n<4) {  stop("not enough sample items")}
  if (!methodF %in% 1:3) methodF<- 1

  m<- length(a)+1
  a[m]<- 1.0
  A<- vector()   # vector of cumulative numbers of sample items with X_i in intervals I_1...I_k
  A[1]<- as.integer(floor(a[1]*n))
  tv<- (A[1]<3)
  if (m>2){
  for (j in 2:(m-1)) {
    A[j]<- as.integer(floor(a[j]*n))
    if (A[j]-A[j-1]<3) {tv<- TRUE}
  }}
  A[m]<- n
  if (A[m]-A[m-1]<3) {tv<- TRUE}   ## a minimum of 3 sample items must be in every interval
  if (tv) {
    stop("vector a not valid")}
    outl<- list()
  u<- rank(x)
  ss1<- ss2<- sf1<- sf2<- sp1<- sp2<- sh1<- sh2<- 0.0
  signi<- 1  # sign of Y in the interval
  for (k in 1:m) {
    ### compute normed difference of ranks
    if (k==1) {
      AA<- 0
      u1<- u[u<=A[1]]
      v1<- y[u<=A[1]]   # sample of Y-values with X_i in I_k
    } else {
      AA<- A[k-1]
      u1<- u[(u>A[k-1])&(u<=A[k])]
      v1<- signi*y[(u>A[k-1])&(u<=A[k])]  # sample of Y-or (-Y)-values with X_i in I_k
    }
    v2<- rank(v1)
    if (methodF==1) {
      uu<- (u1-AA-v2)/(A[k]-AA)  # normed difference of the ranks
      vv<- (u1-A[k]-1+v2)/(A[k]-AA)
    }else{
      if (methodF==2){
        uu<- (u1-AA-v2)/(A[k]-AA+1.0)  # normed difference of the ranks
        vv<- (u1-A[k]-1+v2)/(A[k]-AA+1.0)
      } else {
        uu<- (u1-AA-v2)/sqrt((A[k]-AA)^2-1.0)  # normed difference of the ranks
        vv<- (u1-A[k]-1+v2)/sqrt((A[k]-AA)^2-1.0)
      }}
    ### Spearman coefficient
    if ("Spearman" %in% method) {
      ss1<- ss1+sum(psi(uu,1))
      ss2<- ss2+sum(psi(vv,1))}
    ### Spearman's footrule
    if ("footrule" %in% method) {
      sf1<- sf1+sum(psi(uu,2))
      sf2<- sf2+sum(psi(vv,2))}
    ### power function coefficient
    if ("power" %in% method) {
      sp1<- sp1+sum(psi(uu,3,parp))
      sp2<- sp2+sum(psi(vv,3,parp))}
    ### Huber function coefficient
    if ("Huber" %in% method) {
      for (i in 1:(length(uu))){
        sh1<- sh1+ psi(uu[i],4,parH)
        sh2<- sh2+ psi(vv[i],4,parH)
      }
    }
    signi<- -signi
  }
  ## evaluation of the coefficients
  if ("Spearman" %in% method) {
    ss1<- 1.0-(ss1/n/psif(1))
    ss2<- 1.0-(ss2/n/psif(1))
    outl<- append(outl,list(Spearman=c(ss1,ss2)))}
  if ("footrule" %in% method) {
    sf1<- 1.0-(sf1/n/psif(2))
    sf2<- 1.0-(sf2/n/psif(2))
    outl<- append(outl,list(footrule=c(sf1,sf2)))}
  if ("power" %in% method) {
    sp1<- 1.0-(sp1/n/psif(3,parp))
    sp2<- 1.0-(sp2/n/psif(3,parp))
    outl<- append(outl,list(power=c(sp1,sp2)))}
  if ("Huber" %in% method) {
    sh1<- 1.0-(sh1/n/psif(4,parH))
    sh2<- 1.0-(sh2/n/psif(4,parH))
    outl<- append(outl,list(Huber=c(sh1,sh2)))}
  return(outl)
}

#' Spearman regression coefficient
#'
#' The function spearr evaluates the multivariate Spearman regression coefficient.
#' It describes how well the target variable y can be fit by a function of regressor variables
#' which is increasing w.r.t. some regressors and decreasing w.r.t. the other
#' regressors.
#'
#' @usage spearr(x,y,direction=NULL,out=0)
#' @param x data matrix of regressor variables
#' @param y data vector of the target variable
#' @param out value 1: full output, value 0: reduced output, only coefficients
#' that are largest in absolute value
#' @param direction vector of length d (d is number of regressors),
#' value 1 refers to regressors leading to increasing y whenever this regressor increases,
#' value -1 refers to regressors leading to decreasing y whenever this regressor increases.
#' If direction=NULL, then all coefficients are computed.
#' @return list of Spearman regression coefficients for several directions
#' @references Eckhard Liebscher (2019). A copula-based dependence measure for regression analysis. submitted
#' @examples
#' library(MASS)
#' data <- gilgais
#' spearr(data[,1:3],data[,4],out=1)
#' @export

spearr<-function(x,y,direction=NULL,out=0){
  d<- length(x[1,])
  if (is.null(direction)) { ka<-(2^(d-1)-1)} else {   #check the parameters
    ka<- 0
    if (any(direction==0)){ stop("parameter direction not valid")
  }}
  nn<- length(y)
  if (nn!=length(x[,1])){ stop("number of sample items in x and y are different")}
  if (!is.vector(y)){ stop("data format wrong")}
  if (nn<4) {  stop("not enough sample items")}

  msr<- -10.0 #maximum coefficient
  ret<- list()
  for (k in ka:0){
    if (ka>0){
      direction<- 2*as.numeric(c(1,intToBits(k)[1:(d-1)]))-1   #sign factors as Vector, factor of X1 is plus
      if (k==ka)  {
        u<- apply(x,2,function(cc) (rank(cc)-0.5)/nn)  #columnwise copula-transform

      } else {
        for (j in 1:d){
          if (direction[j]<0) {u[,j]<- (rank(-x[,j])-0.5)/nn  #copula-transform for backwards direction
          }else{
            u[,j]<- (rank(x[,j])-0.5)/nn}}  #copula transform for forwards direction
      }
    } else {  #  case ka=0
      u<- matrix(nrow=nn,ncol=d)
      for (j in 1:d){
        if (direction[j]<0) {u[,j]<- (rank(-x[,j])-0.5)/nn  #copula-transform for backwards direction
        }else{
          u[,j]<- (rank(x[,j])-0.5)/nn}}    #copula-transform for forwards direction
    }

    if (d>2) {z<- apply(u,1,function(rw) (-prod(1.0-rw)+prod(rw))) #negative phi-function (rowwise)
    } else{ z<- apply(u,1,function(rw) sum(rw)-1.0)}

    v<- 2.0*rank(y)-1.0-nn  # transform the y-values
    Bn<- sum(z*v)  #numerator in the estimator
    An<- sum(sort(z)*(2*(1:nn)-1.0-nn)) #denominator in the estimator
    if (out==1){   # full output
      if (An==0.0) {ret[[k+1]]<- list(dcoeff="not defined",dir=direction)
      }else{ret[[k+1]]<- list(dcoeff=Bn/An,dir=direction)}  #An,Bn,
    }else{
      if (An!=0.0){  ##reduced output
        h<- Bn/An
        if (abs(h)>msr) {
          msr<- abs(h)
          ret<- list(dcoeff=h,dir=direction)
        }
      }
    }
  }
  return(ret)
}

#' Spearman regression coefficient for split domains
#'
#' The function spearrs evaluates the multivariate Spearman regression coefficient
#' for two regressors and split regressor region.
#' It describes how well the target variable can be fit in each split region
#' by a function which is increasing w.r.t. some regressors and decreasing
#' w.r.t. the other regressors.
#'
#' @usage spearrs(x,y,splitp=NULL)
#' @param x datamatrix of regressor variables with two columns,
#' @param y data vector of the target variable
#' @param splitp vector of length 2 of the splitting points,
#' If p1 is the first component of this vector, then the point splits the domain of the
#' first regressor into a left region of fraction p1 of data items and a right region
#' of the remaining data items. The same is done for the second regressor. As the
#' result we obtain 4 subregions of the regressor domain. default=c(0.5,0.5)
#' @return list of Kendall regression coefficients for the 4 split regions
#' and the total coefficient together with the corresponding optimal directions.
#'   direction ++ means that y increases whenever both regressors increases
#'   direction +- means that y increases whenever the first regressor increases and the
#' other regressor decreases..etc.
#' @references Eckhard Liebscher (2019). A copula-based dependence measure for regression analysis. submitted
#' @examples
#' library(MASS)
#' data<- gilgais
#' spearrs(data[,1:2],data[,3],splitp=c(0.4,0.6))
#' @export

spearrs<-function(x,y,splitp=NULL){
  nn<- length(y)
  if (length(x[,1])!=nn){ stop("number of sample items in x and y are different")} #check the parameters
  if (!is.vector(y)){ stop("data format wrong")}
  if (nn<4) {  stop("not enough sample items")}

  cg<- 0.0  # total coefficient
  if (length(x[1,])>2) {x<- x[,1:2]} #first columns are selected in the case d>2

  if (is.null(splitp)) {
    splitp<- c(0.5,0.5)}  # default case: Split points in the middle
  if ((min(splitp)<=0.0)|(max(splitp)>=1.0)){ stop("parameter(s) not valid")}

  ret<- list()
  dca<- dcb<- array(0.0,dim=c(2,2,2))
  dc<- vector(length=4)
  u<- apply(x,2,function(cc) (rank(cc)-0.5)/nn)  #columnwise copula transform

  us1<- (u[,1]<=splitp[1])  #logical vector for Split region 1 (left)
  us2<- (u[,2]<=splitp[2])  #logical vector for Split region 2 (right)
  for (k in 1:2){
    k0<- 1
    for (l1 in 1:2){  #Split for u[,1]
      l01<- (l1==1)
      for (l2 in 1:2){  #Split for u[,2]
        l02<- (l2==1)
        u0<- subset(u,(us1==l01)&(us2==l02)) #choose the subset
        y0<- subset(y,(us1==l01)&(us2==l02))
        n1<- length(u0[,1])  #number of sample items in the subset
        if (k==1) {   #k=1 --> "++" coefficient
          u0<- apply(u0,2,function(cc) (rank(cc)-0.5)/n1)
        } else{       #k=2 --> "+-" coefficient
          u0[,1]<- (rank(u0[,1])-0.5)/n1
          u0[,2]<- (rank(-u0[,2])-0.5)/n1
        }
        z<- apply(u0,1,function(rw) sum(rw)-1.0) #apply phi-function rowwise

        v<- 2.0*rank(y0)-1.0-n1
        Bn<- sum(z*v)  #numerator of the estimator
        An<- sum(sort(z)*(2*(1:n1)-1.0-n1)) #denominator
        if (k==1) {
          dc[k0]<- ifelse((An!=0.0),Bn/An,0.0)
        } else{
          if (An!=0.0) {h<- Bn/An}else{h<- 0.0}
          ret[[k0]]<- dc[k0]
          names(ret)[[k0]]<- paste0("dcoeff++",toString(l1),l2)
          ret[[k0+1]]<- h
          names(ret)[[k0+1]]<- paste0("dcoeff+-",toString(l1),l2)
        }
        dca[k,l1,l2]<- ifelse((An>0.0),n1*An,0.0)
        dcb[k,l1,l2]<- n1*Bn
        k0<- k0+2
      }}
  }
  mc<- 0.0
  for (t1 in 1:2){
    for (t2 in 1:2){
      for (t3 in 1:2){
        for (t4 in 1:2){
          h<- (abs(dcb[t1,1,1])+abs(dcb[t2,1,2])+abs(dcb[t3,2,1])+abs(dcb[t4,2,2]))/
            (dca[t1,1,1]+dca[t2,1,2]+dca[t3,2,1]+dca[t4,2,2])
          if (h>mc) {
            mc<- h
            if (t1==1){ ic<- ifelse((dcb[t1,1,1]>0),"++","--")
            } else { ic<- ifelse((dcb[t1,1,1]>0),"+-","-+")}
            if (t2==1){ ic<- paste(ic,ifelse((dcb[t2,1,2]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t2,1,2]>0),"+-","-+"))}
            if (t3==1){ ic<- paste(ic,ifelse((dcb[t3,2,1]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t3,2,1]>0),"+-","-+"))}
            if (t4==1){ ic<- paste(ic,ifelse((dcb[t4,2,2]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t4,2,2]>0),"+-","-+"))}
            #cat(ic,"\n")
          }
        }}}}
  ret[[k0]]<- mc     #output total coefficient
  names(ret)[[k0]]<- "totalcoeff"
  ret[[k0+1]]<- ic
  names(ret)[[k0+1]]<- "directions"
  return(ret)
}

#' Kendall regression coefficient
#'
#' The function kendr evaluates the multivariate Kendall regression coefficient.
#' It describes how well the target variable y can be fit by a function of regressor variables
#' which is increasing w.r.t. some regressors and decreasing w.r.t. the other
#' regressors.
#'
#' @usage kendr(x,y,direction=NULL,out=0)
#' @param x data matrix of regressor variables
#' @param y data vector of the target variable
#' @param out value 1: full output, value 0: reduced output, only coefficients
#' that are largest in absolute value
#' @param direction vector of length d (d is number of regressors),
#' value 1 refers to regressors leading to increasing y whenever this regressor increases,
#' value -1 refers to regressors leading to decreasing y whenever this regressor increases.
#' If direction=NULL, then all coefficients are computed.
#' @return list of Kendall regression coefficients for several directions
#' @references Eckhard Liebscher (2019). Kendall regression Coefficient. submitted
#' @examples
#' library(MASS)
#' data <- gilgais
#' kendr(data[,1:3],data[,4],out=1)
#' @export

kendr<-function(x,y,direction=NULL,out=0){
  nn<- length(y)
  d<- length(x[1,])
  msr<- -10.0 #maximaler Koeffizient
  if (is.null(direction)) { ka<-(2^(d-1)-1)} else {ka<- 0}
  ret<- list()

  for (k in ka:0){
    u<- as.matrix(x)
    if (ka>0){
      direction<- 2*as.numeric(c(1,intToBits(k)[1:(d-1)]))-1   #Vorzeichenfaktoren als Vektor, X1 immer plus
      if (k!=ka)  {
        for (j in 1:d){
          if (direction[j]<0) {u[,j]<- -x[,j]}
        }}
    } else {  #  Fall ka=0
      u<- x
      for (j in 1:d){
        if (direction[j]<0) {u[,j]<- -x[,j] }
      }
    }
    v<- as.matrix(cbind(u,y))
    A<-(sum(copula::F.n(u,u,offset = 0, smoothing = "none"))-1.0)/(nn-1.0) #Nenner
    B<-(sum(copula::F.n(v,v,offset = 0, smoothing = "none"))-1.0)/(nn-1.0) #Zaehler
    if (out==1){   # Ausgabe aller Koeffizienten
      if (A==0.0) {ret[[k+1]]<- list(dcoeff="not defined",dir=direction)
      }else{ret[[k+1]]<- list(dcoeff=2*B/A-1.0,dir=direction)}  #An,Bn,
    }else{
      if (A!=0.0){  ##reduzierte Ausgabe des betragsm. größten Wertes
        h<- 2*B/A-1.0
        if (abs(h)>msr) {
          msr<- abs(h)
          ret[[1]]<- list(dcoeff=h,dir=direction)
        }
      }
    }
  }
  return(ret)
}

#' Kendall regression coefficient for split domains
#'
#' The function kendrs evaluates the multivariate Kendall regression coefficient
#' for two regressors and split regressor region.
#' It describes how well the target variable can be fit in each split region
#' by a function which is increasing w.r.t. some regressors and decreasing
#' w.r.t. the other regressors.
#'
#' @usage kendrs(x,y,splitp=NULL)
#' @param x datamatrix of regressor variables with two columns,
#' @param y data vector of the target variable
#' @param splitp vector of length 2 of the splitting points,
#' If p1 is the first component of this vector, then the point splits the domain of the
#' first regressor into a left region of fraction p1 of data items and a right region
#' of the remaining data items. The same is done for the second regressor. As the
#' result we obtain 4 subregions of the regressor domain. default=c(0.5,0.5)
#' @return list of Kendall regression coefficients for the 4 split regions
#' and the total coefficient together with the corresponding optimal directions.
#' direction ++ means that y increases whenever both regressors increases
#' direction +- means that y increases whenever the first regressor increases and the
#' other regressor decreases..etc.
#' @references Eckhard Liebscher (2019). Kendall regression coefficient. submitted
#' @examples
#' library(MASS)
#' data<- gilgais
#' kendrs(data[,1:2],data[,3],splitp=c(0.4,0.6))
#' @export

kendrs<-function(x,y,splitp=NULL){
  nn<- length(y)
  cg<- 0.0  # total coefficient
  if (length(x[1,])>2) {x<- x[,1:2]} #the first two columns are chosen for d>2

  if (is.null(splitp)) {
    splitp<- c(0.5,0.5)}  # Split points in the middle
  ret<- list()
  dca<- dcb<- array(0.0,dim=c(2,2,2))
  dc<- vector(length=4)
  s1<- stats::quantile(x[,1],splitp[1])
  s2<- stats::quantile(x[,2],splitp[2])
  us1<- (x[,1]<=s1)  #logical Vector for split region 1
  us2<- (x[,2]<=s2)  #logical Vector for split region 2
  for (k in 1:2){
    k0<- 1
    for (l1 in 1:2){  #Split for u[,1]
      l01<- (l1==1)
      for (l2 in 1:2){  #Split for u[,2]
        l02<- (l2==1)
        x0<- as.matrix(subset(x,(us1==l01)&(us2==l02))) #choose the subsets
        y0<- subset(y,(us1==l01)&(us2==l02))
        n1<- length(x0[,1])  #number of items
        if (k==2) {   #k=1 --> "++" coefficient, k=2 --> "+-" coefficient
          x0[,2]<- -x0[,2]
        }
        v<- as.matrix(cbind(x0,y0))
        B<- (sum(copula::F.n(x0,x0,offset = 0, smoothing = "none"))-1.0)/(n1-1.0) #denominator
        A<- (sum(copula::F.n(v,v,offset = 0, smoothing = "none"))-1.0)/(n1-1.0) #numerator
        if (k==1) {
          dc[k0]<- ifelse((B!=0.0),2*A/B-1.0,0.0)
        } else{
          if (B!=0.0) {h<- 2*A/B-1.0}else{h<- 0.0}
          ret[[k0]]<- dc[k0]
          names(ret)[[k0]]<- paste0("dcoeff++",toString(l1),l2)
          ret[[k0+1]]<- h
          names(ret)[[k0+1]]<- paste0("dcoeff+-",toString(l1),l2)
        }
        dca[k,l1,l2]<- n1*A
        dcb[k,l1,l2]<- ifelse((B>0.0),n1*B,0.0)
        k0<- k0+2
      }}
  }
  mc<- 0.0
  for (t1 in 1:2){
    for (t2 in 1:2){
      for (t3 in 1:2){
        for (t4 in 1:2){
          h<- -1.0+2.0*(abs(dca[t1,1,1])+abs(dca[t2,1,2])+abs(dca[t3,2,1])+abs(dca[t4,2,2]))/
            (dcb[t1,1,1]+dcb[t2,1,2]+dcb[t3,2,1]+dcb[t4,2,2])
          if (h>mc) {
            mc<- h
            if (t1==1){ ic<- ifelse((dcb[t1,1,1]>0),"++","--")
            } else { ic<- ifelse((dcb[t1,1,1]>0),"+-","-+")}
            if (t2==1){ ic<- paste(ic,ifelse((dcb[t2,1,2]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t2,1,2]>0),"+-","-+"))}
            if (t3==1){ ic<- paste(ic,ifelse((dcb[t3,2,1]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t3,2,1]>0),"+-","-+"))}
            if (t4==1){ ic<- paste(ic,ifelse((dcb[t4,2,2]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t4,2,2]>0),"+-","-+"))}
          }
        }}}}
  ret[[k0]]<- mc
  names(ret)[[k0]]<- "gencoeff"
  ret[[k0+1]]<- ic
  names(ret)[[k0+1]]<- "directions"
  return(ret)
}
