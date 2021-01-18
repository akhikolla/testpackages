################
# R call for the F procedures
# functional depth computation
# Stanislav Nagy
# nagy at karlin.mff.cuni.cz
# 09/10/2018
#
# Stanislav Nagy, Irene Gijbels & Daniel Hlubinka (2017) Depth-Based Recognition of Shape Outlying Functions, Journal of Computational and Graphical Statistics, 26:4, 883-893, DOI: 10.1080/10618600.2017.1336445
# Stanislav Nagy, Frederic Ferraty (2018) Data Depth for Noisy Random Functions. Under review.
################

if(FALSE){
  # for roxygenize only
  unlink("depth.fd",recursive=TRUE)
  library(roxygen2)
  library(devtools)
  package.skeleton("depth.fd",code_files="depth.fd.R")
  roxygenize("depth.fd")
}

#' @useDynLib depth.fd
#' @export FKS
#' @export shape.fd.analysis
#' @export shape.fd.outliers
#' @export depthf.BD
#' @export depthf.ABD
#' @export depthf.fd1
#' @export depthf.fd2
#' @export depthf.HR
#' @export depthf.hM
#' @export depthf.hM2
#' @export depthf.RP1
#' @export depthf.RP2
#' @export derivatives.est
#' @export L2metric
#' @export Cmetric
#' @export depth.sample
#' @export infimalRank
#' @export dataf2rawfd
#' @export rawfd2dataf

#' @title Transform a \code{dataf} object to raw functional data
#' 
#' @description
#' From a (possibly multivariate) functional data object \code{dataf} constructs an array of the functional values
#' evaluated at an equi-distant grid of points.
#' 
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords functional
#'
#' @param dataf Functions to be transformed, represented by a (possibly multivariate) \code{dataf} object of their arguments
#' and functional values. \code{m} stands for the number of functions. The grid of observation points for the 
#' functions in \code{dataf} may not be the same.
#'
#' @param range The common range of the domain where the fucntions \code{dataf} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{dataf}. If the range is not provided, the smallest interval in which all the arguments from the data functions
#' are contained is chosen as the domain.
#' 
#' @param d Grid size to which all the functional data are transformed. All functional observations are 
#' transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation, see Nagy et al. (2016).
#'
#' @return If the functional data are univariate (scalar-valued), a matrix of size \code{m*d} is given, with each row
#' corresponding to one function. If the functional data are \code{k}-variate with k>1, an array of size \code{m*d*k}
#' of the functional values is given.
#' 
#' @examples
#' ## transform a matrix into a functinal data set and back
#' n = 5
#' d = 21
#' X = matrix(rnorm(n*d),ncol=d)
#' R = rawfd2dataf(X,range=c(0,1))
#' R2 = dataf2rawfd(R,range=c(0,1),d=d)
#' all.equal(X,R2)
#' 
#' ## transform a functional dataset into a raw matrix of functional values
#' dataf = dataf.population()$dataf
#' dataf2rawfd(dataf,range=c(1950,2015),d=66)
#' 
#' ## transform an array into a multivariate functional data set and back
#' k = 3
#' X = array(rnorm(n*d*k),dim=c(n,d,k))
#' R = rawfd2dataf(X,range=c(-1,1))
#' dataf2rawfd(R,range=c(-1,1),d=50)
#'
#' @seealso \code{\link{rawfd2dataf}}
#' @seealso \code{\link{depthf.fd1}}
#' @seealso \code{\link{depthf.fd2}}

# M = dataf2rawfd(dataf.growth()$dataf)

dataf2rawfd = function(dataf, range = NULL, d = 101){
  # transform dataf format for functional data to a raw matrix
  # of functional values evaluated at a grid common to all functions
  # range: range of the common grid, if not specified the range of all the functions
  # d: no of discretized points in the grid, equidistant
  # approximation procedure follows Nagy et al. (2016, JMVA)
  
  # Check "dataf"
  if (!is.list(dataf))
    stop("Argument 'dataf' must be a list")
  
  if(is.vector(dataf[[1]]$vals)){
    mv = FALSE
    for (df in dataf) 
      if (!(is.list(df) && length(df) == 2 &&
            !is.null(df$args) && !is.null(df$vals) &&
            is.vector(df$args) && is.vector(df$vals) &&
            is.numeric(df$args) && is.numeric(df$vals) &&
            length(df$args) == length(df$vals) &&
            is.sorted(df$args)))
        stop("Argument 'dataf' must be a list containing lists (functions) 
             of two vectors of equal length, named 'args' and 'vals': 
             arguments sorted in ascending order and corresponding them 
             values respectively") 
  }
  
  if(is.matrix(dataf[[1]]$vals)){
    mv = TRUE
    for (df in dataf)
      if (!(is.list(df) && length(df) == 2 &&
            !is.null(df$args) && !is.null(df$vals) &&
            is.vector(df$args) && is.matrix(df$vals) &&
            is.numeric(df$args) && is.numeric(df$vals) &&
            length(df$args) == nrow(df$vals) &&
            is.sorted(df$args)))
        stop("Argument 'dataf' must be a list containing lists (functions) 
             of a vector named 'args' and a matrix named 'vals'. The arguments
             of 'args' must be sorted in ascending order. To each element of 'args'
             the corresponding row in 'vals' represents 
             the functional values at this point") 
  }
  
  # range construction
  rng = numeric(0)
  for (df in dataf)	rng = range(c(rng,df$args)) # common range of all the data
  if(!is.null(range)){ 
    if(!(length(range) == 2 && is.numeric(range) && 
         range[1]<=rng[1] && range[2]>=rng[2])) 
      stop("Argument 'range' must be a numeric vector of two components
           that defines the range of the domain of functional data. All
           functional data must have 'args' vales inside this domain.")
  } else range = rng
  if(!(range[1]<range[2])) stop("Argument 'range' must define a non-degenerate interval.")	
  
  t = seq(range[1],range[2],length=d)
  n = length(dataf)
  if(!mv){
    X = matrix(nrow=n,ncol=d)
    # functional data interpolation / extrapolation
    for(i in 1:n){
      ni = length(dataf[[i]]$args)
      X[i,] = approx(dataf[[i]]$args,dataf[[i]]$vals,t)$y
      X[i,t<=dataf[[i]]$args[1]] = dataf[[i]]$vals[1]
      X[i,t>=dataf[[i]]$args[ni]] = dataf[[i]]$vals[ni]
    }
  } else {
    k = ncol(dataf[[1]]$vals)
    X = array(dim=c(n,d,k))
    # functional data interpolation / extrapolation
    for(i in 1:n){
      ni = length(dataf[[i]]$args)
      for(j in 1:k) X[i,,j] = approx(dataf[[i]]$args,dataf[[i]]$vals[,j],t)$y
      X[i,t<=dataf[[i]]$args[1],] = dataf[[i]]$vals[1,]
      X[i,t>=dataf[[i]]$args[ni],] = dataf[[i]]$vals[ni,]
    }
  }
  return(X)
  }

#' @title Transform raw functional data to a \code{dataf} object
#' 
#' @description
#' Constructs a (possibly multivariate) functional data object given by an array of its functional values
#' evaluated at an equi-distant grid of points, and transforms it into a \code{dataf} object more suitable 
#' for work in the \code{ddalpha} package.
#' 
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords functional
#'
#' @param X Either a matrix of size \code{n*d}, or an array of dimension \code{n*d*k} of functional values. Here \code{n}
#' stands for the number of functions, \code{d} is the number of equi-distant points in the domain where the functional
#' values are evaluated, and if applicable, \code{k} is the dimensionality of the (vector-valued) functional data.
#'
#' @param range A vector of size two that represents the endpoints of the common domain of all functions \code{X}.
#'
#' @return A (possibly multivariate) \code{dataf} object corresponding to the functional data \code{X} evaluated at an
#' equi-distant grid of points.
#' 
#' @examples
#' ## transform a matrix into a functinal data set
#' n = 5
#' d = 21
#' X = matrix(rnorm(n*d),ncol=d)
#' rawfd2dataf(X,range=c(0,1))
#' 
#' ## transform an array into a multivariate functional data set
#' k = 3
#' X = array(rnorm(n*d*k),dim=c(n,d,k))
#' rawfd2dataf(X,range=c(-1,1))
#'
#' @seealso \code{\link{dataf2rawfd}}
#' @seealso \code{\link{depthf.fd1}}
#' @seealso \code{\link{depthf.fd2}}

rawfd2dataf = function(X, range){
  # transform a raw array of functional data values 
  # to a dataf format, where the domain is assumed to be
  # an equidistant grid in the interval given by range
  
  if(is.vector(X)) X = matrix(X,nrow=1) # if X is a single vector, it is considered to be a matrix of one scalar function
  
  # Check "rawfd"
  if (!is.array(X))
    stop("Argument 'X' must be an array (multivariate functional data)
         or a matix (univariate functional data)")
  mv = !is.matrix(X)
  
  # range construction
  if(!(range[1]<range[2])) stop("Argument 'range' must define a non-degenerate interval.")	
  n = dim(X)[1]
  d = dim(X)[2]
  if(mv) k = dim(X)[3]
  t = seq(range[1],range[2],length=d)
  dataf = list()
  for(i in 1:n){
    if(mv) df = list(args = t, vals = X[i,,]) else df = list(args = t, vals = X[i,])
    dataf[[i]] = df
  }
  return(dataf)
}

#' @title Fast kernel smoothing
#'
#' @description
#' Produces a kernel smoothed version of a function based on
#' the vectors given in the input. Bandwidth is selected using cross-validation.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords smoothing kernel functional
#'
#' @details
#' A vector of the same length as \code{Tout}
#' corresponding to the values of the
#' function produced using kernel smoothing, is provided. Bandwidth is selected using the
#' \code{K}-fold cross-validation of randomly shuffled input values.
#'
#' @param dataf A set of functional data given by a \code{dataf} object that are to be smoothed.
#'
#' @param Tout vector of values in the domain of the functions at which the
#' resulting smoothed function is evaluated
#' 
#' @param kernel Kernel used for smoothing. Admissible values are \code{uniform}, 
#' \code{triangular}, \code{Epanechnikov}, \code{biweight}, \code{triweight} and \code{Gaussian}.
#' By default, \code{uniform} is used.
#' 
#' @param m Number of points in the grid for choosing the cross-validated bandwidth.
#' 
#' @param K Performs \code{K}-fold cross-validation based on randomly shuffled data.
#'
#' @return A \code{dataf} object corresponding to \code{Tout} of smoothed functional values.
#' 
#' @examples
#' d = 10
#' T = sort(runif(d))
#' X = T^2+ rnorm(d,sd=.1)
#' Tout = seq(0,1,length=101)
#' 
#' plot(T,X)
#' dataf = list(list(args=T,vals=X))
#' data.sm = FKS(dataf,Tout,kernel="Epan")
#' lines(data.sm[[1]]$args,data.sm[[1]]$vals,col=2)
#' 
#' datafs = structure(list(dataf=dataf,labels=1:length(dataf)),class="functional")
#' plot(datafs)
#' points(T,X)
#' data.sms = structure(list(dataf=data.sm,labels=1:length(data.sm)),class="functional")
#' plot(data.sms)
#' 
#' n = 6
#' dataf = list()
#' for(i in 1:n) dataf[[i]] = list(args = T<-sort(runif(d)), vals = T^2 + rnorm(d,sd=.1))
#' data.sm = FKS(dataf,Tout,kernel="triweight")
#' data.sms = structure(list(dataf=data.sm,labels=1:length(data.sm)),class="functional")
#' plot(data.sms)

FKS = function(dataf,Tout,kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"),m=51,K=20){
  # kernel smoothing
  # produces directly the kernel smoother with CV value of h
  # Hmax 	the maximal H for CV
  # m		no of elements of H for CV
  # K		K-fold CV for bandwidth selection
  #
  # Args:
  #   dataf:  list containing lists (functions) of two vectors of equal length, 
  #           named "args" and "vals": arguments sorted in ascending order and 
  #           corresponding them values respectively
  #   labels: output labels of the functinal observations
  
  # Check "dataf"
  if (!is.list(dataf))
    stop("Argument 'dataf' must be a list")
  for (df in dataf)
    if (!(is.list(df) && length(df) == 2 &&
          !is.null(df$args) && !is.null(df$vals) &&
          is.vector(df$args) && is.vector(df$vals) &&
          is.numeric(df$args) && is.numeric(df$vals) &&
          length(df$args) == length(df$vals) &&
          is.sorted(df$args)))
      stop("Argument 'dataf' must be a list containing lists (functions) 
           of two vectors of equal length, named 'args' and 'vals': 
           arguments sorted in ascending order and corresponding them 
           values respectively")  
  
  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel)
  kernI = which(kernel==krn)
  for(df in dataf) if(sum(is.na(df$args))>0) stop("NA values in the functional args vector are not allowed")
  # produce the output structure
  sm.dataf = list()
  Tout = sort(Tout)
  for(j in 1:length(dataf)){
    T = dataf[[j]]$args
    X = dataf[[j]]$vals
    ST = dataf[[j]]$args
    eps = 10^(-6)
    Hmin = max(c(ST[1], (1 - ST[length(ST)]), max(ST[-1] - ST[-length(ST)])/2)) + eps
    Hmax = max(ST[length(ST)]-ST[1],(max(Tout)-min(Tout)))
    if (Hmin >= Hmax) Hmin = Hmax/10 # in the extreme case when Hmin > Hmax take Hmin very small
    H = seq(Hmin,Hmax,length=m)
    #
    n = length(ST)
    nR = max(1,floor(n/K))
    if(nR==1) K = n
    TRE = rep(NA,nR*K)
    XRE = rep(NA,nR*K)
    TNRE = rep(NA,K*(n-nR))
    XNRE = rep(NA,K*(n-nR))
    SH = sample(n,n)			# random shuffle of the points for CV
    for(i in 1:K){
      S = SH[(i-1)*nR+(1:nR)] 
      TRE[(i-1)*nR+(1:nR)] = T[S]
      XRE[(i-1)*nR+(1:nR)] = X[S]
      TNRE[(i-1)*(n-nR)+(1:(n-nR))] = T[-S]
      XNRE[(i-1)*(n-nR)+(1:(n-nR))] = X[-S]
    }
    Res = .Fortran("CVKERNSM",
                   as.double(T),
                   as.double(X),
                   as.double(Tout),
                   as.integer(length(T)),
                   as.integer(length(Tout)),
                   as.double(H),
                   as.integer(length(H)),
                   as.integer(kernI),
                   as.double(TRE),
                   as.double(XRE),
                   as.double(TNRE),
                   as.double(XNRE),
                   as.integer(nR),		
                   as.integer(K),
                   as.double(rep(0,length(Tout)))
    )
    Res[Res[[15]]>10^6] = Inf
    sm.dataf[[j]] = list(args = Tout, vals = Res[[15]])
  }
  return(sm.dataf)
}

#' @title Fast depth computation for univariate and bivariate random samples
#'
#' @description
#' Faster implementation of the halfspace and the simplicial depth. Computes the depth 
#' of a whole random sample of a univariate or a bivariate data in one run.	
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth
#'
#' @details
#' The function returns vectors of sample halfspace and simplicial depth values.
#'
#' @param A Univariate or bivariate points whose depth is computed, represented by a matrix of 
#' size \code{m*2}. \code{m} stands for the number of points, \code{d} is 1 for univariate and 2 
#' for bivariate data.
#'
#' @param B Random sample points with respect to which the depth of \code{A} is computed. 
#' \code{B} is represented by a matrix of size \code{n*2}, where \code{n} is the sample size.
#'
#' @return Vector of length \code{m} of depth halfspace depth values is returned.
#' 
#' @examples
#' n = 100
#' m = 150
#' A = matrix(rnorm(2*n),ncol=2)
#' B = matrix(rnorm(2*m),ncol=2)
#' depth.sample(A,B)
#' system.time(D1<-depth.halfspace(A,B))
#' system.time(D2<-depth.sample(A,B))
#' max(D1-D2$Half)
#'
#' A = rnorm(100)
#' B = rnorm(150)
#' depth.sample(A,B)
#' # depth.halfspace(matrix(A,ncol=1),matrix(B,ncol=1))
#'
#' @seealso \code{\link{depth.halfspace}}
#' @seealso \code{\link{depth.simplicial}}

depth.sample = function(A,B){
  # bivariate halfspace depth
  # A points whose depth I compute, M*2 matrix
  # B points wrt whose the depth is computed, N*2 matrix
  if(is.null(dim(A))){ # for univariate data
    A1 = as.vector(A)
    B1 = as.vector(B)	
    m = length(A)
    n = length(B)
    FD = .Fortran("dpth1",	
                  as.numeric(A1),		#	A1
                  as.numeric(B1),		#	B1
                  as.integer(m),		#	m
                  as.integer(n),		#	n
                  sdep=as.numeric(rep(-1,m)),
                  hdep=as.numeric(rep(-1,m))
    )
    return(list(Simpl = FD$sdep, Half = FD$hdep))
  } else {
    A1 = as.vector(A[,1])
    A2 = as.vector(A[,2])
    B1 = as.vector(B[,1])
    B2 = as.vector(B[,2])
    m = dim(A)[1]
    n = dim(B)[1]
    if((dim(B)[2]!=2)|(dim(A)[2]!=2)) stop("Computation for two dimensions only")
    FD = .Fortran("dpth2",	
                  as.numeric(A1),		#	A1
                  as.numeric(A2),		#	A2
                  as.numeric(B1),		#	B1
                  as.numeric(B2),		#	B2
                  as.integer(m),		#	m
                  as.integer(n),		#	n
                  sdep=as.numeric(rep(-1,m)),
                  hdep=as.numeric(rep(-1,m))
    )
    return(list(Simpl = FD$sdep, Half = FD$hdep))
  }
}

#' @title Univariate integrated and infimal depth for functional data
#'
#' @description
#' Usual, and order extended integrated and infimal depths for real-valued functional data based on the
#' halfspace and simplicial depth.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth functional
#'
#' @details
#' The function returns vectors of sample integrated and infimal depth values.
#'
#' @param datafA Functions whose depth is computed, represented by a \code{dataf} object of their arguments
#'  and functional values. \code{m} stands for the number of functions.
#'
#' @param datafB Random sample functions with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a \code{dataf} object of their arguments
#'  and functional values. \code{n} is the sample size. 
#' The grid of observation points for the 
#' functions \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation, see Nagy et al. (2016).
#'
#' @param order The order of the order extended integrated and infimal depths.
#' By default, this is set to \code{1}, meaning that the usual univariate depths of 
#' the functional values are computed. For \code{order=2} or \code{3}, the second
#' and the third order extended integrated and infimal depths are computed, 
#' respectively. 
#'
#' @param approx Number of approximations used in the computation of the order extended depth
#' for \code{order} greater than \code{1}. For \code{order=2}, the default
#' value is set to \code{0}, meaning that the depth is computed at all possible \code{d^order}
#' combinations of the points in the domain. For \code{order=3}, 
#' the default value is set to \code{101}. When \code{approx} is a positive integer, \code{approx}
#' points are randomly sampled in \code{[0,1]^order} and at these points the \code{order}-variate depths of the
#' corresponding functional values are computed.
#'
#' @return Four vectors of length \code{m} of depth values are returned:
#'	\itemize{
#'	\item \code{Simpl_FD} the integrated depth based on the simplicial depth,
#'	\item \code{Half_FD} the integrated depth based on the halfspace depth,
#'	\item \code{Simpl_ID} the infimal depth based on the simplicial depth,
#'	\item \code{Half_ID} the infimal depth based on the halfspace depth.
#'	}
#' In addition, two vectors of length \code{m} of the relative area of smallest depth values is returned:
#'	\itemize{
#'	\item \code{Simpl_IA} the proportions of points at which the depth \code{Simpl_ID} was attained,
#'	\item \code{Half_IA} the proportions of points at which the depth \code{Half_ID} was attained.
#'	}
#' The values \code{Simpl_IA} and \code{Half_IA} are always in the interval [0,1]. 
#' They introduce ranking also among functions having the same
#' infimal depth value - if two functions have the same infimal depth, the one with larger infimal area
#' \code{IA} is said to be less central. 	
#' For \code{order=2} and \code{m=1}, two additional matrices of pointwise depths are also returned:
#'	\itemize{
#'    \item \code{PSD} the matrix of size \code{d*d} containing the computed 
#'    pointwise bivariate simplicial depths used for the computation of \code{Simpl_FD} and \code{Simpl_ID},
#'    \item \code{PHD} the matrix of size \code{d*d} containing the computed 
#'    pointwise bivariate halfspace depths used for the computation of \code{Half_FD} and \code{Half_ID}.
#'	}
#' For \code{order=3}, only \code{Half_FD} and \code{Half_ID} are provided.
#' 
#' @references Nagy, S., Gijbels, I. and Hlubinka, D. (2016). 
#' Weak convergence of discretely observed functional data with applications. 
#' \emph{Journal of Multivariate Analysis}, \bold{146}, 46--62.
#' 
#' @references Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
#' Depth-based recognition of shape outlying functions. 
#' \emph{Journal of Computational and Graphical Statistics}, \bold{26:4}, 883--893.
#'
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#' depthf.fd1(datafA,datafB)
#' depthf.fd1(datafA,datafB,order=2)
#' depthf.fd1(datafA,datafB,order=3,approx=51)
#'
#' @seealso \code{\link{depthf.fd2}}, \code{\link{infimalRank}}

depthf.fd1 = function(datafA,datafB,range=NULL,d=101,order=1,approx=0){
  # univariate integrated depth
  # A functions whose depth I compute, M*D matrix
  # B functions wrt whose the depth is computed, N*D matrix
  # both 1dimensional, n*d, n nr of functions, d dimensionality
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)
  
  if(order==1){
    A1 = as.vector(A)
    B1 = as.vector(B)
    d = dim(A)[2]
    m = dim(A)[1]
    n = dim(B)[1]
    if(dim(B)[2]!=d) stop("dimension mismatch")
    FD = .Fortran("funD1",	
                  as.numeric(A1),		#	A
                  as.numeric(B1),		#	B
                  as.integer(m),		#	m
                  as.integer(n),		#	n
                  as.integer(d),		#	d
                  funsdep=as.numeric(rep(-1,m)),	
                  funhdep=as.numeric(rep(-1,m)),
                  fIsdep =as.numeric(rep(-1,m)),
                  fIhdep =as.numeric(rep(-1,m)),
                  IAsdep =as.integer(rep(-1,m)),
                  IAhdep =as.integer(rep(-1,m))
    )
    return(list(	Simpl_FD = FD$funsdep, Half_FD = FD$funhdep,
                 Simpl_ID = FD$fIsdep, Half_ID = FD$fIhdep, 
                 Simpl_IA = FD$IAsdep/d, Half_IA = FD$IAhdep/d ))
  }
  if(order==2) return(DiffDepth(A,B,approx))
  if(order==3) return(DiffDepth3D(A,B,approx))
}

#' @title Fast computation of the \eqn{L^2} metric for sets of functional data
#'
#' @description
#' Returns the matrix of \eqn{L^2} distances between two sets of functional data.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords metric functional
#'
#' @details
#' For two sets of functional data of sizes \code{m} and \code{n}
#' represented by matrices of their functional values on the common domain {1,...,d}, 
#' this function returns the symmetric matrix of size \code{m*n} whose entry in the
#' \code{i}-th row and \code{j}-th column is the approximated \eqn{L^2} distance of the 
#' \code{i}-th function from the first set, and the \code{j}-th function from the second set.
#' This function is utilized in the computation of the h-mode depth.
#'
#' @param A Functions of the first set, represented by a matrix of their functional values of 
#' size \code{m*d}. \code{m} stands for the number of functions, \code{d}
#' is the number of the equi-distant points {1,...,d} in the domain of the data [1,d] at which the functional
#' values of the \code{m} functions are evaluated.
#'
#' @param B Functions of the second set, represented by a matrix of their functional values of 
#' size \code{n*d}. \code{n} stands for the number of functions, \code{d}
#' is the number of the equi-distant points {1,...,d} in the domain of the data [1,d] at which the functional
#' values of the \code{n} functions are evaluated. The grid of observation points for the 
#' functions \code{A} and \code{B} must be the same.
#'
#' @return A symmetric matrix of the distances of the functions of size \code{m*n}.
#' 
#' @examples
#' datapop = dataf2rawfd(dataf.population()$dataf,range=c(1950,2015),d=66)
#' A = datapop[1:20,]
#' B = datapop[21:50,]
#' L2metric(A,B)
#'
#' @seealso \code{\link{depthf.hM}}
#' @seealso \code{\link{dataf2rawfd}}

L2metric = function(A,B){
  # computes fast approximation of L2 distance between fctions A and B
  M = .Fortran("metrl2",
               as.numeric(A),
               as.numeric(B),
               as.integer(m<-nrow(A)),
               as.integer(n<-nrow(B)),
               as.integer(d<-length(as.numeric(A))/nrow(A)),
               m = as.numeric(rep(-1,m*n)))$m
  return(M = matrix(M,nrow=m))
}

depthf.M = function(A,B,q=.2){
  # h-mode depth for the L2 metric
  mdist = L2metric(B,B)
  mdist2 = L2metric(A,B)
  hq2 = quantile(mdist[mdist>0],probs=q,type=1)	# probs=.2 as in Cuevas et al. (2007)
  return(rowSums(dnorm(mdist2/hq2)))
}

#' @title Fast computation of the uniform metric for sets of functional data
#'
#' @description
#' Returns the matrix of \eqn{C} (uniform) distances between two sets of functional data.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords metric functional
#'
#' @details
#' For two sets of functional data of sizes \code{m} and \code{n}
#' represented by matrices of their functional values, 
#' this function returns the symmetric matrix of size \code{m*n} whose entry in the
#' \code{i}-th row and \code{j}-th column is the approximated \eqn{C} (uniform) distance of the 
#' \code{i}-th function from the first set, and the \code{j}-th function from the second set.
#' This function is utilized in the computation of the h-mode depth.
#'
#' @param A Functions of the first set, represented by a matrix of their functional values of 
#' size \code{m*d}. \code{m} stands for the number of functions, \code{d}
#' is the number of the equi-distant points in the domain of the data at which the functional
#' values of the \code{m} functions are evaluated.
#'
#' @param B Functions of the second set, represented by a matrix of their functional values of 
#' size \code{n*d}. \code{n} stands for the number of functions, \code{d}
#' is the number of the equi-distant points in the domain of the data at which the functional
#' values of the \code{n} functions are evaluated. The grid of observation points for the 
#' functions \code{A} and \code{B} must be the same.
#'
#' @return A symmetric matrix of the distances of the functions of size \code{m*n}.
#' 
#' @examples
#' datapop = dataf2rawfd(dataf.population()$dataf,range=c(1950,2015),d=66)
#' A = datapop[1:20,]
#' B = datapop[21:50,]
#' Cmetric(A,B)
#'
#' @seealso \code{\link{depthf.hM}}
#' @seealso \code{\link{dataf2rawfd}}

Cmetric = function(A,B){
  # computes fast approximation of \eqn{C} distance between fctions A and B
  M = .Fortran("metrC",
               as.numeric(A),
               as.numeric(B),
               as.integer(m<-nrow(A)),
               as.integer(n<-nrow(B)),
               as.integer(d<-length(as.numeric(A))/nrow(A)),
               m = as.numeric(rep(-1,m*n)))$m
  return(M = matrix(M,nrow=m))
}

depthf.MC = function(A,B,q=.2){
  # h-mode depth for the \eqn{C} metric
  mdist = Cmetric(B,B)
  mdist2 = Cmetric(A,B)
  hq2 = quantile(mdist[mdist>0],probs=q,type=1)	# probs=.2 as in Cuevas et al. (2007)
  return(rowSums(dnorm(mdist2/hq2)))
}

#' @title h-mode depth for functional data
#'
#' @description
#' The h-mode depth of functional real-valued data.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth functional
#'
#' @details
#' The function returns the vectors of the sample h-mode depth values. The kernel 
#' used in the evaluation is the standard Gaussian kernel, the bandwidth value is chosen
#' as a quantile of the non-zero distances between the random sample curves.
#'
#' @param datafA Functions whose depth is computed, represented by a \code{dataf} object of their arguments
#'  and functional values. \code{m} stands for the number of functions.
#'
#' @param datafB Random sample functions with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a \code{dataf} object of their arguments
#'  and functional values. \code{n} is the sample size. 
#' The grid of observation points for the 
#' functions \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @param norm The norm used for the computation of the depth. Two possible 
#' choices are implemented: \code{C} for the uniform norm of continuous functions, 
#' and \code{L2} for the \eqn{L^2} norm of integrable functions.
#'
#' @param q The quantile used to determine the value of the bandwidth \eqn{h}
#' in the computation of the h-mode depth. \eqn{h} is taken as the \code{q}-quantile of
#' all non-zero distances between the functions \code{B}. By default, this value is set
#' to \code{q=0.2}, in accordance with the choice of Cuevas et al. (2007).
#'
#' @return A vector of length \code{m} of the h-mode depth values.
#' 
#' @references Cuevas, A., Febrero, M. and Fraiman, R.  (2007).
#' Robust estimation and classification for functional data via projection-based depth notions. 
#' \emph{Computational Statistics} \bold{22} (3), 481--496.
#' 
#' @references Nagy, S., Gijbels, I. and Hlubinka, D. (2016). 
#' Weak convergence of discretely observed functional data with applications. 
#' \emph{Journal of Multivariate Analysis}, \bold{146}, 46--62.
#'
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#' depthf.hM(datafA,datafB)
#' depthf.hM(datafA,datafB,norm="L2")

depthf.hM = function(datafA,datafB,range=NULL,d=101, norm = c("C","L2"), q=.2){
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)
  
  norm = match.arg(norm)
  switch(norm,
         C = depthf.MC(A,B,q),
         L2 = depthf.M(A,B,q)
  )
}

AdjBDKernL = function(b,v,J=3){
  # v is a matrix m x eval
  m = dim(v)[1]                                   # sample size
  eval = length(b)                                
  com = combn(m,J)
  poc = dim(com)[2]
  s = .Fortran("AdjLP",
               as.integer(eval),
               as.integer(J),
               as.integer(m),
               as.integer(poc),
               as.integer(as.vector(com)),
               as.numeric(b),
               as.numeric(v),
               dj = as.numeric(1))$dj
  return(s)
}

AdjBDL = function(b,v,J=3,K=1){
  eval = length(b)
  m = dim(v)[1]                                   # sample size
  sam = sample.int(m,m)                           # shuffle
  if (K==0) K=1                                   # make sure K!=0
  nk = m%/%K                                      # subsample size
  nlast = m %% K+nk                               # last subsample size
  Dpom = rep(NA,K)	   					# subsample depth
  if (K>1){
    for (k in 1:(K-1)) 
      if (eval>1) Dpom[k] = AdjBDKernL(b,v[sam[((k-1)*nk+1):(k*nk)],],J)
      else  Dpom[k] = AdjBDKernL(b,t(as.matrix(v[sam[((k-1)*nk+1):(k*nk)],])),J)
  }
  {if (nlast>0){	
    if (eval>1) Dpom[K]= AdjBDKernL(b,v[sam[((K-1)*nk+1):m],],J)
    else  Dpom[K]= AdjBDKernL(b,t(as.matrix(v[sam[((K-1)*nk+1):m],])),J)
  }
    else (K=K-1)}
  return(mean(Dpom[1:K]))
}

AdjBDsampleL = function(A,B,J=3,K=1){
  m = dim(A)[1]
  hlb = rep(NA,m)
  for (i in 1:m) hlb[i] = AdjBDL(A[i,],B,J,K)
  return(hlb)
}

AdjBDKernC = function(b,v,J=3){
  # v ma rozmery m x eval
  m = dim(v)[1]						# sample size
  eval = length(b)					
  com = combn(m,J)	
  poc = dim(com)[2]
  s = .Fortran("AdjC",
               as.integer(eval),
               as.integer(J),
               as.integer(m),
               as.integer(poc),
               as.integer(as.vector(com)),
               as.numeric(b),
               as.numeric(v),
               dj = as.numeric(1))$dj
  return(s)
}

AdjBDC = function(b,v,J=3,K=1){
  # v is a matrix m x eval
  eval = length(b)
  m = dim(v)[1]                                   # sample size
  sam = sample.int(m,m)                           # shuffle
  if (K==0) K=1                                   # make sure K!=0
  nk = m%/%K                                      # subsample size
  nlast = m %% K+nk                               # last subsample size
  Dpom = rep(NA,K)	   					# subsample depth
  if (K>1){
    for (k in 1:(K-1)) 
      if (eval>1) Dpom[k] = AdjBDKernC(b,v[sam[((k-1)*nk+1):(k*nk)],],J)
      else  Dpom[k] = AdjBDKernC(b,t(as.matrix(v[sam[((k-1)*nk+1):(k*nk)],])),J)
  }
  {if (nlast>0){
    if (eval>1) Dpom[K]= AdjBDKernC(b,v[sam[((K-1)*nk+1):m],],J)
    else  Dpom[K]= AdjBDKernC(b,t(as.matrix(v[sam[((K-1)*nk+1):m],])),J)
  }
    else (K=K-1)}
  return(mean(Dpom[1:K]))
}

AdjBDsampleC = function(A,B,J=3,K=1){
  m = dim(A)[1]
  hlb = rep(NA,m)
  for (i in 1:m) hlb[i] = AdjBDC(A[i,],B,J,K)
  return(hlb)
}

#' @title Adjusted band depth for functional data
#'
#' @description
#' The adjusted band depth 
#' of functional real-valued data based on either the
#' \eqn{C} (uniform) norm, or on the \eqn{L^2} norm of functions.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth functional
#'
#' @details
#' The function returns the vector of the sample adjusted band depth values. The kernel 
#' used in the evaluation is the function \eqn{K(u) = exp(-u)}.
#'
#' @param datafA Functions whose depth is computed, represented by a \code{dataf} object of their arguments
#'  and functional values. \code{m} stands for the number of functions.
#'
#' @param datafB Random sample functions with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a \code{dataf} object of their arguments
#'  and functional values. \code{n} is the sample size. 
#' The grid of observation points for the 
#' functions \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation, see Nagy et al. (2016).
#'
#' @param norm The norm used for the computation of the depth. Two possible 
#' choices are implemented: \code{C} for the uniform norm of continuous functions, 
#' and \code{L2} for the \eqn{L^2} norm of integrable functions. 
#'
#' @param J The order of the adjusted band depth, that is the maximal number of functions
#' taken in a band. Acceptable values are \code{2}, \code{3},... By default this value is set to \code{2}. 
#' Note that this is NOT the order as
#' defined in the order-extended version of adjusted band depths in Nagy et al. (2016), used
#' for the detection of shape outlying curves.
#' 
#' @param K Number of sub-samples of the functions from \code{B} taken to speed up the
#' computation. By default, sub-sampling is not performed. Values of \code{K} larger than \code{1}
#' result in an approximation of the adjusted band depth.
#'
#' @return A vectors of length \code{m} of the adjusted band depths.
#' 
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#' depthf.ABD(datafA,datafB)
#' depthf.ABD(datafA,datafB,norm="L2")
#'
#' @seealso \code{\link{depthf.BD}}
#'
#' @references Gijbels, I., Nagy, S. (2015).
#' Consistency of non-integrated depths for functional data.
#' \emph{Journal of Multivariate Analysis} \bold{140}, 259--282.
#' 
#' @references Nagy, S., Gijbels, I. and Hlubinka, D. (2016). 
#' Weak convergence of discretely observed functional data with applications. 
#' \emph{Journal of Multivariate Analysis}, \bold{146}, 46--62.
#' 
#' @references Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
#' Depth-based recognition of shape outlying functions. 
#' \emph{Journal of Computational and Graphical Statistics}, \bold{26:4}, 883--893.

depthf.ABD = function(datafA,datafB,range=NULL,d=101, norm = c("C","L2"), J=2, K=1){
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)
  
  norm = match.arg(norm)
  switch(norm,
         C = AdjBDsampleC(A,B,J,K),
         L2 = AdjBDsampleL(A,B,J,K)
  )
}

DiffDepth = function(A,B,approx=0){
  # A functions whose depth I compute
  # B functions wrt whose the depth is computed
  # both 1dimensional, n*d, n nr of functions, d dimensionality
  # approx is number of approximations to be used, 0 for full computation
  
  d = dim(A)[2]
  m = dim(A)[1]
  n = dim(B)[1]
  if(dim(B)[2]!=d) stop("dimension mismatch")
  Av = as.vector(A)
  Bv = as.vector(B)
  if (approx>0){
    C = combn(1:d,2)		# all couples 
    if (approx>=ncol(C)){ 
      ind.sample = 1:ncol(C) 
      approx = ncol(C) } else 
        ind.sample = sample.int(ncol(C),approx)
      RN = matrix(C[,ind.sample],nrow=2)		#	RN = random numbers of size 2*REP from 1:d
  } else RN = 0		
  FD = .Fortran("DiffD",	
                as.numeric(Av),		#	A
                as.numeric(Bv),		#	B
                as.integer(m),		#	m
                as.integer(n),		#	n
                as.integer(d),		#	d
                as.integer(approx),	#	REP
                as.integer(RN),		#	RN
                funsdep=as.numeric(rep(-1,m)),	
                funhdep=as.numeric(rep(-1,m)),
                funsdepm=as.numeric(rep(-1,m)),	
                funhdepm=as.numeric(rep(-1,m)),
                Psdep = as.numeric(rep(-1,d*d)),
                Phdep = as.numeric(rep(-1,d*d)),
                IAsdep =as.integer(rep(-1,m)),
                IAhdep =as.integer(rep(-1,m)))
  if(approx==0){
    S_IA = FD$IAsdep/(d^2)
    H_IA = FD$IAhdep/(d^2)
  } else {
    S_IA = FD$IAsdep/(approx)
    H_IA = FD$IAhdep/(approx)		
  }
  if (m>1) return(list(	Simpl_FD = FD$funsdep, Half_FD = FD$funhdep,
                        Simpl_ID = FD$funsdepm, Half_ID = FD$funhdepm,
                        Simpl_IA = S_IA, Half_IA = H_IA))
  if ((m==1)&(approx==0)){
    PSD = matrix(FD$Psdep,ncol=d)
    PSD = pmax(PSD,0)
    PSD = PSD + t(PSD)
    diag(PSD) = NA
    PHD = matrix(FD$Phdep,ncol=d)
    PHD = pmax(PHD,0)
    PHD = PHD + t(PHD)
    diag(PHD) = NA
    return(list(	Simpl_FD = FD$funsdep, Half_FD = FD$funhdep,
                 Simpl_ID = FD$funsdepm, Half_ID = FD$funhdepm,
                 PSD = PSD, PHD = PHD,
                 Simpl_IA = S_IA, Half_IA = H_IA))
  }
  if ((m==1)&(approx>0)){
    PSD = matrix(FD$Psdep,ncol=d)
    for(i in 1:approx) PSD[RN[2*i-1],RN[2*i]]=PSD[RN[2*i],RN[2*i-1]]	# symmetrize
    PSD[PSD==-1] = NA
    PHD = matrix(FD$Phdep,ncol=d)
    for(i in 1:approx) PHD[RN[2*i-1],RN[2*i]]=PHD[RN[2*i],RN[2*i-1]]	# symmetrize
    PHD[PHD==-1] = NA
    return(list(	Simpl_FD = FD$funsdep, Half_FD = FD$funhdep, 
                 Simpl_ID = FD$funsdepm, Half_ID = FD$funhdepm,
                 PSD = PSD, PHD = PHD,
                 Simpl_IA = S_IA, Half_IA = H_IA))
  }
}

DiffDepth3D = function(A,B,approx=101){
  
  d = dim(A)[2]
  m = dim(A)[1]
  n = dim(B)[1]
  if(dim(B)[2]!=d) stop("dimension mismatch")
  if(approx==0) approx=101
  DT = matrix(nrow=m,ncol=approx)
  for(a in 1:approx){
    I = sample.int(d,3)
    DT[,a] = depth.halfspace(A[,I],B[,I],exact=TRUE)  
  }
  D = apply(DT,1,mean)	# integrated depth
  DI = apply(DT,1,min)	# infimal depth
  IA = apply(DT,1,function(x) sum(x==min(x)))/approx	# infimal area
  return(list(Half_FD=D,Half_ID=DI,Half_IA=IA))
}

#' @title Bivariate h-mode depth for functional data based on the \eqn{L^2} metric
#'
#' @description
#' The h-mode depth 
#' of functional bivariate data (that is, data of the form \eqn{X:[a,b] \to R^2},
#' or \eqn{X:[a,b] \to R} and the derivative of \eqn{X}) based on the
#' \eqn{L^2} metric of functions.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth functional derivatives
#'
#' @details
#' The function returns the vectors of sample h-mode depth values. The kernel 
#' used in the evaluation is the standard Gaussian kernel, the bandwidth value is chosen
#' as a quantile of the non-zero distances between the random sample curves.
#'
#' @param datafA Bivariate functions whose depth is computed, represented by a multivariate \code{dataf} object of 
#' their arguments (vector), and a matrix with two columns of the corresponding bivariate functional values. 
#' \code{m} stands for the number of functions.
#'
#' @param datafB Bivariate random sample functions with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a multivariate \code{dataf} object of their arguments
#'  (vector), and a matrix with two columns of the corresponding bivariate functional values.
#' \code{n} is the sample size. The grid of observation points for the 
#' functions \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @param q The quantile used to determine the value of the bandwidth \eqn{h}
#' in the computation of the h-mode depth. \eqn{h} is taken as the \code{q}-quantile of
#' all non-zero distances between the functions \code{B}. By default, this value is set
#' to \code{q=0.2}, in accordance with the choice of Cuevas et al. (2007).
#'
#' @return Three vectors of length \code{m} of h-mode depth values are returned:
#'	\itemize{
#'	\item \code{hM} the unscaled h-mode depth,
#'	\item \code{hM_norm} the h-mode depth \code{hM} linearly transformed so that its range is [0,1],
#'	\item \code{hM_norm2} the h-mode depth \code{FD} linearly transformed by a transformation such that 
#'    the range of the h-mode depth of \code{B} with respect to \code{B} is [0,1]. This depth may give negative values.
#'	}
#'
#' @references Cuevas, A., Febrero, M. and Fraiman, R.  (2007). 
#' Robust estimation and classification for functional data via projection-based depth notions. 
#' \emph{Computational Statistics} \bold{22} (3), 481--496.
#' 
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#'
#' datafA2 = derivatives.est(datafA,deriv=c(0,1))
#' datafB2 = derivatives.est(datafB,deriv=c(0,1))
#' 
#' depthf.hM2(datafA2,datafB2)
#'
#' depthf.hM2(datafA2,datafB2)$hM
#' # depthf.hM2(cbind(A2[,,1],A2[,,2]),cbind(B2[,,1],B2[,,2]))$hM
#' # the two expressions above should give the same result
#'
#' @seealso \code{\link{depthf.hM}}

depthf.hM2 = function(datafA,datafB,range=NULL,d=101,q=.2){
  # h-Mode depth Cuevas_etal2007
  # A functions whose depth is computed :	either M*D matrix (M functions, D points per function), 
  #							or M*D*2 array (2 derivative levels), 
  # B functions of random sample :		as A
  # q :							quantile, bandwidth value in the resulting kernel estimate
  #							computes the same as depthf.M (with cbind-ed 2 levels if derivatives are involved)
  #							depthf.M(cbind(A[,,1],A[,,2]),cbind(B[,,1],B[,,2]))
  #							depthf.M2(A,B)$FD
  #							depthf.M2(cbind(A[,,1],A[,,2]),cbind(B[,,1],B[,,2]))$FD
  # 							all should give the same result
  
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)
  
  # q=.2 Cuevas et al, 0.15 default v usc.fda
  M = .Fortran(	"funMD",
                as.numeric(A),
                as.numeric(B),
                as.integer(m<-nrow(A)),
                as.integer(n<-nrow(B)),
                as.integer(d<-length(as.numeric(A))/nrow(A)),
                as.numeric(q),
                md = as.numeric(rep(-1,m))
  )$md
  #	because of scaling in depth.mode() in fda.usc
  M.ori = .Fortran(	"funMD",
                    as.numeric(B),
                    as.numeric(B),
                    as.integer(n),
                    as.integer(n),
                    as.integer(d),
                    as.numeric(q),
                    md = as.numeric(rep(-1,n))
  )$md
  Mn = (M-min(M))/(max(M)-min(M))
  Mn2 = (M-min(M.ori))/(max(M.ori)-min(M.ori))
  return(list(hM = M, hM_norm = Mn, hM_norm2 = Mn2))
}

#' @title Univariate random projection depths for functional data
#'
#' @description
#' Random projection depth and random functional depth for functional data.	
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth functional
#'
#' @details
#' The function returns the vectors of sample random projection, and random functional depth values. 
#' The random projection depth described in Cuevas et al. (2007) is based on the average univariate depth
#' of one-dimensional projections of functional data. The projections are taken randomly as a sample of standard
#' normal \code{d}-dimensional random variables, where \code{d} stands for the dimensionality of the discretized
#' functional data. 
#'
#' The random functional depth (also called random Tukey depth, or random halfspace depth) is described in
#' Cuesta-Albertos and Nieto-Reyes (2008). The functional data are projected into the real line in random 
#' directions as for the random projection depths. Afterwards, an approximation of the halfspace (Tukey) depth
#' based on this limited number of univariate projections is assessed. 
#'
#' @param datafA Functions whose depth is computed, represented by a \code{dataf} object of their arguments
#'  and functional values. \code{m} stands for the number of functions.
#'
#' @param datafB Random sample functions with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a \code{dataf} object of their arguments
#'  and functional values. \code{n} is the sample size. 
#' The grid of observation points for the 
#' functions \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @param nproj Number of projections taken in the computation of the random projection depth. By default taken
#' to be \code{51}.
#'
#' @param nproj2 Number of projections taken in the computation of the random functional depth. By default taken
#' to be \code{5}. \code{nproj2} should be much smaller than \code{d}, the dimensionality of the discretized 
#' functional data.
#'
#' @return Three vectors of depth values of length \code{m} are returned:
#'	\itemize{
#'	\item \code{Simpl_FD} the random projection depth based on the univariate simplicial depth,
#'	\item \code{Half_FD} the random projection depth based on the univariate halfspace depth,
#'	\item \code{RHalf_FD} the random halfspace depth.
#'	}
#'
#' @references Cuevas, A., Febrero, M. and Fraiman, R.  (2007).
#' Robust estimation and classification for functional data via projection-based depth notions, 
#' \emph{Computational Statistics} \bold{22} (3), 481--496.
#'
#' @references Cuesta-Albertos, J.A. and Nieto-Reyes, A. (2008).
#'  The random Tukey depth. 
#' \emph{Computational Statistics & Data Analysis} \bold{52} (11), 4979--4988.
#' 
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#'
#' depthf.RP1(datafA,datafB)
#'
#' @seealso \code{\link{depthf.RP2}}

depthf.RP1 = function(datafA,datafB,range=NULL,d=101,nproj=50,nproj2=5){
  # simple 1D projection depth
  # nproj nr of projections taken
  #
  #	SFD : projection \eqn{L^2} -> R -> mean of univariate Simplicial Depth (nproj projections)
  #	HFD : projection \eqn{L^2} -> R -> mean of univariate Halfspace Depth (nproj projections)
  #	RFD : projection \eqn{L^2} -> R -> infimum of univariate Halfspace Depths (nproj2 projections)
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)  
  
  A1 = as.vector(A)
  B1 = as.vector(B)
  d = dim(A)[2]
  m = dim(A)[1]
  n = dim(B)[1]
  V = rnorm(nproj*d)
  if(dim(B)[2]!=d) stop("dimension mismatch")
  FD = .Fortran("funRPD1",	
                as.numeric(A1),		#	A
                as.numeric(B1),		#	B
                as.integer(m),		#	m
                as.integer(n),		#	n
                as.integer(d),		#	d
                as.integer(nproj),	#	nproj
                as.integer(nproj2),	#	nproj2
                as.numeric(V),		#	projections
                funsdep=as.numeric(rep(-1,m)),	
                funhdep=as.numeric(rep(-1,m)),
                funrdep=as.numeric(rep(-1,m)))
  return(list(Simpl_FD = FD$funsdep, Half_FD = FD$funhdep, RHalf_FD = FD$funrdep))
}

#' @title Bivariate integrated and infimal depth for functional data
#'
#' @description
#' Integrated and infimal depths 
#' of functional bivariate data (that is, data of the form \eqn{X:[a,b] \to R^2},
#' or \eqn{X:[a,b] \to R} and the derivative of \eqn{X}) based on the
#' bivariate halfspace and simplicial depths.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth functional derivatives
#'
#' @details
#' The function returns the vectors of sample integrated and infimal depth values.
#'
#' @param datafA Bivariate functions whose depth is computed, represented by a multivariate \code{dataf} object of 
#' their arguments (vector), and a matrix with two columns of the corresponding bivariate functional values. 
#' \code{m} stands for the number of functions.
#'
#' @param datafB Bivariate random sample functions with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a multivariate \code{dataf} object of their arguments
#'  (vector), and a matrix with two columns of the corresponding bivariate functional values.
#' \code{n} is the sample size. The grid of observation points for the 
#' functions \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @return Four vectors of length \code{m} are returned:
#'	\itemize{
#'	\item \code{Simpl_FD} the integrated depth based on the bivariate simplicial depth,
#'	\item \code{Half_FD} the integrated depth based on the bivariate halfspace depth,
#'	\item \code{Simpl_ID} the infimal depth based on the bivariate simplicial depth,
#'	\item \code{Half_ID} the infimal depth based on the bivariate halfspace depth.
#'	}
#' In addition, two vectors of length \code{m} of the relative area of smallest depth values is returned:
#'	\itemize{
#'	\item \code{Simpl_IA} the proportions of points at which the depth \code{Simpl_ID} was attained,
#'	\item \code{Half_IA} the proportions of points at which the depth \code{Half_ID} was attained.
#'	}
#' The values \code{Simpl_IA} and \code{Half_IA} are always in the interval [0,1]. 
#' They introduce ranking also among functions having the same
#' infimal depth value - if two functions have the same infimal depth, the one with larger infimal area
#' \code{IA} is said to be less central. 	
#'
#' @references Hlubinka, D., Gijbels, I., Omelka, M. and Nagy, S. (2015). 
#' Integrated data depth for smooth functions and its application in supervised classification. 
#' \emph{Computational Statistics}, \bold{30} (4), 1011--1031.
#'
#' @references Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
#' Depth-based recognition of shape outlying functions. 
#' \emph{Journal of Computational and Graphical Statistics}, \bold{26:4}, 883--893.
#'  
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#'
#' dataf2A = derivatives.est(datafA,deriv=c(0,1))
#' dataf2B = derivatives.est(datafB,deriv=c(0,1))
#' depthf.fd2(dataf2A,dataf2B)
#'
#' @seealso \code{\link{depthf.fd1}}, \code{\link{infimalRank}}

depthf.fd2 = function(datafA,datafB,range=NULL,d=101){
  # A functions whose depth I compute
  # B functions wrt whose the depth is computed
  # both 2dimensional, n*d*2, n nr of functions, d dimensionality
  # now provides also infimal depth (inf_D2)
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)
  
  A1 = as.vector(A[,,1])
  A2 = as.vector(A[,,2])
  B1 = as.vector(B[,,1])
  B2 = as.vector(B[,,2])
  d = dim(A)[2]
  m = dim(A)[1]
  n = dim(B)[1]
  if(dim(B)[2]!=d) stop("dimension mismatch")
  FD = .Fortran("funD2",	
                as.numeric(A1),		#	A1
                as.numeric(A2),		#	A2
                as.numeric(B1),		#	B1
                as.numeric(B2),		#	B2
                as.integer(m),		#	m
                as.integer(n),		#	n
                as.integer(d),		#	d
                funsdep=as.numeric(rep(-1,m)),	
                funhdep=as.numeric(rep(-1,m)),
                fIsdep =as.numeric(rep(-1,m)),
                fIhdep =as.numeric(rep(-1,m)),
                IAsdep =as.integer(rep(-1,m)),
                IAhdep =as.integer(rep(-1,m))
  )
  return(list(	Simpl_FD = FD$funsdep, Half_FD = FD$funhdep,
               Simpl_ID = FD$fIsdep,	Half_ID = FD$fIhdep,
               Simpl_IA = FD$IAsdep/d, Half_IA = FD$IAhdep/d ))
}

#' @title Bivariate random projection depths for functional data
#'
#' @description
#' Double random projection depths of functional bivariate data (that is, data of the form \eqn{X:[a,b] \to R^2},
#' or \eqn{X:[a,b] \to R} and the derivative of \eqn{X}).	
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth functional derivatives
#'
#' @details
#' The function returns the vectors of sample double random projection depth values. 
#' The double random projection depths are described in Cuevas et al. (2007). They are of two types: RP2 type, and
#' RPD type. Both types of depths are based on bivariate projections of the bivariate functional data. 
#' These projections are taken randomly as a sample of standard
#' normal \code{d}-dimensional random variables, where \code{d} stands for the dimensionality of the internally 
#' represented discretized
#' functional data. For RP2 type depths, the average bivariate depth of the projected quantities is assessed.
#' For RPD type depths, further univariate projections of these bivariate projected quantities are evaluated, and
#' based on these final univariate quantities, the average univariate depth is computed.
#' 
#' @param datafA Bivariate functions whose depth is computed, represented by a multivariate \code{dataf} object of 
#' their arguments (vector), and a matrix with two columns of the corresponding bivariate functional values. 
#' \code{m} stands for the number of functions.
#'
#' @param datafB Bivariate random sample functions with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a multivariate \code{dataf} object of their arguments
#'  (vector), and a matrix with two columns of the corresponding bivariate functional values.
#' \code{n} is the sample size. The grid of observation points for the 
#' functions \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @param nproj Number of projections taken in the computation of the double random projection depth. By default taken
#' to be \code{51}.
#'
#' @return Five vectors of length \code{m} are returned:
#'	\itemize{
#'	\item \code{Simpl_FD} the double random projection depth RP2 based on the bivariate simplicial depth,
#'	\item \code{Half_FD} the double random projection depth RP2 based on the bivariate halfspace depth,
#'	\item \code{hM_FD} the double random projection depth RP2 based on the bivariate h-mode depth,
#'	\item \code{Simpl_DD} the double random projection depth RPD based on the univariate simplicial depth,
#'	\item \code{Half_DD} the random projection depth RPD based on the univariate halfspace depth,
#'	}
#'
#' @references Cuevas, A., Febrero, M. and Fraiman, R.  (2007).
#' Robust estimation and classification for functional data via projection-based depth notions.
#' \emph{Computational Statistics} \bold{22} (3), 481--496.
#' 
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#'
#' dataf2A = derivatives.est(datafA,deriv=c(0,1))
#' dataf2B = derivatives.est(datafB,deriv=c(0,1))
#' depthf.RP2(dataf2A,dataf2B)
#'
#'
#' @seealso \code{\link{depthf.RP1}}

depthf.RP2 = function(datafA,datafB,range=NULL,d=101,nproj=51){
  # double projection depth
  # nproj nr of projections taken
  #
  #	SFD :	(\eqn{L^2})^2 -> R^2 -> 2D Mean Simplicial Depth
  #	HFD : (\eqn{L^2})^2 -> R^2 -> 2D Mean Halfspace Depth
  #	MFD : (\eqn{L^2})^2 -> R^2 -> 2D Mean h-Mode Depth
  #	SDD : (\eqn{L^2})^2 -> R^2 -> R^1 -> Mean Simplicial Depth
  #	HDD : (\eqn{L^2})^2 -> R^2 -> R^1 -> Mean Halfspace Depth
  #
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)
  
  q = .2	# quantile for modal depth
  A1 = as.vector(A[,,1])
  A2 = as.vector(A[,,2])
  B1 = as.vector(B[,,1])
  B2 = as.vector(B[,,2])
  d = dim(A)[2]
  m = dim(A)[1]
  n = dim(B)[1]
  V = rnorm(nproj*d+nproj*2)
  if(dim(B)[2]!=d) stop("dimension mismatch")
  FD = .Fortran("funRPD2",	
                as.numeric(A1),		#	A1
                as.numeric(A2),		#	A2
                as.numeric(B1),		#	B1
                as.numeric(B2),		#	B2
                as.integer(m),		#	m
                as.integer(n),		#	n
                as.integer(d),		#	d
                as.integer(nproj),	#	nproj
                as.numeric(V),		#	projections
                as.numeric(q),		# 	q
                funsdep=as.numeric(rep(-1,m)),	
                funhdep=as.numeric(rep(-1,m)),
                funmdep=as.numeric(rep(-1,m)),
                funsddep=as.numeric(rep(-1,m)),
                funhddep=as.numeric(rep(-1,m)))
  return(list(Simpl_FD = FD$funsdep, Half_FD = FD$funhdep, hM_FD = FD$funmdep, Simpl_DD = FD$funsddep, Half_DD = FD$funhddep))
}

#' @title Half-region depth for functional data
#'
#' @description
#' The half-region depth 
#' for functional real-valued data.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth functional
#'
#' @details
#' The function returns the vector of the sample half-region depth values.
#'
#' @param datafA Functions whose depth is computed, represented by a \code{dataf} object of their arguments
#'  and functional values. \code{m} stands for the number of functions.
#'
#' @param datafB Random sample functions with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a \code{dataf} object of their arguments
#'  and functional values. \code{n} is the sample size. 
#' The grid of observation points for the 
#' functions \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @return A vector of length \code{m} of the half-region depth values.
#'
#' @references Lopez-Pintado, S. and Romo, J. (2011).
#' A half-region depth for functional data.
#' \emph{Computational Statistics & Data Analysis} \bold{55} (4), 1679--1695.
#' 
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#' depthf.HR(datafA,datafB)

depthf.HR = function(datafA,datafB,range=NULL,d=101){
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)
  
  d = dim(A)[2]
  m = dim(A)[1]
  n = dim(B)[1]
  if(dim(B)[2]!=d) stop("dimension mismatch")
  FD = .Fortran("HRD",	
                as.numeric(A),		#	A
                as.numeric(B),		#	B
                as.integer(m),		#	m
                as.integer(n),		#	n
                as.integer(d),		#	d
                FD=as.numeric(rep(-1,m)))
  return(FD$FD)
}

#' @title Band depth for functional data
#'
#' @description
#' The (unadjusted) band depth 
#' for functional real-valued data of order \code{J=2}.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth functional
#'
#' @details
#' The function returns the vector of the sample (unadjusted) band depth values.
#'
#' @param datafA Functions whose depth is computed, represented by a \code{dataf} object of their arguments
#'  and functional values. \code{m} stands for the number of functions.
#'
#' @param datafB Random sample functions with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a \code{dataf} object of their arguments
#'  and functional values. \code{n} is the sample size. 
#' The grid of observation points for the 
#' functions \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @return A vector of length \code{m} of the band depth values.
#' 
#' @references Lopez-Pintado, S. and Romo, J. (2009), On the concept of depth for functional data,
#' \emph{J. Amer. Statist. Assoc.} \bold{104} (486), 718 - 734.
#'
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#' depthf.BD(datafA,datafB)
#'
#' @seealso \code{\link{depthf.ABD}}, \code{\link{depthf.fd1}}

depthf.BD = function(datafA,datafB,range=NULL,d=101){
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)  
  
  d = dim(A)[2]
  m = dim(A)[1]
  n = dim(B)[1]
  if(dim(B)[2]!=d) stop("dimension mismatch")
  FD = .Fortran("BD",	
                as.numeric(A),		#	A
                as.numeric(B),		#	B
                as.integer(m),		#	m
                as.integer(n),		#	n
                as.integer(d),		#	d
                FD=as.numeric(rep(-1,m))
  )
  return(FD$FD)
}

#' @title Adjusted ranking of functional data based on the infimal depth
#'
#' @description
#' Returns a vector of adjusted depth-based ranks for infimal depth for functional data.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords rank depth functional
#'
#' @details
#' Infimal depths for functional data tend to give to many functional observations the same 
#' value of depth. Using this function, the data whose depth is the same is ranked according
#' to the infimal area indicator. This indicator is provided in functions \code{depthf.fd1} along
#' the value of the infimal depth. 
#'
#' @param ID The vector of infimal depths of the curves of length \code{n}.
#'
#' @param IA The vector of the infimal areas corresponding to the infimal depths from \code{ID}
#' of length \code{n}.
#'
#' @param ties.method Parameter for breaking ties in infimal area index. By default \code{max}, see 
#' \code{rank}.
#' 
#' @return A vector of length \code{n}. Low depth values mean high ranks, i.e. potential outlyingess. 
#' If some of the infimal depths are identical, the ranking of these functions is made according to the 
#' values of the infimal area. There, higher infimal area index means higher rank, i.e. non-centrality.
#' 
#' @references Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
#' Depth-based recognition of shape outlying functions. 
#' \emph{Journal of Computational and Graphical Statistics}, \bold{26:4}, 883--893. 
#' 
#' @examples
#' datafA = dataf.population()$dataf[1:20]
#' datafB = dataf.population()$dataf[21:50]
#' D = depthf.fd1(datafA,datafB)
#' infimalRank(D$Half_ID,D$Half_IA) 
#' 
#' ID = c(0,1,0,0,0,1,1)
#' IA = c(2,3,1,0,2,4,1)
#' infimalRank(ID,IA)


infimalRank = function(ID,IA,ties.method="max"){
  # finds the adjusted rank for appropriate for the infimal depth for functional data
  # ID is the vector of infimal depths
  # IA is the vector of the corresponding infimal areas
  # returns a vector of ranks
  n = length(ID)
  if(length(IA)!=n) stop("Lengths of the vectors differ")
  U = sort(unique(ID),decreasing=TRUE)
  R = rep(NA,n)
  cR = 0						# currently assigned rank	
  for(u in U){
    Iu = (ID==u)
    R[Iu] = cR+rank(IA[Iu],ties.method=ties.method)
    cR = sum(!is.na(R))
  }
  return(R)
}

#' @title Estimation of the first two derivatives for functional data
#'
#' @description
#' Returns the estimated values of derivatives of functional data.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords derivatives kernel functional
#' 
#' @details
#' If the input \code{dataf} is a functional random sample of size \code{m}, 
#' the function returns a \code{dataf} object of \code{nd}-dimensional functional data, where 
#' in the elements of the vector-valued functional data represent the estimated values of the 
#' derivatives of \code{dataf}. All derivatives are evaluated at an equi-distant grid of \code{d}
#' points in the domain given by \code{range}. \code{nd} here stands for \code{1}, \code{2} or \code{3}, 
#' depending on how many derivatives of \code{dataf} are
#' requested to be computed. For the estimation, functions \code{D1ss} and \code{D2ss} from the package
#' \code{sfsmisc} are utilized. 
#'
#' @param dataf Functional dataset, represented by a \code{dataf} object of their arguments
#'  and functional values. \code{m} stands for the number of functions.
#'  
#' @param range The common range of the domain where the fucntions \code{dataf} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{dataf}.
#' 
#' @param d Grid size to which all the functional data are transformed. For computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @param spar If provided, this parameter is passed to functions \code{D1ss} and \code{D2ss} from package \code{sfsmisc}
#' as the value of the smoothing spline parameter in order to numerically approximate
#' the derivatives of \code{dataf}.
#' 
#' @param deriv A vector composed of \code{0}, \code{1}, and \code{2} of the demanded 
#' functional values / derivatives of the functions in the rows of \code{dataf}.
#' \code{0} stands for the functional values, \code{1} for the first derivatives, 
#' \code{2} for the second derivatives.
#' 
#' @return A multivariate \code{dataf} object of the functional values and / or the derivatives of \code{dataf}. 
#' The dimensionality of the vector-valued functional data is \code{nd}. The arguments of the data are all equal to 
#' an equi-distant grid of \code{d} points in the domain given by \code{range}. \code{nd} is the demanded number 
#' of derivatives at the output, i.e. the length of the vector \code{deriv}.
#' 
#' @seealso \code{\link[sfsmisc]{D1ss}} in package sfsmisc
#' @seealso \code{\link[sfsmisc]{D2ss}} in package sfsmisc
#' 
#' @examples
#' dataf = dataf.population()$dataf
#' derivatives.est(dataf,deriv=c(0,1,2))

derivatives.est = function(dataf,range=NULL,d=101,spar=NULL,deriv=c(0,1)){
  
  X = dataf2rawfd(dataf, range = range, d = d)
  derd = length(deriv)
  XK = array(dim=c(dim(X),derd))
  # range construction
  rng = numeric(0)
  for (df in dataf)	rng = range(c(rng,df$args)) # common range of all the data
  if(!is.null(range)){ 
    if(!(length(range) == 2 && is.numeric(range) && 
         range[1]<=rng[1] && range[2]>=rng[2])) 
      stop("Argument 'range' must be a numeric vector of two components
           that defines the range of the domain of functional data. All
           functional data must have 'args' vales inside this domain.")
  } else range = rng
  if(!(range[1]<range[2])) stop("Argument 'range' must define a non-degenerate interval.")
  
  t = seq(range[1],range[2],length=d)
  n = nrow(X)	
  pb = txtProgressBar(min = 0, max = n, style = 3)	
  for (nd in 1:derd){
    if (deriv[nd]==0) XK[,,nd] = X else
      for(i in 1:n){
        if (deriv[nd]==1){
          if(is.null(spar))	XK[i,,nd] = D1ss(t,X[i,]) else
            XK[i,,nd] = D1ss(t,X[i,],spl.spar=spar)
        }
        if (deriv[nd]==2){
          if(is.null(spar))	XK[i,,nd] = D2ss(t,X[i,])$y else
            XK[i,,nd] = D2ss(t,X[i,],spl.spar=spar)$y
        }
        setTxtProgressBar(pb, i)
      }
  }
  close(pb)
  return(rawfd2dataf(XK,range=range))
  }

#' @title Diagnostic plot for first and second order integrated and infimal depths
#'
#' @description
#' Produce the diagnostic plot based on the fist or second order extended integrated / infimal depths.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth shape outlier plot functional
#'
#' @details
#' Plots a diagnostic plot of pointwise univariate (or bivariate) depths for all possible points (or couples of points) from the domain of the 
#' functional data. From such a plot it is possible to infer into the first order (or second order) properties of a single function \emph{x} with respect 
#' to the given set of functional data. For \code{order=1}, the integral of the displayed function is the integrated depth of \emph{x}, 
#' the smallest value of the function is the infimal depth of \emph{x}. 
#' For \code{order=2}, the bivariate integral of the displayed surface gives the second order extended 
#' integrated depth of \emph{x}, the infimum of this bivariate function gives the second order infimal depth of \emph{x}. 
#' For details see Nagy et al. (2016) and \code{\link{depthf.fd1}}.
#'
#' @param datafA A single function whose depth is computed, represented by a 
#' \code{dataf} object of arguments and functional values.
#'
#' @param datafB Functional dataset with respect to which the depth of \code{datafA} is computed. 
#' \code{datafB} is represented by a \code{dataf} object of arguments and functional values. 
#' \code{n} stands for the number of functions. The grid of observation points for the 
#' functions in \code{datafA} and \code{datafB} may not be the same.
#' 
#' @param range The common range of the domain where the fucntions \code{datafA} and \code{datafB} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{datafA} and \code{datafB}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @param order The order of the depth to be used in the plot, for \code{order=1} produces
#' the plot of univariate marginal depth of \code{A} and \code{nfun} functions from \code{B} 
#' over the domain of the functions. For \code{order=2} produces the bivariate contour plot 
#' of the bivariate depths of \code{A} at couples of points from the domain.
#'
#' @param method The depth that is used in the diagnostic plot. possible values are \code{halfspace} for 
#' the halfspace depth, or \code{simplicial} for the simplicial depth. 
#'
#' @param approx For \code{order=2}, the number of approximations used in the computation of the order extended depth. By default
#' this is set to \code{0}, meaning that the depth is computed at all possible \code{d^2}
#' combinations of the points in the domain. When set to a positive integer, \code{approx}
#' bivariate points are randomly sampled in unit square, and at these points the bivariate depths of the
#' corresponding functional values are computed.
#'
#' @param title The title of the diagnostic plot.
#'
#' @param nfun For \code{order=1}, the number of functions from \code{B} whose coordinate-wise
#' univariate depths of functional values should be displayes with the depth of \code{A}.
#' The depth of \code{A} is displayed in solid red line, the depths of the functions from \code{B}
#' in dashed black.
#'
#' @param plot Logical: should the function by plotted?
#' 
#' @return For \code{order=1} two depth values, and two vectors of pointwise depths:
#'	\itemize{
#'	\item \code{Simpl_FD} the first order integrated depth based on the simplicial depth,
#'	\item \code{Half_FD} the first order integrated depth based on the halfspace depth,
#'	\item \code{Simpl_ID} the first order infimal depth based on the simplicial depth,
#'	\item \code{Half_ID} the first order infimal depth based on the halfspace depth,
#'    \item \code{PSD} the vector of length \code{d} containing the computed 
#'    pointwise univariate simplicial depths used for the computation of \code{Simpl_FD} and \code{Simpl_ID},
#'    \item \code{PHD} the vector of length \code{d} containing the computed 
#'    pointwise univariate halfspace depths used for the computation of \code{Half_FD} and \code{Half_ID}.
#'	}
#'    In addition, the first order integrated / infimal depth diagnostic plot of the function \code{A} with respect to
#'    the random sample given by the functions corresponding to the rows of the matrix \code{B} is produced.
#'
#'    For \code{order=2} four depth values, and two matrices of pointwise depths:
#'	\itemize{
#'	\item \code{Simpl_FD} the second order integrated depth based on the simplicial depth,
#'	\item \code{Half_FD} the second order integrated depth based on the halfspace depth,
#'	\item \code{Simpl_ID} the second order infimal depth based on the simplicial depth,
#'	\item \code{Half_ID} the second order infimal depth based on the halfspace depth,
#'    \item \code{PSD} the matrix of size \code{d*d} containing the computed 
#'    pointwise bivariate simplicial depths used for the computation of \code{Simpl_FD} and \code{Simpl_ID},
#'    \item \code{PHD} the matrix of size \code{d*d} containing the computed 
#'    pointwise bivariate halfspace depths used for the computation of \code{Half_FD} and \code{Half_ID}.
#'	}
#'    In addition, the second order integrated / infimal depth diagnostic plot of the function \code{A} with respect to
#'    the random sample given by the functions corresponding to the rows of the matrix \code{B} is produced.
#' 
#' @references Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
#' Depth-based recognition of shape outlying functions. 
#' \emph{Journal of Computational and Graphical Statistics}, \bold{26:4}, 883--893.
#'
#' @examples
#' datafA = dataf.population()$dataf[1]
#' dataf = dataf.population()$dataf[2:20]
#' shape.fd.analysis(datafA,dataf,order=1)
#' shape.fd.analysis(datafA,dataf,order=2,approx=0)
#'
#' @seealso \code{\link{depthf.fd1}}

shape.fd.analysis = function(datafA,datafB,range=NULL,d=101,order=1,method=c("halfspace","simplicial"),approx=0,title="",nfun=10,plot=TRUE){
  
  A = dataf2rawfd(datafA, range = range, d = d)
  B = dataf2rawfd(datafB, range = range, d = d)  
  
  if (nrow(A)>1) stop("works only for a single function A (matrix 1*d)")
  method = match.arg(method)
  #
  # diagonal
  d = ncol(A)
  HDP = matrix(nrow=nrow(B)+1,ncol=d)
  SDP = HDP
  for(i in 1:d){
    FL = ecdf(B[,i])
    FU = ecdf(-B[,i])
    HDP[,i] = pmin(FL(rbind(B,A)[,i]),FU(-rbind(B,A)[,i]))
    SDP[,i] = FL(rbind(B,A)[,i])*FU(-rbind(B,A)[,i])*nrow(B)^2/choose(nrow(B),2)
  }
  if(order==2){
    DD = DiffDepth(A,B,approx)
    if (method=="halfspace") diag(DD$PHD) = HDP[nrow(B)+1,] else diag(DD$PSD) = SDP[nrow(B)+1,]
    if (plot){
      if (method=="halfspace") filled.contour(DD$PHD,main=title,color.palette= function(x)rev(heat.colors(x))) else
        filled.contour(DD$PSD,main=title,color.palette = function(x)rev(heat.colors(x)))
    }
    return(DD)
  }
  if(order==1){
    nfun = min(nrow(B),nfun)
    t = seq(0,1,length=ncol(B))
    DD = list(Simpl_FD = mean(SDP[nrow(B)+1,]), Half_FD = mean(HDP[nrow(B)+1,]), 
              Simpl_ID = min(SDP[nrow(B)+1,]), Half_ID = min(HDP[nrow(B)+1,]), 
              PHD = HDP[nrow(B)+1,], PSD = SDP[nrow(B)+1,]) 
    if (plot) { if (method=="halfspace"){
      plot(rep(t,nrow(B)+1),HDP,type="n",main=title,ylim = c(0,max(HDP)),ann=FALSE)
      for(i in 1:nfun) lines(t,HDP[i,],lwd=.75,lty=2)
      lines(t,HDP[nrow(B)+1,],col=2,lwd=3,lty=1)
    } else {
      plot(rep(t,nrow(B)+1),SDP,type="n",main=title,ylim = c(0,max(SDP)),ann=FALSE)
      for(i in 1:nfun) lines(t,SDP[i,],lwd=.75,lty=2)
      lines(t,SDP[nrow(B)+1,],col=2,lwd=3,lty=1)
    }
    }	
    return(DD)
  }
}

#' @title Functional depth-based shape outlier detection
#'
#' @description
#' Detects functional outliers of first three orders, based on the order extended integrated depth for functional data.		
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#' @keywords depth
#' @keywords outlier
#' @keywords functional
#'
#' @details
#' Using the procedure described in Nagy et al. (2016), the function uses the order extended integrated depths for functions, 
#' see \code{\link{depthf.fd1}} and \code{\link{shape.fd.analysis}}, to perform informal functional shape outlier detection. 
#' Outliers of the first order (horizontal shift outliers) are found as the functions with \code{q} \% of smallest (first order)
#' integrated depth values. Second and third order outliers (shape outliers) are found using the extension of the boxplot method
#' for depths as described in the paper Nagy et al. (2016).
#'
#' @param dataf Functional dataset, represented by a \code{dataf} object of their arguments
#'  and functional values. \code{n} stands for the number of functions.
#'  
#' @param range The common range of the domain where the fucntions \code{dataf} are observed.
#' Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
#' \code{dataf}.
#' 
#' @param d Grid size to which all the functional data are transformed. For depth computation, 
#' all functional observations are first transformed into vectors of their functional values of length \code{d}
#' corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
#' points are reconstructed using linear interpolation, and extrapolation.
#'
#' @param q The quantile presenting a threshold for the first order outlier detection. Functions with first order integrated depth
#' smaller than the \code{q} quantile of this sample of depths are flagged as potential outliers. If set to \code{NULL}, the
#' the outliers are detected from the first order integrated depth after the log-transformation, as for higher order outliers.
#'
#' @param method The depth that is used in the diagnostic plot. possible values are \code{halfspace} for 
#' the halfspace depth, or \code{simplicial} for the simplicial depth. 
#'
#' @param approx For the computation of the third order integrated depth,
#' the number of approximations used in the computation of the order extended depth. By default
#' this is set to \code{100}, meaning that \code{100}
#' trivariate points are randomly sampled in unit cube, and at these points the trivariate depths of the
#' corresponding functional values. May be set to \code{0} to compute the depth at all possible \code{d^3}
#' combinations of the points in the domain. This choice may result in very slow computation, see also \code{\link{depthf.fd1}}. 
#'
#' @param print If the rows of \code{X} are named, \code{print=TRUE} enables a graphical output when the names of the outlying curves
#' are displayed.
#'
#' @param plotpairs If set to \code{TRUE}, the scatter plot of the computed depths for orders \code{1}, \code{2} and \code{3} is
#' is displayed. Here, the depths corresponding to the flagged outliers are plotted in colour.
#'
#' @param max.order Maximal order of shape outlyingness to be computed, can be set to \code{1}, \code{2}, or \code{3}.
#'
#' @param exclude.out Logical variable; exclude the detected lower order outliers in the flagging process? By default \code{TRUE}.
#'
#' @param output Output method, can be set to \code{matrix} for a matrix with logical entries (\code{TRUE} for outliers), or \code{list} for 
#' a list of outliers.
#' 
#' @param identifiers A vector of names for the data observations. Facilitates identification of outlyig functions.
#' 
#' @return A matrix of logical values of size \code{n*4}, where \code{n} is the sample size. In the first three rows indicators of outlyingness
#' of the corresponding functions for orders \code{1}, \code{2} and \code{3} are given, in the fourth row the indicator of outlyingness
#' with respect to the comparison of the first, and third order depths is given. That is, the fist row corresponds to the first order outliers, 
#' the second row to the second order outliers, and the last two rows formally to the third order outliers. Please consult Nagy et al. (2016)
#' to interpret the notion of shape outlyingness.
#' 
#' @references Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
#' Depth-based recognition of shape outlying functions. 
#' \emph{Journal of Computational and Graphical Statistics}, \bold{26:4}, 883--893.
#'
#' @examples
#' n = 30
#' dataf = dataf.population()$dataf[1:n]
#' shape.fd.outliers(dataf,print=TRUE,plotpairs=TRUE,
#' identifiers=unlist(dataf.population()$identifier)[1:n])
#'
#' @seealso \code{\link{depthf.fd1}}, \code{\link{shape.fd.analysis}}

shape.fd.outliers = function(dataf,range=NULL,d=101,q=.05,method=c("halfspace","simplicial"),
                             approx=100,print=FALSE,plotpairs=FALSE,max.order=3,
                             exclude.out = TRUE, output=c("matrix","list"), identifiers = NULL){
  
  X = dataf
  
  method = match.arg(method)
  output = match.arg(output)
  if(is.null(identifiers) | length(identifiers)!=length(dataf)){
    print = FALSE
    warning("Inconsistent identifiers, print is set to FALSE")
  }
  if((max.order>3)|(max.order<1)){
    max.order=3
    warning("Maximal order set to 3")
  }
  if (max.order==3){
    if(approx<50){
      approx=50
      warning("Too small approx value, approximation set to 50")
    }
  }
  n = length(X)									# set up the depths
  D = matrix(nrow=max.order,ncol=n)
  if(method=="halfspace"){
    D[1,] = depthf.fd1(X,X)$Half_FD
    if(max.order>1) D[2,] = depthf.fd1(X,X,range=range,d=d,order=2)$Half_FD
    if(max.order>2) D[3,] = depthf.fd1(X,X,range=range,d=d,order=3,approx=approx)$Half_FD
  }
  if(method=="simplicial"){
    D[1,] = depthf.fd1(X,X)$Simpl_FD
    if(max.order>1) D[2,] = depthf.fd1(X,X,range=range,d=d,order=2)$Simpl_FD
    if(max.order>2) D[3,] = depthf.fd1(X,X,range=range,d=d,order=3,approx=approx)$Simpl_FD
  }
  D = apply(D,1:2,function(x) max(0,x-1/n))
  if(max.order==1) out.rows = 1
  if(max.order==2) out.rows = 2
  if(max.order==3) out.rows = 4
  O = matrix(nrow=out.rows,ncol=ncol(D))					# compute the outliers
  if (print){
    colnames(O) = identifiers
    print("first order outliers: ")
  }
  if(!is.null(q)) O[1,]<-D[1,]<quantile(D[1,],q)
  if(is.null(q)){                  # for q null perform the log-transform
    S = -log(D[1,]/1)
    S[D[1,]==0] = Inf
    B = boxplot(S,plot=FALSE)
    O[1,]<-((S>B$stats[5,])|(S==Inf))            
  } 	
  if (print) print(which(O[1,]))
  if(max.order>1) for(i in 1:(max.order-1)){
    if (print){
      if (i==1) print("second order outliers: ")
      if (i==2) print("third order outliers: ")
    }
    S = -log(D[i+1,]/D[i,])
    S[D[i+1,]==0] = Inf
    B = boxplot(S,plot=FALSE)
    if(exclude.out) O[i+1,]<-(((S>B$stats[5,])|(S==Inf))&(O[i,]==FALSE)&(O[1,]==FALSE))
    if(!exclude.out) O[i+1,]<-((S>B$stats[5,])|(S==Inf))
    if (print) print(which(O[i+1,]))
  }
  if(max.order==3){
    if (print) print("1/3 outliers")
    S = -log(D[3,]/D[1,])
    S[D[3,]==0] = Inf
    B = boxplot(S,plot=FALSE)
    if(exclude.out) O[4,]<-(((S>B$stats[5,])|(S==Inf))&(O[3,]==FALSE)&(O[2,]==FALSE)&(O[1,]==FALSE))
    if(!exclude.out) O[4,]<-((S>B$stats[5,])|(S==Inf))
    if (print) print(which(O[4,]))
  }
  if(max.order==1) rownames(O) = "1st"
  if(max.order==2) rownames(O) = c("1st","2nd")
  if(max.order==3){
    rownames(O) = c("1st","2nd","3rd","1/3rd")
    if (plotpairs) DpairsPlot(D,O)
  }
  if(output=="matrix") return(O)
  if(output=="list") return(apply(O,1,which))
}

DpairsPlot = function(DB2,O,sp.index=NULL){
  # plots a 2x2 scatter of all the pairs of order extended integrated depth values
  # for orders 1, 2, and 3 as presented in the Nagy et al. (2016), Figure 5
  cexaxis = 1.3
  cexlab = 1.5
  n = ncol(O)
  col = rep("grey",n)
  colout1 = "olivedrab3"
  colout2 = "navy"
  colout3 = "darkorange"
  colout = c(colout1,colout2,colout3)
  col[O[1,]] = colout1
  col[O[2,]] = colout2
  col[O[3,]] = colout3
  col[O[4,]] = colout3
  pch = rep(16,n)
  pch[O[1,]] = 1
  pch[O[2,]] = 2
  pch[O[3,]] = 18
  pch[O[4,]] = 18
  cx = 1.5
  cex = rep(1,n)
  cex[O[1,]] = cx
  cex[O[2,]] = cx
  cex[O[3,]] = cx
  cex[O[4,]] = cx
  OI = (colSums(O)>0)	# outlier indicator
  
  op<-par(cex.axis=cexaxis,cex.lab=cexlab,mfrow = c(2, 2),     # 2x2 layout
          oma = c(3, 3.5, 0.45, 0.45), # two rows of text at the outer left and bottom margin
          mar = c(1, 1, 0, 0), # space for one row of text at ticks and to separate plots
          mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
          xpd = NA)            # allow content to protrude into outer margin (and beyond)
  plot(c(0,M<-max(DB2)),c(0,max(DB2)),type="n",ann=FALSE,axes=FALSE,frame=TRUE)
  axis(2)
  title(ylab=expression(FD[2]),line=2.5)
  points(DB2[1,],DB2[2,],col=col,pch=pch,cex=cex,lwd=1.35*cex)
  for(i in 1:n) if(OI[i]) points(DB2[1,i],DB2[2,i],col=col[i],pch=pch[i],cex=cex[i],lwd=1.35*cex[i])
  if (!is.null(sp.index)) points(DB2[1,sp.index],DB2[2,sp.index],pch=16,col="orange")
  segments(0,0,M,M,lty=3,lwd=2)
  plot(c(0,1),c(0,1),type="n",ann=FALSE,axes=FALSE)
  legend(0,.7,c("1st order outliers","2nd order outliers","3rd order outliers"),col=colout,pch=c(1,2,18),cex=cx,bty="n")
  plot(c(0,max(DB2)),c(0,max(DB2)),type="n",ann=FALSE)
  axis(2)
  title(ylab=expression(FD[3]^A),line=2.5)
  axis(1)
  title(xlab=expression(FD[1]),line=3)
  points(DB2[1,],DB2[3,],col=col,pch=pch,cex=cex,lwd=1.35*cex)
  for(i in 1:n) if(OI[i]) points(DB2[1,i],DB2[3,i],col=col[i],pch=pch[i],cex=cex[i],lwd=1.35*cex[i])
  if (!is.null(sp.index)) points(DB2[1,sp.index],DB2[3,sp.index],pch=16,col="orange")	
  segments(0,0,M,M,lty=3,lwd=2)
  plot(c(0,max(DB2)),c(0,max(DB2)),type="n",ann=FALSE,axes=FALSE,frame=TRUE)
  axis(1)
  title(xlab=expression(FD[2]),line=3)
  points(DB2[2,],DB2[3,],col=col,pch=pch,cex=cex,lwd=1.35*cex)
  for(i in 1:n) if(OI[i]) points(DB2[2,i],DB2[3,i],col=col[i],pch=pch[i],cex=cex[i],lwd=1.35*cex[i])
  if (!is.null(sp.index)) points(DB2[2,sp.index],DB2[3,sp.index],pch=16,col="orange")
  segments(0,0,M,M,lty=3,lwd=2)
  par(op)
}