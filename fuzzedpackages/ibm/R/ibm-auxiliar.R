# Local integration -------------------------------------------------------

localIntegration = function(x, y, R1, R2=R1) {
  
  if(ncol(x)!=ncol(y)) stop("Dimensions of x and y must agree.")
  if(ncol(x)>3) stop("Dimensions must be less than 3.")
  Nx = nrow(x) 
  Px = nrow(y)
  
  if(all(Nx!=0,Px!=0)) { 
    
    int    = local.int.cpp(posx=x, posy=y, R1=R1, R2=R2)
    intX   = int$x
    intY   = int$y
    
  } else {
    
    intX = rep(0, length=Nx)
    intY = rep(0, length=Px)	
    
  }
  
  return(list(x=intX, y=intY))
}


# Brownian motion ---------------------------------------------------------


brownian1D = function(object, sd, N=NULL, ...) {
  if(N==0) return(object)
  if(any(is.null(N), N>length(object), N<0)) N = length(object)
  object[seq_len(N)] 	= object[seq_len(N)] + rnorm(N, mean=0, sd=sd)
  return(object)
}

brownian2D = function(object, sd, N=NULL, ...) {
  if(N==0) return(object)
  if(any(is.null(N), N>length(object), N<0)) N	= nrow(object)
  angulo	= runif(N)
  dist	= rnorm(N, mean=0, sd=sd)
  object[seq_len(N), ] 	= object[seq_len(N), ] + dist*cbind(cos(pi*angulo),sin(pi*angulo))
  return(object)
}

brownian3D = function(object, sd, N=NULL, ...) {
  if(N==0) return(object)
  if(any(is.null(N), N>length(object), N<0)) N	= nrow(object)
  angulo1	= runif(N)
  angulo2	= runif(N)
  dist	= rnorm(N, mean=0, sd=sd)
  object[seq_len(N), , ] 	= object[seq_len(N), , ] + 
                            dist*cbind(sin(pi*angulo2)*cos(pi*angulo1), 
                                       sin(pi*angulo2)*sin(pi*angulo1), 
                                       cos(pi*angulo2))
  return(object)
}


# sampling ----------------------------------------------------------------

.checkRates = function(rates, n) {
  if(length(rates)==1) rates = rep(rates, length=n)
  if(length(rates)!=n) stop("Rates and population size do not match.") 
  rates = pmin(pmax(rates, 0),1)
  if(any(is.na(rates))) warning("NA rates are being taken as zero.")
  rates[is.na(rates)] = 0
  return(rates)
}


# Misc --------------------------------------------------------------------

.getVector = function(x, size) {
  if(is.null(size)) return(x)
  if(size>length(x)) return(x)
  return(x[seq_len(size)])
}

.getMatrix = function(x, size) {
  if(is.null(size)) return(x)
  if(size>nrow(x)) return(x)
  return(x[seq_len(size), , drop=FALSE])
}

updateMatrixSize = function(x, n, max) {
  ndim = dim(x)
  ndim[1] =  min(max(5*ndim[1], 2*n), max)
  out = array(dim=ndim)
  out[seq_len(dim(x)[1]), ] = x
  return(out)
}


# Outputs -----------------------------------------------------------------


.summary = function(x, alpha, nmax) {
  
  nmax = min(nmax, ncol(x))
  output = list()
  output$mean = rowMeans(x, na.rm=TRUE)
  output$median = apply(x, 1, median, na.rm=TRUE)
  output$ll = apply(x, 1, quantile, prob=(1-alpha)/2, na.rm=TRUE)
  output$ul = apply(x, 1, quantile, prob=1-(1-alpha)/2, na.rm=TRUE)
  output$rep = x[, sample(ncol(x), nmax)]
  output$xlim = c(0, nrow(x))
  output$ylim = quantile(x, prob=c((1-alpha)/4, 1-(1-alpha)/4), na.rm=TRUE)
  
  return(output)
  
}


