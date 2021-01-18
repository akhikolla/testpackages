Design.matrix <- function(u, x, kn, degree)
{
  n = length(u)
  q = kn+degree+1
  u.k = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots = as.numeric(stats::quantile(u, u.k))
  # pi.u = bs(u, knots=Knots, intercept=TRUE, degree=degree,Boundary.knots=c(-0.3,10))[,1:q]
  pi.u = splines::bs(u, knots=Knots, intercept=TRUE, degree=degree)[,1:q]
  pi.u = cbind(1,pi.u[,-1])
  Pi.1 = apply(x, 2, function(x) x*pi.u)

  X1 = matrix(Pi.1[,1], nrow=n)             # intercept
  X2 = Pi.1[c(1:n), -1]                     # constant part
  X3 = matrix(Pi.1[-c(1:n), -1], nrow=n)    # varying part
  X = cbind(X1, X2, X3)

  X4 = matrix(Pi.1[, -1], nrow=n)    # varying part
  Xns = cbind(X1, X4)

  design = list(pi.u = pi.u, X = X, Xns= Xns, Pi.1 = Pi.1)
  return(design)
}


Selection.CI <- function(GS.r, L, level){
  lt = (1-level)/2; ut= 1-lt
  limits = apply(GS.r, 2, stats::quantile, probs = c(lt, ut))
  temp = matrix(abs(sign(limits[1,]) + sign(limits[2,]))==2, nrow = L)
  (apply(temp, 2, sum)>=1)*1
}


.onUnload <- function (libpath) {
  library.dynam.unload("spinBayes", libpath)
}
