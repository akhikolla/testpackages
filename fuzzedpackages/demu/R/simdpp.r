 #    simdpp.r: R functions implementing determinantal point process based
 #              sampling of (approximately) optimal designs for GP regression.
 #    Copyright (C) 2018  Matthew T. Pratola <mpratola@stat.osu.edu>

 #    This program is free software: you can redistribute it and/or modify
 #    it under the terms of the GNU Affero General Public License as published
 #    by the Free Software Foundation, either version 3 of the License, or
 #    (at your option) any later version.

 #    This program is distributed in the hope that it will be useful,
 #    but WITHOUT ANY WARRANTY; without even the implied warranty of
 #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #    GNU Affero General Public License for more details.

 #    You should have received a copy of the GNU Affero General Public License
 #    along with this program.  If not, see <http://www.gnu.org/licenses/>.





########################################################################################
# DPP Functions ########################################################################
########################################################################################


mcmc.modal.dpp<-function(N,l.dez,pts,sw=0.05,burn=0.1*N,corr.model="gaussian")
{
	if(corr.model!="gaussian") stop("only gaussian correlation model is currently supported\n")

	# initialize vector
	draw.rhos=rep(NA,N)
	draw.rhos[1]=0.5 

	# initialize acceptance rate
	accept.rhos=1

	for(i in 2:N)
	{
		draw.rhos[i]=draw.rhos[i-1]
		rho.new=draw.rhos[i]+runif(1,-sw,sw)
		a=-Inf
		if(rho.new>0 && rho.new<1)
			a=log.prob.modal.gaussian.dpp(rho.new,l.dez,pts)-log.prob.modal.gaussian.dpp(draw.rhos[i],l.dez,pts)

		a=min(0,a)
		if(log(runif(1))<a)
		{
			draw.rhos[i]=rho.new
			accept.rhos=accept.rhos+1
		}

		if(i%%10==0)
			cat(i/N*100," percent complete\r")
	}

	draw.rhos=draw.rhos[(burn+1):N]
	accept.rhos=accept.rhos/N

	return(list(rhos=draw.rhos,accept=accept.rhos))
}

log.prob.modal.gaussian.dpp<-function(rho,l.dez,pts)
{
	d=length(l.dez)
	rho=rep(rho,d)
	R=rhomat(l.dez,rho)$R

	return(log(prob.modal.dpp(R,pts)))
}

# prob.modal.dpp assumes probability calculated using modal eigenvectors
# This is only proportional to the probabiliity since we don't compute the normalizing constant.
prob.modal.dpp<-function(R,pts)
{
	if(Matrix::isDiagonal(R)) stop("The matrix R is diagonal!\n")
	n=length(pts)
	N=nrow(R)
	e=eigen(R)
	set=1:n # because modal

	V=e$vectors[,set]
	rm(e)

	prob.points=1
	for(i in 1:n)
	{
		prob=apply(V[pts,,drop=FALSE],1,sumsq)
		ix=which(prob==max(prob))[1]
		prob.points=prob.points*prob[ix]/(n-i+1)

		v.orthto=V[,1]
		V=V[,-1,drop=FALSE]
		V=V-v.orthto%o%V[pts[ix],]/v.orthto[pts[ix]]
		pts=pts[-ix]
		if(i<(n-1))
			V=subspace.3(V)
	}

	return(prob.points)
}

sim.dpp.modal.fast<-function(R,n)
{
	if(n<=0) stop("Invalid n!\n")
	if(n>nrow(R)) stop("Invalid n!\n")
	if(class(R)=="spam") {
		warning("Converted SPAM matrix to type dgCMatrix.\n")
		R=spam::as.dgCMatrix.spam(R)
	}
	if(class(R)!="dgCMatrix") stop("Expecting matrix of class dgCMatrix!\n")
	if(Matrix::isDiagonal(R)) stop("The matrix R is diagonal!\n")

	pts=simDppModal_(R,n)

	# because indexing in C goes from 0...n-1 we have to add 1
	pts=pts+1

	return(pts)
}

sim.dpp.modal.fast.seq<-function(curpts, R, n)
{
    if(Matrix::isDiagonal(R)) stop("The matrix R is diagonal!\n")
    if(length(curpts)>=(nrow(R)-n)) stop("Requested n too large -- available points exhausted!\n")
	curpts=as.vector(curpts)
    ix = 1:nrow(R)
    notcur = ix[-curpts]
    R.curpts = R[curpts, curpts]
    R.notcur = R[notcur, notcur]
    R.curpts.notcur = R[curpts, notcur]
    Rprime = R.notcur - Matrix::t(R.curpts.notcur) %*% Matrix::chol2inv(Matrix::chol(R.curpts)) %*% 
        R.curpts.notcur
    if(!Matrix::isSymmetric(Rprime)) 
        Rprime = (Rprime + Matrix::t(Rprime))/2
    pts.new = sim.dpp.modal.fast(Rprime,n)
    pts.new = notcur[pts.new]
    return(list(pts.new = pts.new, pts.old = curpts, pts.allin = c(curpts, 
        pts.new)))
}

remove.projections<-function(curpts,X)
{
	p=ncol(X)
	projpts=NULL

	for(i in 1:length(curpts))
	for(j in 1:p)
	{
		dvec=abs(X[curpts[i],j]-X[,j])
		newprojpts=which(dvec==0)
		projpts=c(projpts,newprojpts)
	}
	projpts=unique(projpts)
	for(i in 1:length(curpts)) {
		ix=which(projpts==curpts[i])
		if(length(ix>0)) projpts=projpts[-ix]
	}
	allpts=unique(c(curpts,projpts))

	return(list(curpts=curpts,projpts=projpts,allpts=allpts))
}

sim.dpp.modal.seq<-function(curpts,R,n)
{
	if(Matrix::isDiagonal(R)) stop("The matrix R is diagonal!\n")
    if(length(curpts)>=(nrow(R)-n)) stop("Requested n too large -- available points exhausted!\n")
	curpts=as.vector(curpts)
	ix=1:nrow(R)
	notcur=ix[-curpts]
	R.curpts=R[curpts,curpts]
	R.notcur=R[notcur,notcur]
	R.curpts.notcur=R[curpts,notcur]
	Rprime=R.notcur-Matrix::t(R.curpts.notcur)%*%chol2inv(chol(R.curpts))%*%R.curpts.notcur
	if(!isSymmetric(Rprime))
		Rprime=(Rprime+Matrix::t(Rprime))/2 # A hack, Rprime may not come back as numerically symmetric even though it must be theoretically.


	pts.new=sim.dpp.modal(Rprime,n)
	pts.new=notcur[pts.new]

	return(list(pts.new=pts.new,pts.old=curpts,pts.allin=c(curpts,pts.new)))
}

# Draw modal DPP sample using a grid-based Nystrom approximation.
# This has limitations as the number of grid points would need to grow like ngrid^p so won't scale
# to high dimensions.
sim.dpp.modal.nystrom<-function(Xin,rho,n=0,ngrid=NULL,method="Nystrom")
{
	# NULL means we aren't do much of an approXinimation...
	if(is.null(ngrid)) stop("Missing number of grid points to reconstruct over entire space!\n")

	p=ncol(Xin)
	l.d=makedistlist(Xin)

	xgrid=seq(0,1,length=ngrid)
	ll=vector("list",p)
	for(i in 1:p) ll[[i]]=xgrid
	Xapp=as.matrix(expand.grid(ll))
	Xapp=as.matrix(rbind(Xapp,Xin))
	l.all=makedistlist.subset(Xapp,1:(nrow(Xapp)),(ngrid^p+1):(ngrid^p+nrow(Xin)))

	R=rhomat(l.d,rho)$R
	e=eigen(R)
	Rall=rhomat(l.all,rho)$R
	eall=e
	eall$vectors=Rall%*%e$vectors
	for(i in 1:ncol(eall$vectors)) eall$vectors[,i]=eall$vectors[,i]/eall$values[i]

	pts=sim.dpp.modal(NULL,n=n,eigs=eall)

	return(list(pts=pts,X=Xapp,design=Xapp[pts,,drop=FALSE]))
}

# Use kmeans-based Nystrom approximation of Li et al (2010) to draw an (approximate) modal sample
# from an observational dataset.
# Xin is the dataset, n is the number of design points from Xall that we want.
# rho is the parameter vector under the Gaussian correlation model.
# m is the number of landmark points for the kmeans-based Nystrom approximation.
# initializer is how MiniBatchKmeans initializes, recommend "kmeans++" or "random".  The
# default in package ClusterR is actually broken and marked experimental, so not sure why it's the default.
# Additional arguments are passed to MiniBatchKmeans - see the docs in ClusterR for other options.
# Because MiniBatchKmeans does not return indicies (ugh), we return the candidate matrix and the index
# within the candidate matrix, as well as the landmark points.
sim.dpp.modal.nystrom.kmeans<-function(Xin,rho=rep(0.01,ncol(Xin)),n,m=max(ceiling(nrow(Xin)*0.1),n),method="KmeansNystrom",initializer="kmeans++",...)
{
	# m<n is going to be a problem...
	if(n>m) stop("Require number of landmark points to be greater than requested sample size!\n")

	X=as.matrix(MiniBatchKmeans(Xin,m,initializer=initializer,...)$centroids)

	l.d=makedistlist(X)
	Xapp=as.matrix(rbind(Xin,X))
	l.all=makedistlist.subset(Xapp,1:(m+nrow(Xin)),(nrow(Xin)+1):(nrow(Xin)+m))

	R=rhomat(l.d,rho)$R
	e=eigen(R)
	Rall=rhomat(l.all,rho)$R
	eall=e
	eall$vectors=Rall%*%e$vectors
	for(i in 1:ncol(eall$vectors)) eall$vectors[,i]=eall$vectors[,i]/eall$values[i]

	pts=sim.dpp.modal(NULL,n=n,eigs=eall)

	return(list(pts=pts,X=Xapp,X.lm=X,design=Xapp[pts,,drop=FALSE]))
}

# Generate modal sample of size n in p dimensions with parameter vector rho assuming Gaussian kernel.
# m is the number of landmark points to use in forming the kmeans nystrom approximation and N is the 
# overall approximation size (using uniform sampling).
# Uses sim.dpp.modal.nystrom.kmeans() to draw the design.
# Returns the candidates and the index into the candidates as well as the landmark points.
sim.dpp.modal.np<-function(n,p,N,rho,m=max(ceiling(N*0.1),n),...)
{
	Xall=matrix(runif(N*p),ncol=p)

	return(sim.dpp.modal.nystrom.kmeans(Xall,rho,n,m,...))
}

# This surprisingly does not work well, seems the matrix approximation is very rank deficient.
# So much for combining the nystrom approximation with the kronecker product trick :(
sim.dpp.modal.nykron<-function(X,rho,n=0,ngrid=NULL,method="NystromKronecker")
{
	# NULL means we aren't do much of an approximation...
	if(is.null(ngrid)) stop("Missing number or grid points per dimension!\n")

	p=ncol(X)
	l.d=makedistlist(X)

	l.temp=list(l1=list(m1=1))
	Rmats=vector("list",p)
	for(i in 1:p) {
		l.temp$l1$m1=l.d[[i]][[1]]
		Rmats[[i]]=rhomat(l.temp,rho[i])$R
	}
	xgrid=seq(0,1,length=ngrid)
	Xapp=as.matrix(rbind(matrix(xgrid,nrow=ngrid,ncol=p,byrow=FALSE),X))
	l.all=makedistlist(Xapp)

	# Construct Nystrom approximations
	eveclist=vector("list",p)
	evallist=vector("list",p)
	for(i in 1:p) {
		e=eigen(Rmats[[i]])
		l.temp$l1$m1=l.all[[i]][[1]]
		Rtemp=rhomat(l.temp,rho[i])$R
		eveclist[[i]]=Rtemp[,(ngrid+1):(ngrid+nrow(X))]%*%e$vectors
		for(j in 1:ncol(eveclist[[i]])) eveclist[[i]][,j]=eveclist[[i]][,j]/e$values[j]
		evallist[[i]]=e$values
	}

	# Construct overall evecs,evals using Kronecker rule
	eall=list("values","vectors")
	eall$vectors=1
	eall$values=1
	for(i in 1:p)
	{
		eall$vectors=eall$vectors%x%eveclist[[i]]
		eall$values=eall$values%x%evallist[[i]]
	}
	ix=sort(eall$values,decreasing=TRUE,index.return=TRUE)$ix
	eall$values=eall$values[ix]
	eall$vectors=eall$vectors[,ix]

	return(sim.dpp.modal(NULL,n=n,eigs=eall))
}

sim.dpp.modal<-function(R,n=0,eigs=NULL)
{
	if(is.null(R) && is.null(eigs)) stop("Missing R and/or eigs!\n")
	if(!is.null(R) && Matrix::isDiagonal(R)) stop("The matrix R is diagonal!\n")
	if(!is.null(R)) {
		eigs=eigen(R)
	}
	p=eigs$values/(1+eigs$values)

	N=length(p)
	set=NULL

	# if n>0 we are drawing a conditional realization of n points
	# otherwise n is randomly selected for drawing the points
	if(n>0)
	{
		set=1:n
	}
	else {
		stop("n>0 required\n")
	}

	pts=NULL
	V=eigs$vectors[,set,drop=FALSE]

	for(i in 1:n)
	{
		pvec=rep(0,N)
		pvec=apply(V,1,sumsq)
		pvec=pvec/(n-i+1)

		p.dx=which(pvec==max(pvec))[1]
		pts=c(pts,p.dx)

		v.orthto=V[,1]
		V=V[,-1,drop=FALSE]
		V=V-v.orthto%o%V[p.dx,]/v.orthto[p.dx]
		if(i<(n-1))
			V=subspace.3(V)
	}

	return(pts)
}

sim.dpp.modal.2<-function(R,n=0)
{
	if(Matrix::isDiagonal(R)) stop("The matrix R is diagonal!\n")
	e=eigen(R)
	p=e$values/(1+e$values)

	N=length(p)
	set=NULL

	# if n>0 we are drawing a conditional realization of n points
	# otherwise n is randomly selected for drawing the points
	if(n>0)
	{
		set=1:n
	}
	else {
		stop("n>0 required\n")
	}

	pts=NULL
	V=e$vectors[,set,drop=FALSE]

	for(i in 1:n)
	{
		pvec=rep(0,N)
		pvec=apply(V,1,sumsq)
		pvec=pvec/(n-i+1)

		p.dx=which(pvec==max(pvec))[1]
		pts=c(pts,p.dx)

		v.orthto=V[,1]
		V=V[,-1,drop=FALSE]
		V=V-v.orthto%o%V[p.dx,]/v.orthto[p.dx]
		if(i<(n-1))
			V=subspace.3(V)
	}

	return(pts)
}

proj<-function(vec,onto)
{
	vec_onto=t(onto)%*%vec/(t(onto)%*%onto)
	vec_onto=vec_onto*onto
	return(vec_onto)
}

sumsq<-function(vec)
{
	return(sum(vec^2))
}

sim.dpp<-function(R,n=0,method="default")
{
	if(Matrix::isDiagonal(R)) stop("The matrix R is diagonal!\n")
	e=eigen(R)
	p=e$values/(1+e$values)

	N=length(p)
	set=NULL

	# if n>0 we are drawing a conditional realization of n points
	# otherwise n is randomly selected for drawing the points
	if(n>0)
	{
		if(method!="modal") {
			set=select.evals(n,N,e$values)
		}
		else {
			set=1:n
		}
	}
	else {
		for(i in 1:N) {
			if(runif(1)<p[i])
				set=c(set,i)
		}
		n=length(set)
	}

	pts=NULL
	V=e$vectors[,set,drop=FALSE]

	for(i in 1:n)
	{
		pvec=rep(0,N)
		pvec=apply(V,1,sumsq)
		pvec=pvec/(n-i+1)

		p.dx=select(pvec)
		pts=c(pts,p.dx)

		v.orthto=V[,1]
		V=V[,-1,drop=FALSE]
		V=V-v.orthto%o%V[p.dx,]/v.orthto[p.dx]
		if(i<(n-1))
			V=subspace.3(V)
	}

	return(pts)
}

subspace<-function(V,v.orthto)
{
	Vnew=matrix(0,nrow=nrow(V),ncol=ncol(V))
	Vnew[,1]=v.orthto

	# Orthogonal basis for Vnew orthogonal to v.orthto
	Vnew[,1]=v.orthto
	for(i in 2:ncol(Vnew)) {
		Vnew[,i]=V[,i]
		for(j in 1:(i-1))
			Vnew[,i]=Vnew[,i]-(t(Vnew[,j])%*%V[,i]/(t(Vnew[,j])%*%Vnew[,j]))*Vnew[,j]
	}

	# Normalize
	for(i in 1:ncol(Vnew)) Vnew[,i]=Vnew[,i]/sqrt(t(Vnew[,i])%*%Vnew[,i])

	# Chop off the first one
	return(Vnew[,2:ncol(Vnew),drop=FALSE])
 }

# Stabilized version of Gram-Schmidt
subspace.2<-function(V,v.orthto)
{
	Vnew=matrix(0,nrow=nrow(V),ncol=ncol(V))
	Vnew[,1]=v.orthto

	# Orthogonal basis for Vnew orthogonal to v.orthto
	for(i in 2:ncol(Vnew)) {
		Vnew[,i]=V[,i]
		for(j in 1:(i-1))
			Vnew[,i]=Vnew[,i]-(t(Vnew[,i])%*%Vnew[,j]/(t(Vnew[,j])%*%Vnew[,j]))*Vnew[,j]
	}

	# Normalize
	for(i in 1:ncol(Vnew)) Vnew[,i]=Vnew[,i]/sqrt(t(Vnew[,i])%*%Vnew[,i])

	# Chop off the first one
	return(Vnew[,2:ncol(Vnew),drop=FALSE])
 }

# Stabalized version of Gram-Schmidt
subspace.3<-function(V)
{
	return(subspace_(V))
}

# Below R code now superceeded by the C++ call above.
# subspace.3<-function(V)
# {
# 	Vnew=matrix(0,nrow=nrow(V),ncol=ncol(V))
# 	Vnew[,1]=V[,1]

# 	# Orthogonal basis for Vnew orthogonal to v.orthto
# 	for(i in 2:ncol(Vnew)) {
# 		Vnew[,i]=V[,i]
# 		for(j in 1:(i-1))
# 			Vnew[,i]=Vnew[,i]-(t(Vnew[,i])%*%Vnew[,j]/(t(Vnew[,j])%*%Vnew[,j]))*Vnew[,j]
# 	}

# 	# Normalize
# 	for(i in 1:ncol(Vnew)) Vnew[,i]=Vnew[,i]/sqrt(t(Vnew[,i])%*%Vnew[,i])

# 	return(Vnew)
#  }

select<-function(vec)
{
	if(!is.vector(vec)) stop("Argument is not a vector\n")
	i=1
	p=runif(1)
	psum=vec[1]
	while(p>psum)
	{
		i=i+1
		psum=psum+vec[i]
	}
	return(i)
}

# See Chen & Liu, Statistica Sinica, vol 7, pp875-892 (1997).
select.evals<-function(n,N,lambda)
{
	set=NULL
	C=1:N
	r=0

	for(k in 1:N)
	{
		pk=lambda[k]*RcB(n-r-1,C[(k+1):N],lambda)/RcB(n-r,C[k:N],lambda)
cat("pk=",pk,"\n")
		if(runif(1)<pk) {
			set=c(set,k)
		}
		r=length(set)

		#early exit of loop if we're done
		if(r==n) break			 
	}

	return(set)
}

# "R" function for the Conditional Bernoulli distribution
# See Chen & Liu, Statistica Sinica, vol 7, pp875-892 (1997).
RcB<-function(k,C,w)
{
	normC=length(C)

	if(k==0) { return(1) }
	if(k>normC) { return(0) }

	R.kC=0
	for(i in 1:k) {
		R.kmi=RcB(k-i,C,w)
		T.ic=TiC(i,C,w)
		R.kC=R.kC+(-1)^(i+1)*T.ic*R.kmi
	}
	R.kC=R.kC/k

	return(R.kC)
}

TiC<-function(i,C,w)
{
	return(sum(w[C]^i))
}



########################################################################################
# Correlation Functions ################################################################
########################################################################################

# rhomat:
# Calculate the correlation matrix for the power exponential/Gaussian model.
#
# The correlation parameters are rho_1,...,rho_p for a p-dim space, each in [0,1] and
# we will have p distance (l.d) matrices.  We construct these in a list of lists, see
# for example the following:
# l1=list(m1=design.distmat.dim1)
# l2=list(m2=design.distmat.dim2)
# l=list(l1=l1,l2=l2)
#
rhomat<-function(l.d,rho,alpha=2)
{
	rho=as.vector(rho)
	if(!is.vector(rho)) stop("non-vector rho!")
	if(any(rho<0)) stop("rho<0!")
	if(any(rho>1)) stop("rho>1!")
	if(any(alpha<1) || any(alpha>2)) stop("alpha out of bounds!")
	if(!is.list(l.d)) stop("wrong format for distance matrices list!")
	if(length(l.d)!=length(rho)) stop("rho vector doesn't match distance list")

	R=matrix(1,nrow=nrow(l.d$l1$m1),ncol=ncol(l.d$l1$m1))
	for(i in 1:length(rho))
		R=R*rho[i]^(abs(l.d[[i]][[1]])^alpha)

	return(list(R=R))
}

# Matern nu=3/2 model
matern32<-function(l.d,theta)
{
	theta=as.vector(theta)
	if(!is.vector(theta)) stop("non-vector theta!")
	if(any(theta<0)) stop("theta<0!")
	if(!is.list(l.d)) stop("wrong format for distance matrices list!")
	if(length(l.d)!=length(theta)) stop("theta vector doesn't match distance list")

	R=matrix(1,nrow=nrow(l.d$l1$m1),ncol=nrow(l.d$l1$m1))
	for(i in 1:length(theta))
	{
		D=(1+sqrt(3)*l.d[[i]][[1]]/theta[i])*exp(-sqrt(3)*l.d[[i]][[1]]/theta[i])
		R=R*D
	}

	return(list(R=R))
}

# Matern nu=5/2 model
matern52<-function(l.d,theta)
{
	theta=as.vector(theta)
	if(!is.vector(theta)) stop("non-vector theta!")
	if(any(theta<0)) stop("theta<0!")
	if(!is.list(l.d)) stop("wrong format for distance matrices list!")
	if(length(l.d)!=length(theta)) stop("theta vector doesn't match distance list")

	R=matrix(1,nrow=nrow(l.d$l1$m1),ncol=nrow(l.d$l1$m1))
	for(i in 1:length(theta))
	{
		D=(1+sqrt(3)*l.d[[i]][[1]]/theta[i]+5*l.d[[i]][[1]]^2/(5*theta[i]^2))*exp(-sqrt(5)*l.d[[i]][[1]]/theta[i])
		R=R*D
	}

	return(list(R=R))
}

# Wendland 1 (see Furrer et al)
wendland1<-function(l.d,theta)
{
	theta=as.vector(theta)
	if(!is.vector(theta)) stop("non-vector theta!")
	if(any(theta<0)) stop("theta<0!")
	if(!is.list(l.d)) stop("wrong format for distance matrices list!")
	if(length(l.d)!=length(theta)) stop("theta vector doesn't match distance list")

	R=matrix(1,nrow=nrow(l.d$l1$m1),ncol=nrow(l.d$l1$m1))
	for(i in 1:length(theta))
	{
		D=(1-l.d[[i]][[1]]/theta[i])
		D[D<0]=0
		D=D^4*(1+4*l.d[[i]][[1]]/theta[i])
		R=R*D
	}

	return(list(R=R))	
}

# Wendland 2 (see Furrer et al)
wendland2<-function(l.d,theta)
{
	theta=as.vector(theta)
	if(!is.vector(theta)) stop("non-vector theta!")
	if(any(theta<0)) stop("theta<0!")
	if(!is.list(l.d)) stop("wrong format for distance matrices list!")
	if(length(l.d)!=length(theta)) stop("theta vector doesn't match distance list")

	R=matrix(1,nrow=nrow(l.d$l1$m1),ncol=nrow(l.d$l1$m1))
	for(i in 1:length(theta))
	{
		D=(1-l.d[[i]][[1]]/theta[i])
		D[D<0]=0
		D=D^6*(1+6*l.d[[i]][[1]]/theta[i]+35*l.d[[i]][[1]]^2/(3*theta[i]^2))
		R=R*D
	}

	return(list(R=R))	
}

# Requires package fields.
# kap is degree of differentiability desired.
generalized.wendland<-function(l.d,theta,kap)
{
	d=length(l.d)
	mu=(d+1)/2  # strictly speaking we need mu>=(d+1)/2 but we can change kappa also so
				# this is fair to assume.

	if(length(theta)>1) stop("theta is incorrect dimension\n")
	if(length(kap)>1) stop("kappa is incorrect dimensions\n")
	if(mu<((d+1)/2) ) stop("mu does not satisfy constraints\n")
	if(kap<0) stop("kappa > 0 required\n")

	D=matrix(0,nrow=nrow(l.d$l1$m1),ncol=nrow(l.d$l1$m1))
	for(i in 1:length(l.d))
		D=D+(l.d[[i]][[1]]^2)
	D=sqrt(D)

	if(kap==0) {
		# kap=0 is essential the Askey correlation
		D=D/theta
		R=1-D
		R[R<0]=0
		R=R^(mu+kap)
	}
	else
	{
		# library(fields) implements the general case
		R=fields::Wendland(D,theta=theta,dimension=d,k=kap)
	}

	rm(D)
	return(list(R=R))
}

# askey.sparse<-function(X,theta)
# {
# 	p=ncol(X)
# 	mu=(p+1)/2  # strictly speaking we need mu>=(p+1)/2 but we can change kappa also so
# 				# this is fair to assume.

# 	if(length(theta)>1) stop("theta is incorrect dimension\n")
# 	if(mu<((p+1)/2) ) stop("mu does not satisfy constraints\n")

# 	theta2=theta^2
# 	n=nrow(X)
# 	R=Matrix(0,n,n)

# 	for(i in 1:n){
# 		for(j in 1:i) {
# 			d=sum((X[i,]-X[j,])^2)
# 			if(d<=theta2) {
# 				R[i,j]=(1-sqrt(d)/theta)^mu
# 				R[j,i]=R[i,j]
# 			}
# 		}
# 	}

# 	return(R)
# }


# getranges
# Return ranges of design inputs
getranges<-function(design)
{
	r=t(apply(design,2,range))
	r[,1]=floor(r[,1])
	r[,2]=ceiling(r[,2])

	return(r)
}

# unscale
unscalemat<-function(mat,r)
{
	p=ncol(mat)
	for(i in 1:p)
		mat[,i]=(mat[,i]*(r[i,2]-r[i,1]))+r[i,1]

	return(mat)
}

# scaledesign
# Rescale the design to the [0,1] hypercube.
scaledesign<-function(design,r)
{
	p=ncol(design)
	for(i in 1:p)
		design[,i]=(design[,i]-r[i,1])/(r[i,2]-r[i,1])

	return(design)
}


# makedistlist
# Make list of distance matrices
makedistlist<-function(design)
{
	design=as.matrix(design)
	p=ncol(design)
	l.d=vector("list",p)
	for(i in 1:p) {
		l.d[[i]]=list(abs(outer(design[,i],design[,i],"-")))
		names(l.d[[i]])=paste("m",i,sep="")
		names(l.d)[[i]]=paste("l",i,sep="")
	}
		
	return(l.d)
}

# for a subset
makedistlist.subset<-function(design,sub1,sub2)
{
	design=as.matrix(design)
	p=ncol(design)
	l.d=vector("list",p)
	for(i in 1:p) {
		l.d[[i]]=list(abs(outer(design[sub1,i],design[sub2,i],"-")))
		names(l.d[[i]])=paste("m",i,sep="")
		names(l.d)[[i]]=paste("l",i,sep="")
	}
		
	return(l.d)	
}
