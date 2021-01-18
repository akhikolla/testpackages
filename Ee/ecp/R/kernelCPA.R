# srcGetV = '
# 	#define SUM(A) std::accumulate(A.begin(), A.end(), 0.0)
# 	using namespace Rcpp;
# 	NumericMatrix K(K_);
# 	int N = K.nrow();
# 	NumericMatrix V(N,N);
	
# 	for(int i=0;i<N;++i)
# 		for(int j=i;j<N;++j)
# 			V(i,j) = V(j,i) = SUM(diag(K(Range(i,j),Range(i,j))))-SUM(K(Range(i,j),Range(i,j)))/(j-i+1);
# 	return wrap(V);
# '
# getVCpp = inline::cxxfunction(signature(K_='matrix'),srcGetV,plugin='Rcpp')

getVCpp = srcGetV

getV = function(X){
	DD = as.matrix(dist(X))**2
	h = getBandwidth(X)
	K = exp(-DD/h)
	V = getVCpp(K)
	return(V)
}

# srcGetBandwidth = '
# 	using namespace Rcpp;
# 	NumericMatrix X(X_);
# 	NumericVector rows(rows_);
# 	int N = rows.size();
# 	NumericVector u(N*N,(double)0);
# 	for(int i=0;i<N;++i){
# 		for(int j=0;j<N;++j){
# 			u[i*N+j] = sum((X(rows[i],_)-X(rows[j],_))*(X(rows[i],_)-X(rows[j],_)));
# 		}
# 	}
# 	return wrap(u);
# '
# getBandwidthCpp = inline::cxxfunction(signature(X_='matrix',rows_='numeric'),srcGetBandwidth,plugin='Rcpp')

getBandwidthCpp = srcGetBandwidth

getBandwidth = function(X){
	if(nrow(X) > 250)
		rows = sample(0:(nrow(X)-1),250,replace=FALSE)
	else
		rows = 0:(nrow(X)-1)
	#rows = rows - 1
	u = getBandwidthCpp(X,rows)
	return(median(u))
}

getVmax = function(X){
	lower = ceiling(.05*nrow(X))
	upper = floor(.95*nrow(X))
	a = cov(as.matrix(X[1:lower,]))
	b = cov(as.matrix(X[upper:nrow(X),]))
	a = sum(diag(a)); b = sum(diag(b))
	if(lower == 1)#in this case only use the trivial bound of vmax <= 1 bc of kernel
		a = 1
	return(max(a,b))
}

# srcKcpa = '
# 	using namespace Rcpp;
# 	NumericMatrix II(II_);
# 	NumericMatrix V(V_);
# 	IntegerMatrix H(H_);
# 	int N = V.nrow();
# 	int L = H.nrow();

# 	for(int k=1;k<L;++k){
# 		for(int i=k;i<N;++i){
# 			for(int j=k-1;j<i;++j){
# 				double tmp = II(k-1,j) + V(j+1,i);
# 				if(tmp < II(k,i)){
# 					II(k,i) = tmp;
# 					H(k,i) = j+1;//fix indexing differences between R and C++
# 				}
# 			}
# 		}
# 	}
# 	return List::create(II,H);
# '
# kcpaCpp = inline::cxxfunction(signature(II_='matrix',V_='matrix',H_='matrix'),srcKcpa,plugin='Rcpp')

kcpaCpp = srcKcpa

kcpa = function(X,L,C){
	L = L+1 #user input L=max # of cps, algo uses L=max # of segments
	v = getVmax(X)
	V = getV(X)
	N = nrow(X)
	II = matrix(Inf,L,N)
	H = matrix(0,L,N)
	II[1,] = V[1,]
	answer = kcpaCpp(II,V,H)
	
	II = answer[[1]]
	H = answer[[2]]
	S = II[,N]
	for(i in 1:L){
		S[i] = S[i] + (C*v*i/N)*(1+log(N/i))
	}
	k = which.min(S)
	cps = N
	cp = cps
	while(k>0){
		cp = H[k,cp]
		cps = c(cp,cps)
		k = k-1
	}
	cps = cps+1
	return(cps)
}
