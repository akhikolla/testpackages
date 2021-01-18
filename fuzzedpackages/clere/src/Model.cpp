#include "Model.h"
#include "Rmath.h"
double static runif(){
  double r = R::runif(0,1); //((double) rand())/((double) RAND_MAX);
  return r; 
};


double static rexp(double lambda){
  double u = runif();
  return -log(u)/lambda;
};
double static rnorm(double mean, double var){
  double u1 = runif();
  double u2 = runif();
  double w  = sqrt( -2.0*log(u1) ) * cos(2.0 * M_PI * u2);
  return(w*sqrt(var) + mean);
};

void static permute(VectorXi &perms){
  double u;
  int p = perms.size();
  int j1,j2;
  for(j1=0;j1<p-1;j1++){
    u = runif();
    if(u>0.5){
      j2 = perms(j1);
      perms(j1)   = perms(j1+1);
      perms(j1+1) = j2;
    }
  }
};
  
int static rmultinom(VectorXd &P, int g){
  double thresh = runif();
  int k  = 0;
  double cdf = 0;
  while( k < g ){
    cdf += P(k);
    if(thresh<=cdf){
      return(k);
    }
    k++;
  }
  return g-1;
};

double rtnorm(double mean, double var, double status){
  double lower = 0.0;
  double upper = 0.0;
  if(status==1.0){
    upper = +INFINITY; 
  }else{
    lower = -INFINITY;
  }
  double sd = sqrt(var);
  double ret,a,u;
  bool cond;
  
  lower = (lower - mean)/sd;
  upper = (upper - mean)/sd;

  ret = 0.0;
  int i = 0;
  if(lower<=0.0 and upper>=0.0){
    cond = false;
    while(!cond){
      ret = rnorm(0.0,1.0);
      cond = (ret>=lower and ret<=upper);
    }
  }
  if(lower>=0.0){
    cond = false;
    while(!cond){
      a = 0.5*(lower + sqrt(lower*lower + 4.0));
      ret = rexp(a) + lower;
      u    = runif();
      cond = (u<=exp(-0.5*(ret-a)*(ret-a)));
      i++;
      if(i>1000){
	//break;
      }   
    }
  }
  if(upper<=0.0){
    cond = false;
    while(!cond){
      a = 0.5*(-upper + sqrt(upper*upper + 4.0));
      ret = rexp(a) - upper;
      u = runif();
      cond = (u<=exp(-0.5*(ret-a)*(ret-a)));
      ret = -ret;
      i++;
      if(i>1000){
	//break;
      }
    }    
  }
 
  double result = ret*sd + mean;
  return result;
}
double static Mstep(MatrixXd &xz,VectorXd &y, VectorXd &v, VectorXd &e,
		  VectorXd &t,Vector2d &var,
		  MatrixXd &D,MatrixXd &D_,VectorXd &h,
		  int g,int n,bool sparse,string family,
		  double lambda){
  VectorXd R_ = VectorXd::Zero(n);
  int i,k;
  for(i=0;i<n;i++){
   R_(i) = 1.0/(1.0 + v(i)/lambda); 
  }  
  D  = xz.transpose() * ( R_.asDiagonal() * xz );
  D_ = D.inverse();
  h  = xz.transpose() * ( R_.asDiagonal() * y  );
  t  = D_*h;
  double ls2;
  ls2 = 0.0;
  if(sparse){
    for(k=0;k<g;k++){
      ls2 += D_(1,k) * h(k);
    }
    ls2 = -ls2/D_(1,1);
  }
  for(k=0;k<g;k++){
    t(k) += ls2*D_(k,1);
  }
  e  = y - xz*t;
  var(0) = e.squaredNorm()/n;
  var(1) = var(0) / lambda;
  // Calculate complete log-likelihood
  double ll = -0.5*n*log(2*M_PI*var(0))-0.5*n;
  for(i=0;i<n;i++){
   ll += 0.5*log(R_(i));
  }
  return ll;
}

double static EM (MatrixXd &xz,VectorXd &y, VectorXd &v, VectorXd &e,
		  VectorXd &t,Vector2d &var,
		  MatrixXd &D,MatrixXd &D_,VectorXd &h,
		  int g,int n,int maxit, double tol, bool sparse,string family){
  int i,k,it;
  D  = xz.transpose() * xz;
  D_ = D.inverse();
  h  = xz.transpose() * y;
  t  = D_*h;
  
  double ls2;
  ls2 = 0.0;
  if(sparse){
    for(k=0;k<g;k++){
      ls2 += D_(1,k) * h(k);
    }
    ls2 = -ls2/D_(1,1);
  }
  for(k=0;k<g;k++){
    t(k) += ls2*D_(k,1);
  }
  if(sparse){
    t(1) = 0.0;
  }
  if(family=="binomial"){
    var(0) = 1.0;
    t(0)   = 0.0;
  }  
  double eps = 1.0;
  double s0,s1;
  double R_i,tmp_a,tmp_b,tmp_c,tmp_d,tmp_e,tmp_f;
  double ll,lold;
  
  it = 0;
  ll = 0.0;
  lold = 0.0;
  
  // For binomial only
  //   VectorXd L1   = VectorXd::Zero(g);
  //   VectorXd L2   = VectorXd::Zero(g);
  //   VectorXd D_L1 = VectorXd::Zero(g);
  //   VectorXd D_L2 = VectorXd::Zero(g);
  
  //   L1(0) = 1.0; // For intercept
  //   L2(1) = 1.0; // For sparse
  //   Matrix2d A;
  //   Vector2d a;
  //   Vector2d lambda;
  //   double aa;
  
  while(it<maxit and eps>tol){
    e  = y-xz*t;
    tmp_c = 0.0;
    tmp_d = 0.0;
    tmp_e = 0.0;
    tmp_f = 0.0;
   
    ll = -n * log(2.0 * M_PI);
    h = D*t;
    for(i=0;i<n;i++){
      R_i   = 1.0/(var(0) + var(1)*v(i));
      tmp_a = e(i)*R_i;
      tmp_b = tmp_a * tmp_a;
      tmp_c += tmp_b;
      tmp_d += v(i) * tmp_b;
      tmp_e += R_i;
      tmp_f += R_i * v(i);
      ll += log( R_i ) - e(i)*tmp_a;
      for(k=0;k<g;k++){
	h(k) += var(0) * xz(i,k) * tmp_a;  
      }
    }   
    ll = 0.5*ll;
    s0 = var(0)*var(0)*(tmp_c-tmp_e) + var(0)*n;
    s1 = var(1)*var(1)*(tmp_d-tmp_f) + var(1)*n;   

    // M step
    if(family=="binomial"){
      var(0) = 1.0;
    }else{
      var(0) = s0/n;
    }
    var(1) = s1/n;
    t      = D_ * h;
    // sparse constraint only
    if(sparse){// and family!="binomial"){
      ls2    = 0.0;
      for(k=0;k<g;k++){
	ls2 += D_(1,k) * h(k);
      }
      ls2 = -ls2/D_(1,1);
      for(k=0;k<g;k++){
	t(k) += ls2*D_(k,1);      
      }
      t(1) = 0.0;
    }
    // binomial constraint of intercept = 0
//     if(false){//family=="binomial" and !sparse){
//       ls2    = 0.0;
//       for(k=0;k<g;k++){
// 	ls2 += D_(0,k) * h(k);
//       }
//       ls2 = -ls2/D_(0,0);
//       for(k=0;k<g;k++){
// 	t(k) += ls2*D_(k,0);      
//       }
//        t(0) = 0.0;
//     }
//     // Both binomial and sparse constraint
//     if(false){//family=="binomial" and sparse){
//       // -- 
//       D_L1 = D_ * L1;
//       D_L2 = D_ * L2;
//       aa = L1.dot(D_ * L2);
//       A(0,0) = L1.dot(D_L1); A(0,1) = aa;  
//       A(1,1) = L2.dot(D_L2); A(1,0) = aa;
//       // --
//       a(0)   = -1.0 * L1.dot(t);
//       a(1)   = -1.0 * L2.dot(t);
//       // --
//       lambda = A.inverse() * a;
//       t += lambda(0)*D_L1;
//       t += lambda(1)*D_L2;
//       t(0) = 0.0;
//       t(1) = 0.0;
//     }
    

    if(it>0){
      eps  = ll-lold;
    }
    lold = ll;
    it++;
  }
  return ll;
};

void static my_bsort(VectorXd &sortedBeta, int p){
  int i,j;
  double tmp;
  for (i = 0; i < (p - 1); ++i){
    for (j = 0; j < p - 1 - i; ++j ){
      if(sortedBeta(j) > sortedBeta(j+1)){
	tmp = sortedBeta(j+1);
	sortedBeta(j+1) = sortedBeta(j);
	sortedBeta(j)   = tmp;
      }
    }
  }
};

Model::Model(int mp,int mg,int mns){
  this->p   = mp;
  this->g   = mg;
  this->pi.resize(this->g);
  this->b.resize(this->g);
  this->B.resize(this->p);
  this->Z.resize(this->p);
  this->nsample = mns;
  this->Zw.resize(this->p,nsample);
  this->Bw.resize(this->p,nsample);
  this->Lmc.resize(nsample);
  this->P.resize(this->p,this->g);
  this->likelihood = -INFINITY;
  this->entropy    = 0.0;
  this->initialized = false;
};

void Model::init_basic(bool sparse){
  int j,k;
  gamma2    = rnorm(0.0,1.0);
  gamma2    = gamma2 * gamma2;
  
  for(k=0;k<g;k++){
    b(k)  = 2 * k * sqrt(gamma2);
    pi(k) = 1.0/g;
  }
  if(sparse){
    b(0) = 0.0;
  }
  // sigma and gamma square
  sigma2    = rnorm(0.0,1.0);
  sigma2    = sigma2 * sigma2;
  
  intercept = rnorm(0.0,1.0);
  for(j=0;j<p;j++){
    Z(j) = rmultinom(pi,g);
    B(j) = rnorm( b(Z(j)), 1.0);
  }
  initialized = true;
};

void Model::init_kmeans(bool sparse){
  int j,k;
  int og,ng,nchange;
  int j1,j2;
  VectorXd sortedBeta = VectorXd::Zero(p);
  double lmin = INFINITY;
  int kmin = 0;
  // k-means  
  // Initialization
  for(k=0;k<g;k++){
    pi(k) = 1.0/g;
  }
  for(j=0;j<p;j++){
    sortedBeta(j) = B(j);    
  }
  my_bsort(sortedBeta,p);

  for(k=0;k<g;k++){
    j1 = k * (p/g);
    j2 = min(p-1, (k + 1) * (p/g) ); // changed p into p-1 => should pass overrun test!
    b(k) = (sortedBeta(j1) + sortedBeta(j2))/2;
    if( abs( b(k) ) < lmin ){
      kmin = k;
      lmin = abs ( b(k) );
    }
  }
  lmin    = b(0);
  b(0)    = b(kmin);
  b(kmin) = lmin;
  
  // Sparse !
  if(sparse){
    b(0) = 0.0;
  }
  
  for(j=0;j<p;j++){
    ng = 0;
    for(k=0;k<g;k++){
      if( abs(B(j)-b(k)) <= abs(B(j)-b(ng)) ){
	ng = k;
      }
    }
    Z(j)   = ng;
  }
  // Calculate the b's
  // Modified date = 04/12/2013
  // Count
  VectorXi ns = VectorXi::Zero(g);
  nchange = 1;
  while(nchange!=0){
    for(k=0;k<g;k++){
      b(k)  = 0.0;
      ns(k) =   0;
    }
    for(j=0;j<p;j++){
      k = Z(j);
      b(k) += B(j);
      ns(k)++;
    }
    for(k=0;k<g;k++){
      if(ns(k)>0){
	b(k) = b(k) / ns(k);
      }
    }
    if(sparse){
      b(0) = 0.0;
    }
    nchange = 0;
    for(j=0;j<p;j++){
      og = Z(j);
      ng = og;
      if(pi(og)>1){
	for(k=0;k<g;k++){
	  if( abs(B(j)-b(k)) < abs(B(j)-b(ng)) ){
	    ng = k;
	    nchange++;
	  }
	}
	Z(j) = ng;
      }
    }
  }
  initialized = true;
};

void Model::updateZ_GibbsRows(IO *io, MatrixXd &xz, VectorXd &e,
			      VectorXi &ns, VectorXd &pdfRow,
			      VectorXi &perms, int nChanges){
  int n = io->n;
  int i,j,k;
  int oj,og,ng;
  double bk_old,bk_new,exj,nxj,norm;
  double lpmax = 0.0;
  og = 0;
  oj = 0;
  //int kmin = 0;
  permute(perms);

  for(k=0;k<g;k++) { 
    ns(k) = 0;
  }
  for(j=0;j<p;j++){
    k = Z(j);
    ns(k)++; 
  }
  // Gibbs Sampling
  for(oj=0;oj<nChanges;oj++){
    j = perms(oj);
    exj = 0.0;
    nxj = 0.0;
    for(i=0;i<n;i++){
      exj += (e(i) * io->x(i,j))/(sigma2 + gamma2*io->v(i));
      nxj += (io->x(i,j) * io->x(i,j))/(sigma2 + gamma2*io->v(i));
    }
    og = Z(j);
    bk_old = b(og);
    lpmax  = -INFINITY;
    for (k = 0; k < g ; k++){
      bk_new = b(k);
      pdfRow(k) = -exj*(bk_old-bk_new) - 0.5*nxj*(bk_old-bk_new)*(bk_old-bk_new)+log(pi(k))-log(pi(og));
      if( pdfRow(k)>lpmax ){
	      lpmax = pdfRow(k);
      }
    }
    norm  = 0.0;
    for(k=0;k<g;k++){
      pdfRow(k) = exp(pdfRow(k)-lpmax);
      norm     += pdfRow(k);
    }
    for(k=0;k<g;k++){
      pdfRow(k) = pdfRow(k)/norm;      
    }
    ng = rmultinom(pdfRow,g);
    // Modify stuffs
    if(ng != og){
      bk_new = b(ng);
      for(i=0;i<n;i++){
	e(i)      += (bk_old-bk_new) * io->x(i,j);
	xz(i,1+og)+= -io->x(i,j);
	xz(i,1+ng)+= +io->x(i,j); 	
      }
      Z(j) = ng;
      ns(ng)++;
      ns(og)--;
    }
  }
};

void Model::updateY_Gibbs(VectorXd &Ytrue, VectorXd &C, VectorXd &a, MatrixXd &Hm,int n,VectorXi &Perms){
  int oi,i,ip;
  double si,mi,vi;
  permute(Perms);   
  for(oi=0;oi<n;oi++){
    i  = Perms(oi);
    si = 0.0;
    vi = 1.0/Hm(i,i);
    for(ip=0;ip<n;ip++){
      if(ip!=i){
	si += Hm(i,ip) * ( Ytrue(ip) - a(ip) );
      }
    }
    mi = a(i) - si*vi;
    Ytrue(i) = rtnorm( mi, vi, C(i) );
  }  
};

// Modofied the 09/09/2013
void Model::fitSEM(IO *io,MatrixXd &Theta){
  int n = io->n;
  int i,j,k;
  int item,itmc;
  double RSS;
  int nItEM = io->nItEM;
  int nItMC = io->nItMC;
  int nBurn = io->nBurn;
  int gp1   = g+1;
    
  MatrixXd Hm;
  VectorXd Ytrue;
  VectorXd Ytrans;
  VectorXd a;
  VectorXi Perms;
  if(io->family=="binomial"){
    Hm.resize(n,n);
    Ytrue.resize(n);
    Ytrans.resize(n);
    a.resize(n);
    Perms.resize(n);
    for(i=0;i<n;i++){
      Perms(i) = i;
    }      
  }
    
  MatrixXd xz       = MatrixXd::Zero(n,gp1);
  VectorXd e        = VectorXd::Zero(n);
  VectorXi ns       = VectorXi::Zero(g);
  VectorXd pdfRow   = VectorXd::Zero(g);
  VectorXi perms    = VectorXi::Zero(p);    
        

  // Tools
  VectorXd t     = VectorXd::Zero(gp1);
  MatrixXd D     = MatrixXd::Zero(gp1,gp1);
  MatrixXd D_    = MatrixXd::Zero(gp1,gp1);
  VectorXd H     = VectorXd::Zero(gp1);
    
  // Additionnal tools
  VectorXd d_    = VectorXd::Zero(p);
  MatrixXd vD_vT = MatrixXd::Zero(p,p);
  MatrixXd vd_   = MatrixXd::Zero(p,p);
  VectorXd mu    = VectorXd::Zero(p);
  VectorXd nu    = VectorXd::Zero(p);  
  VectorXd xy    = VectorXd::Zero(p);
  Vector2d var;
   
  for(j=0;j<p;j++){
    perms(j) = j;
  }
  if(!io->IsModelInitialized){
    double m_j;
    double v_j;
    double i_j;
    for(j=0;j<p;j++){
      m_j  = (io->xTy(j) - io->sx(j)*io->sy/n )/(io->xTx(j) - io->sx(j)*io->sx(j)/n );
      i_j  = (io->sy-m_j*io->sx(j))/n;
      for(i=0;i<n;i++){
	e(i) =  io->y(i) - i_j - m_j * io->x(i,j);	  
      }
      v_j  = e.norm()/n;
      B(j) = rnorm( m_j , v_j/io->xTx(j) );
    }
    init_kmeans(io->sparse); // <= Modifies Z
    gamma2 = 0.0;
    for(j=0;j<p;j++){      
      k = Z(j);
      gamma2 += (B(j) - b(k)) * (B(j) - b(k));
      ns(k)++;      
    }
    gamma2 = gamma2/p;
    intercept = io->sy/n;
    for(k=0;k<g;k++){
      pi(k) = 1.0/g;
    }
  }else{
    intercept = io->intercept;
    for(k=0;k<g;k++){
      b(k)  = io->b(k);
      pi(k) = io->pi(k);
    }
    sigma2 = io->sigma2;
    gamma2 = io->gamma2;
    for(j=0;j<p;j++){
      Z(j) = io->Z0(j); 
    }
  }
        
  for(k=0;k<g;k++){
    ns(k) = 0; 
  }
  for(j=0;j<p;j++){
    k = Z(j);
    ns(k)++;
  }
  // Calculate XZ = X . Z
  for(i=0;i<n;i++){
    xz(i,0) = io->su(i);
    for(k=0;k<g;k++){
      xz(i,1+k) = 0.0;
    }
    for(j=0;j<p;j++){
      xz(i,1+Z(j)) += io->x(i,j);
    }
  }
  t(0) = intercept;
  for(k=0;k<g;k++){
    t(1+k) = b(k);
  }
  if(io->family=="binomial"){
    Ytrue  = io->y;
    Ytrans = io->U.transpose() * Ytrue;
  }
  if(io->family=="binomial"){
    e = Ytrans-xz*t;
  }else{
    e = io->y-xz*t;
  }
  RSS = e.squaredNorm();
  if(io->family=="binomial"){
    sigma2 = 1.0;
  }else{
    sigma2 = RSS/n;
  }
   
  int nMC     = nItMC;
  double llz  = jointLikelihood(io,e);
  
  for(item=0;item<nItEM;item++){     
    if(io->family=="binomial"){
      for(i=0;i<n;i++){
	a(i) = 1.0/(1.0 + gamma2 * io->v(i));
      }
      Hm = io->U * ( a.asDiagonal() * io->U.transpose() );
    }
    if(item>=nBurn){
      nMC = 1;
    }
    for(itmc=0;itmc<nMC;itmc++){       
      updateZ_GibbsRows(io,xz,e,ns,pdfRow,perms,p);
      if(io->family=="binomial"){
	// Just in case xz = x . z
	for(i=0;i<n;i++){
	  xz(i,0) = io->su(i);
	  for(k=0;k<g;k++){
	    xz(i,1+k) = 0.0;
	  }
	  for(j=0;j<p;j++){
	    xz(i,1+Z(j)) += io->x(i,j);
	  }
	}
	a = (io->U) * (xz * t);
	updateY_Gibbs(Ytrue,io->y,a,Hm,n,Perms);
	//e -= Ytrans;
	Ytrans = io->U.transpose() * Ytrue;
	e = Ytrans-xz*t;
      }
    }
     
    // Re-fill void class     
    if(ns.minCoeff()>0){
      // -----> M-step 
      for(k=0;k<g;k++){
	ns(k) = 0; 
      }
      for(j=0;j<p;j++){
	ns(Z(j))++;
      }
      for(k=0;k<g;k++){
	pi(k) = ((double) ns(k))/p;
      }
      var(0) = sigma2;
      var(1) = gamma2;
      
      if(io->lambda<=0.0){
	if(io->family=="binomial"){
	  llz    = EM(xz,Ytrans,io->v,e,t,var,D,D_,H,gp1,n,io->maxit,io->tol,io->sparse,io->family);
	}else{
	  llz    = EM(xz,io->y,io->v,e,t,var,D,D_,H,gp1,n,io->maxit,io->tol,io->sparse,io->family);
	}
      }else{
	if(io->family=="binomial"){
	  llz    = Mstep(xz,Ytrans,io->v,e,t,var,D,D_,H,gp1,n,io->sparse,io->family,io->lambda);
	}else{
	  llz    = Mstep(xz,io->y,io->v,e,t,var,D,D_,H,gp1,n,io->sparse,io->family,io->lambda);
	}
      }
      intercept = t(0);
      for(k=0;k<g;k++){
	b(k) = t(1+k);
      }
      sigma2 = var(0);
      gamma2 = var(1);
     
      for(j=0;j<p;j++){
	llz += log( pi(Z(j)) );
      }
    }else{
      llz = jointLikelihood(io,e);
    }
    Theta(item,0) = intercept;
    for(k=0;k<g;k++){
      Theta(item,1+k)   = b(k);
      Theta(item,1+g+k) = pi(k);
    }
    Theta(item,2*g+1) = sigma2;
    Theta(item,2*g+2) = gamma2;
    Theta(item,2*g+3) = llz;     
  }
  // Calculate theta Maximum Likelihood Estimates
  if(nItEM>nBurn){
    for(k=0;k<g;k++){
      b(k)  = 0.0;
      pi(k) = 0.0;
    }
    sigma2  = 0.0;
    gamma2   = 0.0;
    intercept = 0.0;
    
    for(item=nBurn;item<nItEM;item++){
      intercept   += Theta(item,0);      
      for(k=0;k<g;k++){
	b(k)      += Theta(item,1+k);
	pi(k)     += Theta(item,1+g+k);
      }
      sigma2      += Theta(item,2*g+1);
      gamma2      += Theta(item,2*g+2);
    }
    for(k=0;k<g;k++){
      b(k)  = b(k)/(nItEM-nBurn);
      pi(k) = pi(k)/(nItEM-nBurn);
    }
    sigma2    = sigma2/(nItEM-nBurn);
    gamma2    = gamma2/(nItEM-nBurn);
    intercept = intercept/(nItEM-nBurn);
  }    
  // Calculate Likelihood
  for(i=0;i<n;i++){
    xz(i,0) = io->su(i);
    for(k=0;k<g;k++){
      xz(i,1+k) = 0.0;
    }
    for(j=0;j<p;j++){
      xz(i,1+Z(j)) += io->x(i,j);
    }
  }
  t(0) = intercept;
  for(k=0;k<g;k++){
    t(1+k) = b(k);
  }
  if(io->family=="binomial"){
    e = Ytrans-xz*t;
  }else{
    e = io->y-xz*t;
  }
  CalculateLikelihood(io,xz,e,ns,pdfRow,perms,a,Hm,Ytrue,Ytrans,t,Perms,d_,vD_vT,vd_,mu,nu,xy);
};


void Model::CalculateLikelihood (IO *io,MatrixXd &xz,VectorXd &e,VectorXi &ns,VectorXd &pdfRow,VectorXi &perms,
				 VectorXd &a, MatrixXd &Hm, VectorXd &Ytrue, VectorXd &Ytrans,VectorXd &t,VectorXi &Perms,
				 VectorXd &d_,MatrixXd &vD_vT,MatrixXd &vd_, VectorXd &mu,VectorXd &nu,VectorXd &xy){
  int i,j,k;
  double lambda = sigma2/gamma2;    
  for(j=0;j<p;j++){
    if(j<io->n){
      d_(j) = 1.0/( io->v(j) + lambda );
    }else{
      d_(j) = 1.0/( lambda );
    }
  }
  vD_vT = (io->V * d_.asDiagonal())*io->V.transpose();
  for(j=0;j<p;j++){
    d_(j) = sqrt(d_(j));
  }
  vd_ = io->V * d_.asDiagonal();

  if(io->family=="binomial"){
    for(i=0;i<io->n;i++){
      //a(i) = io->v(i) / (io->v(i) + 1.0/gamma2);
      a(i) = 1.0/(1.0 + gamma2 * io->v(i));

    }
    Hm = io->U * ( a.asDiagonal() * io->U.transpose() );
  }
  // Sample now and calculate likelihood and entropy
  for(j=0;j<p;j++){
    for(k=0;k<g;k++){
      P(j,k) = 0.0;
    }
  }
  double llz,Lmax = -INFINITY;
  int itmc,index;
  
  xy = io->xTy;
  for(itmc=0;itmc<(io->dp)*(io->nsample);itmc++){   
    updateZ_GibbsRows(io,xz,e,ns,pdfRow,perms,p);
    if(io->family=="binomial"){
      a = io->U * (xz * t);
      updateY_Gibbs(Ytrue,io->y,a,Hm,io->n,Perms);
      e -= Ytrans;
      Ytrans = io->U.transpose() * Ytrue;
      e += Ytrans;
      xy = io->x.transpose() * Ytrans;
    }
    // sample Beta
    for(j=0;j<p;j++){
      nu(j) = xy(j) - intercept*io->sx(j) + b(Z(j)) * lambda;
    }
    mu = vD_vT * nu;
    for(j=0;j<p;j++){
      nu(j) = rnorm(0.0,sigma2);
    }
    B  = mu + vd_ * nu;     
    if(itmc%(io->dp)==0){
      llz = jointLikelihood(io,e);
      
      index = itmc/(io->dp);
      Lmc(index)  = llz;
      if(Lmc(index) > Lmax ){
	Lmax = Lmc(index);
      }
      for(j=0;j<p;j++){
	Bw(j,index) = B(j);
	Zw(j,index) = Z(j);
	P(j,Z(j))  += 1.0;
      }
    }
  }
  P = P/nsample;
  entropy = 0.0;
  for(j=0;j<p;j++){
    for(k=0;k<g;k++){
      if( P(j,k)>0.0){
	entropy += -P(j,k) * log( P(j,k) );
      }
    }
  }
  likelihood = -log((double) io->nsample) + Lmax;
  double tmp = 0.0;
  for(itmc=0;itmc<nsample;itmc++){
    tmp += exp( Lmc(itmc) - Lmax);
  }
  likelihood += log(tmp);   
};


// Added the 15/10/2013 at 02:33
void Model::fitMCEM(IO *io,MatrixXd &Theta){
  int n = io->n;
  int nsample = io->nsample;
  int i,j,k;
  int item,itmc;
  double RSS;
  int nItEM = io->nItEM;
  int nBurn = io->nBurn;
  int dp    = io->dp;
  
  VectorXd e        = VectorXd::Zero(n);
  VectorXi ns       = VectorXi::Zero(g);
  VectorXd pdfRow   = VectorXd::Zero(g);
  
  // Additionnal tools
  VectorXd d_    = VectorXd::Zero(p);
  MatrixXd vD_vT = MatrixXd::Zero(p,p);
  MatrixXd vd_   = MatrixXd::Zero(p,p);
  VectorXd mu    = VectorXd::Zero(p);
  VectorXd nu    = VectorXd::Zero(p);  
  VectorXd Sig   = VectorXd::Zero(io->nsample);
  
  if(!io->IsModelInitialized){
    double m_j;
    double v_j;
    double i_j;
    for(j=0;j<p;j++){
      m_j  = (io->xTy(j) - io->sx(j)*io->sy/n )/(io->xTx(j) - io->sx(j)*io->sx(j)/n );
      i_j  = (io->sy-m_j*io->sx(j))/n;
      for(i=0;i<n;i++){
	e(i) =  io->y(i) - i_j - m_j * io->x(i,j);	  
      }
      v_j  = e.norm()/n;
      B(j) = rnorm( m_j , v_j/io->xTx(j) );
    }
    init_kmeans(io->sparse); // <= Modifies Z
    gamma2 = 0.0;
    for(j=0;j<p;j++){      
      k = Z(j);
      gamma2 += (B(j) - b(k)) * (B(j) - b(k));
      ns(k)++;      
    }
    gamma2 = gamma2/p;
    intercept = io->sy/n;
    for(k=0;k<g;k++){
      pi(k) = 1.0/g;
    }
  }else{
    intercept = io->intercept;
    for(k=0;k<g;k++){
      b(k)  = io->b(k);
      pi(k) = io->pi(k);
    }
    sigma2 = io->sigma2;
    gamma2 = io->gamma2;
    for(j=0;j<p;j++){
      Z(j) = io->Z0(j); 
    }
  }
  
  ns.setZero();
  for(j=0;j<p;j++){
    ns(Z(j))++;
  }
  
  e = io->y-io->x*B;
  for(i=0;i<n;i++){
    e(i) += -intercept; 
  }
  
  RSS = e.squaredNorm();
  sigma2 = RSS/n;
  
  int M = nsample;
  double llz  = jointLikelihood2(io,e);
  double lambda;
  double Lmax = -INFINITY;
  int m;

  for(item=0;item<nItEM;item++){    
    P.setZero();
    lambda = sigma2/gamma2;
    for(j=0;j<p;j++){
      if(j<io->n){
    	d_(j) = 1.0/( io->v(j) + lambda );
      }else{
    	d_(j) = 1.0/( lambda );
      }
    }
    vD_vT = (io->V * d_.asDiagonal())*io->V.transpose();
    for(j=0;j<p;j++){
      d_(j) = sqrt(d_(j));
    }
    vd_ = io->V * d_.asDiagonal();
    ns.setZero();
    
    for(itmc=0;itmc<nBurn+dp*M;itmc++){
      	// Sample Beta
	for(j=0;j<p;j++){
	  nu(j) = io->xTy(j) - intercept*io->sx(j) + b(Z(j)) * lambda;
	}
	mu = vD_vT * nu;
	for(j=0;j<p;j++){
	  nu(j) = rnorm(0.0,sigma2);
	}
	B = mu + vd_ * nu;
	// Sample Z
	double lpmax,norm;
	int ng;
	for(j=0;j<p;j++){	  
	  lpmax = -INFINITY;
	  for(k=0;k<g;k++){
	    pdfRow(k) = -0.5*(B(j)-b(k))*(B(j)-b(k))/gamma2 + log( pi(k) );
	    if( pdfRow(k)>lpmax ){
	      lpmax = pdfRow(k);
	    }
	  }
	  norm  = 0.0;
	  for(k=0;k<g;k++){
	    pdfRow(k) = exp(pdfRow(k)-lpmax);
	    norm     += pdfRow(k);
	  }
	  for(k=0;k<g;k++){
	    pdfRow(k) = pdfRow(k)/norm;
	  }
	  ng  = rmultinom(pdfRow,g);
	  Z(j) = ng;
	}
	// Calculate residuals and log-likelihood
	if(itmc>=nBurn and (itmc-nBurn)%dp==0){ // >= or just > ?
	  e = io->y-io->x*B;
	  for(i=0;i<io->n;i++){
	    e(i) += -intercept; 
    	  }
    	  m = (itmc-nBurn)/dp;
    	  Sig(m) = e.squaredNorm();
    	  llz = jointLikelihood2(io,e);
    	  Lmc(m) = llz;
    	  if( Lmc(m) > Lmax ){
            Lmax = Lmc(m);
    	  }
    	  for(j=0;j<p;j++){	  
	    k = Z(j);
    	    Bw(j,m) = B(j);
    	    Zw(j,m) = k;
    	    P(j,k)     += 1.0;
	    ns(k)++;
    	  }
	}
      }
	  // Re-fill void class     
	  if(ns.minCoeff()>0){
	    // M-step
	    pi.setZero();
	    b.setZero();
	    ns.setZero();
	    
	    for(j=0;j<p;j++){
	      B(j) = 0.0;
	      for(m=0;m<M;m++){
		  k = Zw(j,m);
		  ns(k)++;
		  b(k) += Bw(j,m);
		  B(j) += Bw(j,m);
	       }
	       B(j) = B(j) / M;
	     }
	     intercept = (io->sy - io->sx.dot(B))/n; // here le 16/10 au soir chez Sadia
	     for(k=0;k<g;k++){
		pi(k) = ( (double) ns(k) ) / (p * M);
		b(k)  = b(k) / ns(k);
	      }
	      if(io->sparse){
		b(0) = 0.0;
	      }
	      gamma2 = 0.0;
	      for(j=0;j<p;j++){
		for(m=0;m<M;m++){
		  k = Zw(j,m);
		  gamma2 += (Bw(j,m)-b(k)) * (Bw(j,m)-b(k));
		}
	      }
	      gamma2 = gamma2 / (p * M);
	      sigma2 = 0.0;
	      for(m=0;m<M;m++){
		sigma2 += Sig(m);
	      }
	      sigma2 = sigma2 / (M * n);
	      llz = Lmc.mean();
	    }
	    
	    Theta(item,0) = intercept;
	    for(k=0;k<g;k++){
	      Theta(item,1+k)   = b(k);
	      Theta(item,1+g+k) = pi(k);
	    }
	    Theta(item,2*g+1) = sigma2;
	    Theta(item,2*g+2) = gamma2;
	    Theta(item,2*g+3) = llz;     
	  }
	  
	  // Calculate Likelihood
	  P = P/M;
	  entropy = 0.0;
	  for(j=0;j<p;j++){
	   for(k=0;k<g;k++){
	    if( P(j,k)>0.0 ){
	      entropy += -P(j,k) * log( P(j,k) );
	    }
	  }
	 }
	 likelihood = -log(M) + Lmax;
	 double tmp = 0.0;
	 for(m=0;m<M;m++){
	  tmp += exp( Lmc(m) - Lmax);
	 }
	 likelihood += log(tmp);
};
    
