#include "IO.h"
IO::IO(Rcpp::S4 & obj){
  int i,j;
  this->seed        =     Rcpp::as<int> (obj.slot("seed"));
  // Parameters for estimation
  this->nItEM       =     Rcpp::as<int> (obj.slot("nItEM"));
  this->nItMC       =     Rcpp::as<int> (obj.slot("nItMC"));
  this->nBurn       =     Rcpp::as<int> (obj.slot("nBurn"));
  this->dp          =     Rcpp::as<int> (obj.slot("dp"));
  this->nsample     =     Rcpp::as<int> (obj.slot("nsamp"));
  this->maxit       =     Rcpp::as<int> (obj.slot("maxit"));
  this->tol         =     Rcpp::as<double> (obj.slot("tol"));
  this->lambda      =     -1.0;//Rcpp::as<double> (obj.slot("lambda"));
  
  // Should we print: trace
  // Should allow one group with b=0
  
  this->sparse         =   Rcpp::as<bool> (obj.slot("sparse"));
  this->IsModelInitialized = Rcpp::as<bool> (obj.slot("initialized"));;
  
  // Size of the problem
  this->n   = Rcpp::as<int> (obj.slot("n"));
  this->p   = Rcpp::as<int> (obj.slot("p"));
  this->g   = Rcpp::as<int> (obj.slot("g"));
  
  this->analysis  = Rcpp::as<string> (obj.slot("analysis"));
  this->family    = "gaussian"; //Rcpp::as<string> (obj.slot("family"));
  this->algorithm = Rcpp::as<string> (obj.slot("algorithm"));
  
  if(this->g == 1){
    analysis = "fit";
  }
  this->instantiated = 1;
  if(this->instantiated!=0){
    this->x.resize(n,p);
    this->y.resize(n);

    this->xTx.resize(p);
    this->xTy.resize(p);
    this->sx.resize(p);
    this->v.resize(n);      
    
    // Read the data
    NumericMatrix Rx(SEXP(obj.slot("x")));
    NumericVector Ry(SEXP(obj.slot("y")));
    convertMatrix<NumericMatrix,MatrixXd>(Rx,x);      
    convertVector<NumericVector,VectorXd>(Ry,y);
    
    // Read initial parameters
    if( this->IsModelInitialized ){
      NumericVector Rb(SEXP(obj.slot("b")));
      NumericVector Rpi(SEXP(obj.slot("pi")));
      IntegerVector RZ0(SEXP(obj.slot("Z0")));

      convertVector<NumericVector,VectorXd>(Rb,b);
      convertVector<NumericVector,VectorXd>(Rpi,pi);
      convertVector<IntegerVector,VectorXi>(RZ0,Z0);

      this->intercept = Rcpp::as<double> (obj.slot("intercept"));
      this->sigma2 = Rcpp::as<double> (obj.slot("sigma2"));
      this->gamma2 = Rcpp::as<double> (obj.slot("gamma2"));
    }
    
    this->sy = this->y.sum();
    for(j=0;j<p;j++){
      this->sx(j)  = 0.0;
      this->xTx(j) = 0.0;
      this->xTy(j) = 0.0;
      for(i=0;i<n;i++){
	this->sx(j) += this->x(i,j);
	this->xTx(j)+= this->x(i,j) * this->x(i,j);
	this->xTy(j)+= this->x(i,j) * this->y(i);
      }
    }      
    this->su.resize(n);
    JacobiSVD<MatrixXd> svd(x, ComputeFullU | ComputeFullV );
    U = svd.matrixU();
    V = svd.matrixV();
    VectorXd s = svd.singularValues();
    int rgX = s.rows();
    for(i=0;i<n;i++){
      if(i<rgX){
	v(i) = s(i) * s(i);	  
      }else{
	v(i) = 0.0;  
      }
      su(i) = U.col(i).sum();
    }
    if(this->algorithm=="SEM"){
      x = U.transpose() * x;
      if(this->family!="binomial"){
	y = U.transpose() * y;
      }
    }
  }
};
void IO::updateObj(Rcpp::S4 & obj){
  // Parameters
  obj.slot("intercept")  = intercept;
  obj.slot("sigma2")     = sigma2;
  obj.slot("gamma2")     = gamma2;
  obj.slot("likelihood") = likelihood;
  obj.slot("entropy") = entropy; 
  obj.slot("b")       = outVector<NumericVector,VectorXd>(b);
  obj.slot("pi")      = outVector<NumericVector,VectorXd>(pi);
  obj.slot("P")       = outMatrix<NumericMatrix,MatrixXd>(P);
  obj.slot("theta")   = outMatrix<NumericMatrix,MatrixXd>(theta);
  obj.slot("Zw")      = outMatrix<IntegerMatrix,MatrixXi>(Zw); 
  obj.slot("Bw")      = outMatrix<NumericMatrix,MatrixXd>(Bw); 
};
