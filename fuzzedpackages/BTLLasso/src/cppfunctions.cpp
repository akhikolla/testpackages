//#define ARMA_64BIT_WORD
// #define ARMA_DONT_USE_CXX11

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;


mat matsqrt2(mat A){
vec eigval;
mat eigvec;

  eig_sym(eigval, eigvec, A);

colvec d = (eigval+abs(eigval))*0.5;
colvec d2 = sqrt(d);
mat B = eigvec*diagmat(d2)*trans(eigvec);
  return B;
}

mat betaupdate(mat betaold, mat delta, mat checktheta, int ntheta){
int thetaOK = 0;
double update_fac = 1;
int rep = 1;
mat betareturn;

while(thetaOK == 0){
  thetaOK = 1;
update_fac = pow(0.5,(rep-1));
betareturn = betaold + update_fac*delta;
checktheta(span(0,(ntheta-1)),span::all) = betareturn(span(0,(ntheta-1)),span::all);
for(int xx = 1; xx < ntheta+1; xx += 1){
  if(checktheta(xx,0)-checktheta(xx-1,0)<0){
    Rcout<<"!Step-Halving!"<<endl;
    thetaOK = 0;
  }
}
rep = rep+1;
}
return betareturn;
}

mat mueta(mat eta){
  mat D;
  D = exp(eta)/(1+exp(eta))%(1-exp(eta)/(1+exp(eta)));
  return D;
}

mat L1(mat xi,
       double c){
         mat ret = 1/sqrt(xi%xi+c);
         return ret;
       }


mat group(mat xi,
            double c){
              int dfj = xi.n_rows;
              mat ret = sqrt((double)(dfj))/sqrt(trans(xi)*xi +c);
              return(repmat(ret,dfj,1));
            }

mat A(mat beta, 
      mat acoefs, 
      double lambda, 
      mat weight, 
      vec index,
      double c,
      double gama,
     std::string norm){

        vec nbs = trans(acoefs)*beta;
        uvec fd2 = find(nbs!=0);
        mat fd = zeros(nbs.n_elem,1);
        fd.rows(fd2) = ones(fd2.n_elem);      

      vec unindex = unique(index);
      mat appro(index.n_elem,1);
      std::string grouped = "grouped";
      std::string l1 = "L1";

      if(norm==grouped){

        int start = 0;
        for(double o=0; o<unindex.n_elem;o++){
          uvec sub = find(index==(o+1));
          int end = start + sub.n_elem -1;

          appro.rows(span(start,end)) = fd.rows(sub)%weight.rows(sub)%group(nbs.rows(sub),c)*lambda;
          start = end+1;
        }
      }else {
        if(norm==l1){
        appro = fd%weight%L1(nbs,c)*lambda;
        }else{
           appro = fd%weight*2*lambda;
        }
      }
mat Aret = trans(acoefs)%repmat(sqrt(appro),1,acoefs.n_rows);

return(trans(Aret)*Aret);
}



// [[Rcpp::export]]
List cumfit(NumericMatrix betanew2,
           double epsilon, 
           int maxiter,
           NumericMatrix acoefs2,
           double lambda,
           NumericMatrix weight2,
           List control,
           NumericMatrix design2,
           int N,
           int n,
           double q,
           NumericMatrix resp2,
           NumericVector index,
           double c,
           double gama,
           std::string norm,
           int hatmatrix,
           double lambda2,
           NumericMatrix checktheta2) { 



mat betanew = as<arma::mat>(betanew2);
mat acoefs = as<arma::mat>(acoefs2);
mat weight = as<arma::mat>(weight2);
mat design = as<arma::mat>(design2);
mat resp = as<arma::mat>(resp2);
mat checktheta = as<arma::mat>(checktheta2);


int i = 1;
double diff = 1.1;
int pp = design.n_cols;

//mat Sigmainv = mat(N,N);

mat eta(N,1);
mat Ai;

int ntheta = floor(q/2);

while((diff > epsilon) & (i < maxiter)){
     mat betaold = betanew;


Ai = A(betaold, 
      acoefs, 
      lambda, 
      weight, 
      index,
      c,
      gama,
      norm);


Ai = Ai + diagmat(ones(Ai.n_cols)*lambda2);

      
 eta = design*betanew;
 mat mu = exp(eta)/(1+exp(eta));

 

     mat D = mueta(eta);

    D = repmat(D,1,design.n_cols);
    
    mat designD = trans(design%D);
     mat scorepart = mat(pp,N);
     
//   if(q>1){
     for(int r=1; r<(n+1); r += 1){

       mat sig1 = mu(span((r-1)*q,r*q-1),0)*trans(1-mu(span((r-1)*q,r*q-1),0));

        sig1 = symmatu(sig1);
        sig1 = sig1 + diagmat(ones(sig1.n_cols)*0.0001);

        mat sig1i = inv(sig1);        
//        Sigmainv(span((r-1)*q,r*q-1),span((r-1)*q,r*q-1)) = sig1i;
        
        mat desDsub = designD(span::all,span((r-1)*q,r*q-1));
        
        scorepart(span::all,span((r-1)*q,r*q-1)) = desDsub*sig1i;
     }
//     }else{
//       Sigmainv = diagmat(1/(mu%(1-mu)));
//       mat scorepart = designD*Sigmainv;
//     }
  


     mat score = scorepart*(resp-mu) -Ai*betaold;

     mat Fisher = (scorepart%trans(D))*design +Ai;

     mat delta = solve(Fisher,score);
     
     betanew = betaupdate(betaold, delta, checktheta, ntheta);

  
     diff = sqrt(accu(pow(betaold-betanew,2)))/sqrt(accu(pow(betaold,2)));

     i += 1;
}


double df = 0;

//if(hatmatrix>0){
//
// eta = design*betanew;
//
//mat D = mueta(eta);
//
//     D = repmat(D,1,Sigmainv.n_cols);
//
//mat help = zeros(pp,pp);
//mat help2(N,pp);
//
//for(int r=1; r<(n+1); r += 1){
//mat Dsub = D(span((r-1)*q,r*q-1),span((r-1)*q,r*q-1));
//mat Wsub = Dsub%trans(Sigmainv(span((r-1)*q,r*q-1),span((r-1)*q,r*q-1))%Dsub);
//mat W12sub = matsqrt2(Wsub);
//         mat desSub = design(span((r-1)*q,r*q-1),span::all);
//        
//       help = help + trans(desSub)*Wsub*desSub ;
//        
//        help2(span((r-1)*q,r*q-1),span::all) = W12sub*design(span((r-1)*q,r*q-1),span::all);
//     }
//
//    help = help + Ai;
//
//
//mat H = help2*solve(help,trans(help2));
//
//rowvec diagH2 = sum(H%H,0);
//colvec dd = diagvec(2*H); 
//rowvec d2 = conv_to< rowvec >::from(dd);
//
//
//df = as_scalar(sum(d2-diagH2));
//}

return List::create(Named("beta.new") = betanew,
                      Named("df") = df);
}



// [[Rcpp::export]]
List binfit(NumericMatrix betanew2,
           double epsilon, 
           int maxiter,
           NumericMatrix acoefs2,
           double lambda,
           NumericMatrix weight2,
           List control,
           NumericMatrix design2,
           int N,
           int n,
           int q,
           NumericMatrix resp2,
           NumericVector index,
           double c,
           double gama,
           std::string norm,
           int hatmatrix,
           double lambda2) { 



mat betanew = as<arma::mat>(betanew2);
mat acoefs = as<arma::mat>(acoefs2);
mat weight = as<arma::mat>(weight2);
mat design = as<arma::mat>(design2);
mat resp = as<arma::mat>(resp2);


int i = 1;
double diff = 1.1;


mat D = mat(N,N);
mat Fisherinv;
mat eta(N,1);
mat Ai;



while((diff > epsilon) & (i < maxiter)){
     mat betaold = betanew;

Ai = A(betaold, 
      acoefs, 
      lambda, 
      weight, 
      index,
      c,
      gama,
      norm);


Ai = Ai + diagmat(ones(Ai.n_cols)*lambda2);

      
 eta = design*betanew;
 mat mu = exp(eta)/(1+exp(eta));
 
// D = diagmat(mueta(eta));

  mat D = mueta(eta);

    D = repmat(D,1,design.n_cols);

     mat score = trans(design)*(resp-mu) -Ai*betaold;

     mat Fisher = trans(design%D)*design +Ai;

     mat delta = solve(Fisher,score);
     betanew = betaold + delta;
     diff = sqrt(accu(pow(betaold-betanew,2)))/sqrt(accu(pow(betaold,2)));

     i += 1;
}

double df = 0;

if(hatmatrix>0){

 eta = design*betanew;

  mat D = mueta(eta);

    D = repmat(D,1,design.n_cols);

     mat W12 = diagmat(sqrt(mueta(eta)));

mat help = inv(trans(design%D)*design+Ai);

mat help2 = W12*design;


mat H = help2*help*trans(help2);

rowvec diagH2 = sum(H%H,0);
colvec dd = diagvec(2*H); 
rowvec d2 = conv_to< rowvec >::from(dd);


df = as_scalar(sum(d2-diagH2));
}

return List::create(Named("beta.new") = betanew,
                      Named("df") = df);
}

