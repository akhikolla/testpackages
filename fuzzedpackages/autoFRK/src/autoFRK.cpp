#include <Rcpp.h>
#include <RcppEigen.h>
#include <SymEigs.h>
#include <iostream>
#include <math.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
typedef Map<MatrixXd> MapMatd;
using Eigen::VectorXd;                  
using Eigen::SelfAdjointEigenSolver;   

using namespace Eigen;
Eigen::VectorXd getASCeigenValues(const Eigen::Map<Eigen::MatrixXd> A) {
    SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    return es.eigenvalues();
}

using namespace Eigen;
// [[Rcpp::export]]
Rcpp::List getASCeigens(const Eigen::Map<Eigen::MatrixXd> A) {
    SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
	Eigen::MatrixXd V;
	Eigen::VectorXd lambda;
	lambda = es.eigenvalues();
	V = es.eigenvectors();
      return Rcpp::List::create(Rcpp::Named("value") = lambda,
                            Rcpp::Named("vector") = V);
}


using namespace Spectra;
using namespace Eigen;
using namespace Rcpp;

void mrtseigencpp(const Eigen::MatrixXd & M, const int ncv, const int k, Eigen::VectorXd &rho, Eigen::MatrixXd &gamma){

  DenseSymMatProd<double> op(M);

  // Construct eigen solver object, requesting the largest three eigenvalues
  SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, k, ncv);

  eigs.init();
  eigs.compute(1000, 1e-10);

  rho = eigs.eigenvalues();
  gamma.noalias() = eigs.eigenvectors();
}


using namespace Rcpp;
using namespace std;
using namespace Eigen;
void tpm2(const MatrixXd P, MatrixXd& L, int p, int d){

  double r;

  if(d == 1){
    for(unsigned int i = 0; i < p; i++){
      for(unsigned int j = i + 1; j < p; ++j){
        r  = abs(P(i, 0) - P(j, 0));
        L(i,j) = pow(r, 3)/12;
      }
    }
  }
  else if(d == 2){
    for(unsigned int i = 0; i < p; i++){
      for(unsigned int j = i + 1; j < p; ++j){
        r  = sqrt(pow(P(i, 0) - P(j, 0),2) +
          (pow(P(i, 1) - P(j, 1), 2)));
        if(r != 0)
          L(i, j) = r*r*log(r)/(8.0*M_PI);
      }
    }
  }
  else if(d ==3){
    for(unsigned int i = 0; i < p; i++){
      for(unsigned int j = i + 1; j < p; ++j){
        r = sqrt(pow(P(i, 0) - P(j, 0), 2) +
          pow(P(i, 1) - P(j, 1), 2) +
          pow(P(i, 2) - P(j, 2), 2));
        L(i, j) = -r/8;
      }
    }
  }

}


using namespace Rcpp;
using namespace std;
using namespace Eigen;
void tpm_predict(const MatrixXd P_new, const MatrixXd P, MatrixXd& L, int d){

  double r;
  int p1 = P_new.rows();
  int p2 = P.rows();

  if(d == 1){
    for(unsigned int i = 0; i < p1; i++){
      for(unsigned int j = 0; j < p2; ++j){
        r  = abs(P_new(i, 0) - P(j, 0));
        L(i,j) = pow(r, 3)/12;
      }
    }
  }
  else if(d == 2){
    for(unsigned int i = 0; i < p1; i++){
      for(unsigned int j = 0; j < p2; ++j){
        r  = sqrt(pow(P_new(i, 0) - P(j, 0),2) +
          (pow(P_new(i, 1) - P(j, 1), 2)));
        if(r != 0)
          L(i, j) = r*r*log(r)/(8.0*M_PI);
      }
    }
  }
  else if(d ==3){
    for(unsigned int i = 0; i < p1; i++){
      for(unsigned int j = 0; j < p2; ++j){
        r = sqrt(pow(P_new(i, 0) - P(j, 0), 2) +
          pow(P_new(i, 1) - P(j, 1), 2) +
          pow(P_new(i, 2) - P(j, 2), 2));
        L(i, j) = -r/8;
      }
    }
  }

}



using namespace Eigen;
using namespace Rcpp;

void mrts(const Eigen::MatrixXd &Xu,
          const Eigen::MatrixXd &xobs_diag,
          const int k,
          const int n,
          const int d,
          Eigen::MatrixXd &H,
          Eigen::MatrixXd &X,
          Eigen::MatrixXd &UZ,
          Eigen::MatrixXd &BBB,
          Eigen::VectorXd &nconst){
  double root = sqrt(n);
  Eigen::MatrixXd B, gammas, X2_temp, gamma;
  Eigen::VectorXd rho;

  H = MatrixXd::Zero(n, n);
  B = MatrixXd::Ones(n, d + 1);
  X =  MatrixXd::Ones(n, k + d + 1);

  tpm2(Xu, H, n, d);

  H += H.transpose().eval();
  B.rightCols(d) = Xu;
  const Eigen::MatrixXd Bt = B.transpose();
  Eigen::MatrixXd BtB(MatrixXd(d+1, d+1).setZero().
                        selfadjointView<Lower>().rankUpdate(Bt));

  const Eigen::LLT<MatrixXd> llt(BtB);
  BBB = llt.solve(Bt);
  const Eigen::MatrixXd AH = H - (H * B) * BBB;
  const Eigen::MatrixXd AHA = AH - BBB.transpose() * (Bt * AH);

  int ncv = min(n, max(2 * k + 1, 20));

  mrtseigencpp(AHA, ncv, k, rho, gamma);


  gammas = (gamma - B * (BBB * gamma)).array().rowwise() / rho.transpose().array() * root;
  Eigen::MatrixXd X_temp = Xu.rowwise() - Xu.colwise().mean();
  nconst = X_temp.colwise().norm();

  X.block(0, 1, n, d) = X_temp.array().rowwise() * (root / nconst.transpose().array());
  X.block(0, d + 1, n, k) = gamma * root;

  UZ = MatrixXd::Zero(n + d + 1, k + d + 1);
  UZ.block(0, 0, n, k) = gammas;
  UZ(n , k) = 1;
  UZ.block(n + 1, k + 1, d, d) = xobs_diag;
  nconst /= root;
}



using namespace Eigen;
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List mrtsrcpp(const Eigen::Map<Eigen::MatrixXd> Xu,
                    const Eigen::Map<Eigen::MatrixXd> xobs_diag,
                    const int k){

  int n(Xu.rows()), d(Xu.cols());
  Eigen::MatrixXd  H, X, UZ, BBB;
  Eigen::VectorXd nconst;

  mrts(Xu, xobs_diag, k, n, d, H, X, UZ, BBB, nconst);


  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("UZ") = UZ,
                            Rcpp::Named("BBBH") = BBB * H,
                            Rcpp::Named("nconst") = nconst);
}


using namespace Eigen;
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List mrtsrcpp_predict0(const Eigen::Map<Eigen::MatrixXd> Xu,
                            const Eigen::Map<Eigen::MatrixXd> xobs_diag,
                            const Eigen::Map<Eigen::MatrixXd> xnew,
                            const int k){


  int n(Xu.rows()), d(Xu.cols()), n2(xnew.rows());
  Eigen::MatrixXd  H, X, UZ, BBB, Hnew;
  Eigen::VectorXd nconst;

  mrts(Xu, xobs_diag, k, n, d, H, X, UZ, BBB, nconst);

  Hnew = MatrixXd::Zero(n2, n);

  tpm_predict(xnew, Xu, Hnew, d);


  Eigen::MatrixXd X1 = Hnew * UZ.block(0, 0, n, k);
  Eigen::MatrixXd B = MatrixXd::Ones(n2, d + 1);
  B.rightCols(d) = xnew;

  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("UZ") = UZ,
                            Rcpp::Named("BBBH") = BBB * H,
                            Rcpp::Named("nconst") = nconst,
                            Rcpp::Named("X1") = X1 - B*((BBB* H) *UZ.block(0, 0, n, k)));
                            
}


using namespace Eigen;
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List mrtsrcpp_predict(const Eigen::Map<Eigen::MatrixXd> Xu,
                            const Eigen::Map<Eigen::MatrixXd> xobs_diag,
                            const Eigen::Map<Eigen::MatrixXd> xnew,
							const Eigen::Map<Eigen::MatrixXd> BBBH,
							const Eigen::Map<Eigen::MatrixXd> UZ,
							const Eigen::Map<Eigen::VectorXd> nconst,
                            const int k){


  int n(Xu.rows()), d(Xu.cols()), n2(xnew.rows());
  Eigen::MatrixXd  Hnew;
  
  Hnew = MatrixXd::Zero(n2, n);

  tpm_predict(xnew, Xu, Hnew, d);


  Eigen::MatrixXd X1 = Hnew * UZ.block(0, 0, n, k);
  Eigen::MatrixXd B = MatrixXd::Ones(n2, d + 1);
  B.rightCols(d) = xnew;

  return Rcpp::List::create(Rcpp::Named("X") = Xu,
                            Rcpp::Named("UZ") = UZ,
                            Rcpp::Named("BBBH") = BBBH,
                            Rcpp::Named("nconst") = nconst,
                            Rcpp::Named("X1") = X1 - B*((BBBH) *UZ.block(0, 0, n, k)));
                            
}
