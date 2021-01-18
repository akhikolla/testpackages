#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericVector determinant_c( int S, int M, Eigen::Map<Eigen::MatrixXi> A, Eigen::Map<Eigen::MatrixXi> B, Eigen::Map<Eigen::VectorXd> c, Eigen::Map<Eigen::MatrixXd> pattern, Eigen::Map<Eigen::MatrixXd> Ineq, int s, int cc){


// for sampling random data
bool found = false;
Eigen::VectorXd thetaTMP;
Eigen::RowVectorXd theta;
Eigen::VectorXd IneqT;
bool check1;
bool check2;

// other vars
Eigen::MatrixXd Theta;
Eigen::VectorXd p;
Eigen::MatrixXd V;
Eigen::VectorXd pc;
Eigen::MatrixXd delta0;
double I;

// for output
NumericVector out(3);

out[0] = 0;
out[1] = 0;
out[2] = 0;

for (int i = 0; i < s; i++) {
  found = false;
  while (!found) {
    check1 = true;
    check2 = true;
    thetaTMP = Rcpp::as<Eigen::VectorXd>(rbeta(S, 0.5, 0.5));
    theta = thetaTMP.transpose();
	
	// added check to avoid multiplication with Ineq if only of size (1,1) (i.e., no inequality restrictions)
	if (Ineq.rows() == 1 && Ineq.cols() == 1 && Ineq(0,0) == 0)  {
		IneqT.resize(1);
		IneqT(0) = 0;
	} else
		IneqT = (theta*Ineq.transpose());
  
    for(int i = 0; i< IneqT.size(); i++){
        if(IneqT(i) < 0)
            check1 = false;
    }
    for(int i = 0; i<S;i++){
        if(theta(i) == 1 || theta(i) == 0)
            check2 = false;
    }
    if (check1 && check2) {
      found = true;
    }
    out[0] ++;
  }
  
  Theta = Eigen::VectorXd::Ones(M)*theta;
  // Rcout << "Theta:" << Theta << "\n";
  //Rcout << "IneqT \\< 0:" << (IneqT < 0) << "\n";
  //Rcout << ((theta*Ineq.transpose())<0) << "\n";
  //Rcout << "thetaTMP:\n"<<thetaTMP<<"\n";
  //Rcout << "theta:\n"<<theta<<"\n";
  
  Eigen::MatrixXd TA_tmp(Theta.rows(), Theta.cols());
  Eigen::MatrixXd TB_tmp(Theta.rows(), Theta.cols());
  
  for (int r = 0; r < Theta.rows(); r++) {
      for (int c = 0; c < Theta.cols(); c++) {
          //Rcout << "Theta: " << Theta(r, c) << "   A: " << A(r, c) << "\n";
          TA_tmp(r, c) = std::pow(Theta(r, c), A(r, c));
          TB_tmp(r, c) = std::pow((1 - Theta(r, c)), B(r, c));
      }
  }
  
  // VectorXd TA_tmp2 = TA_tmp.rowwise().prod();
  // VectorXd TB_tmp2 = TB_tmp.rowwise().prod();
  p = TA_tmp.rowwise().prod().array() * TB_tmp.rowwise().prod().array() * c.array();
  
  //Rcout << "p diagonal:\n" << p.asDiagonal() << "\n";
  
  // MatrixXd P_tmp = p.asDiagonal();
  // MatrixXd V = pattern * P_tmp;
  V = pattern * p.asDiagonal();
  //Rcout << "V:\n"<< V <<"\n";
  //MatrixXd out = Theta.array().pow(A.array());
  
  pc = V.rowwise().sum();
  Eigen::VectorXd D(pc.size());
  for (int c = 0; c < D.size(); c++) {
      if (pc(c) == 0)
          D(c) = 0;
      else D(c) = 1 / pc(c);
  }
  
  // MatrixXd thetaD = theta.asDiagonal();
  //Rcout << "thetaD:\n"<<thetaD<<"\n";
  // MatrixXd ABA = (A-(A+B)).cast<double>();
  // Rcout << "ABA:\n"<<ABA<<"\n";
  delta0 = V*(A.cast<double>()-(A+B).cast<double>()*theta.asDiagonal());
  
  //Rcout << "delta0:\n" << delta0 << "\n";
  
  I = fabs(((delta0.transpose() * D.asDiagonal() * delta0 * (theta.array()*(1 - theta.array())).matrix().asDiagonal().inverse() * M_PI * M_PI)*cc).determinant());
  // double detI = I.determinant();

  // Rcout << "log(i): " << II << "\tI: " << I << "\n";
  
  out[1] += I;
  out[2] += std::sqrt(I);
  
}
return out;


}



