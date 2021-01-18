#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


Eigen::SparseMatrix<double> ComputeEstGraph(Eigen::VectorXd alpha, Eigen::MatrixXd H, SEXP p, SEXP m, SEXP q_H){

using Eigen::Map;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::MappedSparseMatrix;

typedef Eigen::Triplet<double> T;

int p_s = Rcpp::as<int>(p);
int q_H_s = Rcpp::as<int>(q_H);
int m_s = Rcpp::as<int>(m);

VectorXd V_hat_ind = VectorXd(q_H_s).setZero();

int ind_i = 0;
int ind_j = 0;
double alpha_norm_l = 0;
int edge_count = 0;

for(int l = 0; l < q_H_s; l++){

  alpha_norm_l = alpha.segment(l*p_s, p_s).norm();

   if(alpha_norm_l == 0){
    V_hat_ind(l) = 1;
     edge_count += 1;
                       }
}


std::vector<T> tripletList;
tripletList.reserve(2*edge_count);

for(int l = 0; l < q_H_s; l++){

  ind_i = int(H(l, 0)) - 1;
  ind_j = int(H(l, 1)) - 1;

  if(V_hat_ind(l) == 1){

    tripletList.push_back(T(ind_i, ind_j, 1));
    tripletList.push_back(T(ind_j, ind_i, 1));
  }
                             }



Eigen::SparseMatrix<double> est_graph_adj(m_s, m_s);

est_graph_adj.setFromTriplets(tripletList.begin(), tripletList.end());


return est_graph_adj;
                           }
