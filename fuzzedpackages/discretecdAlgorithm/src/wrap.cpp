#include <Rcpp.h>
#include "dBasic.h"
#include "dCD.h"
#include "type.h"
#include <iostream>
using namespace Rcpp;

// temporary solution to "Found no calls to: R_registerRoutines, R_useDynamicSymbols"
// https://github.com/RcppCore/Rcpp/issues/636#issuecomment-280985661
void R_init_ccdrAlgorithm(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
// using Rcpp::as;
typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;
typedef Map<VectorXd> MapVecd;
typedef Map<VectorXi> MapVeci;


// [[Rcpp::export]]
List CD(            int node,
                    int dataSize,
                    Eigen::Map<Eigen::MatrixXi> data,
                    Eigen::Map<Eigen::VectorXi> nlevels,
                    List obsIndex_R,
                    int eor_nr,
                    Eigen::Map<Eigen::MatrixXi> eor,
                    Eigen::Map<Eigen::VectorXd> lambda_seq,
                    int nlam,
                    double eps,
                    double convLb,
                    double qtol,
                    Eigen::Map<Eigen::MatrixXd> weights,
                    double gamma,
                    double upperbound,
                    int threshold
) {

  // convert Rcpp type to C++ type
  // construct obsIndex
 VectorXVXi obsIndex(node);
 for (int i=0; i<node; i++) {
   obsIndex(i) = obsIndex_R[i];
 }

  // construct levelIndex
  MatrixXVXi levelIndex(node, node);
  for (int i=0; i<node; i++) {
    for (int j=0; j<node; j++) {
      if (j != i) {
        levelIndex(j, i).resize(nlevels[j] - 1);
        for (int k=0; k<(nlevels[j]-1); k++) {
          levelIndex(j, i)(k) = k;
        }
      }
      else {
        levelIndex(j, i).resize(nlevels[j]);
        for (int k=0; k<nlevels[j]; k++) {
          levelIndex(j, i)(k) = k;
        }
      }
    }
  }

  MatrixXMXd betaM(node+1, node);
  MatrixXMXd betaN(nlam*(node+1), node);
  Eigen::MatrixXi estimateG = Eigen::MatrixXi::Zero(node*nlam, node);
  Eigen::VectorXd log_like(nlam), dur(nlam);
  for (int i = 0; i < nlam; i++) {
    log_like(i) = 0.0;
    dur(i) = 0.0;
  }

  // Eigen::MatrixXi t_data(Rcpp::as< MapMati >(data));
  Eigen::MatrixXi t_data(dataSize, node);
  for (int i=0; i<dataSize; i++) {
    for (int j=0; j<node; j++) {
      t_data(i, j) = data(i, j);
    }
  }
  Eigen::MatrixXd t_weights(node, node);
  for (int i=0; i<node; i++) {
    for (int j=0; j<node; j++) {
      t_weights(i, j) = weights(i, j);
    }
  }

  Eigen::MatrixXi t_eor(eor_nr, 2);
  for (int i=0; i<eor_nr; i++) {
    for (int j=0; j <2; j++) {
      t_eor(i, j) = eor(i, j);
    }
  }

  Eigen::VectorXi t_nlevels(node);
  for (int i=0; i<node; i++) {
    t_nlevels(i) = nlevels(i);
  }

  Eigen::VectorXd lambdaSeq(nlam);
  for (int i=0; i<nlam; i++) {
    lambdaSeq(i) = lambda_seq(i);
  }

  // Run CDAlgorithm

  CDAlgo(node, dataSize, t_data, t_nlevels, obsIndex, levelIndex, eor_nr, t_eor,
         nlam, eps, convLb, qtol, lambdaSeq, log_like, dur, betaM, betaN,
         estimateG, t_weights, gamma, upperbound, threshold);

  NumericMatrix adaptWeights(node, node);
  for (int i=0; i<node; i++) {
    for (int j=0; j<node; j++) {
      adaptWeights(i, j) = betaM(i, j).norm();
    }
  }

  // IntegerMatrix outputG(node*nlam, node);
  // for (int i=0; i<(node*nlam); i++) {
  //   for (int j=0; j<node; j++) {
  //     outputG(i, j) = estimateG(i, j);
  //   }
  // }
  // should return lambdaSeq, time.
  // return estimateG, adaptive weight matrix, and timing, ;
  return List::create(_["estimateG"] = wrap(estimateG), _["time"] = wrap(dur), _["adaptive_weights"] = wrap(adaptWeights));
}

// [[Rcpp::export]]
double lambdaMax( int node,
         int dataSize,
         Eigen::Map<Eigen::MatrixXi> data,
         Eigen::Map<Eigen::VectorXi> nlevels,
         List obsIndex_R,
         Eigen::Map<Eigen::MatrixXd> weights,
         double gamma,
         double upperbound
) {

  // convert Rcpp type to C++ type
  // construct obsIndex
  VectorXVXi obsIndex(node);
  for (int i=0; i<node; i++) {
    obsIndex(i) = obsIndex_R[i];
  }

  // construct levelIndex
  MatrixXVXi levelIndex(node, node);
  for (int i=0; i<node; i++) {
    for (int j=0; j<node; j++) {
      if (j != i) {
        levelIndex(j, i).resize(nlevels[j] - 1);
        for (int k=0; k<(nlevels[j]-1); k++) {
          levelIndex(j, i)(k) = k;
        }
      }
      else {
        levelIndex(j, i).resize(nlevels[j]);
        for (int k=0; k<nlevels[j]; k++) {
          levelIndex(j, i)(k) = k;
        }
      }
    }
  }

  MatrixXMXd betaM(node+1, node);
  // Eigen::MatrixXi t_data(Rcpp::as< MapMati >(data));
  Eigen::MatrixXi t_data(dataSize, node);
  for (int i=0; i<dataSize; i++) {
    for (int j=0; j<node; j++) {
      t_data(i, j) = data(i, j);
    }
  }
  Eigen::MatrixXd t_weights(node, node);
  for (int i=0; i<node; i++) {
    for (int j=0; j<node; j++) {
      t_weights(i, j) = weights(i, j);
    }
  }

  Eigen::VectorXi t_nlevels(node);
  for (int i=0; i<node; i++) {
    t_nlevels(i) = nlevels(i);
  }

  double lambda = 0.0;

  // Run CDAlgorithm

  maxLambda(node, dataSize, t_data, t_nlevels, obsIndex, levelIndex, betaM, t_weights, lambda, gamma, upperbound);

  // should return lambdaSeq, time.
  // return estimateG;
  return lambda;
}

// [[Rcpp::export]]
IntegerMatrix DatGen(int maxdeg,
                     int node,
                     Eigen::Map<Eigen::MatrixXi> ordex,
                     IntegerVector ts,
                     int dataSize,
                     List ivn,
                     List ivn_vals,
                     bool ivn_rand,
                     IntegerVector coef_length,
                     Eigen::Map<Eigen::VectorXi> nlevels,
                     List coef)
{
  Eigen::MatrixXi t_ordex(maxdeg, node);
  for (int i=0; i<maxdeg; i++) {
    for (int j=0; j<node; j++) {
      t_ordex(i, j) = ordex(i, j);
    }
  }

  std::vector<int> t_ts(node);
  for (int i=0; i<node; i++) {
    t_ts[i] = ts[i];
  }

  VectorXVXi temp_ivn(dataSize);
  for (int i=0; i<dataSize; i++) {
    temp_ivn(i) = ivn[i];
  }
  std::vector< std::vector<int> > t_ivn(dataSize);
  for (int i=0; i<dataSize; i++) {
    for (int j=0; j<temp_ivn(i).size(); j++) {
      t_ivn[i].push_back(temp_ivn(i)(j));
    }
  }

  VectorXVXi temp_ivn_vals(dataSize);
  for (int i=0; i<dataSize; i++) {
    temp_ivn_vals(i) = ivn_vals[i];
  }
  std::vector< std::vector<int> > t_ivn_vals(dataSize);
  for (int i=0; i<dataSize; i++) {
    for (int j=0; j<temp_ivn(i).size(); j++) {
      t_ivn_vals[i].push_back(temp_ivn_vals(i)(j));
    }
  }

  Eigen::VectorXi t_nlevels(node);
  for (int i=0; i<node; i++) {
    t_nlevels(i) = nlevels(i);
  }

  Eigen::MatrixXi data = Eigen::MatrixXi::Zero(dataSize, node);

  std::vector<VectorXMXd> t_coef(node);
  for (int i=0; i<node; i++) {
    List single_nodeList = coef[i];
    if (coef_length[i]) {
      t_coef[i].resize(coef_length[i]);
      for (int l_index=0; l_index < coef_length[i]; l_index++) {
        t_coef[i](l_index) = single_nodeList[l_index];
      }
    }
  }

  // test for input

  // Rcout << "ordex is: \n";
  // Rcout << t_ordex << std::endl;
  // Rcout << std::endl;
  //
  // Rcout << "ts is: \n";
  // Rcout << t_ts[0] << " " << t_ts[1] << " " << t_ts[2] << " " << t_ts[3] <<" " << t_ts[4] << std::endl;
  // Rcout << std::endl;
  //
  // Rcout << "ivn is: \n";
  // for (int i=0; i<ivn.size(); i++) {
  //   for (int j=0; j<1; j++) {
  //     Rcout << t_ivn[i][j] << " ";
  //   }
  //   Rcout << std::endl;
  // }
  // Rcout << std::endl;
  //
  // Rcout << "nlevel is: \n";
  // Rcout << t_nlevels << std::endl;
  // Rcout << std::endl;
  //
  // Rcout << "data is: \n";
  // Rcout << data << std::endl;
  // Rcout << std::endl;
  //
  // Rcout << "coef is : \n";
  // for (int i=0; i<node; i++) {
  //   for (int j=0; j<(coef_length[i]); j++) {
  //     Rcpp::Rcout << t_coef[i](j) << std::endl;
  //   }
  //   Rprintf("\n");
  // }


  // run data generating function
  DatGen(t_ordex, t_ts, t_ivn, t_ivn_vals, ivn_rand, t_nlevels, data, t_coef);

  return wrap(data);
}
