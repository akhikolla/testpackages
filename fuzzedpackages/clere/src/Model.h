#ifndef MODEL_H
#define MODEL_H
#define N_STEP_MAX   100
//#include <vector>
//#include <algorithm>
#include <cstdio>
#include <cmath>
#include <cstdlib>
//#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "IO.h"

using namespace std;
using namespace Eigen;

class Model{
  
 private:
  int  p,g;          // Parameters dimension
  VectorXd b;        // Centers
  VectorXd pi;       // Proportion
  double sigma2;     // Residual Variance
  double gamma2;     // Intra-Component Variance
  double intercept;  // Intercept
  int nsample;
  
  VectorXd B;         // Toy Beta vector
  VectorXi Z;         // Toy Z vector
  
  MatrixXd P;         // Assignment probabilities
  MatrixXi Zw;        // Sample from p(Z|y,x;theta)
  MatrixXd Bw;        // Sample from p(B|y,x;theta)
  VectorXd Lmc;       // Sampled joint log-likelihood
  
  double likelihood;
  double entropy;
  string logMessage;
  bool initialized; 
  
  double jointLikelihood(IO *io,VectorXd &e){
    int i,j;
    double tmp;
    double L = -(io->n)*log(2.0 * M_PI);
    for(i=0;i<io->n;i++){
      tmp = sigma2 + io->v(i) * gamma2;
      L  += -log(tmp) - e(i)*e(i)/tmp;
    }
    L = 0.5*L;
    for(j=0;j<(io->p);j++){
      L += log( pi(Z(j)) );
    }
    return L;
  };
  double jointLikelihood2(IO *io,VectorXd &e){
    int j,k;
    double tmp;
    double L = -0.5*(io->n)*log(2.0 * M_PI * sigma2) - 0.5*e.squaredNorm()/sigma2;
    for(j=0;j<(io->p);j++){
      k = Z(j);
      tmp = B(j)-b(k);
      L += log( pi(k) ) - 0.5*tmp*tmp/gamma2;
    }
    return L;
  };  
  
  
 public:
  Model(){};
  Model(int mp, int mg,int mns);
  ~Model(){};
  void  init_basic(bool sparse);
  void  init_kmeans(bool sparse);

  void   updateZ_GibbsRows(IO *io,MatrixXd &xz,VectorXd &e,VectorXi &ns,VectorXd &pdfRow,VectorXi &perms, int nChanges);
  
  
  void   fitSEM(IO *io, MatrixXd &Theta);
  void   fitMCEM(IO *io, MatrixXd &Theta);
  void   CalculateLikelihood(IO *io,MatrixXd &xz,VectorXd &e,VectorXi &ns,VectorXd &pdfRow,VectorXi &perms,
			     VectorXd &a, MatrixXd &Hm, VectorXd &Ytrue, VectorXd &Ytrans,VectorXd &t,VectorXi &Perms,
			     VectorXd &d_,MatrixXd &vD_vT,MatrixXd &vd_, VectorXd &mu,VectorXd &nu,VectorXd &xy);
  void   updateY_Gibbs(VectorXd &Ytrue,  VectorXd &C, VectorXd &a, MatrixXd &Hm,int n,VectorXi &Perms);
  
  double get_likelihood(){
    return likelihood;
  };
  double get_entropy(){
    return entropy;
  };
  double get_gamma2(){
    return gamma2;
  };
  double get_intercept(){
    return intercept;
  };
  double get_sigma2(){
    return sigma2;
  };
  VectorXd get_b(){
    return b;
  }
  VectorXd get_pi(){
    return pi;
  }
  MatrixXd get_P(){
    return P;
  }  
  MatrixXi get_Zw(){
    return Zw;
  }  
  MatrixXd get_Bw(){
    return Bw;
  }
};
#endif
