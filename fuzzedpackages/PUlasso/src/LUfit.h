#ifndef LUfit_h
#define LUfit_h
#include "groupLasso.h"

using namespace Eigen;
using namespace std;
template <class TX>
class LUfit : public groupLassoFit<TX>
{
protected:
  using groupLassoFit<TX>::X;// with intercept, N by p matrix, p = 1+k1+..+k(J-1)
  using groupLassoFit<TX>::y;// size N
  using groupLassoFit<TX>::beta;// size p
  using groupLassoFit<TX>::gsize;// size J, first group = intercept
  using groupLassoFit<TX>::pen; // size J, first element = 0;
  using groupLassoFit<TX>::lambdaseq; //size k, default 100
  using groupLassoFit<TX>::isUserLambdaseq;
  using groupLassoFit<TX>::pathLength;
  using groupLassoFit<TX>::lambdaMinRatio;
  using groupLassoFit<TX>::maxit;
  using groupLassoFit<TX>::wei;
  using groupLassoFit<TX>::weiOption;
  using groupLassoFit<TX>::tol;
  using groupLassoFit<TX>::verbose;
  using groupLassoFit<TX>::trace;
  
  //Definition Inside
  using groupLassoFit<TX>::centerFlag;
  using groupLassoFit<TX>::default_lambdaseq;
  using groupLassoFit<TX>::resid;
  using groupLassoFit<TX>::grpSIdx;//size J
  using groupLassoFit<TX>::iters;
  using groupLassoFit<TX>::Rinvs;
  using groupLassoFit<TX>::coefficients; //size p*k
  using groupLassoFit<TX>::std_coefficients; //size p*k
  using groupLassoFit<TX>::Xcenter; 
  using groupLassoFit<TX>::iter;// current iterations
  //using groupLassoFit<TX>::intercept_set;
  using groupLassoFit<TX>::activeSet;
  using groupLassoFit<TX>::inactiveSet;
  using groupLassoFit<TX>::inactiveSet1;
  using groupLassoFit<TX>::inactiveSet2;
  using groupLassoFit<TX>::g;
  using groupLassoFit<TX>::convFlag;
  
  //Dimension Information
  using groupLassoFit<TX>::N;
  using groupLassoFit<TX>::J;
  using groupLassoFit<TX>::p;
  using groupLassoFit<TX>::K;
  using groupLassoFit<TX>::linpred;
  using groupLassoFit<TX>::linpred_update;
  using groupLassoFit<TX>::coordinateDescent_0;
  using groupLassoFit<TX>::quadraticBCD;
  using groupLassoFit<TX>::gr;
  //using groupLassoFit<TX>::setupinactiveSets;
  
  ///////////////////////////////////////////
  
  const double t; //Hessian bound
  VectorXd lresp;//latent response
  double pi;
  int nUpdate;
  int max_nUpdates;
  double inner_tol;
  bool useStrongSet;
  
  //
  int nl;
  int nu;
  double wnl;
  double wnu;
  double bias;
  double c0;
  ArrayXi nUpdates;
  
  VectorXd lpred_old;
  VectorXd mu;// mu = 1/(1+exp(-Qbeta))
  VectorXd resid_old;
  VectorXd Deviances;
  double nullDev;
  VectorXd fVals;//function values
  MatrixXd subgrads;//generalized gradients
  MatrixXd fVals_all;//all function values
  MatrixXd beta_all;
  
  //private functions
  VectorXd convert_mu(const VectorXd & beta);
  VectorXd convert_mustar(const VectorXd & beta);
  void compute_mu_mustar(const VectorXd & lpred, VectorXd & mu, VectorXd & mustar);
  
  void updateObjFunc(VectorXd & lpred); //update mu, lresp, resid at new Qbeta
  void setupinactiveSets(int k, const VectorXd & resid, double lam_max,
                         const ArrayXd & lambdaseq,bool useStrongSet);
  double evalDev(const VectorXd & lpred);
  double evalObjective(const VectorXd & lpred, const VectorXd & beta, const ArrayXd & lambda);
  VectorXd  evalObjectiveGrad(const VectorXd & lpred);
  ArrayXd lambda_b(int k, const ArrayXd & pen);
  
public:
  LUfit(TX & X_, VectorXd & z_, VectorXd & icoef_, ArrayXd & gsize_,ArrayXd & pen_,
        ArrayXd & lambdaseq_, bool isUserLambdaseq_,int pathLength_,double lambdaMinRatio_,
        double pi_, int max_nUpdates_, int maxit_, VectorXd & wei_, bool weiOption_,
        double tol_, double inner_tol_,bool useStrongSet_,bool verbose_, int trace_);
  
  void LUfit_main();
  using groupLassoFit<TX>::computeLambdaSequence;
  using groupLassoFit<TX>::getCoefficients;
  using groupLassoFit<TX>::getStdCoefficients;
  using groupLassoFit<TX>::getIters;
  using groupLassoFit<TX>::getconvFlag;
  ArrayXi getnUpdates();
  double getnullDev();
  VectorXd getDeviances();
  VectorXd getfVals();
  MatrixXd getSubGradients();
  SparseMatrix<double> getfVals_all();
  SparseMatrix<double> getbeta_all();
  using groupLassoFit<TX>::back_to_org;
  using groupLassoFit<TX>::org_to_std;
  using groupLassoFit<TX>::decenterX;
  
  
};

template <class TX>
double evalDeviance(const TX & X, const VectorXd & z, const double pi, const VectorXd & coef, const VectorXd & wei, bool weiOption);


#endif /* LUfit_h */
