#ifndef pgGroupLasso_h
#define pgGroupLasso_h

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <algorithm>
#include <Rcpp.h>
// #include <iostream>

using namespace Eigen;
using namespace std;
template <class TX>
class pgGroupLassoFit
{
//protected:
public:
    //Input needed
    TX & X;// with intercept, N by p matrix, p = 1+k1+..+k(J-1)
    VectorXd & y;// size N
    double  pi;
    ArrayXd & gsize;// size J, first group = intercept
    ArrayXd & pen; // size J, first element = 0;
    ArrayXd & lambdaseq;//size K, default 100
    bool isUserLambdaseq;
    int pathLength;
    double lambdaMinRatio;
    int maxit;
    double tol;
    bool verbose;
    int trace;
    
    //
    int iter;// current iterations
    
    //Dimension Information
    int N;
    int nl;
    int nu;
    int J;
    int p;
    int K;

    //Definition Inside
    ArrayXi grpSIdx;//size J
    ArrayXi iters;
    MatrixXd coefficients; //size p*k
    MatrixXd std_coefficients;
    VectorXd Xcenter;
    std::vector<MatrixXd> Rinvs;
    VectorXd beta;// size p
    ArrayXd default_lambdaseq;
    ArrayXi convFlag;

    //These constructors will be called only from derived classes
    pgGroupLassoFit(TX & X_, VectorXd & y_, double pi_, VectorXd & icoef_, ArrayXd & gsize_,ArrayXd & pen_,ArrayXd & lambdaseq_,bool isUserLambdaseq_,  int pathLength_,
                  double lambdaMinRatio_,int maxit_, double tol_, bool verbose_,int trace_);
    ~pgGroupLassoFit(){
        destandardizeX();
    }
    void Rinvs_X();
    VectorXd linpred(bool intercept, const VectorXd & beta, const ArrayXi & ridx);
    double evalObjective(const VectorXd & beta, const ArrayXd & lambda);
    VectorXd gradient(const VectorXd & beta, const ArrayXi & ridx);
    VectorXd subgradient(VectorXd & gradient, VectorXd & beta, ArrayXd lambdaj);
    VectorXd SoftThreshold(const VectorXd & beta, const ArrayXd & thresh);
    VectorXd q(int i);
//public:
    //getters
    ArrayXi getIters();
    MatrixXd getCoefficients();
    MatrixXd getStdCoefficients();
    ArrayXd getLambdaSequence();
    ArrayXi getconvFlag();
    
    //Misc functions
    VectorXd back_to_org(const VectorXd & beta);
    VectorXd org_to_std(const VectorXd & coef);
    ArrayXd computeLambdaSequence(const VectorXd & resp);
    void checkDesignMatrix(const TX & X);
    void destandardizeX();
    void standardizeX();
};

#endif /* pgGroupLasso_h */
