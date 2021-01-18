#ifndef pgLUfit_h
#define pgLUfit_h
#include "pgGroupLasso.h"
#include "sgd.h"

// #include <random>
using namespace Eigen;
using namespace std;
template <class TX>
class pgLUfit : public pgGroupLassoFit<TX>
{
protected:
    using pgGroupLassoFit<TX>::X;// with intercept, N by p matrix, p = 1+k1+..+k(J-1)
    using pgGroupLassoFit<TX>::y;// size N
    using pgGroupLassoFit<TX>::beta;// size p
    using pgGroupLassoFit<TX>::pi;
    using pgGroupLassoFit<TX>::gsize;// size J, first group = intercept
    using pgGroupLassoFit<TX>::pen; // size J, first element = 0;
    using pgGroupLassoFit<TX>::lambdaseq; //size k, default 100
    using pgGroupLassoFit<TX>::isUserLambdaseq;
    using pgGroupLassoFit<TX>::pathLength;
    using pgGroupLassoFit<TX>::lambdaMinRatio;
    using pgGroupLassoFit<TX>::maxit;
    using pgGroupLassoFit<TX>::tol;
    using pgGroupLassoFit<TX>::verbose;
    using pgGroupLassoFit<TX>::trace;

//
    //Definition Inside
    VectorXd stepSizeSeq;
    using pgGroupLassoFit<TX>::default_lambdaseq;
    using pgGroupLassoFit<TX>::grpSIdx;//size J
    using pgGroupLassoFit<TX>::iters;
//    using pgGroupLassoFit<TX>::Rinvs;
    using pgGroupLassoFit<TX>::coefficients; //size p*k
    using pgGroupLassoFit<TX>::std_coefficients; //size p*k
//    using pgGroupLassoFit<TX>::Xcenter;
    using pgGroupLassoFit<TX>::iter;// current iterations
    using pgGroupLassoFit<TX>::convFlag;
//
    //Dimension Information
    using pgGroupLassoFit<TX>::N;
    using pgGroupLassoFit<TX>::J;
    using pgGroupLassoFit<TX>::p;
    using pgGroupLassoFit<TX>::K;
//
//    //function
//    using pgGroupLassoFit<TX>::linpred;
//
//    ///////////////////////////////////////////
//
    
    double stepSize;
    double stepSizeAdj;
    int batchSize;
    int updateFreq;
    bool useLipschitz;
    std::vector<double> samplingProbabilities;
    std::string method;
    VectorXd Deviances;
    double nullDev;
    VectorXd fVals;
    MatrixXd subgrads;
    MatrixXd fVals_all;
    MatrixXd beta_all;
    VectorXd L;
//

    using pgGroupLassoFit<TX>::q;
    double evalDev(const VectorXd & lpred);
    
public:
    pgLUfit(TX & X_, VectorXd & z_, VectorXd & icoef_, ArrayXd & gsize_,ArrayXd & pen_,ArrayXd & lambdaseq_,bool isUserLambdaseq_,int pathLength_,double lambdaMinRatio_, double pi_, int maxit_, double tol_,bool verbose_, double stepSize_, double stepSizeAdj_, int batchSize_, int updateFreq_, std::vector<double> samplingProbabilities_,bool useLipschitz_,std::string method_,int trace_);
   
    void pgLUfit_main();
    using pgGroupLassoFit<TX>::computeLambdaSequence;
    using pgGroupLassoFit<TX>::getCoefficients;
    using pgGroupLassoFit<TX>::getStdCoefficients;
    using pgGroupLassoFit<TX>::getIters;
    using pgGroupLassoFit<TX>::getconvFlag;
    double getnullDev();
    VectorXd getDeviances();
    VectorXd getfVals();
    SparseMatrix<double> getfVals_all();
    SparseMatrix<double> getbeta_all();
    MatrixXd getSubGradients();
    double getStepSize();
    std::vector<double> getSamplingProbabilities();
    using pgGroupLassoFit<TX>::back_to_org;
    using pgGroupLassoFit<TX>::org_to_std;
    using pgGroupLassoFit<TX>::evalObjective;
    using pgGroupLassoFit<TX>::standardizeX;
    using pgGroupLassoFit<TX>::destandardizeX;
    
//    using pgGroupLassoFit<TX>::decenterX;
    
};

#endif /* pgLUfit_h */

