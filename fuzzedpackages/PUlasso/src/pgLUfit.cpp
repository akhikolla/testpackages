#include "pgLUfit.h"
using namespace Eigen;
using namespace std::placeholders;

template <class TX>
pgLUfit<TX>::pgLUfit(TX & X_, VectorXd & z_, VectorXd & icoef_, ArrayXd & gsize_,ArrayXd & pen_,ArrayXd & lambdaseq_,bool isUserLambdaseq_,int pathLength_,double lambdaMinRatio_, double pi_, int maxit_, double tol_,bool verbose_, double stepSize_, double stepSizeAdj_, int batchSize_, int updateFreq_, std::vector<double> samplingProbabilities_,bool useLipschitz_,std::string method_,int trace_):
pgGroupLassoFit<TX>(X_,z_,pi_,icoef_,gsize_,pen_,lambdaseq_,isUserLambdaseq_,pathLength_,lambdaMinRatio_,maxit_,tol_,verbose_,trace_),stepSize(stepSize_),stepSizeAdj(stepSizeAdj_),batchSize(batchSize_), updateFreq(updateFreq_), useLipschitz(useLipschitz_),samplingProbabilities(samplingProbabilities_),method(method_)
{
    VectorXd qi;
    L.resize(N);L.setZero();
    double L0(0);
    bool is_user_sampProb(true);
    if(samplingProbabilities_.size()!=N){
        is_user_sampProb = false;
        samplingProbabilities.resize(N);
    }
    
    //Either to calculate L0 or samplingProb
    if(useLipschitz){
        for (int i=0;i<N;i++){
            qi = q(i);
            L(i)=qi.adjoint()*qi;
            if(!is_user_sampProb){samplingProbabilities[i]= L(i);}
        }
        if(!is_user_sampProb){
            double temp_sum; // calculate sum Li, the normalize samplingProbabilities
            temp_sum = std::accumulate(samplingProbabilities.begin(),samplingProbabilities.end(),0);
            for (int i=0;i<N;i++){samplingProbabilities[i]/=temp_sum;}
        }
        if(stepSize==0){L0 = (L*3/(4*N)).sum();stepSize=stepSizeAdj/L0;}
    }
    stepSizeSeq.resize(maxit);
    Deviances = VectorXd::Zero(K);
    fVals = VectorXd::Zero(K);
    subgrads = MatrixXd::Zero(p,K);

    //trace
    //trace == 0 , no trace; trace == 1 , trace beta; trace==2, trace fVal, trace==3, trace all.
    if(trace>=1){
        fVals_all.resize(maxit+1, K);fVals_all.setZero();
        beta_all.resize(p*K, maxit+1);beta_all.setZero();
    }
    
    VectorXd lpred0(N),beta0(p);
    lpred0 = VectorXd::Ones(N)*std::log(pi/(1-pi));
    beta0 << std::log(pi/(1-pi)),VectorXd::Zero(p-1); //I-only, no need to standardize
    nullDev = evalDev(lpred0);
    
    //Initialize beta
    beta = org_to_std(icoef_);
    if(!beta.segment(1,(p-1)).any())
    {
        beta = beta0 ;
    }
    
    VectorXd resp0(N);
    VectorXd z0(N);
    z0 =1.0-z_.array();
    resp0 = z_+(z0*pi);
    
    default_lambdaseq = computeLambdaSequence(resp0);
    if(!isUserLambdaseq){lambdaseq = default_lambdaseq;}

};

//Getters
template <typename TX>
double pgLUfit<TX>::getnullDev(){return nullDev;}
template <typename TX>
VectorXd pgLUfit<TX>::getDeviances(){return Deviances;}
template <typename TX>
VectorXd pgLUfit<TX>::getfVals(){return fVals;}
template <typename TX>
SparseMatrix<double> pgLUfit<TX>::getfVals_all(){return fVals_all.sparseView();}
template <typename TX>
SparseMatrix<double> pgLUfit<TX>::getbeta_all(){return beta_all.sparseView();}
template <typename TX>
MatrixXd pgLUfit<TX>::getSubGradients(){return subgrads;}
template <typename TX>
double pgLUfit<TX>::getStepSize(){return stepSize;}
template <typename TX>
std::vector<double> pgLUfit<TX>::getSamplingProbabilities(){return samplingProbabilities;}
//calculate deviance using precalculated lpred
template <class TX>
double pgLUfit<TX>::evalDev(const VectorXd & lpred)
{
    int nl(y.sum());
    int nu = N-nl;
    const double c = std::log(nl/(pi*nu));
    VectorXd pred, logExpLpred, logExpPred,obslogL;
    logExpLpred = (lpred.array().exp().array()+1).array().log();
    pred = c+lpred.array()-logExpLpred.array();
    logExpPred = (1+pred.array().exp()).array().log();
    obslogL = (y.array()*pred.array()-logExpPred.array());//response z_lu = y
    return -2*obslogL.sum();
}

template <class TX>
void pgLUfit<TX>::pgLUfit_main(){
    
    std::function<double(VectorXd,ArrayXd)> f=std::bind(&pgGroupLassoFit<TX>::evalObjective,this,_1,_2);
    std::function<VectorXd(VectorXd,const ArrayXi &)> g=std::bind(&pgGroupLassoFit<TX>::gradient,this,_1,_2);
    std::function<VectorXd(int)> q = std::bind(&pgGroupLassoFit<TX>::q,this,_1);
    std::function<VectorXd(const VectorXd &, const ArrayXd &)> ST = std::bind(&pgGroupLassoFit<TX>::SoftThreshold,this,_1,_2);
    std::function<VectorXd(VectorXd &, VectorXd &, ArrayXd &)> subgradient = std::bind(&pgGroupLassoFit<TX>::subgradient,this,_1,_2,_3);
    ArrayXd lambda_k(K);
    
    int method_int(0);
    if(method=="GD"){
        method_int=1;
    }else if(method=="SGD"){
        method_int=2;
    }else if(method=="SVRG"){
        method_int=3;
    }else if(method=="SAG"){
        method_int=4;
    }else{
        method_int=9;
    }

    for(int k=0; k<K; k++)
    {
        lambda_k = lambdaseq(k)* pen;
        VectorXd subgrad_k(p);
        double fVal_k(0);
        // Define without reserving memory
        VectorXd fVal_all_k;
        MatrixXd betaMat;
        if(trace>=1){
            fVal_all_k.resize(maxit);    fVal_all_k.setZero();
            betaMat.resize(p,(maxit+1)); betaMat.setZero();
        }
        
        bool converged(false);
        if(verbose){Rcpp::Rcout<<"Fitting "<<k<<"th lambda\n";}
        Rcpp::checkUserInterrupt();
        switch(method_int){
            case 1:
                std::tie(beta,fVal_k,subgrad_k,fVal_all_k,betaMat,iter,converged) = GD(f,g,N,q,beta,stepSize,maxit,ST,subgradient,lambda_k,tol,trace);
                break;
            case 2:
                for(int i=0;i<maxit;i++){
                  // stepSizeSeq(i) = stepSize/(1.0+stepSize*lambdaseq(k)*i);}
                stepSizeSeq(i)=stepSize;} // constant step size
                std::tie(beta,fVal_k,subgrad_k,fVal_all_k,betaMat,iter,converged) = SGD(f,g,q,beta,samplingProbabilities,stepSizeSeq,batchSize,maxit,ST,subgradient,lambda_k,tol,trace);
                break;
            case 3:
                 std::tie(beta,fVal_k,subgrad_k,fVal_all_k,betaMat,iter,converged) = SVRG(f,g,q,beta,samplingProbabilities,stepSize,updateFreq,batchSize,maxit,ST,subgradient,lambda_k,tol,trace);
                break;
            case 4:
                std::tie(beta,fVal_k,subgrad_k,fVal_all_k,betaMat,iter,converged)= SAG(f,g,q,beta,samplingProbabilities,stepSize,batchSize,maxit,ST,subgradient,lambda_k,true,tol,trace);
                break;
            default:
                std::tie(beta,fVal_k,subgrad_k,fVal_all_k,betaMat,iter,converged) = GD(f,g,N,q,beta,stepSize,maxit,ST,subgradient,lambda_k,tol,trace);
                break;}
        
        Map<MatrixXd> coefficients_k(&coefficients.coeffRef(0, k),p,1);
        Map<MatrixXd> std_coefficients_k(&std_coefficients.coeffRef(0, k),p,1);
        coefficients_k = back_to_org(beta);
        std_coefficients_k = beta;
        iters(k)=iter;
        fVals(k)=fVal_k;
        double penVal(0);
        VectorXd bj;
        for (int j=0;j<J;j++){
            bj=beta.segment(grpSIdx(j)+1,gsize(j));
            penVal+=lambda_k(j)*bj.lpNorm<2>();
        }
        Deviances(k) = (fVal_k-N*penVal)*2;
        subgrads.col(k)=subgrad_k;
        if(trace>=1){
            fVals_all.col(k)=fVal_all_k;
            beta_all.block(k*p,0,p,(iter+1))=betaMat.block(0,0,p,(iter+1));
        }
        
        if(!converged){convFlag(k)=1;}
       if(verbose&&convFlag(k)==0){Rcpp::Rcout<<"converged at "<<iter<<"th iterations\n";}
    }
}

////The explicit instantiation part
template class pgLUfit<MatrixXd>;
template class pgLUfit<SparseMatrix<double> >;
template class pgLUfit<Map<MatrixXd> >;
