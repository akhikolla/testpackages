#include "sgd.h"
using namespace Eigen;

std::tuple<VectorXd,double,VectorXd,VectorXd,MatrixXd,int,bool> SGD(std::function<double(const VectorXd &, const ArrayXd &)> objective,std::function<VectorXd(const VectorXd &,const ArrayXi &)> gradient,std::function<VectorXd(int)> x, const VectorXd & ibeta, std::vector<double> samplingProbabilities,const VectorXd stepSize, int batchsize, int nIters, std::function<VectorXd(const VectorXd &, const ArrayXd &)> SoftThreshold,std::function<VectorXd(VectorXd &, VectorXd &, ArrayXd &)> subgradient,  ArrayXd lambdaj, double eps, int trace)
{
    int N(samplingProbabilities.size()), m(batchsize);
    VectorXd beta_n(ibeta), beta_n1(ibeta);
    MatrixXd betaMat;
    VectorXd subgrad(ibeta.size()); subgrad.setZero();//generalized gradient
    double fVal(0);
    VectorXd fVal_all;fVal_all.setZero();
    VectorXd s_gradient_full(N),grad_full(ibeta.size()); //Full gradient
    VectorXd s_gradient_m(m), grad_n(ibeta.size()),xs;
    ArrayXi ridx_f(N),ridx(m);
    for (int s=0;s<N;s++){ridx_f(s)=s;}
    double gs(0);//scalar gradient
    s_gradient_m.setZero(); grad_n.setZero();xs.setZero();
    int iter(0);
    double betadiff(1);
    bool converged(false);
    
    if(trace>=1){
        betaMat.resize(ibeta.size(),(nIters+1));betaMat.setZero();
        fVal_all.resize(nIters+1); fVal_all.setZero();
    }
    // C++ version random number generator
     // std::discrete_distribution<> dist(samplingProbabilities.begin(),samplingProbabilities.end());
     // std::random_device r;
     // std::default_random_engine generator(r());
    
    // R
   Rcpp::NumericVector ridx_Rf(N), ridx_R(m);
   for(int i=0;i<N;i++){ridx_Rf[i]=i;}
   Rcpp::NumericVector prob(samplingProbabilities.begin(),samplingProbabilities.end());
    
    while(iter<nIters&&!converged){
        
        switch(trace){
            case 1:
                betaMat.col(iter) = beta_n; break;
            case 2:
                fVal_all(iter)= objective(beta_n,lambdaj); break;
            case 3:
                betaMat.col(iter) = beta_n;
                fVal_all(iter)= objective(beta_n,lambdaj); break;
            default:
                break;
        }
        
        // Select random m samples
        // c++
         // for(int s=0;s<m;s++){
         //     ridx(s) = dist(generator);
         // }
        // R
       ridx_R = Rcpp::sample(ridx_Rf,m,true,prob);
       std::vector<int> ret_vec = Rcpp::as<std::vector<int> >(ridx_R);
       ridx = Eigen::Map<Eigen::ArrayXi>(ret_vec.data(),m);
        
        // Calculate scalar part of gradients of m points at beta_n
        // grad_i(beta_n) = gi*(1,xi)^T
        s_gradient_m = gradient(beta_n,ridx);// m by 1
        
        // Average over m points. grad_n = weighted mean(grad_i)
        grad_n.setZero();
        for(int s=0;s<m;s++){
            // Rcpp::Rcout<<"sample: "<<ridx[s]<<endl;
            xs = x(ridx(s));
            //            cout<<"s_gradient_m(s)"<<s_gradient_m(s)<<endl;
            gs = s_gradient_m(s)/samplingProbabilities[ridx(s)];
            grad_n += gs*xs;
        }
        grad_n/=(m*N);
        beta_n -= stepSize(iter)*grad_n;
        beta_n = SoftThreshold(beta_n, lambdaj*stepSize(iter));
        betadiff=(beta_n-beta_n1).array().abs().maxCoeff();
        beta_n1 = beta_n;
        if(betadiff<eps){converged=true;}
        iter++;
    }
    // At the end, calculate fVal, subGradient
    fVal = objective(beta_n,lambdaj);
    
    switch(trace){
        case 1:
            betaMat.col(iter) = beta_n; break;
        case 2:
            fVal_all(iter)= objective(beta_n,lambdaj); break;
        case 3:
            betaMat.col(iter) = beta_n;
            fVal_all(iter)= objective(beta_n,lambdaj); break;
        default:
            break;
    }
    
    grad_full.setZero();
    s_gradient_full = gradient(beta_n,ridx_f);
    for (int s=0; s<N; s++){
        xs = x(s);
        gs = s_gradient_full(s);
        grad_full += gs*xs;
    }
    grad_full/=N;
    subgrad = subgradient(grad_full,beta_n,lambdaj);
    return std::make_tuple(beta_n,fVal,subgrad,fVal_all,betaMat,iter,converged);
}


