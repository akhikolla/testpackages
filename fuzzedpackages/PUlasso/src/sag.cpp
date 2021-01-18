#include "sgd.h"
using namespace Eigen;

std::tuple<VectorXd,double,VectorXd,VectorXd,MatrixXd,int,bool> SAG(std::function<double(const VectorXd &, const ArrayXd &)> objective,std::function<VectorXd(const VectorXd &,const ArrayXi &)> gradient, std::function<VectorXd(int)> x, const VectorXd & ibeta, std::vector<double> samplingProbabilities, double stepSize, int batchsize, int nIters, std::function<VectorXd(const VectorXd &, const ArrayXd &)> SoftThreshold, std::function<VectorXd(VectorXd &, VectorXd &, ArrayXd &)> subgradient, ArrayXd lambdaj, bool sampleSizeAdjustment, double eps, int trace)
{
    int N(samplingProbabilities.size()), m(batchsize);
    VectorXd beta_n=ibeta;
    MatrixXd betaMat;
    VectorXd subgrad(ibeta.size()); subgrad.setZero();
    double fVal(0);
    VectorXd fVal_all;
    VectorXd s_gradient_p(N), sum_grad_p(ibeta.size());//gradient evaluated at different theta
    VectorXd s_gradient_full(N),grad_full(ibeta.size()); //Full gradient
    VectorXd s_gradient_m(m), grad_update(ibeta.size());
    VectorXd xs;
    ArrayXi ridx_f(N),ridx(m);
    for (int s=0;s<N;s++){ridx_f(s)=s;}
    ArrayXi sampleIndicator(N); int sampleSize(0);
    double gs(0);
    s_gradient_p.setZero();
    s_gradient_m.setZero();
    sum_grad_p.setZero();
    xs.setZero();
    sampleIndicator.setZero();
    int iter(0);
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
//        if(trace){fVal_all(iter)= objective(beta_n,lambdaj);}
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
        // for(int s=0;s<m;s++){ridx(s) = dist(generator);}
        // R
               ridx_R = Rcpp::sample(ridx_Rf,m,true,prob);
               std::vector<int> ret_vec = Rcpp::as<std::vector<int> >(ridx_R);
               ridx = Eigen::Map<Eigen::ArrayXi>(ret_vec.data(),m);
            
        s_gradient_m = gradient(beta_n,ridx);// m by 1
        grad_update.setZero();
        
        // gradients at m points, adjusted by previous gradients at m points
        for(int s=0;s<m;s++){
            //            cout<<"sample: "<<ridx[s]<<endl;
            xs = x(ridx(s));
            gs = (s_gradient_m(s)-s_gradient_p(ridx(s)))/samplingProbabilities[ridx(s)];
            //            cout<<"gs: "<<gs<<endl;
            grad_update += gs*xs;
            s_gradient_p(ridx(s)) = s_gradient_m(s);
            if(sampleSizeAdjustment&&sampleSize<=N){
                sampleIndicator(ridx(s))=1;
                sampleSize = sampleIndicator.sum();
            }
        }
        grad_update/=(m*N);
        //        cout<<"grad_update:\n"<<grad_update<<endl;
        sum_grad_p +=grad_update;//update sum_grad_p
        //        cout<<"mean_grad:\n"<<(sum_grad_p/N)<<endl;
        if(sampleSizeAdjustment){
            beta_n -= stepSize*(sum_grad_p/sampleSize);
        }else{
            beta_n -= stepSize*(sum_grad_p/N);
        }
        
        beta_n = SoftThreshold(beta_n, lambdaj*stepSize);
        //        cout<<"grad_n:\n"<<grad_n<<endl;
        //        cout<<"stepSize:"<<stepSize<<endl;
        //        cout<<"beta_n:\n"<<beta_n<<endl;
        if(iter%N==0&& iter!=0){
            grad_full.setZero();
            s_gradient_full = gradient(beta_n,ridx_f);
            for (int s=0; s<N; s++){
                xs = x(s);
                gs = s_gradient_full(s);
                grad_full += gs*xs;
            }
            grad_full/=N;
            
            subgrad = grad_full;
            subgrad = subgradient(grad_full,beta_n,lambdaj);
            if(subgrad.array().abs().maxCoeff()<eps){converged=true;}
        }
        iter++;
    }
    fVal = objective(beta_n,lambdaj);
//    if(trace){fVal_all(iter)= objective(beta_n,lambdaj);}
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



