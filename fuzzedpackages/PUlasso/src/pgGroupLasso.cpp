#include "pgGroupLasso.h"
using namespace Eigen;

//Constructor
template <class TX>
pgGroupLassoFit<TX>::pgGroupLassoFit(TX & X_, VectorXd & y_, double pi_, VectorXd & icoef_, ArrayXd & gsize_,ArrayXd & pen_,ArrayXd & lambdaseq_, bool isUserLambdaseq_,int pathLength_,double lambdaMinRatio_,int maxit_, double tol_, bool verbose_, int trace_)
:X(X_),y(y_), pi(pi_), gsize(gsize_), pen(pen_),lambdaseq(lambdaseq_), isUserLambdaseq(isUserLambdaseq_),pathLength(pathLength_),lambdaMinRatio(lambdaMinRatio_),maxit(maxit_), tol(tol_),verbose(verbose_),trace(trace_),iter(0)
{
    checkDesignMatrix(X);
    N = static_cast<int>(X.rows());
    nl = static_cast<int>(y.sum());
    nu = N-nl;
    p = static_cast<int>(X.cols())+1;
    J = static_cast<int>(gsize.size());
    K = isUserLambdaseq?(static_cast<int>(lambdaseq.size())):(pathLength);
    
    grpSIdx=ArrayXi::Zero(J);
    
    for(int ii=2;ii<J;++ii)
    {
        grpSIdx(ii)=grpSIdx(ii-1)+gsize(ii-1);
    }
    
    iters = ArrayXi::Zero(K);
    coefficients = MatrixXd::Zero(p, K);
    std_coefficients = MatrixXd::Zero(p, K);
    
    //Calculate Xcenter, Rinvs
    //For a dense class X, X = P0X, Sparse or Map class, no change in X
    
    Xcenter = VectorXd::Ones(p-1);
    Rinvs.resize(J);
    //if(verbose){Rcpp::Rcout<<"QR decompositions\n";}
    Rinvs_X();
    
    //Initialize beta
    beta =org_to_std(icoef_);
    
    convFlag.resize(K);
    convFlag.setZero();
    
    //if(verbose){Rcpp::Rcout<<"end of construction\n";}
}

//Getters
template <typename TX>
MatrixXd pgGroupLassoFit<TX>::getCoefficients(){return coefficients;}

template <typename TX>
MatrixXd pgGroupLassoFit<TX>::getStdCoefficients(){return std_coefficients;}

template <typename TX>
ArrayXi pgGroupLassoFit<TX>::getIters(){return iters;}

template <typename TX>
ArrayXd pgGroupLassoFit<TX>::getLambdaSequence(){return lambdaseq;}

template <typename TX>
ArrayXi pgGroupLassoFit<TX>::getconvFlag(){return convFlag;}

//Misc functions
template <typename TX>
VectorXd pgGroupLassoFit<TX>::back_to_org(const VectorXd & beta)
{
    VectorXd gamma(beta);
//    cout<<"back_to_org"<<endl;
    for(int j=1;j<J;++j)
    {
        gamma.segment(grpSIdx(j)+1,gsize(j))= Rinvs[j]*beta.segment(grpSIdx(j)+1,gsize(j));
//        cout<<"Rinvj:\n"<<Rinvs[j]<<endl;
//        cout<<"beta_seg:\n"<<beta.segment(grpSIdx(j)+1,gsize(j))<<endl;
//        cout<<"gamma_seg:\n"<<gamma.segment(grpSIdx(j)+1,gsize(j))<<endl;
    }
//    cout<<"gamma:\n"<<gamma<<endl;
    gamma(0) = beta(0)-gamma.segment(1,p-1).adjoint()*Xcenter;
//    cout<<"gamma:\n"<<gamma<<endl;
    return gamma;
}

template <typename TX>
VectorXd pgGroupLassoFit<TX>::org_to_std(const VectorXd & gamma)
{
    VectorXd beta(gamma);
    
    for(int j=1;j<J;++j)
    {
        beta.segment(grpSIdx(j)+1,gsize(j)) =
        Rinvs[j].inverse()*gamma.segment(grpSIdx(j)+1,gsize(j));
    }
    
    beta(0)= gamma(0)+ gamma.segment(1,p-1).adjoint()*Xcenter ;
    
    return beta;
}

//linpred = beta0+Q1beta1+...+Qpbetap where Qj = P0*Xj*Rinvj
//If intercept= false, beta0 = 0
//If dense matrix, X is already centered
template <typename TX>
VectorXd pgGroupLassoFit<TX>::linpred(bool intercept, const VectorXd & beta, const ArrayXi & ridx)
{
    VectorXd lpred(ridx.size());
    if(intercept){lpred.setConstant(beta(0));}else{lpred.setZero();}
    VectorXd tmp;

    for (int i=0;i<ridx.size();i++){
        for (int j=1; j<J; ++j){
            tmp=X.block(ridx(i),grpSIdx(j),1,gsize(j))*beta.segment(grpSIdx(j)+1,gsize(j));
            lpred(i)+=tmp.coeffRef(0, 0);
        }
    }
    
    return lpred;
}

template <>
VectorXd pgGroupLassoFit<SparseMatrix<double> >::linpred(bool intercept, const VectorXd & beta, const ArrayXi & ridx)
{
    VectorXd lpred(ridx.size());
    lpred.setZero();
    VectorXd tmp;
    std::vector<VectorXd> RinvBeta;
    RinvBeta.resize(J);
    VectorXd muRinvBeta(J);
    RinvBeta[0]=MatrixXd::Zero(1,1);muRinvBeta(0)=0;
    for (int j=1; j<J; ++j){
        RinvBeta[j] = Rinvs[j]*beta.segment(grpSIdx(j)+1,gsize(j));
        muRinvBeta(j) = Xcenter.segment(grpSIdx(j), gsize(j)).transpose()*RinvBeta[j];
//        cout<<"beta:\n"<<beta<<endl;
//        cout<<"RinvBeta["<<j<<"]\n"<<RinvBeta[j]<<endl;
//        cout<<"muRinvBeta["<<j<<"]\n"<<muRinvBeta(j)<<endl;
    }
    
    for (int i=0;i<ridx.size();i++){
        for (int j=1; j<J; ++j){
            tmp=X.block(ridx(i),grpSIdx(j),1,gsize(j))*RinvBeta[j];
            lpred(i)+=tmp.coeffRef(0, 0);
        }
    }
    if(intercept){
        lpred= lpred.array()-(muRinvBeta.sum()-beta(0));
    }else{
        lpred =lpred.array()-(muRinvBeta.sum());
    }
    return lpred;
}

template <class TX>
VectorXd pgGroupLassoFit<TX>::SoftThreshold(const VectorXd & beta, const ArrayXd & thresh)
{
    double bjnorm;
    VectorXd STbeta(beta);
    //No threshold for beta(0)
    for(int j=1; j<J;j++){
//        Map<VectorXd> bj(&beta.coeffRef(grpSIdx(j)+1),gsize(j));
        Map<VectorXd> STbj(&STbeta.coeffRef(grpSIdx(j)+1),gsize(j));
        bjnorm = beta.segment(grpSIdx(j)+1,gsize(j)).norm();
//        bjnorm = bj.norm();
        STbj = ((bjnorm>thresh(j))?(1-(thresh(j)/bjnorm)):0)*beta.segment(grpSIdx(j)+1,gsize(j));
    }
    return STbeta;
}


template <class TX>
ArrayXd pgGroupLassoFit<TX>::computeLambdaSequence(const VectorXd & y)
{
    ArrayXd lambda_path(pathLength);
//    std::vector<VectorXd> g;
//    g.resize(J);
    VectorXd gradnorm(J);
    VectorXd gj;
    double lammax(1);
    double TOLERANCE(1e-08);
    
    gradnorm.setZero();
    
    
    for (int j=1; j<J;++j)
    {
        int sind = grpSIdx(j);
        gj = (X.block(0,sind,N,gsize(j)).adjoint()*y)/N;
        gradnorm(j) = gj.norm()/pen(j);
    }
    
    lammax = gradnorm.maxCoeff()+TOLERANCE;
    
    double logDiff=std::log(lammax)-std::log(lambdaMinRatio*lammax);
    double ratio=std::exp(-logDiff/(pathLength-1));
    
    lambda_path(0)=lammax;
    for (int i=1; i<pathLength;++i)
    {
        lambda_path(i)=lambda_path(i-1)*ratio;
    }
    
    return lambda_path;
}

template <>
ArrayXd pgGroupLassoFit<SparseMatrix<double> >::computeLambdaSequence(const VectorXd & y)
{
    ArrayXd lambda_path(pathLength);
    //    std::vector<VectorXd> g;
    //    g.resize(J);
    VectorXd gradnorm(J);
    VectorXd gj;
    double lammax(1);
    double TOLERANCE(1e-08);
    
    gradnorm.setZero();
    VectorXd ycentered = y.array()-y.mean();
    
    for (int j=1; j<J;++j)
    {
        int sind = grpSIdx(j);
        gj = (Rinvs[j].adjoint()*(X.block(0,sind,N,gsize(j)).adjoint()*ycentered))/N;
        gradnorm(j) = gj.norm()/pen(j);
        
    }
    
    lammax = gradnorm.maxCoeff()+TOLERANCE;
    
    double logDiff=std::log(lammax)-std::log(lambdaMinRatio*lammax);
    double ratio=std::exp(-logDiff/(pathLength-1));
    
    lambda_path(0)=lammax;
    for (int i=1; i<pathLength;++i)
    {
        lambda_path(i)=lambda_path(i-1)*ratio;
    }
    
    return lambda_path;
}

//Rinvs
//X scaled and centered after
template <typename TX>
void pgGroupLassoFit<TX>::Rinvs_X()
{
    
    for (int l=0;l<(p-1);++l)
    {
        Xcenter(l) = X.col(l).mean();
        X.col(l) = X.col(l).array()-Xcenter(l);
    }
    
    for(int j=1;j<J;++j){
        int sind = grpSIdx(j);
        if(gsize(j)>1)
        {
            //Do QR decomposition
            ColPivHouseholderQR<MatrixXd> qr(X.block(0,sind,N,gsize(j)));
            
            if(qr.rank() < gsize(j)){throw std::invalid_argument("X(j) does not have full column rank");}
            
            MatrixXd R = qr.matrixR().topLeftCorner(qr.rank(), qr.rank()).triangularView<Upper>();
            MatrixXd P =qr.colsPermutation();
            R=R*P.inverse()/std::sqrt(N);
            
            Rinvs.at(j)= R.inverse();// QtQ = NIn. R' = R/sqrt(N)
        }
        else
        {
            Rinvs.at(j) = X.block(0,sind,N,gsize(j)).adjoint()*X.block(0,sind,N,gsize(j))/N;
            Rinvs.at(j) = Rinvs.at(j).array().sqrt().inverse();
        }
        X.block(0,sind,N,gsize(j)) = X.block(0,sind,N,gsize(j))*Rinvs[j];
    }
}

template <class TX>
void pgGroupLassoFit<TX>::destandardizeX()
{
    for (int j=1;j<J;j++){
        X.block(0,grpSIdx(j),N,gsize(j))*=Rinvs[j].inverse();
        X.block(0,grpSIdx(j),N,gsize(j)).rowwise()+= Xcenter.segment(grpSIdx(j), gsize(j)).transpose();
    }
    
}

template <class TX>
void pgGroupLassoFit<TX>::standardizeX()
{
    for (int j=1;j<J;j++){
        X.block(0,grpSIdx(j),N,gsize(j)).rowwise()-= Xcenter.segment(grpSIdx(j), gsize(j)).transpose();
        X.block(0,grpSIdx(j),N,gsize(j))*=Rinvs[j];
    }
}
template <class TX>
double pgGroupLassoFit<TX>::evalObjective(const VectorXd & beta, const ArrayXd & lambdaj)
{
    const double c = std::log(nl/(pi*nu));
    double l12norm(0);
    VectorXd lpred, pred, logExpLpred, logExpPred,obslogL;
    VectorXd bj;
    ArrayXi ridx(N); for(int i=0;i<N;i++){ridx(i)=i;}
    
    lpred = linpred(true,beta,ridx);
    logExpLpred = (lpred.array().exp().array()+1).array().log();
    pred = c+lpred.array()-logExpLpred.array();
    logExpPred = (1+pred.array().exp()).array().log();
    obslogL = (y.array()*pred.array()-logExpPred.array());//response z_lu = y
    
    for (int j=0;j<J;j++){
        bj=beta.segment(grpSIdx(j)+1,gsize(j));
        l12norm+=lambdaj(j)*bj.lpNorm<2>();
    }
    return -obslogL.sum()+N*l12norm;
}
template <class TX>
VectorXd pgGroupLassoFit<TX>::gradient(const VectorXd & beta, const ArrayXi & ridx){
    VectorXd eta, h_eta;
    VectorXd exp_eta;//exp(eta)
    VectorXd p_h; // 1/(1+exp(-h))
    VectorXd p1_eta; //1/(1+exp(eta))
    VectorXd gr(ridx.size()); gr.setZero();
    
    eta = linpred(true,beta,ridx);
//    cout<<"eta:"<<eta.block(0,0,1,1)<<endl;
    exp_eta= eta.array().exp();

    h_eta = std::log(nl/(pi*nu))+eta.array()-(exp_eta.array()+1).array().log().array();//bias+eta-log(1+exp(eta))
    p_h = (-h_eta).array().exp();
    p_h = 1/(1+p_h.array());
    p1_eta = 1/(1+exp_eta.array());
//    cout<<"Gradient: \n p_h:\n"<<p_h<<",\n p1_eta:\n"<<p1_eta<<endl;
    for(int i=0; i<ridx.size();i++){
//        cout<<"y("<<ridx(i)<<"):"<<y(ridx(i))<<endl;
        gr(i)= -(y(ridx(i))-p_h(i))*p1_eta(i);
    }
    return gr;
}

template <class TX>
VectorXd pgGroupLassoFit<TX>::subgradient(VectorXd & gradient, VectorXd & beta, ArrayXd lambdaj){
    double gjnorm, bjnorm;
    VectorXd subgrad(gradient);
    subgrad.setConstant(1);
    subgrad(0) = gradient(0);
    for (int j=1;j<J;j++){
        Map<VectorXd> gj(&gradient.coeffRef(grpSIdx(j)+1),gsize(j));
        Map<VectorXd> bj(&beta.coeffRef(grpSIdx(j)+1),gsize(j));
        bjnorm=bj.norm();
        gjnorm=gj.norm();
        if(bjnorm>0){
            if(lambdaj(j)>1e-10&&gjnorm>1e-10){
                subgrad.segment(grpSIdx(j)+1,gsize(j))=(1-(lambdaj(j)/gjnorm))*gj;
            }else{
                subgrad.segment(grpSIdx(j)+1,gsize(j))=gj;
            }
        }else{
            subgrad.segment(grpSIdx(j)+1,gsize(j))= ((gjnorm>lambdaj(j))?(1-(lambdaj(j)/gjnorm)):0)*gj;
        }
    }
    return subgrad;
}




template <class TX>
VectorXd pgGroupLassoFit<TX>::q(int i){
    VectorXd q1i(p), qi(X.row(i));
    q1i<<1,qi;
//    cout<<"q["<<i<<",]:\n"<<q1i<<endl;
    return q1i;
}

template <>
VectorXd pgGroupLassoFit<SparseMatrix<double> >::q(int i){
    MatrixXd qi(X.row(i));
    VectorXd q1i(p);
//    cout<<"q["<<i<<",]:\n"<<qi<<endl;
    for (int j=1;j<J;j++){
        qi.block(0,grpSIdx(j),1,gsize(j)).rowwise()-= Xcenter.segment(grpSIdx(j), gsize(j)).transpose();
        qi.block(0,grpSIdx(j),1,gsize(j))*=Rinvs[j];
    }
    q1i(0) = 1; q1i.segment(1, p-1) =qi.transpose();
    return q1i;
}

////////////////////////////////////////////////////////
//MatrixXd Specialization
///////////////////////////////////////////////////////



///////////////////////////////////////////////////////
//Specialization Sparse
///////////////////////////////////////////////////////
template <>
void pgGroupLassoFit<SparseMatrix<double> >::Rinvs_X()
{
    MatrixXd Xcentered;
    MatrixXd Xdl;
    for (int l=0;l<(p-1);++l)
    {
        Xdl = X.col(l);
        Xcenter(l) = Xdl.mean();
    }
    
    for(int j=1;j<J;++j)
    {
        int sind = grpSIdx(j);
        Xcentered = X.block(0,sind,N,gsize(j));
        
        int k(0);
        for(int l=sind; l<(sind+gsize(j)); ++l)
        {
            Xdl = X.col(l);
            k = l-sind;
            Xcentered.col(k) = Xdl.array()-Xcenter(l);
        }
        if(gsize(j)>1)
        {
            //Do QR decomposition
            ColPivHouseholderQR<MatrixXd> qr(Xcentered);
            
            if(qr.rank() < gsize(j)){throw std::invalid_argument("X(j) does not have full column rank");}
            
            MatrixXd R = qr.matrixR().topLeftCorner(qr.rank(), qr.rank()).triangularView<Upper>();
            MatrixXd P =qr.colsPermutation();
            R=R*P.inverse()/std::sqrt(N);
            
            Rinvs.at(j)= R.inverse();// QtQ = NIn. R' = R/sqrt(N)
        }
        else
        {
            Rinvs.at(j) = Xcentered.adjoint()*Xcentered/N;
            Rinvs.at(j) = Rinvs.at(j).array().sqrt().inverse();
        }
        
    }
    
}
template <>
void pgGroupLassoFit<SparseMatrix<double> >::destandardizeX()
{
    //do nothing
}

template <>
void pgGroupLassoFit<SparseMatrix<double> >::standardizeX()
{
    //do nothing
}


////////////////////////////////////////////////////////////////////////////////
template<class TX>
void pgGroupLassoFit<TX>::checkDesignMatrix(const TX & X)
{
    for(int j=0;j<X.cols();j++)
    {
        if((X.col(j).array()==0).all()){throw std::invalid_argument("each column should have at least one non-zero element");}
    }
}


template<>
void pgGroupLassoFit<SparseMatrix<double> >::checkDesignMatrix(const SparseMatrix<double> & X)
{
    for(int j=0;j<X.cols();j++)
    {
        if(X.col(j).nonZeros()==0){throw std::invalid_argument("each column should have at least one non-zero element");}
    }
}

//Explicit Instantiation
template class pgGroupLassoFit<MatrixXd>;
template class pgGroupLassoFit<SparseMatrix<double> >;
template class pgGroupLassoFit<Map<MatrixXd> >;
