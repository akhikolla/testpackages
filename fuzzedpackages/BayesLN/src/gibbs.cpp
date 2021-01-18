
extern "C"
{
#include "gig_generator.h"
}
#include <RcppEigen.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
void rmn_mu_S2(Eigen::VectorXd &rnd_vec,const int p,
               const Eigen::VectorXd mu,
               const Eigen::MatrixXd S){
  Eigen::VectorXd z=Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(p));
  rnd_vec=mu+S.llt().matrixL()*z;
}

// [[Rcpp::export]]
Rcpp::List MCMC_alg(const Eigen::VectorXd y,
                            const Eigen::MatrixXd X,
                            const List Z_list,
                            const List K_gamma_list,
                            const Eigen::MatrixXd S_beta_pri,
                            double l_s,Eigen::VectorXd l_t,
                            double d_s,Eigen::VectorXd d_t,
                            double g_s,Eigen::VectorXd g_t,
                            int s,
                            int nsamp, int verbose,
                            Eigen::VectorXd beta_init,
                            double sigma2_init,
                            Eigen::VectorXd tau2_init){

  //Fixed quantites
  const Eigen::MatrixXd
  XtX=X.transpose()*X, Xt=X.transpose(),
    Q_beta_pri=S_beta_pri.inverse();
  const int
    p=X.cols(), n=X.rows();

  //Creation arrays of dimension s
  std::vector<Eigen::SparseMatrix<double> >
    Z(s), ZtZ(s), K_gamma(s);
  std::vector<Eigen::VectorXd>
    gamma_init(s);
  std::vector<int>
    K(s);

  //Initializations arrays of dimension s
  for(int j=0;j<s;j++){
    Z[j]=Rcpp::as<Eigen::SparseMatrix<double> >(Z_list[j]);
    ZtZ[j]=Z[j].transpose()*Z[j];
    K[j]=Z[j].cols();
    K_gamma[j]=Rcpp::as<Eigen::SparseMatrix<double> >(K_gamma_list[j]);
    gamma_init[j]=Eigen::VectorXd::Zero(K[j]);
  }

  //Final outputs - storage
  Eigen::VectorXd
    s2eps_mc(nsamp);
  Eigen::MatrixXd
    beta_mc(nsamp,p),s2gamma_mc(nsamp,s);
  std::vector<Eigen::MatrixXd>
    gamma_mc(s);
  for(int j=0; j<s; j++){
    gamma_mc[j]=Eigen::MatrixXd::Zero(nsamp,K[j]);
  }

  //Current values of the parameters and functions
  double
    s2eps_curr=sigma2_init;
  Eigen::VectorXd
    Zgamma=Eigen::VectorXd::Zero(n), beta_curr=beta_init,
      s2gamma_curr=tau2_init;
  std::vector<Eigen::VectorXd>
    gamma_curr(s);
  for(int j=0; j<s; j++){
    gamma_curr[j]=gamma_init[j];
    Zgamma+=Z[j]*gamma_curr[j];
  }

  // Auxiliary quantites
  Eigen::VectorXd
    res_beta(n),res_eps(n), mu_beta_fc_curr(p), Xbeta_curr(n);
  Eigen::MatrixXd
    Q_beta_fc_curr(p,p),S_beta_fc_curr(p,p);
  std::vector<Eigen::SparseMatrix<double> >
    Q_gamma_fc_curr(s);
  std::vector<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower> >
    Chol_Q_gamma_fc_curr(s);
  std::vector<Eigen::VectorXd>
    e(s),res_gamma(s), c(s), rnd_z_gamma(s),
    b_gamma_fc_curr(s),mu_gamma_fc_curr(s);
  std::vector<int>
    nc(s);

  // Initialization auxiliary quantites
  for(int j=0; j<s; j++){
    res_gamma[j]=Eigen::VectorXd::Zero(n);
    b_gamma_fc_curr[j]=Eigen::VectorXd::Zero(K[j]);
    mu_gamma_fc_curr[j]=Eigen::VectorXd::Zero(K[j]);
    rnd_z_gamma[j]=Eigen::VectorXd::Zero(K[j]);
    Q_gamma_fc_curr[j]=(ZtZ[j]/s2eps_curr+K_gamma[j]/s2gamma_curr[j]);
    Chol_Q_gamma_fc_curr[j].analyzePattern(Q_gamma_fc_curr[j]);
    Chol_Q_gamma_fc_curr[j].factorize(Q_gamma_fc_curr[j]);
  }

  // counter storage
  int it=0;
  //time
  time_t now;

  //start sampling
  for(int k=0;k<nsamp;k++){
    /////////////////
    //Sampling beta
    res_beta=y-Zgamma;
    //Par full conditional
    Q_beta_fc_curr=(XtX/s2eps_curr+Q_beta_pri);
    S_beta_fc_curr=Q_beta_fc_curr.inverse();
    mu_beta_fc_curr=S_beta_fc_curr*((1.0/s2eps_curr)*Xt*res_beta);
    //Generate

    //beta_curr=rmn_mu_S(mu_beta_fc_curr,S_beta_fc_curr,p);
    rmn_mu_S2(beta_curr,p,mu_beta_fc_curr,S_beta_fc_curr);
    //update quantities
    Xbeta_curr=X*beta_curr;

    /////////////////
    //Sampling gamma
    Zgamma=Zgamma*0.0;
    for(int j=0; j<s; j++){
      res_gamma[j]=res_gamma[j]*0.0;
      res_gamma[j]=y-Xbeta_curr;
      for(int i=0; i<s; i++){
        if(i!=j){
          res_gamma[j]-=Z[i]*gamma_curr[i];
        }
      }
      // Full conditionals parameters
      Q_gamma_fc_curr[j]=(ZtZ[j]/s2eps_curr+K_gamma[j]/s2gamma_curr[j]);
      Chol_Q_gamma_fc_curr[j].factorize(Q_gamma_fc_curr[j]);
      b_gamma_fc_curr[j]=(1.0/s2eps_curr)*Z[j].transpose()*res_gamma[j];
      mu_gamma_fc_curr[j]=Chol_Q_gamma_fc_curr[j].solve(b_gamma_fc_curr[j]);
      //Generate
      rnd_z_gamma[j]=Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(K[j]));
      gamma_curr[j]=Chol_Q_gamma_fc_curr[j].permutationPinv() *
        Chol_Q_gamma_fc_curr[j].matrixU().solve(rnd_z_gamma[j]);
      gamma_curr[j]+=mu_gamma_fc_curr[j];
      //Update quantity
      Zgamma+=Z[j]*gamma_curr[j];
    }

    /////////////////
    //Sampling s2eps
    res_eps=y-Xbeta_curr-Zgamma;
    s2eps_curr=rgig(l_s - n / 2.0, d_s * d_s + res_eps.squaredNorm(), g_s * g_s);

    ///////////////////
    //Sampling s2gamma
    for(int j=0; j<s; j++){

      s2gamma_curr[j]= rgig(l_t[j] - K[j] / 2.0, d_t[j] *d_t[j] + gamma_curr[j].squaredNorm(), g_t[j] * g_t[j]);

    }


    s2eps_mc(k)= s2eps_curr;
    s2gamma_mc.row(k)= s2gamma_curr;
    beta_mc.row(k)=beta_curr;
    for(int j=0;j<s;j++){
      gamma_mc[j].row(k)=gamma_curr[j];
    }


    //Progress bar
    double Perc=100*(k+1)/(nsamp);
    if(((k+1)%(nsamp/100))==0 && (verbose==1)){
      Rprintf("-");
      R_FlushConsole();
      //R_ProcessEvents(); for windows
      if(((k+1)%(nsamp/10))==0){
        Rprintf(":%.1f%%\n",Perc);
        R_FlushConsole();
        //R_ProcessEvents(); for windows
      }
    }
    if (k % 1000 == 0){
      Rcpp::checkUserInterrupt();
    }

  }





  return List::create(Named("beta") = beta_mc,
                      Named("u") = gamma_mc, Named("tau2") = s2gamma_mc, Named("sigma2") = s2eps_mc);

}



// [[Rcpp::export]]
Eigen::MatrixXd post_pred(List output, Eigen::MatrixXd Xrep, List Zrep_list, int s, int nsamp){
  int n_rep = Xrep.rows();
  Eigen::VectorXd sigma2_mc=Rcpp::as<Eigen::VectorXd>(output["sigma2"]);
  Eigen::MatrixXd beta_mc=Rcpp::as<Eigen::MatrixXd>(output["beta"]),
    tau2_mc=Rcpp::as<Eigen::MatrixXd>(output["tau2"]);
  std::vector<Eigen::MatrixXd> gamma_mc= Rcpp::as<std::vector<Eigen::MatrixXd> >(output["u"]);
  std::vector<Eigen::SparseMatrix<double> > Zrep(s);

  Eigen::MatrixXd Zgamma_mc = Eigen::MatrixXd::Zero(n_rep, nsamp);
  for(int j=0;j<s;j++){
    Zrep[j]=Rcpp::as<Eigen::SparseMatrix<double> >(Zrep_list[j]);
    Zgamma_mc += Zrep[j] * gamma_mc[j].transpose();
  }
  Eigen::MatrixXd Xbeta_mc = Xrep * beta_mc.transpose();
  Eigen::MatrixXd y_rep(nsamp, n_rep);
  for(int i=0; i<nsamp; i++){
    y_rep.row(i) = Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(n_rep))* pow(sigma2_mc[i], 0.5);
    y_rep.row(i) += Xbeta_mc.col(i);
    y_rep.row(i) += Zgamma_mc.col(i);
  }


  return(y_rep);
}



