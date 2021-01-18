#include <RcppEigen.h>
#include <numeric>

//[[Rcpp::depends(RcppEigen)]]

//' @title Geodesic distance
//' @description Evaluate geodesic distance (shortest path) between all pairs of nodes in the network.
//'
//' @param M Input adjacency matrix
//'
//' @return Matrix containing all the pairwise geodesic distances
//' @examples dst(example_adjacency_matrix)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd dst(const Eigen::Map<Eigen::MatrixXd> M){

  int g=M.rows();
  Eigen::MatrixXd  Yr=M;
  Eigen::MatrixXd Dst=M;
  for (int i(0); i<g;i++){
    for (int j(0); j<g;j++){
      if (M(i,j)==1){
        Dst(i,j)=1;
      }
      else {
        Dst(i,j)=g;
      }
    }
  }
  for (int r(2);r<g;r++){
    Yr=Yr*M;
    for (int i(0); i<g;i++){
      for (int j(0); j<g;j++){
        if (Yr(i,j)>0 && Dst(i,j)==g){
          Dst(i,j)=r;
        }
      }
    }
    
  }
  return Dst;
}

 
//' @title Distance between latent positions
//' @description Compute the square root of the Euclidean distances between latent positions and 
//' return them with a negative sign.
//'
//' @param Z Latent positions matrix. The matrix size must be \code{(n,k)}, where \code{n} and \code{k} denote respectively 
//' the number of nodes in the network and the latent space dimensionality.
//' @return Matrix containing the negative square root of the Euclidean distances between latent positions
//' @examples pos = matrix(rnorm(20), ncol=2)
//' lpz_dist(pos)
//' @export
// [[Rcpp::export]] 
Eigen::MatrixXd lpz_dist(Eigen::MatrixXd Z){

  Eigen::MatrixXd ZtZ=Z*(Z.adjoint());
  int k=Z.rows();
  Eigen::MatrixXd temp;
  temp.setOnes(1,k);
  Eigen::MatrixXd temp2=ZtZ.diagonal()*temp;
  Eigen::MatrixXd mg=temp2;
  mg=mg+temp2.adjoint()-ZtZ*2;
  for (int i(0);i<k;i++){
    for (int j(0);j<k;j++){ 
      mg(i,j)=-sqrt(mg(i,j));
    }
  }
  return mg;
}


//' @title Network log-likelihood
//' @description Compute the log-likelihood of the whole observed network based on the
//' latent positions estimates and the model assumptions. See \link[BLSM]{BLSM} for more information. 
//' 
//' @param Y Adjacency matrix of the observed network
//' @param lpz Matrix containing the negative square root of the Euclidean distances between latent positions 
//' (output of \link[BLSM]{lpz_dist})
//' @param alpha Model variable \eqn{\alpha}
//' @param W BLSM Weights matrix of the observed network
//' 
//' @return Log-likelihood of the observed network
// [[Rcpp::export]] 
double lpY (Eigen::MatrixXd Y, Eigen::MatrixXd lpz, double alpha, Eigen::MatrixXd W){

  double val(0);
  double lpg;
  int k=lpz.rows();
  
  for (int i(0);i<k;i++){
    for (int j(0);j<k;j++){
      if (i!=j){
        lpg=W(i,j)*lpz(i,j)+alpha;
        val=val+Y(i,j)*lpg-log(1+exp(lpg));
      }
    }
  }
  return val;
}


//' @title Network (positive) log-likelihood 
//' @description Compute the (positive) log-likelihood of the whole observed network based on the
//' latent positions estimates and the model assumptions. The inputs are slightly different from those of \link[BLSM]{lpY},
//' so the function basically applies some preprocessing before calling \link[BLSM]{lpY} and returning its value with the opposite sign. 
//' 
//' @param avZ Vector containing the \eqn{\alpha} value and the latent positions 
//' @param Y Adjacency matrix of the observed network
//' @param W BLSM Weights matrix of the observed network
//' 
//' @return Log-likelihood of the observed network
// [[Rcpp::export]] 
double mlpY (Eigen::VectorXd avZ, Eigen::MatrixXd Y, Eigen::MatrixXd W){

  int k=Y.rows();
  double val(0);
  double alpha=avZ(0);
  
  Eigen::MatrixXd Z=avZ.tail(avZ.size()-1);
  Z.resize(k,2);
  Eigen::MatrixXd lpz=lpz_dist(Z);
  
  val = lpY (Y, lpz, alpha, W);
  return -val;
}


//' @title lpz_dist optimized for individual updates
//' @description Compute the square root of the Euclidean distances between a specific coordinate in the latent space
//' and all the others. The function follows almost the same approach as \link[BLSM]{lpz_dist}, but it is
//' more suitable for the individual updates occurring during the simulation.
//' 
//' @param Z Latent positions matrix
//' @param node Specific node in the network corresponding to the latent coordinate which will be used as reference
//' @param diag Diagonal from \code{t(Z)\%*\%Z} matrix, passed to speed up the process.
//' @return Vector containing the negative square root of the Euclidean distances between latent positions
// [[Rcpp::export]] 
Eigen::VectorXd lpz_distNODE(Eigen::MatrixXd Z, int node, Eigen::VectorXd diag){

  int k=Z.rows();
  Eigen::VectorXd ZtZ=Z.row(node)*(Z.adjoint());
  Eigen::VectorXd temp(k);
  temp.fill(diag(node));
  Eigen::VectorXd mg;
  mg=temp+diag-ZtZ*2;
  for (int i(0);i<k;i++){
    mg(i)=-sqrt(mg(i));
  }
  return mg;
}


//' @title Network log-likelihood for individual updates
//' @description Compute the log-likelihood of the whole observed network based on the
//' latent positions estimates and the model assumptions. The function follows almost the same approach as \link[BLSM]{lpY}, but it is
//' more suitable for the individual updates occurring during the simulation.
//' @param Y Adjacency matrix of the observed network
//' @param Z Latent positions matrix
//' @param alpha Model variable \eqn{\alpha}
//' @param node Specific node in the network corresponding to the latent coordinate which will be used as reference
//' @param diag Diagonal from \code{t(Z)\%*\%Z} matrix, passed to speed up the process.
//' @param W BLSM Weights matrix of the observed network
//' 
//' @return Log-likelihood of the observed network
// [[Rcpp::export]] 
double lpYNODE (Eigen::MatrixXd Y, Eigen::MatrixXd Z, double alpha, int node, Eigen::VectorXd diag, Eigen::MatrixXd W){
  
  int k=Z.rows();
  double val(0);
  Eigen::VectorXd lpz=lpz_distNODE(Z, node, diag);
  double lpg;  
  for (int i(0);i<k;i++){
    if (i!=node){
    lpg=W(i,node)*lpz(i)+alpha;
      val+=Y(i,node)*lpg-log(1+exp(lpg));
    }
  }
  return val;
}


//' @title Update step for the latent positions
//' @description Accept/reject the proposals for the latent positions
//'  
//' @param Y Adjacency matrix of the observed network
//' @param Z Latent positions matrix
//' @param W BLSM Weights matrix of the observed network
//' @param alpha Model variable \eqn{\alpha}
//' @param zdelta Standard deviation of the Gaussian proposal for latent positions
//' @param mu_z Mean of the Gaussian prior distribution for latent positions 
//' @param sd_z Standard deviation of the Gaussian prior distribution for latent positions
//' 
//' @return Updated latent positions matrix
// [[Rcpp::export]] 
Eigen::MatrixXd Z_up (Eigen::MatrixXd Y, Eigen::MatrixXd Z, Eigen::MatrixXd W, double alpha, double zdelta, double mu_z, double sd_z){
  Rcpp::RNGScope scope;
  int k=Z.rows();
  int h=Z.cols();
  Eigen::MatrixXd Znew = Z;
  double lnew, lold, hr ;
  Rcpp::NumericVector UP(h), D1(h), D2(h), OLD(h);
  Eigen::VectorXd ZtZ=(Z*(Z.adjoint())).diagonal();
  for (int i(0);i<k;i++){
    for (int j(0);j<h;j++){
      OLD[j] = Z(i,j);
      UP[j] = Rcpp::rnorm(1,OLD[j],zdelta)[0];
      Znew(i,j) = UP[j];
    }
    D1 = Rcpp::dnorm(UP, mu_z, sd_z, TRUE);
    D2 = Rcpp::dnorm(OLD, mu_z, sd_z, TRUE);
    lold=lpYNODE(Y,Z,alpha,i, ZtZ, W);
    ZtZ(i) = Znew.row(i)*Znew.row(i).adjoint();
    lnew=lpYNODE(Y,Znew,alpha, i, ZtZ, W);
    hr = 2*(lnew-lold) + std::accumulate(D1.begin(),D1.end(),0) - std::accumulate(D2.begin(),D2.end(),0);
    if (Rcpp::runif(1).at(0)>exp(hr)) {
      Znew.row(i)=Z.row(i);
      ZtZ(i) = Znew.row(i)*Znew.row(i).adjoint();
    }
    else {
      Z.row(i)=Znew.row(i);
    }
  }
  return Znew;
}


//' @title Update step for the \eqn{\alpha} variable
//' @description Accept/reject the proposal for the \eqn{\alpha} model variable 
//'  
//' @param Y Adjacency matrix of the observed network
//' @param lpz Matrix containing the negative square root of the Euclidean distances between latent positions 
//' @param W BLSM Weights matrix of the observed network
//' @param alpha Model variable \eqn{\alpha}
//' @param adelta The uniform proposal for \eqn{\alpha} is defined on the \eqn{[-adelta,+adelta]} interval
//' @param a_a Shape parameter of the Gamma prior distribution for \eqn{\alpha}. The value is usually set to 1, so the prior is actually an exponential distribution.
//' @param a_b Rate parameter of the Gamma prior distribution for \eqn{\alpha}. 
//' 
//' @return Updated value of the \eqn{\alpha} variable
// [[Rcpp::export]] 
double alpha_up (Eigen::MatrixXd Y, Eigen::MatrixXd lpz, Eigen::MatrixXd W, double alpha, double adelta, double a_a, double a_b){
  Rcpp::RNGScope scope;
  double lnew, lold, hr;
  Rcpp::NumericVector alphanew(1),alphaV(1);
  alphanew[0] = std::abs(alpha+Rcpp::runif(1,-adelta,adelta)[0]);
  alphaV[0] = alpha;
  lnew=lpY(Y,lpz,alphanew[0], W);
  lold=lpY(Y,lpz,alpha, W);
  hr=lnew-lold + Rcpp::dgamma(alphanew,a_a,a_b,TRUE)[0] - Rcpp::dgamma(alphaV,a_a,a_b,TRUE)[0];
    
  if(Rcpp::runif(1).at(0)>exp(hr)){ 
    alphanew[0]=alpha;
  }
    
  return alphanew[0];
}
