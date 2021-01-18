// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//' @importFrom Rcpp sourceCpp
//' @useDynLib BoltzMM

#include "RcppArmadillo.h"

#define NDEBUG
#define BOOST_DISABLE_ASSERTS
#include <boost/assert.hpp>
#include <boost/dynamic_bitset.hpp>

//for returning a vector within a list
Rcpp::NumericVector export_vec(arma::vec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}


//creating binary vectors
arma::vec bin_vec(int y,  int n)
{
  const boost::dynamic_bitset<> b(n, y);
  arma::vec x = arma::zeros(n);
  for(int i=0; i<n; i++){
    x(i) = 2*(b[i]-0.5) ;
  }
  return(x);
}

//'Probability mass function of a fully-visible Boltzmann machine evaluated for an individual vector.
//'@description Compute the probability of a string of n>1 binary spin variables (i.e. each element is -1 or 1) arising from a fully-visible Boltzmann machine with some specified bias vector and interaction matrix.
//'@param xval Vector of length n containing binary spin variables.
//'@param bvec Vector of length n containing real valued bias parameters.
//'@param Mmat Symmetric n by n matrix, with zeros along the diagonal, containing the interaction parameters.
//'@return The probability of the random string \code{xval} under a fully-visible Boltzmann machine with bias vector \code{bvec} and interaction matrix \code{Mmat}.
//'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
//'@author Andrew T. Jones and Hien D. Nguyen
//'@examples # Compute the probability of the vector xval=(-1,1,-1), under bvec and Mmat.
//'xval <- c(-1,1,-1)
//'bvec <- c(0,0.5,0.25)
//'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
//'pfvbm(xval,bvec,Mmat)
//'@export
// [[Rcpp::export]]
double pfvbm(arma::vec xval, arma::vec bvec, arma::mat Mmat) {
    int n = bvec.n_elem;
    double prob = 0.0;
    double norm = 0.0;
    int count = std::pow(2,n);

    if(xval.n_elem!= n || Mmat.n_rows!=n || Mmat.n_cols!=n || Mmat.n_rows!=Mmat.n_cols ){
      Rcpp::Rcerr << "Input variable dimensions do not match";
    }else{

        for(int i=0; i<count; i++){
            arma::vec zeta_i = bin_vec(i,n);
            norm += as_scalar(arma::exp(0.5*zeta_i.t()*Mmat*zeta_i+arma::dot(bvec,zeta_i)));
        }

        prob = as_scalar(arma::exp(0.5*xval.t()*Mmat*xval+ arma::dot(bvec,xval)));
        prob /= norm;
    }
    return prob;
}

//'Probability mass function of a fully-visible Boltzmann machine evaluated for all possible vectors.
//'@description Compute the probability of all 2^n strings of n>1 binary spin variables (i.e. each element is -1 or 1) arising from a fully-visible Boltzmann machine with some specified bias vector and interaction matrix.
//'@param bvec Vector of length n containing real valued bias parameters.
//'@param Mmat Symmetric n by n matrix, with zeros along the diagonal, containing the interaction parameters.
//'@return A vector of the probabilities of all 2^n binary spin vectors under a fully-visible Boltzmann machine with bias vector \code{bvec} and interaction matrix \code{Mmat}. Probabilities are reported in ascending order of the binary strings; i.e for n=2 the reporting order is (-1,1), (-1,1), (1,-1), and (1,1).
//'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
//'@author Andrew T. Jones and Hien D. Nguyen
//'@examples # Compute the probability of every length n=3 binary spin vector under bvec and Mmat.
//'bvec <- c(0,0.5,0.25)
//'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
//'allpfvbm(bvec,Mmat)
//'@export
// [[Rcpp::export]]
arma::rowvec allpfvbm(arma::vec bvec, arma::mat Mmat) {
    int n = bvec.n_elem;
    double norm = 0.0;
    int count = std::pow(2,n);
    arma::rowvec probvec = arma::zeros(count).t();

    if(Mmat.n_rows!=n || Mmat.n_cols!=n || Mmat.n_rows!=Mmat.n_cols ){
      Rcpp::Rcerr << "Input variable dimensions do not match";
    }else{
        for(int i=0; i<count; i++){
          arma::vec zeta_i = bin_vec(i,n);
          double prob = as_scalar(arma::exp(0.5*zeta_i.t()*Mmat*zeta_i+arma::dot(bvec, zeta_i)));
          probvec(i) = prob;
          norm +=  prob;
        }

        probvec /= norm;
    }
    return(probvec);
}

//'Random data generation from a fully-visible Boltzmann machine.
//'@description Generate N random strings of n>1 binary spin variables (i.e. each element is -1 or 1) arising from a fully-visible Boltzmann machine with some specified bias vector and interaction matrix.
//'@param num Number N of random strings to be generated.
//'@param bvec Vector of length n containing real valued bias parameters.
//'@param Mmat Symmetric n by n matrix, with zeros along the diagonal, containing the interaction parameters.
//'@return An N by n matrix, where each row contains a random spin variable string from a fully-visible Boltzmann machine with bias vector \code{bvec} and interaction matrix \code{Mmat}.
//'@note The function \code{allpfvbm} must be called each time this function is run. Thus, it is much more efficient to generate N strings all at once, than to generate strings one at a time.
//'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
//'@author Andrew T. Jones and Hien D. Nguyen
//'@examples # Generate num=10 random strings of n=3 binary spin variables under bvec and Mmat.
//'num <- 10
//'bvec <- c(0,0.5,0.25)
//'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
//'rfvbm(num,bvec,Mmat)
//'@export
// [[Rcpp::export]]
arma::mat rfvbm(int num, arma::vec bvec, arma::mat Mmat) {
  int n = bvec.n_elem;
  arma::mat returnmat   = arma::zeros(num,n);

  if(Mmat.n_rows!=n || Mmat.n_cols!=n || Mmat.n_rows!=Mmat.n_cols ){
    Rcpp::Rcerr << "Input variable dimensions do not match";
  }else{

      arma::rowvec cumprob = arma::cumsum(allpfvbm(bvec,Mmat));
      //need to fix this, cant call system RNG if CRAN
      arma::vec random_nums = Rcpp::runif(num);
      for(int i=0; i<num; i++){
        int j = arma::as_scalar(find(cumprob > random_nums(i), 1, "first"));
        arma::vec zeta_j = bin_vec(j,n);
        returnmat.row(i) = zeta_j.t();
      }
  }
  return(returnmat);
}


//'Maximum pseudolikelihood estimation of a fully-visible Boltzmann machine.
//'@description Estimates the bias vector and interaction matrix of a fully-visible Boltzmann machine via maximum pseudolikelihood estimation using an MM algorithm.
//'@param data An N by n matrix, where each of the N rows contains a length n string of spin variables  (i.e. each element is -1 or 1).
//'@param bvec Initial estimate for a vector of length n containing real valued bias parameters.
//'@param Mmat Initial estimate for a symmetric n by n matrix, with zeros along the diagonal, containing the interaction parameters.
//'@param delta_crit Real threshold value for the convergence criterion, based on the relative change in the Euclidean distance of parameter estimates from consecutive iterations.
//'@param max_it Integer value indicating the maximum number of iterations that the algorithm is to run for.
//'@return A list containing 4 objects: the final log-pseudolikelihood value \code{pll}, a vector containing the estimate of the bias parameters \code{bvec}, a matrix containing the estimate of the interaction parameters \code{Mmat}, and the number of algorithm iterations \code{itt}.
//'@references H.D. Nguyen and I.A. Wood (2016), A block successive lower-bound maximization algorithm for the maximum pseudolikelihood estimation of fully visible Boltzmann machines, Neural Computation, vol 28, pp. 485-492
//'@author Andrew T. Jones and Hien D. Nguyen
//'@examples # Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.
//'num <- 1000
//'bvec <- c(0,0.5,0.25)
//'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
//'data <- rfvbm(num,bvec,Mmat)
//'# Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.
//'fitfvbm(data,bvec,Mmat)
//'@export
// [[Rcpp::export]]
Rcpp::List fitfvbm(arma::mat data, arma::vec bvec, arma::mat Mmat, double delta_crit = 0.001, int max_it = 1000){
    //New parameters transfer into old parameters
    int N = data.n_rows;
    int D = bvec.n_elem;
    int itt = 0;

    if(Mmat.n_rows!=D || Mmat.n_cols!=D || data.n_cols !=D || data.n_cols!=Mmat.n_cols || data.n_cols!=Mmat.n_rows || Mmat.n_rows!=Mmat.n_cols){
        Rcpp::Rcerr << "Input variable dimensions do not match";
        return(Rcpp::List::create());
    }

    if(arma::any(arma::any(arma::abs(data)!=1))){
      Rcpp::Rcerr << "Input data must consist of only 1 or -1 values.";
      return(Rcpp::List::create());
    }

    arma::mat MM = Mmat;
    arma::mat dataj = data;
    arma::mat datak = data;
    arma::vec temp = arma::zeros(N);
    double delta = delta_crit+10.0;
    arma::mat par  = arma::join_rows(bvec,MM);
    arma::mat old_par = par;

    double DERIV = 0.0;
    double LIKE = 0.0;
    double sumj = 0.0;
    double sumk = 0.0;

    while ((delta > delta_crit) && (itt<max_it))
    {
        itt++;
        old_par = par;

        for(int j=0; j<D; j++)
        {
            DERIV = arma::sum(data.col(j));
            dataj = data.each_row()%MM.col(j).t();
            DERIV -= arma::as_scalar(arma::sum(arma::tanh(sum(dataj,1) + bvec(j))));
            bvec(j) +=  DERIV/N;
        }



        for(int j=0; j<D; j++)
        {
            for(int k=(j+1); k<D; k++)
            {
                DERIV = 2.0*arma::dot(data.col(j),data.col(k));

                dataj = data.each_row()%MM.col(j).t();
                datak = data.each_row()%MM.col(k).t();

                sumk = arma::as_scalar(arma::dot(data.col(k),(arma::tanh(sum(dataj,1) + bvec(j)))));
                sumj = arma::as_scalar(arma::dot(data.col(j),(arma::tanh(sum(datak,1) + bvec(k)))));

                DERIV-=(sumk+sumj);

                MM(j,k) += DERIV/(2.0*N);
                MM(k,j) = MM(j,k);
            }
        }

        par = arma::join_rows(bvec,MM);
        delta = std::sqrt(arma::accu(arma::pow((par-old_par),2.0)))/std::max(std::sqrt(std::pow(arma::accu(old_par),2.0)),1.0);
    }

    LIKE = 0.0;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<D; j++)
        {
            LIKE += data(i,j)*arma::dot(MM.col(j),data.row(i)) + bvec(j)*data(i,j) - std::log(std::cosh(arma::dot(MM.col(j),data.row(i))+bvec(j))) - std::log(2.0);
        }
    }

    Rcpp::List retList = Rcpp::List::create(
        Rcpp::Named("pll")= LIKE,
        Rcpp::Named("bvec")= export_vec(bvec),
        Rcpp::Named("Mmat")= MM,
        Rcpp::Named("itt")= itt
    );

    return(retList);
}

//'Partial derivatives of the log-pseudolikelihood function for a fitted fully-visible Boltzmann machine.
//'@description Computes the partial derivatives for all unique parameter elements of the bias vector and interaction matrix of a fully-visible Boltzmann machine, for some random length n string of spin variables (i.e. each element is -1 or 1) and some fitted parameter values.
//'@param data Vector of length n containing binary spin variables.
//'@param model List generated from \code{fitfvbm}.
//'@return A list containing 2 objects: a vector containing the partial derivatives corresponding to the bias parameters \code{bvec}, and a matrix containing the partial derivatives corresponding to the interaction parameters \code{Mmat}.
//'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
//'@author Andrew T. Jones and Hien D. Nguyen
//'@examples # Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.
//'num <- 1000
//'bvec <- c(0,0.5,0.25)
//'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
//'data <- rfvbm(num,bvec,Mmat)
//'# Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.
//'model <- fitfvbm(data,bvec,Mmat)
//'# Compute the partial derivatives evaluated at the first observation of data.
//'fvbmpartiald(data,model)
//'@export
//[[Rcpp::export]]
Rcpp::List fvbmpartiald(arma::mat data, Rcpp::List model){
    //int N = data.n_rows;
    int D = data.n_cols;

    arma::mat Mmat=Rcpp::as<arma::mat>(model(2));
    arma::vec bvec=Rcpp::as<arma::vec>(model(1));


    arma::mat partiald = arma::zeros(D,D+1);

    for(int j = 0; j<D; j++) {

      arma::mat dataj = data.each_row()%Mmat.col(j).t();
      arma::vec tempj = arma::tanh(sum(dataj,1) + bvec(j));
      double sumj = arma::as_scalar(arma::sum(tempj));
      partiald(j,0) = arma::sum(data.col(j)) - sumj;
      arma::mat temp2 = data.col(j)-tempj;
      partiald(j,arma::span(1,D)) = arma::sum(temp2.t()*data,0);

    }

    arma::vec bvecpartial = arma::zeros(D);
    arma::mat Mmatpartial = arma::zeros(D,D);


    bvecpartial = partiald.col(0);
    Mmatpartial = partiald(arma::span(0,D-1),arma::span(1,D));
    Mmatpartial = Mmatpartial + Mmatpartial.t();
    Mmatpartial.diag().zeros();


    Rcpp::List retList = Rcpp::List::create(
        Rcpp::Named("bvec")= export_vec(bvecpartial),
        Rcpp::Named("Mmat")= Mmatpartial
    );

    return(retList);

}




//'Sandwich estimator of the covariance matrix for a fitted fully-visible Boltzmann machine.
//'@description Computes the sandwich estimator of the covariance matrix for a maximum pseudolikelihood estimated fully-visible Boltzmann machine.
//'@param data An N by n matrix, where each of the N rows contains a length n string of spin variables  (i.e. each element is -1 or 1).
//'@param model List generated from \code{fitfvbm}.
//'@param fvbmHess A function that computes the Hessian of the parameter elements. Currently, the only implemented method is the default \code{fvbmHess} function.
//'@return The n+choose(n,2) by n+choose(n,2) sandwich covariance matrix, estimated using \code{data} and evaluated at the fitted parameter values provided in \code{model}. Each row (column) is a unique element of the bias vector and interaction matrix. The rows are arranged in lexicographical order with the bias elements first, followed by the interaction elements. For example, if n=3, the order would be bias[1], bias[2] bias[3], interaction[1,2], interaction[1,3], and interaction[2,3].
//'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
//'@author Andrew T. Jones and Hien D. Nguyen
//'@examples # Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.
//'num <- 1000
//'bvec <- c(0,0.5,0.25)
//'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
//'data <- rfvbm(num,bvec,Mmat)
//'# Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.
//'model <- fitfvbm(data,bvec,Mmat)
//'# Compute the sandwich covariance matrix using the data and the model.
//'fvbmcov(data,model,fvbmHess)
//'@export
//[[Rcpp::export]]
arma::mat fvbmcov(arma::mat data, Rcpp::List model, Rcpp::Function fvbmHess){
    int N = data.n_rows;
    int D = data.n_cols;

    arma::mat Mmat=Rcpp::as<arma::mat>(model(2));
    int hessDim = D+D*(D-1)/2;
    arma::mat I_1 = -(1.0/N)*Rcpp::as<arma::mat>(fvbmHess(data,model));
    arma::mat I_2 = arma::zeros(hessDim,hessDim);

    //remove N loop? need to vectorise fvbmpartiald
    for(int i =0; i< N; i++) {
      Rcpp::List Partial_res = fvbmpartiald(arma::mat(data.row(i)),model);
      arma::vec Extract = arma::zeros(hessDim);
      Extract(arma::span(0,D-1)) =  Rcpp::as<arma::vec>(Partial_res(0));
      Extract(arma::span(D,hessDim-1)) = arma::nonzeros(arma::trimatl(Rcpp::as<arma::mat>(Partial_res(1)), -1));
      I_2 += Extract*Extract.t();
     }

    I_2 = (1.0/N)*I_2;
    arma::mat I_1_s = arma::pinv(I_1);
    arma::mat Covar = I_1_s*I_2*I_1_s;

    return(Covar);

}


