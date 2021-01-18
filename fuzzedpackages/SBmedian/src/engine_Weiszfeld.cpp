#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

// auxiliary functions
// 1. maximum of three values
double triple_max(double a, double b, double c){
  arma::vec tmp(3,fill::zeros);
  tmp(0) = a;  tmp(1) = b;  tmp(2) = c;
  double output = tmp.max();
  return(output);
}
// 2. minimum of three values
double triple_min(double a, double b, double c){
  arma::vec tmp(3,fill::zeros);
  tmp(0) = a;  tmp(1) = b;  tmp(2) = c;
  double output = tmp.min();
  return(output);
}
// 3. log RBF kernel matrix 
//    x (m x m), y (n x n), output (m x n)
arma::mat compute_kernel (double sigma, arma::mat x, arma::mat y) {
  unsigned int xsize = x.n_rows;
  unsigned int ysize = y.n_rows;
  
  // definition
  arma::mat xprod(xsize, ysize, fill::zeros);
  arma::mat yprod(ysize, xsize, fill::zeros);
  arma::mat xyprod(xsize, ysize, fill::zeros);
  arma::mat kernmat(xsize, ysize, fill::zeros);
  
  // compute xprod
  arma::vec vxprod = arma::diagvec(x*arma::trans(x));
  for (int i=0;i<ysize;i++){
    xprod.col(i) = vxprod;
  }
  // compute yprod
  arma::vec vyprod = arma::diagvec(y*arma::trans(y));
  for (int i=0;i<xsize;i++){
    yprod.col(i) = vyprod;
  }
  // compute xyprod and kernel matrix
  xyprod     = (x*arma::trans(y));
  kernmat    = (2.0*sigma*xyprod)-(sigma*xprod)-(sigma*arma::trans(yprod));
  
  // return kernel matrix
  return(kernmat);
}
// 4. squared distance of subset posteriors from M-posterior
arma::colvec compute_wweights (double sigma, arma::field<arma::mat> subsetAtomsList, arma::mat medianAtoms, arma::field<arma::mat> subsetProbsList, arma::mat medianProbs) {
  unsigned int ii;
  unsigned int nsubset = subsetAtomsList.n_elem; // # subset posteriors
  unsigned int natom = medianAtoms.n_rows;       // # M-posterior atoms
  arma::colvec norms(nsubset, fill::zeros);      // sq-distance between M-posterior and subset posteriors
  const double SMALL = 1e-15;

  double omax = 0.0;  
  
  // matrices for storing 1) distances between M-posterior and subset posterior atoms and 2) probabilities of M-posterior and subset posterior atoms
  mat kernmat1, kernmat2(natom, natom), kernmat12, wtmat1, wtmat2(natom, natom), wtmat12, subWt1, medNorm2, subMed12, subsetAtoms, subsetProbs, subsetProbsMat; 
  
  kernmat2 = compute_kernel(sigma, medianAtoms, medianAtoms); 
  wtmat2   = log(medianProbs * arma::trans(medianProbs));
  medNorm2 = kernmat2 + wtmat2;
  
  // one round of Weiszfeld updates
  for (ii = 0; ii < nsubset; ++ii) {    
    subsetAtoms = subsetAtomsList(ii);
    subsetProbs = subsetProbsList(ii);
    
    kernmat1  = compute_kernel(sigma, subsetAtoms, subsetAtoms);
    wtmat1    = log(subsetProbs * arma::trans(subsetProbs));
    
    kernmat12 = compute_kernel(sigma, subsetAtoms, medianAtoms);
    wtmat12   = log(subsetProbs * arma::trans(medianProbs));
    
    subWt1   = kernmat1 + wtmat1;
    subMed12 = kernmat12 + wtmat12;
    
    omax = triple_max(subWt1.max(), medNorm2.max(), subMed12.max());
    norms(ii) = std::exp(static_cast<float>(omax)) * (accu(exp(subWt1 - omax)) + arma::accu(exp(medNorm2 - omax)) - 2.0 * arma::accu(exp(subMed12 - omax)));
    if (norms(ii) < SMALL) {
      // Rcpp::Rcout << "norm " << (ii + 1) << " is small; truncating to 1e-15." << std::endl;      
      norms(ii) = SMALL;
    }
    
    subsetAtoms.reset(); 
    subsetProbs.reset(); 
    subsetProbsMat.reset();
    kernmat1.reset();  wtmat1.reset();
    kernmat12.reset(); wtmat12.reset();
  }
  
  return(norms); 
}

// Main Computation
// [[Rcpp::export]]
Rcpp::List engine_main(Rcpp::List subsetAtomsList, double sigma, unsigned int maxit, double tol, bool showprog = false, std::string myfname="mpost.euc") {
  unsigned int nsubset = subsetAtomsList.size(); // # subset posteriors
  unsigned int ii, jj, kk, ndim, idx;
  
  arma::field<mat> subsetPosteriorSamplesList(nsubset);
  arma::mat subsetAtoms;
  SEXP subAtoms;
  for (ii = 0; ii < nsubset; ++ii) {
    subAtoms = subsetAtomsList[ii];
    subsetPosteriorSamplesList(ii) = Rcpp::as<arma::mat>(subAtoms);
  }
  
  ndim = subsetPosteriorSamplesList(0).n_cols; // # dimensions
  
  // vectors for storing sq-dist between M-posterior and subset posteriors
  arma::colvec snorms(nsubset,fill::zeros);
  arma::colvec norms(nsubset, fill::zeros);
  arma::colvec normsOld(nsubset, fill::zeros);
  
  // vectors for storing 1) # atoms in and 2) Weiszfeld wts of subset posteriors
  arma::colvec natoms(nsubset, fill::zeros);
  arma::colvec weiszfeldWts(nsubset, fill::zeros);

  arma::field<arma::mat> empiricalMeasureProbList(nsubset); // store the prob of atoms in each subset posterior 
  
  // initialize empirical subset posterior measure 
  for (ii = 0; ii < nsubset; ++ii) {
    empiricalMeasureProbList(ii) = ones(subsetPosteriorSamplesList(ii).n_rows, 1) / subsetPosteriorSamplesList(ii).n_rows;
    natoms(ii) = subsetPosteriorSamplesList(ii).n_rows;
  }
  
  arma::mat histWts(nsubset, maxit, fill::zeros);
  
  // initialize M-posterior atoms and their measures
  int sum_natoms = sum(natoms);
  arma::mat medianEmpiricalMeasureAtoms(sum_natoms, ndim, fill::zeros);
  arma::mat medianEmpiricalMeasureProbs(sum_natoms, 1,    fill::zeros);
  idx = 0;
  for (ii = 0; ii < nsubset; ++ii) {
    for (jj = 0; jj < natoms(ii); ++jj) {
      medianEmpiricalMeasureProbs(idx, 0) = empiricalMeasureProbList(ii)(jj, 0);
      for (kk = 0; kk < ndim; ++kk) {
        medianEmpiricalMeasureAtoms(idx, kk) = subsetPosteriorSamplesList(ii)(jj, kk);
      }
      idx = idx + 1;
    }
  }
  medianEmpiricalMeasureProbs = medianEmpiricalMeasureProbs / arma::accu(medianEmpiricalMeasureProbs);
  
  for (jj = 0; jj < maxit; ++jj) {
    if (((jj + 1) % 5 == 0)&&(showprog)){
      Rcpp::Rcout << " * " << myfname << " : iteration " << jj + 1 << "/" << maxit << " complete.." << std::endl;  
    } 
    
    norms = compute_wweights(sigma, subsetPosteriorSamplesList, medianEmpiricalMeasureAtoms, empiricalMeasureProbList, medianEmpiricalMeasureProbs);
    snorms = sqrt(norms);
    weiszfeldWts = (ones(nsubset) / snorms) / sum(ones(nsubset) / snorms);
    
    idx = 0;
    for (ii = 0; ii < nsubset; ++ii) {
      histWts(ii, jj) =  weiszfeldWts(ii);
      for (kk = 0; kk < natoms(ii); ++kk) {
        medianEmpiricalMeasureProbs(idx, 0) = weiszfeldWts(ii) / natoms(ii);
        idx = idx + 1;
      }
    }
    
    if (((sum(abs(norms - normsOld)) / nsubset) < tol) & (jj > 10)) {
      if (showprog){
        Rcpp::Rcout << " * " << myfname << " : iteration " << jj + 1 << "/" << maxit << " reached convergence." << std::endl;
      }
      break;
    }
    normsOld = norms;
  }  
  
  if (jj < maxit){
    histWts = histWts.cols(0, jj);
  } 
  
  Rcpp::List result = List::create(Named("natoms", natoms),
                                   Named("weiszfeldWts", weiszfeldWts),
                                   Named("historyWeiszfeldWts", histWts),
                                   Named("medianAtoms", medianEmpiricalMeasureAtoms)
  );  
  
  return(result);
}

