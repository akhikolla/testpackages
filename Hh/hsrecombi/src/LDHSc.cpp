#include <Rcpp.h>
using namespace Rcpp;


//' Count genotype combinations at 2 SNPs
//'
//' @name countNumbers
//' @param X numeric matrix of genotypes
//' @return \code{count} vector of counts of 9 possible genotypes at SNP pair
//'
IntegerVector countNumbers(NumericMatrix X){
  int co = 0;
  IntegerVector count (9, 0);

  for(int k = 2; k >= 0; --k){
    for(int l = 2; l >= 0; --l){
      for(int i = 0; i < X.nrow(); ++i) {
        if((X(i, 0) == k) & (X(i, 1) == l)) count(co) += 1;
      }
      co += 1;
    }
  }

  return count;
}


//' Calculate log-likelihood function
//'
//' @name loglikfun
//'
//' @param counts integer vector of observed 2-locus genotype
//' @param fAA frequency of maternal haplotype 1-1
//' @param fAB frequency of maternal haplotype 1-0
//' @param fBA frequency of maternal haplotype 0-1
//' @param fBB frequency of maternal haplotype 0-0
//' @param theta paternal recombination rate
//' @return \code{lik} value of log likelihood at parameter estimates
//'
double loglikfun(IntegerVector counts, double fAA, double fAB, double fBA, double fBB, double theta){
  NumericVector fak (9);
  double lik = 0;

  fak(8) = (1 - theta) / 2 * fBB;
  fak(7) = (1 - theta) / 2 * fBA + theta / 2 * fBB;
  fak(6) = theta / 2 * fBA;
  fak(5) = (1 - theta) / 2 * fAB + theta / 2 * fBB;
  fak(4) = (1 - theta) / 2 * (fBB + fAA) + theta / 2 * (fAB + fBA);
  fak(3) = (1 - theta) / 2 * fBA + theta / 2 * fAA;
  fak(2) = theta / 2 * fAB;
  fak(1) = (1 - theta) / 2 * fAB + theta / 2 * fAA;
  fak(0) = (1 - theta) / 2 * fAA;

  for(int l = 0; l < 9; ++l){
    if(fak(l) > 1e-6) lik += counts(l) * std::log(fak(l));
  }

  return lik;
}


//' Expectation Maximisation (EM) algorithm
//'
//' @name LDHScpp
//'
//' @param XGF1 numeric matrix of progeny genotypes in genomic family 1
//' @param XGF2 numeric matrix of progeny genotypes in genomic family 2
//' @param fAA frequency of maternal haplotype 1-1
//' @param fAB frequency of maternal haplotype 1-0
//' @param fBA frequency of maternal haplotype 0-1
//' @param theta paternal recombination rate
//' @param display logical for displaying additional information
//' @param threshold convergence criterion
//' @return list of parameter estimates
//' \describe{
//'  \item{\code{D}}{maternal LD}
//'  \item{\code{fAA}}{frequency of maternal haplotype 1-1}
//'  \item{\code{fAB}}{frequency of maternal haplotype 1-0}
//'  \item{\code{fBA}}{frequency of maternal haplotype 0-1}
//'  \item{\code{fBB}}{frequency of maternal haplotype 0-0}
//'  \item{\code{p1}}{Maternal allele frequency (allele 1)}
//'  \item{\code{p2}}{Maternal allele frequency (allele 0)}
//'  \item{\code{nfam1}}{size of genomic family 1}
//'  \item{\code{nfam2}}{size of genomic family 2}
//'  \item{\code{error}}{0 if computations were without error; 1 if EM algorithm
//'   did not converge}
//'  \item{\code{iteration}}{number of EM iterations}
//'  \item{\code{theta}}{paternal recombination rate}
//'  \item{\code{r2}}{\eqn{r^2} of maternal LD}
//'  \item{\code{logL}}{value of log likelihood function}
//' }
//' @export
// [[Rcpp::export]]
List LDHScpp(Nullable<NumericMatrix> XGF1,
               Nullable<NumericMatrix> XGF2,
               double fAA, double fAB, double fBA, double theta, bool display, double threshold){

  int nAAAA, nAAAB, nAABB, nABAA, nABAB, nABBB, nBBAA, nBBAB, nBBBB, N, iter;
  int nfam1 = 0, nfam2 = 0, maxit = 10000, errorlevel = 1;
  IntegerVector n1 (9), n2 (9);
  double fAA1, fAB1, fBA1, fBB1, theta1, fAA2, fAB2, fBA2, fBB2, theta2;
  double fAANew, fABNew, fBANew, fBBNew, thetaNew, den, LD, r2, p1, p2, loglik = 0;
  double fBB = 1 - fAA - fAB - fBA;


  if (XGF1.isNotNull()){
    NumericMatrix fam1 = as<NumericMatrix>(XGF1);
    nfam1 = fam1.nrow();
    n1 = countNumbers(fam1);
    if(display) Rcout << n1 << "\n";
    if(fam1.ncol() > 2) Rcout << "WARNING: only the first two columns of *genomic family 1* will be used \n";
  }
  if (XGF2.isNotNull()){
    NumericMatrix fam2 = as<NumericMatrix>(XGF2);
    nfam2 = fam2.nrow();
    n2 = countNumbers(fam2);
    if(display) Rcout << n2 << "\n";
    if(fam2.ncol() > 2) Rcout << "WARNING: only the first two columns of *genomic family 2* will be used \n";
  }

  if((nfam1 + nfam2 > 0) & (fBB >= 0)){

    for(iter = 0; iter <= maxit; ++iter){

      if(nfam1 > 0){
        nAAAA = n1(0); nAAAB = n1(1); nAABB = n1(2);
        nABAA = n1(3); nABAB = n1(4); nABBB = n1(5);
        nBBAA = n1(6); nBBAB = n1(7); nBBBB = n1(8); N = nfam1;
        den = theta * (fAB + fBA) + (1 - theta) * (fAA + fBB);

        fAA1 = (nAAAA + theta  * fAA * nAAAB / (theta * fAA + (1 - theta) * fAB) +
          theta * fAA * nABAA / (theta * fAA + (1 - theta) * fBA) + (1 - theta) * fAA * nABAB / den) / N;
        fAA1 = std::max(fAA1, threshold);

        fAB1 = (nAABB + (1 - theta) * fAB *nAAAB / (theta * fAA + (1 - theta) * fAB) +
          (1 - theta) * fAB * nABBB / (theta * fBB + (1 - theta) * fAB) + theta * fAB * nABAB / den) / N;
        fAB1 = std::max(fAB1, threshold);

        fBA1 = (nBBAA + (1 - theta) * fBA * nABAA / (theta * fAA + (1 - theta) * fBA) +
          (1 - theta) * fBA * nBBAB / (theta * fBB + (1 - theta) * fBA) + theta * fBA * nABAB / den) / N;
        fBA1 = std::max(fBA1, threshold);

        fBB1 = (nBBBB + theta * fBB * nABBB / (theta * fBB + (1 - theta) * fAB) +
          theta * fBB * nBBAB / (theta * fBB + (1 - theta) * fBA) + (1 - theta) * fBB * nABAB / den) / N;
        fBB1 = std::max(fBB1, threshold);

        /* corrected version (pers. comm. GR 25/09/15) */
        theta1 = (nAABB + nBBAA + theta * fAA * nAAAB / (theta * fAA + (1 - theta) * fAB) +
          theta * fAA * nABAA / (theta * fAA + (1 - theta) * fBA) + theta * (fAB + fBA) * nABAB / den +
          theta * fBB * nABBB / (theta * fBB + (1 - theta) * fAB) + theta * fBB * nBBAB / (theta * fBB + (1 - theta) * fBA)) / N;
        theta1 = std::max(theta1, threshold);
        theta1 = std::min(theta1, 1.0);
      } else fAA1 = fAB1 = fBA1 = fBB1 = theta1 = 0;

      if(nfam2 > 0){
        nAAAA = n2(0); nAAAB = n2(1); nAABB = n2(2);
        nABAA = n2(3); nABAB = n2(4); nABBB = n2(5);
        nBBAA = n2(6); nBBAB = n2(7); nBBBB = n2(8); N = nfam2;
        den = (1 - theta) * (fAB + fBA) + theta * (fAA + fBB);

        fAA2 = (nAAAA + (1 - theta) * fAA * nAAAB / ((1 - theta) * fAA + theta * fAB) +
          (1 - theta) * fAA * nABAA / ((1 - theta) * fAA + theta * fBA) + theta * fAA * nABAB / den) / N;
        fAA2 = std::max(fAA2, threshold);

        fAB2 = (nAABB + theta * fAB * nAAAB / ((1 - theta) * fAA + theta * fAB) +
          theta * fAB * nABBB / ((1 - theta) * fBB + theta * fAB) + (1 - theta) * fAB * nABAB / den) / N;
        fAB2 = std::max(fAB2, threshold);

        fBA2 = (nBBAA + theta * fBA * nABAA / ((1 - theta) * fAA + theta * fBA) +
          theta * fBA * nBBAB / ((1 - theta) * fBB + theta * fBA) + (1 - theta) * fBA * nABAB / den) / N;
        fBA2 = std::max(fBA2, threshold);

        fBB2 = (nBBBB + (1 - theta) * fBB * nABBB / ((1 - theta) * fBB + theta * fAB) +
          (1 - theta) * fBB * nBBAB / ((1 - theta) * fBB + theta * fBA) + theta * fBB * nABAB / den) / N;
        fBB2 = std::max(fBB2, threshold);

        /* corrected version (pers. comm. GR 25/09/15) */
        theta2 = (nAAAA + nBBBB + theta * fAB * nAAAB / (theta * fAB + (1 - theta) * fAA) +
          theta * fBA * nABAA / (theta * fBA + (1 - theta) * fAA) + theta * (fAA + fBB) * nABAB / den +
          theta * fAB * nABBB / (theta * fAB + (1 - theta) * fBB) + theta * fBA * nBBAB / (theta * fBA + (1 - theta) * fBB)) / N;
        theta2 = std::max(theta2, threshold);
        theta2 = std::min(theta2, 1.0);
      } else fAA2 = fAB2 = fBA2 = fBB2 = theta2 = 0;

      fAANew = (nfam1 * fAA1 + nfam2 * fAA2) / (nfam1 + nfam2);
      fABNew = (nfam1 * fAB1 + nfam2 * fAB2) / (nfam1 + nfam2);
      fBANew = (nfam1 * fBA1 + nfam2 * fBA2) / (nfam1 + nfam2);
      fBBNew = (nfam1 * fBB1 + nfam2 * fBB2) / (nfam1 + nfam2);
      thetaNew = (nfam1 * theta1 + nfam2 * theta2) / (nfam1 + nfam2);

      if((std::abs(fAANew - fAA) < threshold) & (std::abs(fBANew - fBA) < threshold) &
           (std::abs(fABNew - fAB) < threshold) & (std::abs(fBBNew - fBB) < threshold) &
           (std::abs(thetaNew - theta) < threshold)) {
        if(display) Rcout << "Convergence after " << iter << " iterations\n";
        break;
      }

      fAA = fAANew;
      fBA = fBANew;
      fAB = fABNew;
      fBB = fBBNew;
      theta = thetaNew;

      if(iter == maxit){
        Rcout << "Warning: No convergence \n";
      } else errorlevel = 0;
    }

    LD = fAA * fBB - fAB * fBA;
    p1 = fAA + fAB;
    p2 = fAA + fBA;
    r2 = LD * LD / (p1 * (1 - p1) * p2 * (1 - p2));

    if(nfam1 > 0){
      loglik += loglikfun(n1, fAA, fAB, fBA, fBB, theta);
    }
    if(nfam2 > 0){
      loglik += loglikfun(n2, fAA, fAB, fBA, fBB, 1 - theta);
    }

  } else if(fBB < 0) Rcerr << "ERROR: allele frequencies (disregard marker pair)\n";

  return List::create(_["D"] = LD, _["fAA"] = fAA, _["fAB"] = fAB, _["fBA"] = fBA, _["fBB"] = fBB,
                      _["p1"] = p1, _["p2"] = p2, _["nfam1"] = nfam1, _["nfam2"] = nfam2,
                      _["error"] = errorlevel, _["iteration"] = iter, _["theta"] = theta, _["r2"] = r2, _["logL"] = loglik);
}



