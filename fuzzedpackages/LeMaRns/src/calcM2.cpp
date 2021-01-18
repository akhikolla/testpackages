 #include <RcppArmadillo.h>
#include <Rinternals.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @title Calculate predation mortality
//'
//' @description Calculates the predation mortality for each species in each length class.
//' @param N A matrix with dimensions \code{nsc} and \code{nfish} representing the number of individuals in each length class for the current time step.
//' @param ration A matrix with dimensions \code{nsc} and \code{nfish} representing the amount of food required for fish of a given species and length class to grow according to the von Bertalanffy growth curve in a time step.
//' @param wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class.
//' @param nfish A numeric value representing the number of species in the model.
//' @param nsc A numeric value representing the number of length classes in the model.
//' @param other A numeric value representing the amount of other food (g) available from prey that is not explicitly represented in the model.
//' @param sc_Linf A numeric vector of length \code{nsc} representing the length class at which each species reaches its asymptotic length.
//' @param suit_M2 A list object of length \code{nfish}. Each element in the list is an array of dimensions \code{nsc}, \code{nsc} and \code{nfish} containing a value between zero and 1 representing prey preference and prey suitability for each species and length class.
//' @details The predation mortality of the \code{i}th species in the \code{j}th length class is
//' @details \code{sum_m(sum_n(I[j,i]*N[j,i]*suit_M2[[m]][n,j,i]/}
//' @details \code{(sum_k(sum_l(suit_M2[[m]][n,l,k]wgt[l,k]N[l,k]))+other)))}
//' @details where \code{sum_m} represents the sum over all \code{m}, \code{sum_n} represents the sum over all \code{n}, \code{sum_l} represents the sum over all \code{l} and \code{sum_k} represents the sum over all \code{k}. This equation corresponds to a Holling type-II functional response. See equation 8 of Hall et al. (2006) for more details.
//' @return A matrix with dimensions \code{nsc} and \code{nfish} representing the the predation mortality for each species in each length class.
//' @references Hall, S. J., Collie, J. S., Duplisea, D. E., Jennings, S., Bravington, M., & Link, J. (2006). A length-based multispecies model for evaluating community responses to fishing. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 63(6):1344-1359.
//' @examples
//' # Set up the inputs to the function - species-independent parameters
//' nfish <- nrow(NS_par)
//' nsc <- 32
//' maxsize <- max(NS_par$Linf)*1.01 # the biggest size is 1% bigger than the largest Linf
//' l_bound <- seq(0, maxsize, maxsize/nsc); l_bound <- l_bound[-length(l_bound)]
//' u_bound <- seq(maxsize/nsc, maxsize, maxsize/nsc)
//' mid <- l_bound+(u_bound-l_bound)/2
//'
//' # Set up the inputs to the function - species-specific parameters
//' Linf <- NS_par$Linf # the von-Bertalanffy asymptotic length of each species (cm).
//' W_a <- NS_par$W_a # length-weight conversion parameter.
//' W_b <- NS_par$W_b # length-weight conversion parameter.
//' k <- NS_par$k # the von-Bertalnaffy growth parameter.
//' Lmat <- NS_par$Lmat # the length at which 50\% of individuals are mature (cm).
//'
//' # Get phi_min
//' tmp <- calc_phi(k, Linf, nsc, nfish, u_bound, l_bound, calc_phi_min=FALSE,
//'                   phi_min=0.1) # fixed phi_min
//' phi <- tmp$phi
//' phi_min <- tmp$phi_min
//'
//' # Calculate growth increments
//' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
//' ration <- tmp$ration
//' sc_Linf <- tmp$sc_Linf
//' wgt <- tmp$wgt
//' g_eff <- tmp$g_eff
//'
//' # Calculate predator-prey size preferences
//' prefs <- calc_prefs(pred_mu=-2.25, pred_sigma=0.5, wgt, sc_Linf)
//'
//' # Calculate prey preference and prey suitability
//' suit_M2 <- calc_suit_vect(nsc, nfish, sc_Linf, prefs, NS_tau)
//'
//'
//' # Get an initial population
//' N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
//'
//' # Calculate the predation mortality
//' M2 <- calc_M2(N0, ration, wgt, nfish, nsc, other=1e12, sc_Linf, suit_M2)
//' @export
// [[Rcpp::export]]
arma::mat calc_M2(arma::mat N,
                     arma::mat ration,
                     arma::mat wgt,
                     double nfish,
                     double nsc,
                     double other,
                     Rcpp::NumericVector sc_Linf,
                     List suit_M2) {

  // Allocate the variables
	arma::mat M2(nsc, nfish);
  M2.fill(0);
	int j = 0;
	int i = 0;
	double denom = 0;

  // Calculations outside loop
  arma::mat wgtN = wgt % N;

// Calculations inside loop
  // Outer loop
  for(j = 0; j < nfish; j++) {

  // We extract the matrix for each fish in a loop
  arma::cube suite = suit_M2[j];

  // Local variable declaration
  int noi = sc_Linf(j);

  // Inner loop
  for(i = 0; i < noi; i++) {

// temp is now an element of the list of matrices
	  arma::mat temp = suite.row(i);
	  arma::mat tempmult = (temp % wgtN);
	  denom = accu(tempmult) + other;

	  if(denom > 0) {
	  M2 += ration(i,j) * N(i,j) * (temp/denom);
	  }

  }



  }
  return M2;

}
