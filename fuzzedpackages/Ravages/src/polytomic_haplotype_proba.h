#include <Rcpp.h>
#include <cmath>
#ifndef Polytomic_haplo_probs
#define Polytomic_haplo_probs

using namespace Rcpp;

// Affreuse fonction bijective qui permet de passer de n à (i,j) avec i <= j
// la réciproque est bien sûr (i,j) -> i*(i+1)/2 + j
template<typename T1, typename T2>
inline void one_to_pair(T1 n, T2 & i, T2 & j) {
  i = (T2) ( 0.5*(sqrt(8.0*n +1.0) - 1.0) ); // horreur de passer dans R
  j = n - (T2) ((i*(i+1))/2);
}



class haplo_probs {
  int n;
  std::vector<double> burden;
  double sd;
  double th1;
  double th2;
  std::vector<double> probas;
  const bool probs;
  public:
  // un vecteur de fardeaux, l'écart-type sd_E de la loi normale à ajouter aux fardeaux
  // (le fardeau d'un individu hi hj est bi + bj + E)
  // les seuils sur la loi normale
  // tous les haplos sont équiprobables a priori
  haplo_probs(NumericVector burden_, double sd_E, double low_threshold, double high_threshold) :
    n(burden_.size()), burden(as<std::vector<double> >(burden_)), sd(sd_E), 
    th1(low_threshold), th2(high_threshold), probas(0), probs(false) { 
  }

  // idem mais on donne des probas haplotypiques
  haplo_probs(NumericVector burden_, NumericVector probas_, double sd_, double low_threshold, double high_threshold) :
      n(burden_.size()), burden(as<std::vector<double> >(burden_)), sd(sd_), 
      th1(low_threshold), th2(high_threshold), probas(as<std::vector<double> >(probas_)), 
      probs(true) {
    if(probas.size() != burden.size())
      stop("Dimensions mismatch");
  }

  // renvoie un nombre proportionnel à p( i,j | liability \in (th1, th2) )
  // NOTE c'est juste Bayes : p(liability | i,j ) * p(i,j) / p(liability \in (th1, th2) )
  // on a p(  liability \in (th1, th2) | i,j ) =  p(  burden(i) + burden(j) + E \in (th1, th2) )
  // on ignore le dénominateur qui est constant
  // !!!! il faut 0 <= i <= j <= n, ça n'est pas testé pour plus de rapidité !!!!
  inline double operator()(size_t i, size_t j) {
    double b = burden[i] + burden[j];
    // P(th1 < E < th2) pour E d'espérance b et d'écart type sd
    double p = R::pnorm(th2, b, sd, 1, 0) - R::pnorm(th1, b, sd, 1, 0);
    if(probs) {
      if(i != j) 
        return 2*probas[i]*probas[j]*p;
      else
        return probas[i]*probas[i]*p;
    } else {
      if(i != j) 
        return 2*p;
      else
        return p;
    }
  }

  // idem avec la paire d'haplotypes indexée par k entre 0 en n(n+1)/2
  inline double operator()(size_t k) {
    size_t i, j;
    one_to_pair(k, i, j);
    return this->operator()(i,j);
  }
};
#endif
