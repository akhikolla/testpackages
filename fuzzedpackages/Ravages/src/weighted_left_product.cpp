/*******************************************
 * LOURDEMENT SIMILAIRE A m4_prod_ms.cpp   *
 *                                         *
 * ce produit n'est pas centré             *
 * mais il est pondéré...                  *
 *                                         *
 *******************************************/
#include "weighted_left_product.h"

using namespace Rcpp;
using namespace RcppParallel;

// WLP = weighted left product
struct paraWLP : public Worker {
  const uint8_t ** data;
  const double * p; // le vecteur de fréquences allélique (freq A2)
  const size_t nrow;
  const size_t ncol;
  const size_t true_ncol;
  const std::vector<double> we;
  const size_t r; // nb cols résultats
  double * v;

  //output
  double * Av;

  //constructeur
  paraWLP(const uint8_t ** data, const double * p, size_t nrow, size_t ncol, size_t true_ncol, std::vector<double> we, 
          size_t r, double * v, double * Av)
          : data(data), p(p), nrow(nrow), ncol(ncol), true_ncol(true_ncol), we(we), r(r), v(v), Av(Av) { }


  void operator()(size_t beg, size_t end) {
    double gg[4];
    gg[0] = 0;
    for(size_t i = beg; i < end; i++) {
// std::cout << "i = " << i << "\n";
      gg[1] = we[i];
      gg[2] = 2*we[i];
      gg[3] = 2*p[i]; // imputation par le "génotype moyen"
      for(size_t c = 0; c < r; c++) {
        size_t k = c*ncol;
        for(size_t j = 0; j < true_ncol; j++) {
          uint8_t x = data[i][j];
          for(int ss = 0; ss < 4 && (4*j + ss < ncol); ss++) {
// std::cout << "v[k]*gg[x&3] = " << v[k]*gg[x&3] << "\n";
            Av[nrow*c+i] += v[k++]*gg[x&3];
            x >>= 2;
          }
        }
      }
    }
  }
};

//[[Rcpp::export]]
NumericMatrix WLP(XPtr<matrix4> pA, NumericVector p, const std::vector<double> & we, NumericMatrix & v) {
  return WLP(const_cast<const uint8_t **>(pA->data), const_cast<const double *>(&p[0]), pA->nrow, pA->ncol, pA->true_ncol, we, v);
}

// calcule R = v' GW avec G = génotype 0 1 2 (donné par pA), W = matrice diagonale des poids (donnés par we)
// le résultat a dimensions nb_snps x v.ncol() ...
NumericMatrix WLP(const uint8_t ** A_data, const double * p, size_t A_nrow, size_t A_ncol, size_t A_true_ncol, 
                  const std::vector<double> & we, NumericMatrix & v) {
  int n = A_nrow; // nb snps
  int m = A_ncol; // nb inds
  // Rcout << "m = " << m << " v.nrow = " << v.nrow() << "\n";
  // Rcout << "n = " << n << " we.size = " << we.size() << "\n";
  if(m != v.nrow() || n != we.size()) stop("Dimensions mismatch");
  int r = v.ncol();

  NumericMatrix R(n,r);
  paraWLP X(A_data, p, A_nrow, A_ncol, A_true_ncol, we, r, v.begin(), R.begin());

  parallelFor(0, n, X);
  return R;
}


