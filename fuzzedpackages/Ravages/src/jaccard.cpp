// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include "gaston/matrix4.h"

using namespace Rcpp;
using namespace RcppParallel;


// ************* calcule deux matrices
// A = matrice du nombre de fois où les individus sont tous les deux à 1 (ou 2)
// B = matrice du nombre de fois où un des deux est à 1 (ou 2)
// le vecteur inverse permet d'inverser les genotypes si le génotype fréquent n'est pas 0
// les génotypes manquants sont imputés par "fréquent"

struct jaccard_para : public Worker {
  // input and others
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol;
  std::vector<bool> inverse;
  int size;

  // output
  int * A;
  int * B;
  
  // constructeurs
  jaccard_para(uint8_t ** data, const size_t ncol, const size_t true_ncol, std::vector<bool> inverse) :
          data(data), ncol(ncol), true_ncol(true_ncol), inverse(inverse), size((4*true_ncol)*(4*true_ncol+1)/2) { 
          A = new int [size];  // padded to a multiple of 4...
          B = new int [size];  
          std::fill(A, A+size, 0);
          std::fill(B, B+size, 0);
        }
  jaccard_para(jaccard_para & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), inverse(Q.inverse), size(Q.size) {
          A = new int [size]; 
          B = new int [size];  
          std::fill(A, A+size, 0);
          std::fill(B, B+size, 0);
        }

  // destructeur
  ~jaccard_para() { 
          delete [] A; 
          delete [] B; 
  }

  // worker !
  void operator()(size_t beg, size_t end) {
    for(size_t i = beg; i < end; i++) {
      uint8_t * dd = data[i];
      size_t k = 0; 
      if(inverse[i]) { // inversion du genotype
         // 2 (et 3 : manquant) sont le génotype fréquent
         // 1 (et 0) le génotype rare
        for(size_t j1 = 0; j1 < true_ncol; j1++) {
          uint8_t x1 = dd[j1];
          for(unsigned int ss1 = 0; (ss1 < 4); ss1++) {
            if((x1&3) == 2 || (x1&3) == 3) { // genotype fréquent
              // dans A on n'aura que des 0 à ajouter -> ràfaire
              // dans B on ajoute 1 quand le génotype de l'individu est 1 (ou 0)
              for(size_t j2 = 0; j2 < j1; j2++) {
                uint8_t x2 = ~(dd[j2]);
                B[k++] += ((x2&2)>>1);
                B[k++] += ((x2&8)>>3);
                B[k++] += ((x2&32)>>5);
                B[k++] += ((x2&128)>>7);
              }
              uint8_t x2 = ~(dd[j1]);
              for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
                B[k++] += ((x2&2)>>1);
                x2 >>= 2;
              } 
            } else { // génotype rare
              // dans B on ajoute des 1 pour tout le monde
              for(size_t jj = 0; jj <= 4*j1+ss1; jj++) B[k+jj]++;
              // dans A on ajoute 1 quand le génotype de l'individu est 1 (ou 0)
              for(size_t j2 = 0; j2 < j1; j2++) {
                uint8_t x2 = ~(dd[j2]);
                A[k++] += ((x2&2)>>1);
                A[k++] += ((x2&8)>>3);
                A[k++] += ((x2&32)>>5);
                A[k++] += ((x2&128)>>7);
              }
              uint8_t x2 = ~(dd[j1]);
              for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
                A[k++] += ((x2&2)>>1);
                x2 >>= 2;
              } 
            }
            x1 >>= 2; // next
          }
        }
      } else {
         // 0 (et 3 : manquant) sont le génotype fréquent
         // 1 (et 2) le génotype rare
        for(size_t j1 = 0; j1 < true_ncol; j1++) {
          uint8_t x1 = dd[j1];
          for(unsigned int ss1 = 0; (ss1 < 4); ss1++) {
            if((x1&3) == 0 || (x1&3) == 3) { // genotype fréquent
              // dans A on n'aura que des 0 à ajouter -> ràfaire
              // dans B on ajoute 1 quand le génotype de l'individu est 1 ou 2
              for(size_t j2 = 0; j2 < j1; j2++) {
                uint8_t x2 = dd[j2];
                x2 ^= (x2>>1);
                B[k++] += (x2&1);
                B[k++] += ((x2&4)>>2);
                B[k++] += ((x2&16)>>4);
                B[k++] += ((x2&64)>>6);
              }
              uint8_t x2 = dd[j1];
              x2 ^= (x2>>1);
              for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
                B[k++] += (x2&1);
                x2 >>= 2;
              } 
            } else { // génotype rare
              // dans B on ajoute des 1 pour tout le monde
              for(size_t jj = 0; jj <= 4*j1+ss1; jj++) B[k+jj]++;
              // dans A on ajoute 1 quand le génotype de l'individu est 1 (ou 2)
              for(size_t j2 = 0; j2 < j1; j2++) {
                uint8_t x2 = dd[j2];
                x2 ^= (x2>>1);
                A[k++] += (x2&1);
                A[k++] += ((x2&4)>>2);
                A[k++] += ((x2&16)>>4);
                A[k++] += ((x2&64)>>6);
              }
              uint8_t x2 = dd[j1];
              x2 ^= (x2>>1);
              for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
                A[k++] += (x2&1);
                x2 >>= 2;
              } 
            }
            x1 >>= 2; // next
          }
        }
      } 
    }
  }

  // recoller
  void join(const jaccard_para & Q) {
    std::transform(A, A + size, Q.A, A, std::plus<int>());
    std::transform(B, B + size, Q.B, B, std::plus<int>());
  }

};

//[[Rcpp::export]]
List Jaccard(XPtr<matrix4> p_A, LogicalVector which_snps, LogicalVector inverse) {
  int nb_snps = sum(which_snps);

  if(which_snps.length() != p_A->nrow)
    stop("Dimensions mismatch");

  uint8_t ** data = new uint8_t * [nb_snps];
  std::vector<bool> inv(nb_snps);
  size_t k = 0;
  for(size_t i = 0; i < p_A->nrow; i++) {
    if(which_snps[i]) {
      inv[k] = inverse[i];
      data[k++] = p_A->data[i];
    }
  }
  jaccard_para X(data, p_A->ncol, p_A->true_ncol, inv);
  parallelReduce(0, nb_snps, X);

  delete [] data;

  NumericMatrix A(p_A->ncol,p_A->ncol);
  NumericMatrix B(p_A->ncol,p_A->ncol);
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      A(j,i) = (double) X.A[k++];
    }
  }
  // symmetriser
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      A(i,j) = (double) X.A[k++]; // ou A(j,i)
    }
  }

  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      B(j,i) = (double) X.B[k++];
    }
  }
  // symmetriser
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      B(i,j) = (double) X.B[k++]; // ou A(j,i)
    }
  }

  List L;
  L["A"] = A; L["B"] = B;
  return L;
}

RcppExport SEXP oz_jaccard(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type inverse(inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(Jaccard(p_A, which_snps, inverse));
    return rcpp_result_gen;
END_RCPP
}

