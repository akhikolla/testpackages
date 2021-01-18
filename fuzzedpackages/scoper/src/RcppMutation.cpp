#include <Rcpp.h>
#include <iostream>
#include <map>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// pairwiseMutMatrix
// [[Rcpp::export]]
List pairwiseMutMatrixRcpp(NumericVector informative_pos,
                           StringMatrix mutMtx,
                           NumericMatrix motifMtx) {
    int i, j, k;
    double nu, shmut, allmut, nonmutab;
    int n = mutMtx.nrow();
    int m = mutMtx.ncol();
    // allocate  matrices
    NumericMatrix sh_mtx(n,n);
    NumericMatrix tot_mtx(n,n);
    NumericMatrix mutab_mtx(n,n);
    StringVector row_i(m);
    StringVector row_j(m);
    
    // initialize  matrices
    std::fill(sh_mtx.begin(), sh_mtx.end(), 0);
    std::fill(tot_mtx.begin(), tot_mtx.end(), 0);
    std::fill(mutab_mtx.begin(), mutab_mtx.end(), 1);
    
    // fill matrices
    for (i = 0; i < n-1; i++) {
        row_i = mutMtx.row(i);
        for (j = i+1; j < n; j++) {
            row_j = mutMtx.row(j);
            shmut = 0; allmut = 0; nonmutab = 1.;
            for (k = 0; k < m; k++) {
                if ((row_i(k) == row_j(k)) && (row_i(k) != "NA") && (row_j(k) != "NA")) {
                    shmut++;
                    nonmutab *= (1. - (motifMtx(i,k) + motifMtx(j,k))/2.);
                }
                if (row_i(k) != "NA") allmut++;
                if (row_j(k) != "NA") allmut++;
            }
            nu = 2./(informative_pos(i) + informative_pos(j));
            sh_mtx(i,j) = 2.*nu*shmut;
            tot_mtx(i,j) = nu*allmut;
            mutab_mtx(i,j) = nonmutab;
        }
    }
    
    List results;
    results["sh_mtx"] = sh_mtx;
    results["tot_mtx"] = tot_mtx;
    results["mutab_mtx"] = mutab_mtx;
    return results;
}
