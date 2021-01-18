// funchisq.cpp
//
// Created: April 30, 2016. Extracted from ExactFunctionalTest.cpp

#include "define.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>

mydouble funchisq(const std::vector<std::vector<int> > & O, mydouble & estimate,
                  const std::string index_kind){
    mydouble fc = 0.0;
    if (O.size() == 0) {
        return fc;
    } else if(O[0].size() == 0) {
        return fc;
    }
    
    int n = 0;
    std::vector <int> colsums ((int)O[0].size(), 0);
    std::vector <int> rowsums ((int)O.size(), 0);
    
    for (size_t i=0; i<O.size(); i++) {
        for (size_t j=0; j<O[i].size(); j++) {
            n += O[i][j];
            colsums[j] += O[i][j];
            rowsums[i] += O[i][j];
        }
    }
    
    if(n == 0)return fc;
    
    size_t nrows = O.size();  // number of rows
    size_t ncols = O[0].size();  // number of columns
    
    mydouble ej = n / (mydouble) ncols;
    mydouble col_chisq = 0.0;
    
    if(ej>0){
        for (size_t j=0; j<ncols; ++j) {
            col_chisq += (colsums[j] - ej) * (colsums[j] - ej) / ej;
        }
    }
    fc -= col_chisq;
    
    for (size_t i=0; i<nrows; ++i) {
        // Expected cound for cell (i,j):
        mydouble eij = rowsums[i] / (mydouble) ncols;
        if (eij > 0) {
            for (size_t j=0; j<ncols; ++j) {
                fc += (O[i][j] - eij) * (O[i][j] - eij) / eij;
            }
        }
    }
    
    mydouble maxfc = -1;
    
    // Version 1: Current. Added June 11, 2017. Bound is maximum reachable
    if (index_kind == "conditional") { // conditional on column marginal sums
        
        // max.fun.chisq <- sum(sort(col.sum, decreasing=TRUE)[seq(min(dim(x)))]) *
        //  ncol(x) - n - col.chisq
        
        // std::sort(&colsums[0], &colsums[ncols], std::greater<int>());
        // maxfc = std::accumulate(&colsums[0], &colsums[std::min(nrows, ncols)], 0);
        // maxfc = maxfc * ncols - n - col_chisq;
        
        maxfc = n * ncols - n - col_chisq;
        
    } else if(index_kind == "unconditional"){
        maxfc = n * ncols * (1 - 1.0 / std::min(nrows, ncols));
    }
    
    //  Version 0: bound is not reachable when nrow(x) < ncol(x)
    else if (index_kind == "conditional-version-0") {
        maxfc = n * (ncols - 1) - col_chisq;
    }else if(index_kind == "unconditional-version-0"){
        maxfc =  n * (ncols - 1);
    }
    
    if(maxfc > 0) {
        estimate = std::sqrt(std::abs(fc) / maxfc);
    } else {
        estimate = 0;
    }
    
    return fc;
}

mydouble funchisq0(const std::vector<std::vector<int> > & O, const std::vector<int> & rowsums,
                  const std::vector<int> & colsums, int n)
{
    mydouble fc = 0.0;
    
    if (n == 0 || O.size() == 0) {
        return fc;
    } else if(O[0].size() == 0) {
        return fc;
    }
    
    size_t nrows = O.size();  // number of rows
    size_t ncols = O[0].size();  // number of columns
    
    mydouble ej = n / (mydouble) ncols;
    if(ej>0){
        for (size_t j=0; j<ncols; ++j) {
            fc -= (colsums[j] - ej) * (colsums[j] - ej) / ej;
        }
    }
    
    for (size_t i=0; i<nrows; ++i) {
        // Expected cound for cell (i,j):
        mydouble eij = rowsums[i] / (mydouble) ncols;
        if (eij > 0) {
            for (size_t j=0; j<ncols; ++j) {
                fc += (O[i][j] - eij) * (O[i][j] - eij) / eij;
            }
        }
    }
    
    return fc;
}


mydouble funchisq(const std::vector<std::vector<int> > & O, const std::vector<int> & rowsums,
                   const std::vector<int> & colsums, int n)
{
    mydouble fc = 0.0;
    
    if (n == 0 || O.size() == 0) {
        return fc;
    } else if(O[0].size() == 0) {
        return fc;
    }
    
    size_t nrows = O.size();  // number of rows
    size_t ncols = O[0].size();  // number of columns
    
    for (size_t j=0; j<ncols; ++j) {
        fc -= colsums[j] * colsums[j];
    }
    fc = fc * ncols / n;
    
    for (size_t i=0; i<nrows; ++i) {
        if (rowsums[i] > 0) {
            mydouble fc_row_i=0;
            for (size_t j=0; j<ncols; ++j) {
                fc_row_i += O[i][j] * O[i][j];
            }
            fc_row_i = fc_row_i * ncols / rowsums[i];
            fc += fc_row_i;
        }
    }
    
    return fc;
}
