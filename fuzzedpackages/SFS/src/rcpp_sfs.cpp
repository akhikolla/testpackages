/* -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

   rcpp_sfs.cpp -- R wrappers for SFS class
   
   This file is part of SFS.
   
   SFS is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   SFS is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with SFS.  If not, see <http://www.gnu.org/licenses/>.
   
   (C) 2016 Utz-Uwe Haus, Cray EMEA Research Lab, Cray Inc.
*/


#include "RcppArmadillo.h"

#include "SFSMatrix.h"

#include <iostream>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

template <typename T>
static std::ostream&
operator<< (std::ostream& os, const std::vector<T>& v) 
{
  bool first=true;
  os << "[";
  for (typename std::vector<T>::const_iterator ii = v.begin();
       ii != v.end();
       ++ii,first=false)
    os << (first? "" : " ") << *ii;
  os << "]";
  return os;
}

static
Rcpp::DataFrame
sfs__matrix_to_DataFrame(SEXP E, const double zero_eps) {
    
    //create a numeric matrix to read SEXP and create 3-columns Data Frame
    Rcpp::NumericMatrix matrix(E);
    Rcpp::IntegerVector rowidx,colidx;
    Rcpp::NumericVector value;
    for(size_t row=0; row<(unsigned)matrix.nrow(); row++) {
        for(size_t col=0; col<(unsigned)matrix.ncol(); col++) {
            double v = matrix(row,col);
            if(row!= col && fabs(v)>zero_eps) {
                rowidx.push_back(row+1);
                colidx.push_back(col+1);
                value.push_back(v);
            }
        }
    }
    return Rcpp::DataFrame::create(Rcpp::Named("V1")=rowidx,Rcpp::Named("V2")=colidx,Rcpp::Named("V3")=value);
}

// from seriation/src/lt.h: col/row 1-based
#ifndef LT_POS
#define LT_POS(n, col, row)					\
    (col)==(row) ? 0                                            \
          : (col)<(row)                                         \
                  ? n*((col)-1) - (col)*((col)-1)/2 + (row)-(col) -1    \
                  : n*((row)-1) - (row)*((row)-1)/2 + (col)-(row) -1
#endif

static
Rcpp::DataFrame
sfs__dist_to_DataFrame(SEXP R_dist, const double zero_eps) {
    // dist objects as used in the seriation package: lower-triangular,
    // without diagonal elements. Typically rather dense; access via
    // LT_POS(size,col_idx,row_idx) (1-based)
    size_t length = LENGTH(R_dist);
    size_t n = (size_t)(1.0 + sqrt((float)(2*length)));
    if(length != n*(n-1)/2) {
        throw std::runtime_error("dist object has invalid length");
    }
    const double *dist = REAL(R_dist);

    Rcpp::IntegerVector rowidx,colidx;
    Rcpp::NumericVector value;
    for(size_t col=0; col<n; col++) {
        for(size_t row=0; row<col; row++) {
            double v = dist[LT_POS(n,col+1,row+1)];
            if(row!= col && fabs(v)>zero_eps) {
                rowidx.push_back(row+1);
                colidx.push_back(col+1);
                value.push_back(v);
                // add symmetric edges
                rowidx.push_back(col+1);
                colidx.push_back(row+1);
                value.push_back(v);
            }
        }
    }

    return Rcpp::DataFrame::create(Rcpp::Named("V1")=rowidx,Rcpp::Named("V2")=colidx,Rcpp::Named("V3")=value);
}

static
Rcpp::DataFrame
sfs__dataframe_to_DataFrame(SEXP E, const double zero_eps, bool symmetric, bool identical_val) {
    Rcpp::DataFrame D = Rcpp::as<Rcpp::DataFrame>(E);
    size_t num_cols = D.size();
    size_t num_rows = D.nrows();
    Rcpp::IntegerVector rowidx,colidx;
    Rcpp::NumericVector value;
    
    if(num_cols==3)  //3 columns assume they're (row, col, val) format
    {
        Rcpp::IntegerVector df_rowidx = D[0];
        Rcpp::IntegerVector df_colidx = D[1];
        Rcpp::NumericVector df_value = D[2];
        for(size_t k=0; k<num_rows; k++)
        {
            if (df_rowidx[k]!= df_colidx[k] && fabs(df_value[k]) > zero_eps)
            {
                rowidx.push_back(df_rowidx[k]);
                colidx.push_back(df_colidx[k]);
                value.push_back(df_value[k]);
                if (symmetric && !identical_val) //then add symmetric edges
                {
                    rowidx.push_back(df_colidx[k]);
                    colidx.push_back(df_rowidx[k]);
                    value.push_back(df_value[k]);
                }
                
                //eliminate identical values (otherwise problems with spMat)
//                if (identical_val) //then add symmetric edges
//                {
//                    int i = 0;
//                    int j = 0;
//                    while (i < rowidx.size() -1)
//                    {
//                        j = i + 1;
//                        while (j < rowidx.size())
//                        {
//                            if (rowidx[i] == rowidx[j] && colidx[i] == colidx[j])
//                            {
//                                rowidx.erase(rowidx.begin()+j);
//                                colidx.erase(colidx.begin()+j);
//                            }
//                            else
//                            {
//                                j++;
//                            }
//                        }
//                        i++;
//                    }
//                }
            }
        }
    }
    else //more than 3 columns assume that is a matrix format
    {
        for(size_t i=0; i<num_rows; i++)
        {
            for(size_t j=0; j<num_rows; j++)
            {
                Rcpp::NumericVector df_rowidx = D[i];
                if (i != j && fabs(df_rowidx[j]) > zero_eps)
                {
                    rowidx.push_back(i+1);
                    colidx.push_back(j+1);
                    value.push_back(df_rowidx[j]);
                }
            }
        }
    }
    return Rcpp::DataFrame::create(Rcpp::Named("V1")=rowidx,Rcpp::Named("V2")=colidx,Rcpp::Named("V3")=value);
}


// exported functionality: Read different input data formats and conver it to 3-columns (row,col,val) 'data frame'
// We support 'matrix', the dense numeric format
//            'dataframe', which can be also used to read from file
//            'Matrix', sparse matrices from the Matrix package
//            'dist', lower-triangular distance matrices as used in the seriation package,
// Optional arguments:
// data -- a representation of the (similarity or dissimilarity) between pairs of objects.
// zero_epsilon -- a numeric value that determines which entries are considered 0.
//              (default: 1.0e-200)
// symmetric -- a boolean value equal to TRUE is the input matrix is symmetric.
//              (default: TRUE)
// identical_val -- a boolean value equal to TRUE if the data is a 3-columns data frame and symmetric edges are included.
//              (default = TRUE)


// [[Rcpp::export]]
Rcpp::DataFrame
read(SEXP data,
    double zero_epsilon = 1.0e-200,
    bool symmetric = true,
    bool identical_val = false) {
    Rcpp::DataFrame M;
    
    //if data is not symmetric, I do not have to add symmetric edges
    if(symmetric == false)
    {
        if (identical_val == true)
        {
            SFSout << "warning: identical_val is set to TRUE because the matrix is not symmetric." << std::endl;
        }
        identical_val = true;
    }
    
    //read data and create SpMat object
    if(Rf_isMatrix(data)) {
        M = sfs__matrix_to_DataFrame(data,zero_epsilon);
    } else if(Rf_inherits(data, "data.frame")) {
        M = sfs__dataframe_to_DataFrame(data, zero_epsilon, symmetric, identical_val);
    }else {
        // some object
        Rcpp::RObject o(data);
        if (o.inherits("dist")) {
            M = sfs__dist_to_DataFrame(data,zero_epsilon);
        } else {
            std::string c = o.hasAttribute("class")
                    ? "(class " + Rcpp::as<std::string>(o.attr("class")) +")"
                    : "(no class)";
            std::string errmsg = "sfs(): unknown object " + c + " in dissimilarity argument";
            throw std::runtime_error(errmsg);
        }
    }
    
    return M;
}

// exported functionality: solve SFS given a 3-columns (row,col,val) 'data frame'
// Optional arguments:
// matrix -- a 3-columns data frame with symmetric entries, representing the list of all similarities (or dissimilarities) between the pairs of objects to reorder.
// zero_sfs -- a numeric value that determines the threshold for entries of being considered the same.
//              (default: 1.0e-200)
// dissimilarity -- a boolean value equal to TRUE is the input data is a dissimilarity.
//              (default: FALSE)
// Robinsonian -- a boolean value equal to TRUE is one wants to recognize a Robinsonian matrix
//              (default: FALSE)
// num_sweeps -- an integer value that determines how many iterations of SFS shall be repeated.
//              (default: 4)


// [[Rcpp::export]]
arma::Row<int>
sfs(SEXP matrix,
    double sfs_epsilon = 0,
    bool dissimilarity = false,
    bool Robinsonian = false,
    int num_sweeps = 4) {
    
    //check parameters bounds
    if (num_sweeps <= 0)
    {
        throw std::runtime_error("number of sweeps must be strictly bigger than zero");
    }
    if (sfs_epsilon < 0)
    {
        throw std::runtime_error("the SFS-epsilon cannot be negative");
    }
    
    //check if input matrix is a "data frame" with three columns
    if (!Rf_inherits(matrix, "data.frame"))
    {
        throw std::runtime_error("The input must be a data frame.");
    }
    Rcpp::DataFrame M = Rcpp::as<Rcpp::DataFrame>(matrix);
    if (M.size() != 3)
    {
        throw std::runtime_error("The data frame must have 3 columns. Preprocess the data throught the read() function first.");
    }
    
    // prepare batch insertion format for armadillo
    std::vector<int> rowidx = M[0];
    std::vector<int> colidx = M[1];
    std::vector<double> val = M[2];
    arma::umat locations(2,rowidx.size());
    arma::vec  values(rowidx.size());
    for(size_t k=0; k<rowidx.size(); k++) {
        locations(0,k) = rowidx[k];
        locations(1,k) = colidx[k];
        values(k) = val[k];
    }
    
    SFSMatrix::SpMat A = SFSMatrix::SpMat(locations,values,true);
    SFSMatrix S(A,sfs_epsilon, dissimilarity, Robinsonian, num_sweeps); //create objects SFSMatrix
    arma::Row<int> permutation(S.solve()); //run SFS
    
    return permutation;
}
