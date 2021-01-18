//
//  EFT_QP.h
//  eft-mhg
//
//  Created by Joe Song on 4/18/15.
//  Copyright (c) 2015 New Mexico State University. All rights reserved.
//
//  Revision history:
//  2019-02-25 (Hien Nguyen): the file name changed from
//     ExactFunctionalTest.h to EFT_QP.h to distinguish from other
//     implementations of the exact functional test.


//#include <iostream>
#include "define.h"
#include <vector>
#include <cmath>
#include "boost/math/special_functions/factorials.hpp"

using namespace::std;
using namespace::boost::math;

enum LBOUND { LBON = 1, LBOFF = 0 };
enum UBOUND { UB_BY_ROW = 1<<1, UB_BY_ELE = 1<<2, UBOFF = 0 };
enum PVAL {PVAL, ONE_MINUS};

mydouble traverse
(vector<vector<int> > &A, size_t i, size_t j,
 mydouble A_running_stat, mydouble A_running_prob,
 vector<vector<int> > & A_running_rowsums,
 vector<vector<int> > & A_running_colsums,
 const vector<int> & O_rowsums, const vector<int> & O_colsums,
 const mydouble O_stat,
 enum LBOUND lb_method, enum UBOUND ub_method
 );

mydouble funchisq(const vector<vector<int> > & O, mydouble & estimate,
                  const string index_kind);

mydouble funchisq(const vector<vector<int> > & O, const vector<int> & rowsums,
                const vector<int> & colsums, int n);

mydouble upper_bound(const vector<vector<int> > &A, size_t i,
                        const mydouble A_running_stat,
                        const vector<vector<int> > & A_running_rowsums,
                        const vector<vector<int> > & A_running_colsums,
                        const vector<int> & O_rowsums,
                        const vector<int> & O_colsums,
                        mydouble O_stat);

mydouble upper_bound(const vector<vector<int> > &A, size_t i, size_t j,
                   const mydouble A_running_stat,
                   const vector<vector<int> > & A_running_rowsums,
                   const vector<vector<int> > & A_running_colsums,
                   const vector<int> & O_rowsums,
                   const vector<int> & O_colsums,
                   mydouble O_stat);

mydouble lower_bound(const vector<vector<int> > &A, size_t i, size_t j,
                        const mydouble A_running_stat,
                        const vector<vector<int> > & A_running_rowsums,
                        const vector<vector<int> > & A_running_colsums,
                        const vector<int> & O_rowsums,
                        const vector<int> & O_colsums,
                        mydouble O_stat);

mydouble prob_entire_branch(vector<vector<int> > &A,
                               size_t i, size_t j,
                               mydouble A_running_prob,
                               const vector<vector<int> > & A_running_rowsums,
                               const vector<vector<int> > & A_running_colsums,
                               const vector<int> & O_rowsums,
                               const vector<int> & O_colsums);

mydouble enumerate_next (vector<vector<int> > &A,
                            size_t i, size_t j,
                            mydouble A_running_stat,
                            mydouble A_running_prob,
                            vector<vector<int> > & A_running_rowsums,
                            vector<vector<int> > & A_running_colsums,
                            const vector<int> & O_rowsums,
                            const vector<int> & O_colsums,
                            const mydouble O_stat,
                            enum LBOUND lb_method,
                            enum UBOUND ub_method,
                            mydouble (*traverse)
                            (vector<vector<int> > &A,
                             size_t i, size_t j,
                             mydouble A_running_stat,
                             mydouble A_running_prob,
                             vector<vector<int> > & A_running_rowsums,
                             vector<vector<int> > & A_running_colsums,
                             const vector<int> & O_rowsums,
                             const vector<int> & O_colsums,
                             const mydouble O_stat,
                             enum LBOUND lb_method,
                             enum UBOUND ub_method)
                            );

mydouble traverse_ge_observed_stat
                            (vector<vector<int> > &A,
                             size_t i, size_t j,
                             mydouble A_running_stat,
                             mydouble A_running_prob,
                             vector<vector<int> > & A_running_rowsums,
                             vector<vector<int> > & A_running_colsums,
                             const vector<int> & O_rowsums,
                             const vector<int> & O_colsums,
                             const mydouble O_stat,
                             enum LBOUND lb_method,
                             enum UBOUND ub_method
                             );

mydouble traverse_lt_observed_stat
                            (vector<vector<int> > &A,
                             size_t i, size_t j,
                             mydouble A_running_stat,
                             mydouble A_running_prob,
                             vector<vector<int> > & A_running_rowsums,
                             vector<vector<int> > & A_running_colsums,
                             const vector<int> & O_rowsums,
                             const vector<int> & O_colsums,
                             const mydouble O_stat,
                             enum LBOUND lb_method,
                             enum UBOUND ub_method
                             );

mydouble exact_func_test_multi_hypergeometric
                            (const vector<vector<int> > &O, mydouble & fc,
                             enum LBOUND lb_method, enum UBOUND ub_method,
                             enum PVAL pval_method);


