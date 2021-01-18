//
//  EFT_QP.cpp
//    Exact functional test implementaiton using quadratic programming
//
//  eft-mhg
//
//  Created by Joe Song on 4/18/15.
//  Copyright (c) 2015 New Mexico State University. All rights reserved.
//
//  Revision history:
//  2019-02-25 (Hien Nguyen): the file name changed from
//     ExactFunctionalTest.cpp to EFT_QP.cpp to distinguish from other
//     implementations of the exact functional test.
//
//  6/28/2015 (Hua Zhong): Increased float comparison precision
//     MS 6/16/2018. Improved numerical precision by using an alternative
//       form of FunChisq definition.

//#include "ExactFunctionalTest.h"
#include "EFT_QP.h"
// #include <Rcpp.h>//for deleted
// using namespace Rcpp;//for deleted

//Hua added, Jun 28,2015.
//To make float comparison more precise
bool is_close(const mydouble & a, const mydouble & b, mydouble epsilon=1e-6){ // a==b?
    //Rcpp::Rcout<<a<<'\t'<<b<<'\t'<<fabs((a - b)/a)<<'\t'<<fabs((a - b)/b)<<'\t'<<((fabs((a - b)/a) < epsilon && fabs((a - b)/b) < epsilon) ? true : false)<<endl;
    return (fabs((a - b)/a) < epsilon && fabs((a - b)/b) < epsilon) ? true : false;
}

bool ge (const mydouble & a, const mydouble & b){ // a>=b?
    if((a-b) >= 0 || is_close(a,b)) return true;
    return false;
}

bool le (const mydouble & a, const mydouble & b){ // a<=b?
    if((a-b) <= 0 || is_close(a,b)) return true;
    return false;
}

bool gg (const mydouble & a, const mydouble & b){ // a>b?
    if((a-b) > 0 && !is_close(a,b)) return true;
    return false;
}

bool ll (const mydouble & a, const mydouble & b){ // a<b?
    if((a-b) < 0 && !is_close(a,b)) return true;
    return false;
}
////

mydouble upper_bound_EXPERIMENTAL(const vector<vector<int> > &A, size_t i,
                                  const mydouble A_running_stat,
                                  const vector<vector<int> > & A_running_rowsums,
                                  const vector<vector<int> > & A_running_colsums,
                                  const vector<int> & O_rowsums,
                                  const vector<int> & O_colsums,
                                  mydouble O_stat)
{
    // Version: EXPERIMENTAL

    // return 10e10;
    mydouble upper_bound = A_running_stat;

    vector<int> U(O_colsums);
    size_t ncols = A[0].size(); // Number of colmns
    size_t nrows = A.size(); // Number of rows

    if (i > 0) {

        int maxU = 0;

        for (size_t q=0; q<ncols; ++q) {
            U[q] = O_colsums[q] - A_running_colsums[i-1][q];
            if(U[q] > maxU) {
                maxU = U[q];
            }
        }
        if(1) {
            int ub = 0;
            for (size_t l=i; l<nrows; ++l) {
                if(O_rowsums[l] > maxU) {
                    ub += maxU;
                } else {
                    ub += O_rowsums[l];
                }
            }
            ub *= ncols;
            if(ll(upper_bound + ub, O_stat)) {
                return upper_bound + ub; // Not a tight upperbound. This will speed up 30\% on 5x4 matrix
            }
        }
    } /* else {
       for (size_t q=0; q<ncols; ++q) {
       if(O_colsums[q] > maxU) {
       maxU = O_colsums[q];
       }
       }
       } */

    vector<size_t> order( ncols );
    for (size_t q=0; q<ncols; ++q) {
        order[q] = q;
    }

    // sort U in decreasing order
    sort(order.begin(), order.end(),
         [&U](size_t i1, size_t i2) {return U[i1] > U[i2];});

    if(1) {
        if(nrows - i < ncols) {
            int ub = 0;
            for (size_t k=0; k<nrows - i; ++k) {
                ub += U[order[k]];
            }
            ub *= ncols;
            if(ll(upper_bound + ub, O_stat)) {
                return upper_bound + ub; // Not a tight upperbound.
            }
        }
    }

    for (size_t l=i; l<nrows; ++l) {
        // find lower bound for row l
        int runsum = 0;
        mydouble el = O_rowsums[l] / (mydouble) ncols;
        for (size_t k=0; k<ncols; ++k) {
            // accumulate the lower bound
            int xmax = O_rowsums[l] - runsum;
            if (U[order[k]] < xmax) {
                mydouble d = U[order[k]] - el;
                if(el>0) upper_bound += d * d / el;
                runsum += U[order[k]];
            } else if (xmax != 0) {
                mydouble d = xmax - el;
                if(el>0) upper_bound += d * d / el;
                runsum += xmax;
            } else {
                upper_bound += el * (ncols - k);
                break;
            }
            if (gg(upper_bound, O_stat)) {
                return upper_bound;
            }
        }
    }

    return upper_bound;

}

mydouble upper_bound_0(const vector<vector<int> > &A, size_t i,
                     const mydouble A_running_stat,
                     const vector<vector<int> > & A_running_rowsums,
                     const vector<vector<int> > & A_running_colsums,
                     const vector<int> & O_rowsums,
                     const vector<int> & O_colsums,
                     mydouble O_stat)
{   // Version 0:

    // return 10e10;
    mydouble upper_bound = A_running_stat;

    vector<int> U(O_colsums);

    size_t ncols = A[0].size(); // Number of colmns
    size_t nrows = A.size(); // Number of rows
    if (i > 0) {
        for (size_t q=0; q<ncols; ++q) {
            U[q] = O_colsums[q] - A_running_colsums[i-1][q];
        }
    }

    vector<size_t> order( ncols );
    for (size_t q=0; q<ncols; ++q) {
        order[q] = q;
    }

    // sort U in decreasing order
    sort(order.begin(), order.end(),
         [&U](size_t i1, size_t i2) {return U[i1] > U[i2];});

    for (size_t l=i; l<nrows; ++l) {
        // find lower bound for row l
        int runsum = 0;

        mydouble el = O_rowsums[l] / (mydouble) ncols;
        for (size_t k=0; k<ncols; ++k) {
            // accumulate the lower bound
            int xmax = O_rowsums[l] - runsum;
            if (U[order[k]] < xmax) {
                mydouble d = U[order[k]] - el;
                if(el>0) upper_bound += d * d / el;
                runsum += U[order[k]];
            } else if (xmax != 0) {
                mydouble d = xmax - el;
                if(el>0) upper_bound += d * d / el;
                runsum += xmax;
            } else {
                upper_bound += el * (ncols - k);
                break;
            }
            if (gg(upper_bound, O_stat)) {
                return upper_bound;
            }
        }
    }

    return upper_bound;

}

mydouble upper_bound(const vector<vector<int> > &A, size_t i,
                       const mydouble A_running_stat,
                       const vector<vector<int> > & A_running_rowsums,
                       const vector<vector<int> > & A_running_colsums,
                       const vector<int> & O_rowsums,
                       const vector<int> & O_colsums,
                       mydouble O_stat)
{   // Version 1: MS 6/16/2018. Improved numerical precision

    // return 10e10;
    mydouble upper_bound = A_running_stat;

    vector<int> U(O_colsums);

    size_t ncols = A[0].size(); // Number of colmns
    size_t nrows = A.size(); // Number of rows
    if (i > 0) {
        for (size_t q=0; q<ncols; ++q) {
            U[q] = O_colsums[q] - A_running_colsums[i-1][q];
        }
    }

    vector<size_t> order( ncols );
    for (size_t q=0; q<ncols; ++q) {
        order[q] = q;
    }

    // sort U in decreasing order
    sort(order.begin(), order.end(),
         [&U](size_t i1, size_t i2) {return U[i1] > U[i2];});

    for (size_t l=i; l<nrows; ++l) {
        // find upper bound for row l
        if(O_rowsums[l] > 0) {
            int runsum = 0;
            for (size_t k=0; k<ncols; ++k) {
                // accumulate the upper bound
                int xmax = O_rowsums[l] - runsum;
                if (U[order[k]] < xmax) {
                    upper_bound += ncols * U[order[k]] * U[order[k]]
                        / (mydouble) O_rowsums[l];
                    runsum += U[order[k]];
                } else if (xmax != 0) {
                    upper_bound += ncols * xmax * xmax / (mydouble) O_rowsums[l];
                    runsum += xmax;
                } else {
                    break;
                }
                if (gg(upper_bound, O_stat)) {
                    return upper_bound;
                }
            }
        }
    }

    return upper_bound;

}


mydouble upper_bound(const vector<vector<int> > &A, size_t i, size_t j,
                     const mydouble A_running_stat,
                     const vector<vector<int> > & A_running_rowsums,
                     const vector<vector<int> > & A_running_colsums,
                     const vector<int> & O_rowsums,
                     const vector<int> & O_colsums,
                     mydouble O_stat)
{
    // return 10e10;
    mydouble upper_bound = A_running_stat;

    vector<int> U(O_colsums);

    size_t ncols = A[0].size(); // Number of colmns
    size_t nrows = A.size(); // Number of rows
    if (i > 0) {
        for (size_t q=0; q<ncols; ++q) {
            if (q < j) {
                U[q] = O_colsums[q] - A_running_colsums[i][q];
            } else {
                U[q] = O_colsums[q] - A_running_colsums[i-1][q];
            }
        }

        /*
         for (size_t q=0; q<ncols; ++q) {
         U[q] = O_colsums[q] - A_running_colsums[i-1][q];
         } */
    }

    vector<size_t> order( ncols );
    for (size_t q=0; q<ncols; ++q) {
        order[q] = q;
    }

    // sort U in decreasing order
    sort(order.begin(), order.end(),
         [&U](size_t i1, size_t i2) {return U[i1] > U[i2];});

    for (size_t l=i; l<nrows; ++l) {
        // find lower bound for row l

        if (l == i && j > 0) {

            vector<size_t> ord( ncols );
            for (size_t q=0; q<ncols; ++q) {
                ord[q] = q;
            }

            // sort U in decreasing order
            sort(ord.begin()+j, ord.end(),
                 [&U](size_t i1, size_t i2) {return U[i1] > U[i2];});

            // find upper bound for row l
            int runsum = A_running_rowsums[i][j-1];
            mydouble el = O_rowsums[l] / (mydouble) ncols;
            for (size_t k=j; k<ncols; ++k) {
                // accumulate the upper bound
                int xmax = O_rowsums[l] - runsum;
                if (U[order[k]] < xmax) {
                    mydouble d = U[order[k]] - el;
                    if(el>0) upper_bound += d * d / el;
                    runsum += U[order[k]];
                } else if (xmax != 0) {
                    mydouble d = xmax - el;
                    if(el>0) upper_bound += d * d / el;
                    runsum += xmax;
                } else {
                    upper_bound += el * (ncols - k);
                    break;
                }
                if (ge(upper_bound, O_stat)) {
                    return upper_bound;
                }
            }
            continue;
        }

        int runsum = 0;
        mydouble el = O_rowsums[l] / (mydouble) ncols;
        for (size_t k=0; k<ncols; ++k) {
            // accumulate the lower bound
            int xmax = O_rowsums[l] - runsum;
            if (U[order[k]] < xmax) {
                mydouble d = U[order[k]] - el;
                if(el>0) upper_bound += d * d / el;
                runsum += U[order[k]];
            } else if (xmax != 0) {
                mydouble d = xmax - el;
                if(el>0) upper_bound += d * d / el;
                runsum += xmax;
            } else {
                upper_bound += el * (ncols - k);
                break;
            }
            if (ge(upper_bound, O_stat)) {
                return upper_bound;
            }
        }
    }

    return upper_bound;

}

mydouble lower_bound_0(const vector<vector<int> > &A, size_t i, size_t j,
                     const mydouble A_running_stat,
                     const vector<vector<int> > & A_running_rowsums,
                     const vector<vector<int> > & A_running_colsums,
                     const vector<int> & O_rowsums,
                     const vector<int> & O_colsums,
                     mydouble O_stat)
{
    mydouble lower_bound = A_running_stat;

    vector<int> U(O_colsums);

    size_t nrows = A.size();
    size_t ncols = A[0].size();

    if (i > 0) {
        for (size_t q=0; q<ncols; ++q) {
            if (q < j) {
                U[q] = O_colsums[q] - A_running_colsums[i][q];
            } else {
                U[q] = O_colsums[q] - A_running_colsums[i-1][q];
            }
        }
    }

    vector<size_t> order( ncols );
    for (size_t q=0; q<ncols; ++q) {
        order[q] = q;
    }

    // sort U in increasing order
    sort(order.begin(), order.end(),
         [&U](size_t i1, size_t i2) {return U[i1] < U[i2];});

    for (size_t l=i; l<nrows; ++l) {
        // find lower bound for row l
        int runsum = 0;
        mydouble e = O_rowsums[l] / (mydouble) ncols;
        for (size_t k=0; k<ncols; ++k) {
            // accumulate the lower bound
            mydouble xavg = (O_rowsums[l]-runsum)/(mydouble)(ncols-k);
            if (U[order[k]] < xavg) {
                if(e>0) lower_bound += (U[order[k]] - e) * (U[order[k]] - e) / e;
                runsum += U[order[k]];
            } else {
                if(e>0) lower_bound += ((xavg - e) * (xavg - e) / e) * (ncols-k);
                break;
            }
            if (ge(lower_bound, O_stat)) {
                return lower_bound;
            }
        }
    }

    if(lower_bound < 0.0) lower_bound = 0.0;

    return lower_bound;
}

mydouble lower_bound(const vector<vector<int> > &A, size_t i, size_t j,
                     const mydouble A_running_stat,
                     const vector<vector<int> > & A_running_rowsums,
                     const vector<vector<int> > & A_running_colsums,
                     const vector<int> & O_rowsums,
                     const vector<int> & O_colsums,
                     mydouble O_stat)
{   // version 1: MS 6/16/2018. Improved numerical precision

    mydouble lower_bound = A_running_stat;

    vector<int> U(O_colsums);

    size_t nrows = A.size();
    size_t ncols = A[0].size();

    if (i > 0) {
        for (size_t q=0; q<ncols; ++q) {
            if (q < j) {
                U[q] = O_colsums[q] - A_running_colsums[i][q];
            } else {
                U[q] = O_colsums[q] - A_running_colsums[i-1][q];
            }
        }
    }

    vector<size_t> order( ncols );
    for (size_t q=0; q<ncols; ++q) {
        order[q] = q;
    }

    // sort U in increasing order
    sort(order.begin(), order.end(),
         [&U](size_t i1, size_t i2) {return U[i1] < U[i2];});

    for (size_t l=i; l<nrows; ++l) {
        // find lower bound for row l
        if(O_rowsums[l] > 0) {
            int runsum = 0;
            for (size_t k=0; k<ncols; ++k) {
                // accumulate the lower bound
                mydouble xavg = (O_rowsums[l]-runsum) / (mydouble)(ncols-k);
                if (U[order[k]] < xavg) {
                    lower_bound += ncols * U[order[k]] * U[order[k]] / (mydouble) O_rowsums[l];
                    runsum += U[order[k]];
                } else {
                    lower_bound += (ncols-k) * ncols * xavg * xavg / (mydouble) O_rowsums[l];
                    break;
                }
                if (ge(lower_bound, O_stat)) {
                    return lower_bound;
                }
            }
        }
    }

    return lower_bound;
}
mydouble prob_entire_branch(vector<vector<int> > &A,
                            size_t i, size_t j,
                            mydouble A_running_prob,
                            const vector<vector<int> > & A_running_rowsums,
                            const vector<vector<int> > & A_running_colsums,
                            const vector<int> & O_rowsums,
                            const vector<int> & O_colsums)
{
    mydouble prob = A_running_prob;

    if(j != 0) throw "ERROR: can only compute whole rows";

    if (i == 0) {
        prob = 1.0;
    } else {
        int remaing_rowsum=0;
        size_t nrows = A.size();
        size_t ncols = A[0].size();

        for (size_t l=i; l<nrows; ++l) {
            prob /= factorial<mydouble>(O_rowsums[l]);
            remaing_rowsum += O_rowsums[l];
        }
        prob *= factorial<mydouble>(remaing_rowsum);

        for(size_t q=0; q<ncols; q++) {
            prob /= factorial<mydouble>(O_colsums[q] - A_running_colsums[i-1][q]);
        }
    }
    return prob;
}

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
                         )
{
    mydouble prob=0.0;

    size_t nrows = A.size();
    size_t ncols = A[0].size();

    size_t j_next = j+1;
    size_t i_next = i;
    if(j_next == ncols) {
        j_next = 0;
        i_next += 1;
    }

    int Lij=0;
    //int Uij;

    if (i == nrows - 1) { // last row

        Lij = O_colsums[j] - A_running_colsums[i-1][j];

    } else if (j == ncols - 1) { // last column

        Lij = O_rowsums[i] - A_running_rowsums[i][j-1];

    } /*else {

       Lij = 0;

       }*/

    int Uij = min(O_rowsums[i] - (j > 0 ? A_running_rowsums[i][j-1]:0),
                  O_colsums[j] - (i > 0 ? A_running_colsums[i-1][j]:0));


    // mydouble eij = O_rowsums[i] / (mydouble) ncols;

    // mydouble A_running_stat_before = A_running_stat;

    for(int x=Lij; x <= Uij; x++) {

        A[i][j] = x;

        if (A[i][j] == Lij) {
            A_running_prob /= factorial<mydouble>( Lij );
        } else {
            A_running_prob /= A[i][j];
        }

        // update running statistics
        A_running_rowsums[i][j] = (j > 0 ? A_running_rowsums[i][j-1] : 0) + A[i][j];
        A_running_colsums[i][j] = (i > 0 ? A_running_colsums[i-1][j] : 0) + A[i][j];

        mydouble stat_ij = 0;

        /* if(0) {
            mydouble d = A[i][j] - eij;
            if(eij>0) stat_ij = d * d / eij;
        } else */
        {
            if(O_rowsums[i] > 0) {
                stat_ij = ncols * A[i][j] * A[i][j] / (mydouble) O_rowsums[i];
            }
        }

        // A_running_stat += stat_ij;

        prob += (*traverse)
        (A, i_next, j_next, A_running_stat + stat_ij, A_running_prob,
         A_running_rowsums, A_running_colsums,
         O_rowsums, O_colsums, O_stat, lb_method, ub_method);

        // BEGIN ***
        // NUMERICALLY VERY CRITICAL. Must use option 1 instead of 2
        // Option 1:
        //   A_running_stat = A_running_stat_before;
        // Option 2: A_running_stat -= stat_ij;
        // Option 3: No need to change A_running_stat

        // END ***

        // A_running_prob *= fac_ij;

        // restore running statistics
        A_running_rowsums[i][j] = 0; // -= A[i][j];
        A_running_colsums[i][j] = 0; // -= A[i][j];
        A[i][j] = 0;

    }

    return prob;
}

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
 )
{
    mydouble prob=0.0;

    size_t nrows = A.size();

    if(i >= nrows) { // calculate probability of A
        prob = ge(A_running_stat, O_stat) ? A_running_prob : 0.0;
        //prob = (A_running_stat >= O_stat || is_close(A_running_stat, O_stat)) ? A_running_prob : 0.0;
        //Rcpp::Rcout<<A_running_stat<<'\t'<<O_stat<<'\t'<<A_running_prob<<'\t'<<prob<<'\t'<<(A_running_stat >= O_stat)<<'\t'<<is_close(A_running_stat, O_stat)<<std::endl;
        /*    } else if (ub_method == UB_BY_ELE // && j == 0
         && upper_bound(A, i, j, A_running_stat,
         A_running_rowsums, A_running_colsums,
         O_rowsums, O_colsums, O_stat)
         < O_stat)
         { // check the upper bound on stat

         // cout << "Skip entire branch" << endl;
         //cout << "0";
         prob = 0.0;
         */
    } else if (ub_method == UB_BY_ROW && j == 0
               && ll(upper_bound(A, i, A_running_stat,
                                 A_running_rowsums, A_running_colsums,
                                 O_rowsums, O_colsums, O_stat)
                     , O_stat))
    { // check the upper bound on stat

        // cout << "Skip entire branch" << endl;
        //cout << "0";
        prob = 0.0;


    } else if (lb_method == LBON && j == 0 && ge(A_running_stat, O_stat)) {

        // cout << "Keep entire branch" << endl;
        prob = prob_entire_branch(A, i, j, A_running_prob, A_running_rowsums,
                                  A_running_colsums, O_rowsums, O_colsums);

        //cout << "e";
    } else if (lb_method == LBON && j == 0 &&
               ge(lower_bound(A, i, j, A_running_stat,
                              A_running_rowsums, A_running_colsums,
                              O_rowsums, O_colsums, O_stat)
                  , O_stat))
    { // check lower bound on stat

        // cout << "Keep entire branch" << endl;
        prob = prob_entire_branch(A, i, j, A_running_prob, A_running_rowsums,
                                  A_running_colsums, O_rowsums, O_colsums);

        //cout << "+";

    } else {

        prob = enumerate_next(A, i, j, A_running_stat, A_running_prob,
                              A_running_rowsums, A_running_colsums,
                              O_rowsums, O_colsums, O_stat,
                              lb_method, ub_method,
                              & traverse_ge_observed_stat);
    }
    return prob;
}

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
 )
{
    mydouble prob=0.0;

    size_t nrows = A.size();

    if(i >= nrows) { // calculate probability of A

        prob = ll(A_running_stat, O_stat) ? A_running_prob : 0.0;

        /*    } else if (ub_method == UB_BY_ELE // && j == 0
         && upper_bound(A, i, j, A_running_stat,
         A_running_rowsums, A_running_colsums,
         O_rowsums, O_colsums, O_stat)
         < O_stat) { // check the upper bound on stat

         // cout << "Keep entire branch" << endl;
         prob = prob_entire_branch(A, i, j, A_running_prob, A_running_rowsums,
         A_running_colsums, O_rowsums, O_colsums);
         */
    } else if (ub_method == UB_BY_ROW && j == 0
               && ll(upper_bound(A, i, A_running_stat,
                                 A_running_rowsums, A_running_colsums,
                                 O_rowsums, O_colsums, O_stat)
                     , O_stat)) { // check the upper bound on stat

                   // cout << "Keep entire branch" << endl;
                   // prob = 0.0;
                   prob = prob_entire_branch(A, i, j, A_running_prob, A_running_rowsums,
                                             A_running_colsums, O_rowsums, O_colsums);

               } else if (lb_method == LBON && j == 0 && ge(A_running_stat, O_stat)) {

                   // cout << "Skip entire branch" << endl;
                   prob = 0.0;

                   //cout << "e";
               } else if (lb_method == LBON && j == 0 &&
                          ge(lower_bound(A, i, j, A_running_stat,
                                         A_running_rowsums, A_running_colsums,
                                         O_rowsums, O_colsums, O_stat)
                             , O_stat)) { // check lower bound on stat

                              // cout << "Skip entire branch" << endl;
                              prob = 0.0;

                          } else {

                              prob = enumerate_next(A, i, j, A_running_stat, A_running_prob,
                                                    A_running_rowsums, A_running_colsums,
                                                    O_rowsums, O_colsums, O_stat,
                                                    lb_method, ub_method,
                                                    & traverse_lt_observed_stat);
                          }

    return prob;
}


mydouble exact_func_test_multi_hypergeometric
(const vector<vector<int> > &O, mydouble & fc,
 enum LBOUND lb_method, enum UBOUND ub_method,
 enum PVAL pval_method)
{
    mydouble pval=1.0;
    fc = 0;

    size_t nrows = O.size();
    size_t ncols = O[0].size();

    if (nrows < 2 || ncols < 2) {
        return pval;
    }

    int n=0;

    vector<vector<int> > A(nrows, vector<int>(ncols, 0));
    mydouble A_running_stat = 0.0;
    mydouble A_running_prob = 1.0;
    vector<vector<int> > A_running_rowsums(A);
    vector<vector<int> > A_running_colsums(A);

    // Calculate row and col sums, sample size, and initial probability
    vector<int> O_rowsums(nrows, 0);
    vector<int> O_colsums(ncols, 0);

    for (size_t i=0; i<nrows; ++i) {
        for (size_t j=0; j<ncols; ++j) {
            O_rowsums[i] += O[i][j];
        }
        A_running_prob *= factorial<mydouble>(O_rowsums[i]);
        n += O_rowsums[i];
    }

    // mydouble ej = n / (mydouble) ncols;
    
    for (size_t j=0; j<ncols; ++j) {
        for (size_t i=0; i<nrows; ++i) {
            O_colsums[j] += O[i][j];
        }
        A_running_prob *= factorial<mydouble>(O_colsums[j]);

        /* if(0) {
            if(ej>0) A_running_stat -= (O_colsums[j]-ej) * (O_colsums[j]-ej) / ej;
        } else */ {
            if(n>0) A_running_stat -= ncols * O_colsums[j] * O_colsums[j] / (mydouble) n;
        }
    }

    A_running_prob /= factorial<mydouble>(n);

    fc = funchisq(O, O_rowsums, O_colsums, n);

    switch (pval_method) {
        case PVAL:
            pval = traverse_ge_observed_stat
            (A, 0, 0, A_running_stat, A_running_prob,
             A_running_rowsums, A_running_colsums,
             O_rowsums, O_colsums, fc,
             lb_method, ub_method);

            break;

        case ONE_MINUS:
            pval = 1 - traverse_lt_observed_stat
            (A, 0, 0, A_running_stat, A_running_prob,
             A_running_rowsums, A_running_colsums,
             O_rowsums, O_colsums, fc,
             lb_method, ub_method);

            break;

        default:
            break;
    }

    return pval;
}
