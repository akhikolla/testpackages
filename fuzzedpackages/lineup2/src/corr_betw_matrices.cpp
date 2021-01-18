// correlations between columns of two matrices
#include "fscale.h"
#include "corr_betw_matrices.h"
#include <Rcpp.h>

using namespace Rcpp;

// calculate correlation between columns of x and corresponding columns of y
//
// [[Rcpp::export(".corr_betw_matrices_paired")]]
NumericVector corr_betw_matrices_paired(const NumericMatrix& x, const NumericMatrix& y)
{
    const int n_row = x.rows();
    const int n_col = x.cols();
    if(n_row != y.rows() || n_col != y.cols())
        throw std::invalid_argument("dim(x) != dim(y)");

    NumericVector result(n_col);

    for(int j=0; j<n_col; j++) {
        checkUserInterrupt();  // check for ^C from user

        double sum=0.0;
        int count=0;

        // delicacy regarding scaling... need to omit x values where y is NA and vice versa
        NumericVector xs = fscale(x(_,j), y(_,j));
        NumericVector ys = fscale(y(_,j), x(_,j));

        for(int i=0; i<n_row; i++) {
            if(R_finite(xs[i]) && R_finite(ys[i])) {
                sum += (xs[i] * ys[i]);
                count++;
            }
        }
        if(count > 1) result[j] = sum/(double)(count-1);
        else result[j] = NA_REAL;
    }

    return result;
}

// for each column of left matrix, find the column in the right matrix
// with the highest correlation
//
// [[Rcpp::export(".corr_betw_matrices_unpaired_bestright")]]
List corr_betw_matrices_unpaired_bestright(const NumericMatrix& x,
                                           const NumericMatrix& y)
{
    const int n_row = x.rows();
    if(y.rows() != n_row)
        throw std::invalid_argument("nrow(x) != nrow(y)");
    const int n_col_x = x.cols();
    const int n_col_y = y.cols();

    // to contain the results
    NumericVector corr(n_col_x);
    IntegerVector yindex(n_col_x);

    for(int xcol=0; xcol < n_col_x; xcol++) {
        checkUserInterrupt();  // check for ^C from user

        double best_corr = -2;
        int best_index = NA_INTEGER;

        for(int ycol=0; ycol < n_col_y; ycol++) {
            // delicacy regarding scaling... need to omit x values where y is NA and vice versa
            NumericVector xs = fscale(x(_,xcol), y(_,ycol));
            NumericVector ys = fscale(y(_,ycol), x(_,xcol));

            double sum=0.0;
            int count=0;
            for(int i=0; i<n_row; i++) {
                if(R_finite(xs[i]) && R_finite(ys[i])) {
                    sum += (xs[i] * ys[i]);
                    count++;
                }
            }
            if(count > 1) {
                sum /= (double)(count-1);
                if(sum > best_corr) {
                    best_corr = sum;
                    best_index = ycol+1;
                }
            }

        } /* end loop over col of y */
        if(best_corr < -1.0) {
            corr[xcol] = NA_REAL;
            yindex[xcol] = NA_INTEGER;
        }
        else {
            corr[xcol] = best_corr;
            yindex[xcol] = best_index;
        }

    } /* end loop over col of x */

    // combine the results into a list
    return List::create(Named("corr")=corr,
                        Named("yindex")=yindex);

}




// return correlations between column of left matrix and column of right matrix
// that exceed corr_threshold
//
// [[Rcpp::export(".corr_betw_matrices_unpaired_bestpairs")]]
List corr_betw_matrices_unpaired_bestpairs(const NumericMatrix& x,
                                           const NumericMatrix& y,
                                           const double corr_threshold)
{
    const int n_row = x.rows();
    if(y.rows() != n_row)
        throw std::invalid_argument("nrow(x) != nrow(y)");
    const int n_col_x = x.cols();
    const int n_col_y = y.cols();

    // to contain the results
    std::vector<double> corr;
    std::vector<int> xindex;
    std::vector<int> yindex;

    for(int xcol=0; xcol < n_col_x; xcol++) {
        checkUserInterrupt();  // check for ^C from user
        for(int ycol=0; ycol < n_col_y; ycol++) {
            // delicacy regarding scaling... need to omit x values where y is NA and vice versa
            NumericVector xs = fscale(x(_,xcol), y(_,ycol));
            NumericVector ys = fscale(y(_,ycol), x(_,xcol));

            double sum = 0.0;
            int count = 0;

            for(int i=0; i<n_row; i++) {
                if(R_finite(xs[i]) && R_finite(ys[i])) {
                    sum += (xs[i] * ys[i]);
                    count++;
                }
            }
            if(count > 1) {
                sum /= (double)(count-1);
                if(sum >= corr_threshold) {
                    corr.push_back(sum);
                    xindex.push_back(xcol+1);
                    yindex.push_back(ycol+1);
                }
            }

        } /* end loop over col of y */
    } /* end loop over col of x */

    // combine the results into a list
    return List::create(Named("corr")=corr,
                        Named("xindex")=xindex,
                        Named("yindex")=yindex);
}

// calculate full set of correlations between columns of x and columns of y
//
// [[Rcpp::export(".corr_betw_matrices_unpaired_all")]]
NumericMatrix corr_betw_matrices_unpaired_all(const NumericMatrix& x,
                                              const NumericMatrix& y)
{
    const int n_row = x.rows();
    if(y.rows() != n_row)
        throw std::invalid_argument("nrow(x) != nrow(y)");
    const int n_col_x = x.cols();
    const int n_col_y = y.cols();

    NumericMatrix result(n_col_x, n_col_y);

    for(int ycol=0; ycol<n_col_y; ycol++) {
        checkUserInterrupt();  // check for ^C from user

        for(int xcol=0; xcol<n_col_x; xcol++) {
            // delicacy regarding scaling... need to omit x values where y is NA and vice versa
            NumericVector xs = fscale(x(_,xcol), y(_,ycol));
            NumericVector ys = fscale(y(_,ycol), x(_,xcol));

            double sum=0.0;
            int count = 0;

            for(int i=0; i<n_row; i++) {
                if(R_finite(xs[i]) && R_finite(ys[i])) {
                    sum += (xs[i] * ys[i]);
                    count++;
                }
            }

            if(count > 1) result(xcol, ycol) = sum/(double)(count-1);
            else result(xcol, ycol) = NA_REAL;

        } /* end loop over col of y */
    } /* end loop over col of x */

    return result;
}
