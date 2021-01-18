#include <Rcpp.h>
#include <Eigen/Dense>
#include <Eigen/LU>

// [[Rcpp::depends(BH)]]

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace boost::multiprecision;
using namespace Rcpp;
using namespace Eigen;
using boost::math::constants::pi;

typedef number<cpp_dec_float<100>> mp_float;

// Declare matrix and vector types with multi-precision scalar type
// typedef Eigen::Matrix<mp_float,Dynamic,Dynamic>  MatrixXmp;
// typedef Eigen::Matrix<mp_float,Dynamic,1>        VectorXmp;


template <class T>
T doubleWishart_pfaffian(T xx, int s, T mm, T nn) {
    // Initialize vector
    Eigen::Matrix<T,Dynamic,1> b(s);
    Eigen::Matrix<T,Dynamic,1> tmp(s);
    int d = s + (s % 2);

    // Initialize matrix with zeroes on diagonal
    Eigen::Matrix<T,Dynamic,Dynamic> A(d, d);
    for (int i = 0; i < d; i++) {
        A(i, i) = 0;
    }
    if (s != d) {
        // Fill in extra column
        for (int i = 0; i < s; i++) {
            tmp(i) = boost::math::beta(mm + i + 1, nn + 1, xx);
            A(i, s) = tmp(i);
            A(s, i) = - A(i, s);
        }
    } else {
        for (int i = 0; i < s; i++) {
            tmp(i) = boost::math::beta(mm + i + 1, nn + 1, xx);
        }
    }
    if (s != 1) {
        for (int i = 0; i < s; i++) {
            b(i) = 0.5 * tmp(i) * tmp(i);

            for (int j = i; j < (s-1); j++) {
                b(j+1) = ((mm+j+1)*b(j) - boost::math::beta(2*mm+i+j+2, 2*nn+2, xx))/(mm+j+nn+2);
                A(i,j+1) = tmp(i) * tmp(j+1) -
                    2*b(j+1);
                A(j+1, i) = - A(i, j+1);
            }
            Rcpp::checkUserInterrupt();
        }
    }
    T det;
    if(d > 4) {
        det = A.fullPivLu().determinant();
    } else {
        det = A.determinant();
    }
    return sqrt(det);
}

template <class T>
T doubleWishart_const(int s, T mm, T nn) {
    // Compute scaling constant
    T C1 = 0;
    for (int i = 1; i <= s; i++) {
        C1 += boost::math::lgamma(0.5*(i+2*mm+2*nn+s+2)) - boost::math::lgamma(0.5*i) -
            boost::math::lgamma(0.5*(i+2*mm+1)) - boost::math::lgamma(0.5*(i+2*nn+1));
        Rcpp::checkUserInterrupt();
    }
    T C = pow(pi<T>(), 0.5 * s) * exp(C1);

    return C;
}

// [[Rcpp::export]]
NumericVector doubleWishart_raw(NumericVector x, int s, double m, double n, bool mp) {
    int size = x.size();
    NumericVector result(size);
    if (mp) {
        mp_float mm(m);
        mp_float nn(n);
        mp_float value, xx;
        mp_float constant = doubleWishart_const(s, mm, nn);
        for(int i = 0; i < size; i++) {
            xx = mp_float(x[i]);
            value = constant * doubleWishart_pfaffian(xx, s, mm, nn);
            result[i] = value.convert_to<double>();
            Rcpp::checkUserInterrupt();
        }
    } else {
        double constant = doubleWishart_const(s, m, n);
        for(int i = 0; i < size; i++) {
            result[i] = constant * doubleWishart_pfaffian(x[i], s, m, n);
            Rcpp::checkUserInterrupt();
        }
    }

    return result;
}

