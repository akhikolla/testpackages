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

template <class T>
T mgamma_C(T, int, bool);

template <class T>
T singleWishart_pfaffian(T xx, int n_min, int n_max) {
    // Declare vector types with multi-precision scalar type
    typedef Eigen::Matrix<T,Dynamic,1> VectorXmp;

    T b;
    int n_mat = n_min + (n_min % 2);
    T alpha = 0.5*(n_max - n_min - 1);

    // Initialize matrix with zeroes on diagonal
    Eigen::Matrix<T,Dynamic,Dynamic> A(n_mat, n_mat);
    for (int i = 0; i < n_mat; i++) {
        A(i, i) = 0;
    }

    // Pre-compute some values
    VectorXmp p(n_min);
    VectorXmp g(n_min);
    VectorXmp q(2*n_min - 2);

    if (n_min != n_mat) {
        // Fill in extra column
        VectorXmp alpha_vec(n_min);
        T alpha_last = 0.5*(n_max + n_min + 1);
        for (int i = 0; i < n_min; i++) {
            alpha_vec(i) = alpha + i + 1;
            A(i, n_min) = pow(2, -alpha_last)*boost::math::gamma_p(alpha_vec(i), 0.5*xx)/boost::math::tgamma(alpha_last);
            A(n_min, i) = - A(i, n_min);
        }
    }

    for (int i = 0; i < n_min; i++) {
        p(i) = boost::math::gamma_p(alpha + i + 1, 0.5*xx);
        g(i) = boost::math::tgamma(alpha + i + 1);
    }
    for (int i = 0; i < (2*n_min - 2); i++) {
        q(i) = pow(0.5, 2*alpha + i + 2) * boost::math::tgamma(2*alpha+i+2)*boost::math::gamma_p(2*alpha + i + 2, xx);
    }

    for (int i = 0; i < n_min; i++) {
        b = 0.5*p(i)*p(i);
        for(int j = i; j < (n_min - 1); j++) {
            b -= q(i+j)/(g(i)*g(j+1));
            A(i, j+1) = p(i)*p(j+1) - 2*b;
            A(j+1, i) = -A(i, j+1);
        }
        Rcpp::checkUserInterrupt();
    }

    // Compute Pfaffian
    T det;
    if(n_mat > 4) {
        det = A.fullPivLu().determinant();
    } else {
        det = A.determinant();
    }
    return sqrt(det);
}

double singleWishart_constDP(int n_min, int n_max) {
    int n_mat = n_min + (n_min % 2);
    double alpha = 0.5*(n_max - n_min - 1);
    // Compute constant
    double K1 = pow(M_PI, 0.5*n_min*n_min);
    K1 /= pow(2, 0.5*n_min*n_max)*mgamma_C(0.5*n_max, n_min, false)*mgamma_C(0.5*n_min, n_min, false);
    double K2 = pow(2, alpha*n_mat+0.5*n_mat*(n_mat+1));
    for (int k = 0; k < n_mat; k++) {
        K2 *= boost::math::tgamma(alpha + k + 1);
    }
    return K1*K2;
}

mp_float singleWishart_constMP(int n_min, int n_max) {
    int n_mat = n_min + (n_min % 2);
    mp_float alpha = 0.5*(n_max - n_min - 1);
    // Compute constant
    mp_float K1 = pow(pi<mp_float>(), 0.5*n_min*n_min);
    K1 /= pow(2, mp_float(0.5*n_min*n_max))*mgamma_C(mp_float(0.5*n_max), n_min, false)*mgamma_C(mp_float(0.5*n_min), n_min, false);
    mp_float K2 = pow(2, alpha*n_mat+0.5*n_mat*(n_mat+1));
    for (int k = 0; k < n_mat; k++) {
        K2 *= boost::math::tgamma(alpha + k + 1);
    }
    return K1*K2;
}

// [[Rcpp::export]]
NumericVector singleWishart_raw(NumericVector x, int n_min, int n_max, bool mp) {
    int n = x.size();
    NumericVector result(n);
    if (mp) {
        mp_float constant = singleWishart_constMP(n_min, n_max);
        mp_float value, xx;
        for (int i = 0; i < n; i++) {
            xx = mp_float(x[i]);
            value = constant * singleWishart_pfaffian(xx, n_min, n_max);
            result[i] = value.convert_to<double>();
            Rcpp::checkUserInterrupt();
        }
    } else {
        double constant = singleWishart_constDP(n_min, n_max);
        for(int i = 0; i < n; i++) {
            result[i] = constant * singleWishart_pfaffian(x[i], n_min, n_max);
            Rcpp::checkUserInterrupt();
        }
    }

    return result;
}

template <class T>
T mgamma_C(T x, int m, bool logar) {
    T res;
    if (logar) {
        res = 0.25*m*(m-1) * log(pi<T>());
        for(int i = 0; i < m; i++) {
            res += boost::math::lgamma(x - 0.5*i);
        }
    } else {
        res = pow(pi<T>(), 0.25*m*(m-1));
        for(int i = 0; i < m; i++) {
            res *= boost::math::tgamma(x - 0.5*i);
        }
    }
    return res;
}
