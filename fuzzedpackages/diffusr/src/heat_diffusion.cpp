/** diffusr: network diffusion algorithms in R
 *
 * Copyright (C) 2016 Simon Dirmeier
 * @author Simon Dirmeier
 * @email simon.dirmeier@bsse.ethz.ch
 *
 * This file is part of diffusr.
 *
 * diffusr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * diffusr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with diffusr. If not, see <http://www.gnu.org/licenses/>.
 */

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::plugins(cpp11)]]
#include <cmath>
#include <Eigen/Eigenvalues>

//' Do graph diffusion using an heat diffusion on a Laplacian.
//'
//' @noRd
//' @param v0  matrix of starting heat
//' @param W  the normalized Laplacian of the matrix
//' @param t  time for which heat is measured
//' @return  returns the distribution of heat at time t
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
Eigen::MatrixXd heat_diffusion_(const Eigen::MatrixXd& v0,
                                const Eigen::MatrixXd& W,
                                const double           t)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(W);
    Eigen::MatrixXd                                V  = es.eigenvectors();
    Eigen::VectorXd                                D  = es.eigenvalues();
    Eigen::MatrixXd                                co = V.transpose() * v0;

    // solution to the heat equation at time t
    for (int i = 0; i < co.rows(); ++i)
    {
        for (int j = 0; j < co.cols(); ++j)
        {
            if (j % 25 == 0)
                Rcpp::checkUserInterrupt();
            co(i, j) *= std::exp(-D(i) * t);
        }
    }
    co = V * co;

    return co;
}
