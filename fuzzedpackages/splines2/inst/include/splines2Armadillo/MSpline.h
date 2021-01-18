//
// R package splines2 by Wenjie Wang and Jun Yan
// Copyright (C) 2016-2020
//
// This file is part of the R package splines2.
//
// The R package splines2 is free software: You can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any later
// version (at your option). See the GNU General Public License at
// <https://www.gnu.org/licenses/> for details.
//
// The R package splines2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//

#ifndef SPLINES2_MSPLINE_H
#define SPLINES2_MSPLINE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "SplineBase.h"

namespace splines2 {

    // define a class for M-splines
    class MSpline : public SplineBase
    {
        // inherits constructors
        using SplineBase::SplineBase;

    public:
        // function members

        //! Compute M-spline basis
        //!
        //! @param complete_basis A `bool` value indicating whether to return a
        //! complete spline basis
        //!
        //! @return arma::mat
        inline virtual rmat basis(const bool complete_basis = true)
        {
            // early exit if latest
            if (is_basis_latest_) {
                if (complete_basis) {
                    return spline_basis_;
                }
                // else
                return mat_wo_col1(spline_basis_);
            }
            // else do the generation
            update_x_index();
            update_knot_sequence();
            // define output matrix
            rmat b_mat {
                arma::zeros(x_.n_elem, spline_df_)
            };
            // generate bases of degree 0
            for (size_t i {0}; i < x_.n_elem; ++i) {
                unsigned int j { x_index_(i) };
                double denom { knot_sequence_(j + order_) -
                               knot_sequence_(j + degree_) };
                b_mat(i, j) = 1 / denom;
            }
            // main loop
            for (unsigned int k {1}; k <= degree_; ++k) {
                double dk { static_cast<double>(k) };
                double dk1 { (1.0 + 1.0 / dk) };
                const unsigned int k_offset { degree_ - k };
                // use the Cox-de Boor recursive formula
                for (size_t i {0}; i < x_.n_elem; ++i) {
                    double saved { 0 };
                    // for each x, at most "order" columns are not zero
                    // basis(x) is not zero from t_ii to t_{ii+k+1}
                    // where ii is index of x in terms of bases
                    for (size_t j {0}; j < k; ++j) {
                        size_t j_index { x_index_(i) + j };
                        size_t i1 { j_index + k_offset };
                        size_t i2 { j_index + order_ };
                        double den { knot_sequence_(i2) - knot_sequence_(i1) };
                        double term { dk1 * b_mat(i, j_index) };
                        b_mat(i, j_index) = saved +
                            (knot_sequence_(i2) - x_(i)) * term / den;
                        double den2 {
                            knot_sequence_(i2 + 1) - knot_sequence_(i1 + 1)
                        };
                        saved = (x_(i) - knot_sequence_(i1 + 1)) * term / den2;
                    }
                    b_mat(i, x_index_(i) + k) = saved;
                }
            }
            // about to return
            spline_basis_ = b_mat;
            is_basis_latest_ = true;
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

        // derivatives of M-splines
        inline virtual rmat derivative(
            const unsigned int derivs = 1,
            const bool complete_basis = true
            )
        {
            if (derivs == 0) {
                throw std::range_error(
                    "'derivs' has to be a positive integer.");
            }
            // early exit if derivs is large enough
            unsigned int old_df { spline_df_ };
            if (degree_ < derivs) {
                if (complete_basis) {
                    return arma::zeros(x_.n_elem, old_df);
                }
                if (old_df == 1) {
                    throw std::range_error("No column left in the matrix.");
                }
                return arma::zeros(x_.n_elem, old_df - 1);
            }
            // back up current results if necessary
            bool backup_basis { is_basis_latest_ };
            bool backup_knot_sequence { is_knot_sequence_latest_ };
            rmat old_basis;
            rvec old_knot_sequence;
            if (backup_basis) {
                old_basis = spline_basis_;
            }
            if (backup_knot_sequence) {
                old_knot_sequence = knot_sequence_;
            }
            // get basis matrix for (degree - derivs)
            set_degree(degree_ - derivs);
            rmat d_mat { basis(true) };
            // restore
            is_basis_latest_ = backup_basis;
            is_knot_sequence_latest_ = backup_knot_sequence;
            set_degree(degree_ + derivs);
            if (backup_basis) {
                spline_basis_ = old_basis;
            }
            if (backup_knot_sequence) {
                knot_sequence_ = old_knot_sequence;
            } else {
                update_knot_sequence();
            }
            // add zero columns
            d_mat = add_zero_cols(d_mat, old_df - d_mat.n_cols);
            // main loop
            for (unsigned int k {1}; k <= derivs; ++k) {
                const unsigned int k_offset { derivs - k };
                const size_t numer { degree_ - k_offset };
                for (size_t i {0}; i < x_.n_elem; ++i) {
                    double saved { 0 };
                    for (size_t j {0}; j < numer; ++j) {
                        size_t j_index { x_index_(i) + j };
                        size_t i1 { j_index + k_offset };
                        size_t i2 { j_index + order_ };
                        double den { knot_sequence_(i2) - knot_sequence_(i1) };
                        double term { (numer + 1) * d_mat(i, j_index) };
                        d_mat(i, j_index) = saved - term / den;
                        double den2 {
                            knot_sequence_(i2 + 1) - knot_sequence_(i1 + 1)
                        };
                        saved = term / den2;
                    }
                    d_mat(i, x_index_(i) + numer) = saved;
                }
            }
            // remove the first column if needed
            if (complete_basis) {
                return d_mat;
            }
            // else
            return mat_wo_col1(d_mat);
        }

        // integral of M-splines (I-splines)
        inline virtual rmat integral(const bool complete_basis = true)
        {
            // back up current results
            bool backup_basis { is_basis_latest_ };
            bool backup_knot_sequence { is_knot_sequence_latest_ };
            rmat old_basis;
            rvec old_knot_sequence;
            if (backup_basis) {
                old_basis = spline_basis_;
            }
            if (backup_knot_sequence) {
                old_knot_sequence = knot_sequence_;
            }
            // get basis matrix for (degree - derivs)
            set_degree(degree_ + 1);
            rmat i_mat { basis(false) };
            rvec knot_sequence_ord { knot_sequence_ };
            // restore
            set_degree(degree_ - 1);
            is_basis_latest_ = backup_basis;
            is_knot_sequence_latest_ = backup_knot_sequence;
            if (backup_basis) {
                spline_basis_ = old_basis;
            }
            if (backup_knot_sequence) {
                knot_sequence_ = old_knot_sequence;
            }
            // compute t_{i+k+1} - t_{i}
            arma::rowvec numer1 { arma::zeros<arma::rowvec>(i_mat.n_cols) };
            for (size_t j { 0 }; j < numer1.n_elem; ++j) {
                numer1(j) = knot_sequence_ord(j + order_ + 2) -
                    knot_sequence_ord(j + 1);
            }
            // for each row of i_mat
            for (size_t i {0}; i < x_.n_elem; ++i) {
                size_t k1 { x_index_(i) }, k2 { k1 + degree_ };
                arma::rowvec numer2 { i_mat(i, arma::span(k1, k2)) };
                arma::rowvec numer { numer1.cols(k1, k2) % numer2 };
                numer = rev_cum_sum(numer);
                for (size_t j {0}; j < i_mat.n_cols; ++j) {
                    if (j > k2) {
                        i_mat(i, j) = 0.0;
                    } else if (j >= k1) {
                        i_mat(i, j) = numer(j - k1) / (order_ + 1);
                    } else {
                        i_mat(i, j) = 1.0;
                    }
                }
            }
            // remove the first column if needed
            if (complete_basis) {
                return i_mat;
            }
            return mat_wo_col1(i_mat);
        }

    };
}  // splines2


#endif /* SPLINES2_MSPLINE_H */
