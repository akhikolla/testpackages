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

#ifndef SPLINE2_SPLINE_BASE_H
#define SPLINE2_SPLINE_BASE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"

namespace splines2 {

    // define base class for some regression splines
    class SplineBase {
    protected:
        // set and get
        rvec x_;
        rvec internal_knots_;
        rvec boundary_knots_;
        unsigned int degree_ = 3;
        unsigned int order_ = 4;
        // degree of freedom of complete spline basis
        // notice that this argument is different with the df in splines::bs()
        unsigned int spline_df_ = 4;

        // knot sequence
        rvec knot_sequence_;
        bool is_knot_sequence_latest_ = false;

        // index of x relative to internal knots
        uvec x_index_;
        bool is_x_index_latest_ = false;

        // complete spline matrix
        rmat spline_basis_;
        bool is_basis_latest_ = false;


        // pre-process some inputs
        // check knots, and do assignment
        inline void clean_knots(const rvec& internal_knots = rvec(),
                                const rvec& boundary_knots = rvec())
        {
            // check before assigment
            if (boundary_knots.has_nan()) {
                throw std::range_error("Boundary knots cannot contain NA.");
            }
            // for unspecified boundary knots
            // 1. do generation if no boundary knots have been set
            // 2. or skip checks if a set of boundary knots have been set
            if (boundary_knots.n_elem == 0) {
                if (boundary_knots_.n_elem != 2 && x_.n_elem > 0) {
                    // set boundary knots to be min(x) and max(x)
                    double left { arma::min(x_) };
                    double right { arma::max(x_) };
                    // check if boundary knots are different
                    if (left == right) {
                        throw std::range_error(
                            "Cannot set boundary knots from x."
                            );
                    }
                    boundary_knots_ = arma::zeros(2);
                    boundary_knots_(0) = left;
                    boundary_knots_(1) = right;
                }
            } else {
                // specified boundary knots
                rvec uni_boundary_knots { arma::unique(boundary_knots) };
                if (uni_boundary_knots.n_elem != 2) {
                    throw std::length_error(
                        "Need two distinct boundary knots.");
                }
                boundary_knots_ = uni_boundary_knots;
            }
            if (internal_knots.has_nan()) {
                throw std::range_error("Internal knots cannot contain NA.");
            }
            // for non-empty internal knots
            // 1. get unique values
            // 2. checks if outside of boundary
            if (internal_knots.n_elem > 0) {
                rvec uni_internal_knots { arma::unique(internal_knots) };
                // check internal knots are inside of boundary knots
                double min_int_knots { arma::min(internal_knots) };
                double max_int_knots { arma::max(internal_knots) };
                if (boundary_knots_[0] >= min_int_knots ||
                    boundary_knots_[1] <= max_int_knots) {
                    throw std::range_error(
                        "Internal knots must be set inside of boundary knots."
                        );
                }
                internal_knots_ = uni_internal_knots;
            } else {
                internal_knots_ = internal_knots;
            }
        }

        // compute spline df
        inline void update_spline_df()
        {
            spline_df_ = internal_knots_.n_elem + order_;
        }

        inline rvec get_knot_sequence(const unsigned int order = 1)
        {
            rvec out { arma::zeros(internal_knots_.n_elem + 2 * order) };
            for (size_t i {0}; i < out.n_elem; ++i) {
                if (i < order) {
                    out(i) = boundary_knots_(0);
                } else if (i < out.n_elem - order) {
                    out(i) = internal_knots_(i - order);
                } else {
                    out(i) = boundary_knots_(1);
                }
            }
            return out;
        }
        inline void update_knot_sequence()
        {
            if (! is_knot_sequence_latest_ || knot_sequence_.n_elem == 0) {
                knot_sequence_ = get_knot_sequence(order_);
                is_knot_sequence_latest_ = true;
            }
        }
        // allow x outside of boundary
        // let degree 0 basis take 1 outside boundary
        inline void update_x_index()
        {
            if (! is_x_index_latest_ || x_index_.n_elem == 0) {
                x_index_ = arma::zeros<arma::uvec>(x_.n_elem);
                for (size_t i {0}; i < x_.n_elem; ++i) {
                    size_t left_index {0};
                    size_t right_index { internal_knots_.n_elem };
                    while (right_index > left_index) {
                        size_t cur_index { (left_index + right_index) / 2 };
                        if (x_(i) < internal_knots_(cur_index)) {
                            right_index = cur_index;
                        } else {
                            left_index = cur_index + 1;
                        }
                    }
                    x_index_(i) = left_index;
                }
                is_x_index_latest_ = true;
            }
        }

    public:
        // the default constructor
        SplineBase() {}
        virtual ~SplineBase() {}

        // explicit constructor
        explicit SplineBase(const SplineBase* pSplineBase) :
            x_ { pSplineBase->x_ },
            internal_knots_ { pSplineBase->internal_knots_ },
            boundary_knots_ { pSplineBase->boundary_knots_ },
            degree_ { pSplineBase->degree_  },
            order_ { pSplineBase->order_ },
            is_knot_sequence_latest_ { pSplineBase->is_knot_sequence_latest_ },
            is_x_index_latest_ { pSplineBase->is_x_index_latest_ },
            is_basis_latest_ { false }
        {
            update_spline_df();
            update_knot_sequence();
            update_x_index();
        }

        // constructor with specificied internal_knots
        SplineBase(const rvec& x,
                   const rvec& internal_knots,
                   const unsigned int degree = 3,
                   const rvec& boundary_knots = rvec()) :
            x_ { x },
            degree_ { degree }
        {
            clean_knots(internal_knots, boundary_knots);
            order_ = degree_ + 1;
            update_spline_df();
        }

        // constructor with specified df
        SplineBase(const rvec& x,
                   const unsigned int spline_df,
                   const unsigned int degree = 3,
                   const rvec& boundary_knots = rvec()) :
            x_ { x },
            degree_ { degree }
        {
            order_ = degree_ + 1;
            if (spline_df < order_) {
                throw std::range_error("The specified df was too small.");
            }
            spline_df_ = spline_df;
            // determine internal knots by spline_df and x
            unsigned int n_internal_knots { spline_df_ - order_ };
            if (n_internal_knots == 0) {
                clean_knots(rvec(), boundary_knots);
            } else {
                rvec prob_vec { arma::linspace(0, 1, n_internal_knots + 2) };
                prob_vec = prob_vec.subvec(1, n_internal_knots);
                clean_knots(rvec(), boundary_knots);
                // get quantiles of x within boundary only
                rvec x_inside { get_inside_x(x, boundary_knots_) };
                rvec internal_knots { arma_quantile(x_inside, prob_vec) };
                clean_knots(internal_knots);
            }
        }

        // function members
        // "setter" functions
        inline SplineBase* set_x(const rvec& x)
        {
            x_ = x;
            is_x_index_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline SplineBase* set_x(const double x)
        {
            x_ = num2vec(x);
            is_x_index_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline SplineBase* set_internal_knots(
            const rvec& internal_knots
            )
        {
            clean_knots(internal_knots);
            // update spline df
            update_spline_df();
            is_knot_sequence_latest_ = false;
            is_x_index_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline SplineBase* set_boundary_knots(
            const rvec& boundary_knots
            )
        {
            clean_knots(internal_knots_, boundary_knots);
            is_knot_sequence_latest_ = false;
            is_x_index_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline SplineBase* set_degree(
            const unsigned int degree
            )
        {
            degree_ = degree;
            order_ = degree + 1;
            // update spline df
            update_spline_df();
            is_knot_sequence_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline SplineBase* set_order(const unsigned int order)
        {
            if (order > 0) {
                set_degree(order - 1);
            } else {
                throw std::range_error("The 'order' must be at least 1.");
            }
            return this;
        }

        // "getter" functions
        inline rvec get_x() const
        {
            return x_;
        }
        inline rvec get_internal_knots() const
        {
            return internal_knots_;
        }
        inline rvec get_boundary_knots() const
        {
            return boundary_knots_;
        }
        inline unsigned int get_degree() const
        {
            return degree_;
        }
        inline unsigned int get_order() const
        {
            return order_;
        }
        inline unsigned int get_spline_df() const
        {
            return spline_df_;
        }

        // define pure virtual functions
        inline virtual rmat basis(
            const bool complete_basis = true
            ) = 0;
        inline virtual rmat derivative(
            const unsigned int derivs = 1,
            const bool complete_basis = true
            ) = 0;
        inline virtual rmat integral(
            const bool complete_basis = true
            ) = 0;

    };

}

#endif
