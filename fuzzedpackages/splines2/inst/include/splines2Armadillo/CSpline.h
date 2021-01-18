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

#ifndef SPLINES2_CSPLINE_H
#define SPLINES2_CSPLINE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "SplineBase.h"
#include "MSpline.h"
#include "ISpline.h"

namespace splines2 {

    // define a class for M-splines
    class CSpline : public SplineBase
    {
        // inherits constructors
        using SplineBase::SplineBase;

    private:
        // hide pure virtual function for integral here
        inline rmat integral(const bool complete_basis = true)
        {
            if (complete_basis) {
                // do nothing
            }
            return rmat();
        }

    protected:
        arma::rowvec scales_;

        // compute scales
        inline void compute_scales()
        {
            ISpline isp_obj { this };
            scales_ = mat2rowvec(
                isp_obj.set_x(boundary_knots_(1))->integral(true)
                );
        }
        inline rmat apply_scales(const rmat& x)
        {
            return x.each_row() / scales_;
        }

    public:
        // function members
        inline arma::rowvec get_scales()
        {
            return scales_;
        }

        //! Compute I-spline basis
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
                return mat_wo_col1(spline_basis_);
            }
            // else do generation
            ISpline isp_obj { this };
            spline_basis_ = isp_obj.integral(true);
            is_basis_latest_ = true;
            // compute the scale on the right boundary knot
            scales_ = mat2rowvec(
                isp_obj.set_x(boundary_knots_(1))->integral(true)
                );
            // rescale each column
            spline_basis_.each_row() /= scales_;
            if (complete_basis) {
                return spline_basis_;
            }
            return mat_wo_col1(spline_basis_);
        }

        inline virtual rmat derivative(const unsigned int derivs = 1,
                                       const bool complete_basis = true)
        {
            if (derivs == 0) {
                throw std::range_error(
                    "'derivs' has to be a positive integer."
                    );
            }
            // compute scales
            compute_scales();
            if (derivs == 1) {
                // I-spline
                ISpline isp_obj { this };
                return apply_scales(
                    isp_obj.basis(complete_basis)
                    );
            }
            // else derivs >= 2
            MSpline msp_obj { this };
            if (derivs == 2) {
                return apply_scales(
                    msp_obj.basis(complete_basis)
                    );
            }
            // else derivs >= 3
            return apply_scales(
                msp_obj.derivative(derivs - 2, complete_basis)
                );
        }

    };



}  // spline2


#endif /* SPLINES2_CSPLINE_H */
