/*
 * alpha_k_cpp : Compute Krippendorff's Alpha
 * Copyright (C) 2017 Alexander Staudt
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.r-project.org/Licenses/GPL-2
 */

#include <Rcpp.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

using namespace Rcpp;

#include "alpha.h"
#include "bootstrap_alpha.h"
#include "RngStream.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List alpha_k_cpp(
        NumericMatrix data,
        int metric,
        bool bootstrap,
        bool bootnp,
        int nboot,
        int nnp,
        NumericVector cmrg_seed,
        int n_threads)
{
    // ========================================================================
    // inputs
    // ========================================================================
    // number of coders, number of units
    int nC = data.rows();
    int nU = data.cols();

    // convert R column-major matrix into row-major std::vector
    std::vector<double> reliability_data = as<std::vector<double>>(transpose(data));

    // seed for l'Ecuyer-CMRG random number generator
    std::vector<double> signed_seed = as<std::vector<double>>(cmrg_seed);
    std::vector<unsigned long> seed (6, 0);
    for (int i = 0; i < 6; i++) { // convert to unsigned integer
        seed[i] = (unsigned long) std::floor(signed_seed[i]);
    }

    // check number of threads
    #ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    #else
    int max_threads = 1;
    #endif
    if (n_threads > max_threads) {
        n_threads = max_threads;
        Rprintf(" Note: Number of specified cores higher than number of available cores.\n");
        Rprintf(" Resetting to available number of %d cores.\n", max_threads);
    } else if (n_threads < 1) {
        n_threads = 1;
        Rprintf(" Note: Number of specified cores lower than 1.\n");
        Rprintf(" Resetting to %d core.\n", n_threads);
    }

    // ========================================================================
    // compute alpha
    // ========================================================================
    resultsAlpha res_alpha;
    int rv = get_alpha(reliability_data,
                       nC,
                       nU,
                       metric,
                       res_alpha);

    // =========================================================================
    // bootstrap alpha using Krippendorff's algorithm
    // =========================================================================
    std::vector<double> b_alphas(1, NAN);
    int rv_boot_ka = 0;
    if (bootstrap == true) {
        b_alphas.resize(nboot, NAN);
        rv_boot_ka = bootstrap_alpha(res_alpha.D_e,
                                   reliability_data,
                                   res_alpha.n_coders,
                                   res_alpha.n_units,
                                   res_alpha.coding_values,
                                   res_alpha.contributions,
                                   res_alpha.metric,
                                   nboot,
                                   seed.data(),
                                   n_threads,
                                   b_alphas);
    }

    // =========================================================================
    // nonparametric bootstrap
    // =========================================================================
    std::vector<double> b_alphasNP(1, NAN);
    int rv_boot_np = 0;
    if (bootnp == true) {
        b_alphasNP.resize(nnp, NAN);
        rv_boot_np = bootstrap_alpha_nonparametric(reliability_data,
                                                   nC,
                                                   nU,
                                                   metric,
                                                   nnp,
                                                   seed.data(),
                                                   n_threads,
                                                   b_alphasNP);
    }

    // =========================================================================
    // return results
    // =========================================================================
    int nV = res_alpha.coding_values.size();
    double D_e = res_alpha.D_e;
    double D_o = res_alpha.D_o;
    double alpha = res_alpha.alpha;

    // corresponding Rcpp-objects
    NumericVector codingValues(nV);
    NumericMatrix cMat(nV, nV);
    NumericMatrix dMat(nV, nV);

    // copy values into Rcpp-objects and create result list
    if (rv == 0) {

        for (int i = 0; i < nV; i++) {
            for (int j = 0; j < nV; j++) {
                cMat[i + nV * j] = res_alpha.coincidence_matrix[i * nV + j];
                dMat[i + nV * j] = res_alpha.delta_matrix[i * nV + j];
            }
            codingValues[i] = res_alpha.coding_values[i];
        }

        List results = List::create(
            Named("alpha") = alpha,
            Named("metric") = metric,
            Named("n_coders") = nC,
            Named("n_units") = nU,
            Named("n_values") = nV,
            Named("coding_values") = codingValues,
            Named("coincidence_matrix") = cMat,
            Named("delta_matrix") = dMat,
            Named("D_e") = D_e,
            Named("D_o") = D_o,
            Named("bootstrap") = bootstrap,
            Named("bootnp") = bootnp,
            Named("nboot") = nboot,
            Named("nnp") = nnp,
            Named("bootstraps") = wrap(b_alphas),
            Named("bootstrapsNP") = wrap(b_alphasNP));

        // additional messages regarding bootstrap
        if (rv_boot_ka == -1 && rv_boot_np != -1) {

            Rprintf("Warning: Memory allocation failed in Krippendorff-bootstrap routine.\n");
            Rprintf("Other results remain valid nonetheless.\n");

        } else if (rv_boot_ka != -1 && rv_boot_np == 1) {

            Rprintf("Warning: Memory allocation failed in non-parametric bootstrap routine.\n");
            Rprintf("Other results remain valid nonetheless.\n");

        } else if (rv_boot_ka == -1 && rv_boot_np == -1) {

            Rprintf("Warning: Memory allocation failed in bootstrap routines.\n");
            Rprintf("Other results remain valid nonetheless.\n");

        }

        return results;

    } else if (rv == -2) {

        stop("The provided metric does not exist");
        return R_NilValue;

    } else {

        stop("Unknown error occurred");
        return R_NilValue;

    }
}
