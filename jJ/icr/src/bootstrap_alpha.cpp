/*
 * Copyright (C) 2017-2019  Alexander Staudt
 * 
 * This file is part of icr.
 *
 * icr is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * icr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with icr; if not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

#include <cstdio>
#include <csignal>

#include "RngStream.h"
#include "alpha.h"
#include "bootstrap_alpha.h"

#ifdef _OPENMP
#include <omp.h>
#endif


// =============================================================================
// signal handling
// =============================================================================
namespace
{
    volatile std::sig_atomic_t user_interrupt = 0;
}

// signal handler function: set user-interrupt flag to 1
void sig_handler(int signo)
{
    if (signo == SIGINT) {
        user_interrupt = 1;
    }
}

// =============================================================================
// get number of non-missing codings for each unit
// =============================================================================
std::vector<int> get_m_u(
    const std::vector<double> &data,
    const int nC,
    const int nU)
{
    int i = 0;
    int u = 0;

    // get nonmissing codings
    std::vector<int> nonmissingCodings (nC * nU, 0);
    for (i = 0; i < nC; i++) {
        for (u = 0; u < nU; u++) {
            if (!std::isnan(data[i * nU + u])) {
                nonmissingCodings[i * nU + u] = 1;
            } else {
                nonmissingCodings[i * nU + u] = 0;
            }
        }
    }

    // get number of nonmissing codings for each unit
    std::vector<int> m_u (nU, 0);
    for (u = 0; u < nU; u++) {
        int sum = 0;
        for (i = 0; i < nC; i++) {
            sum = sum + nonmissingCodings[i * nU + u];
        }
        m_u[u] = sum;
    }
    return m_u;
}

// =============================================================================
// get number of pairable observations (n_dotdot)
// =============================================================================
int get_n_pairable_obs(
    const std::vector<int> &m_u)
{
    unsigned int nU = m_u.size();
    unsigned int i = 0;
    int n_dotdot = 0;
    for (i = 0; i < nU; i++) {
        if (m_u[i] > 1) {
            n_dotdot = n_dotdot + m_u[i];
        }
    }
    return n_dotdot;
}

// =============================================================================
// get possible coder-value pairs
// =============================================================================
std::vector<double> get_possible_pairs(
    const std::vector<double> &data,
    const int nC,
    const int nU,
    const std::vector<int> &m_u)
{
    // get total number of unique coder-value pairs contributing to alpha
    int nucvp = 0;
    for (int u = 0; u < nU; u++) {
        nucvp = nucvp + (m_u[u] * (m_u[u] - 1) / 2);
    }

    // get possible coder value pairs
    std::vector<double> possiblePairs(nucvp * 2, 0);
    int hasnan = 0;
    int i_p = 0;
    for (int u = 0; u < nU; u++) {
        for (int i = 0; i < nC; i++) {
            for (int j = 0; j < nC; j++) {
                if (i < j) {
                    hasnan = std::isnan(data[i * nU + u]) || std::isnan(data[j * nU + u]);
                    if (hasnan == 0) {
                        possiblePairs[i_p * 2 + 0] = data[i * nU + u];
                        possiblePairs[i_p * 2 + 1] = data[j * nU + u];
                        i_p = i_p + 1;
                    }
                }
            }
        }
    }
    return possiblePairs;
}

// =============================================================================
// bootstrap alpha (Krippendorff's algorithm)
// =============================================================================
int bootstrap_alpha(
    const double D_e,
    const std::vector<double> &data,
    const int nC,
    const int nU,
    const std::vector<double> &coding_values,
    const std::vector<double> &contributions,
    const int metric,
    const int bootstraps,
    const unsigned long seed[6],
    const int n_threads,
    std::vector<double> &alphas)
{
    // compute number of coded values per unit
    std::vector<int> m_u = get_m_u(data, nC, nU);

    // compute total number of pairable observations
    int n_dotdot = get_n_pairable_obs(m_u);

    // get number of unique pairs per unit
    std::vector<int> uniquePairs (nU);
    for (int u = 0; u < nU; u++) {
        uniquePairs[u] = m_u[u] * (m_u[u] - 1) / 2;
    }

    // total number of unique pairs in reliability data
    int N_o = std::accumulate(uniquePairs.begin(), uniquePairs.end(), 0);

    // get (valid) coder value pairs
    std::vector<double> possiblePairs = get_possible_pairs(data, nC, nU, m_u);

    // cumulative sum (zum Indizieren der coder-value-pairs)
    std::vector<int> cvp_indices (nU + 1, 0); // cumulative sum + index 0
    int sum = 0;
    for (int u = 0; u < nU; u++) {
        sum = sum + uniquePairs[u];
        cvp_indices[u + 1] = sum;
    }

    RngStream::SetPackageSeed(seed);
    std::vector<RngStream> generators(n_threads);

    // compute alpha
    #ifdef _OPENMP
    #pragma omp parallel num_threads(n_threads)
    {
        int thread_id = omp_get_thread_num();
    #pragma omp for
    #else
        int thread_id = 0;
    #endif
        for (int x = 0; x < bootstraps; x++) {

            double alpha = 1;

            // compute deviations per bootstrap
            double summed_deviations = 0;
            for (int u = 0; u < nU; u++) {

                // compute E(r) per unit
                double deviation_u = 0;
                for (int i = cvp_indices[u]; i < cvp_indices[u + 1]; i++) {

                    int r = generators[thread_id].RandInt(0, N_o - 1);
                    double cvp_1 = possiblePairs[r * 2 + 0];
                    double cvp_2 = possiblePairs[r * 2 + 1];

                    // compute delta squared for different metrics
                    double deltaSquared = get_deltasquared(cvp_1, cvp_2, coding_values, contributions, metric);
                    double E_r = 2 * deltaSquared/(n_dotdot * D_e);
                    double d = E_r/(m_u[u] - 1);
                    deviation_u = deviation_u + d;
                }
                summed_deviations = summed_deviations + deviation_u;
            }
            alpha = alpha - summed_deviations;
            if (alpha >= -1) {
                alphas[x] = alpha;
            }
            // check user-interrupt
            signal(SIGINT, sig_handler);
            if ((x % 1000 == 0) && user_interrupt == 1) {
                std::fill(alphas.begin(), alphas.end(), NAN); // invalidate incomplete bootstraps
                x = bootstraps; // jump to end of loop
            }
        }
    #ifdef _OPENMP
    }
    #endif

    // reset interrupt-flag
    user_interrupt = 0;

    // return
    return 0;
}

// =============================================================================
// bootstrap alpha (non-parametric bootstrap)
// =============================================================================
int bootstrap_alpha_nonparametric(
    const std::vector<double> &data,
    const int nC,
    const int nU,
    const int metric,
    const int bootstraps,
    const unsigned long seed[6],
    const int n_threads,
    std::vector<double> &alphas)
{
    std::vector<int> column_indices (bootstraps * nU);
    try {
        column_indices.reserve(bootstraps * nU);
    } catch (const std::bad_alloc &e) {
        return -1;
    }

    RngStream::SetPackageSeed(seed);
    RngStream generator;
    for (int i = 0; i < bootstraps * nU; i++) {
        column_indices[i] = generator.RandInt(0, nU - 1);
    }

    // compute alpha on resampled reliability data
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads)
    #endif
    for (int i = 0; i < bootstraps; i++) {
        std::vector<double> resampled_data (nC * nU, 0);
        for (int j = 0; j < nC; j++) {
            for (int k = 0; k < nU; k++) {
                resampled_data[j * nU + k] = data[j * nU + column_indices[i * nU + k]];
            }
        }
        resultsAlpha results_alpha;
        int rv = get_alpha(resampled_data, nC, nU, metric, results_alpha);
        if (rv == 0) {
            alphas[i] = results_alpha.alpha;
        } else {
            alphas[i] = NAN;
        }
        // check user-interrupt
        signal(SIGINT, sig_handler);
        if ((i % 1000 == 0) && user_interrupt == 1) {
            std::fill(alphas.begin(), alphas.end(), NAN); // invalidate incomplete bootstraps
            i = bootstraps; // jump to end of loop
        }
    }

    // reset interrupt-flag
    user_interrupt = 0;

    // return
    return 0;
}
