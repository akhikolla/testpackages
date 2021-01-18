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

/* language: C++11 */
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include "alpha.h"

// =============================================================================
// indicate missing values in data (vector, matrix, etc.)
// =============================================================================
std::vector<int> check_nonmissing(
    const std::vector<double> &data,
    const int nC,
    const int nU)
{
    std::vector<int> is_nonmissing (nC * nU, 0);
    for (int i = 0; i < nC; i++) {
        for (int u = 0; u < nU; u++) {
            if (!std::isnan(data[i * nU + u])) {
                is_nonmissing[i * nU + u] = 1;
            } else {
                is_nonmissing[i * nU + u] = 0;
            }
        }
    }
    return is_nonmissing;
}

// =============================================================================
// obtain colsums from matrix
// =============================================================================
// vectors/matrices of type double
std::vector<double> get_colsum_double(
    std::vector<double> &data,
    const int nC,
    const int nU)
{
    std::vector<double> colsum (nU, 0);
    for (int u = 0; u < nU; u++) {
        double sum = 0;
        for (int i = 0; i < nC; i++) {
            sum = sum + data[i * nU + u];
        }
        colsum[u] = sum;
    }
    return colsum;
}

// vectors/matrices of type int
std::vector<int> get_colsum_int(
    std::vector<int> &data,
    const int nC,
    const int nU)
{
    std::vector<int> colsum (nU, 0);
    for (int u = 0; u < nU; u++) {
        int sum = 0;
        for (int i = 0; i < nC; i++) {
            sum = sum + data[i * nU + u];
        }
        colsum[u] = sum;
    }
    return colsum;
}


// =============================================================================
// obtain nonmissing values from vector
// =============================================================================
std::vector<double> get_nonmissing_values(
    const std::vector<double> &data)
{
    std::vector<double> nonmissing_values (data.size(), 0);
    int j = 0;
    for (unsigned int i = 0; i < data.size(); i++) {
        if (!std::isnan(data[i])) {
            nonmissing_values[j] = data[i];
            j++;
        }
    }
    nonmissing_values.resize(j);
    return nonmissing_values;
}

// =============================================================================
// obtain unique values from vector
// =============================================================================
std::vector<double> get_unique_values(
    const std::vector<double> &data)
{
    std::vector<double> unique_values = get_nonmissing_values(data);
    std::sort(unique_values.begin(), unique_values.end());
    std::vector<double>::iterator iter;
    iter = std::unique(unique_values.begin(), unique_values.end());
    unique_values.resize(std::distance(unique_values.begin(), iter));

    return unique_values;
}

// =============================================================================
// obtain coder-value pairs
// =============================================================================
std::vector<double> get_coder_value_pairs(
    const std::vector<double> &data,
    const int nC,
    const int nU)
{
    int pairs = nC * (nC - 1);
    std::vector<double> coder_value_pairs (pairs * 2 * nU, 0);
    for (int u = 0; u < nU; u++) {
        int k = 0;
        for (int i = 0; i < nC; i++) {
            for (int j = 0; j < nC; j++) {
                if (i != j) {
                    coder_value_pairs[u * pairs * 2 + k * 2 + 0] = data[i * nU + u];
                    coder_value_pairs[u * pairs * 2 + k * 2 + 1] = data[j * nU + u];
                    k++;
                }
            }
        }
    }
    return coder_value_pairs;
}

// =============================================================================
//  obtain coincidence matrix
// =============================================================================
std::vector<double> get_coincidence_matrix(
    const std::vector<double> coder_value_pairs,
    const std::vector<double> unique_values,
    const std::vector<int> m_u,
    const int nC,
    const int nU)
{
    int pairs = nC * (nC - 1);
    int nV = unique_values.size();
    int occurrence = 0;
    double value_1 = 0;
    double value_2 = 0;
    double contribution = 0;
    std::vector<double> coincidence_matrix(nV * nV, 0);

    for (int i = 0; i < nV; i++) {
        for (int j = 0; j < nV; j++) {
            contribution = 0;
            for (int u = 0; u < nU; u++) {
                if (m_u[u] - 1 > 0) {
                    occurrence = 0;
                    for (int p = 0; p < pairs; p++) {
                        value_1 = coder_value_pairs[u * pairs * 2 + p * 2 + 0];
                        value_2 = coder_value_pairs[u * pairs * 2 + p * 2 + 1];
                        if (value_1 == unique_values[i] && value_2 == unique_values[j]) {
                            occurrence++;
                        }
                    }
                    contribution = contribution + (double) occurrence / (double) (m_u[u] - 1);
                }
            }
            coincidence_matrix[i * nV + j] = contribution;
        }
    }
    return coincidence_matrix;
}


// =============================================================================
// compute sum over consecutive array elements
// =============================================================================
double get_array_sum(
    const std::vector<double> &vector,
    const int first,
    const int last)
{
    double sum = 0;
    int i = 0;
    for (i = first; i <= last; i++) {
        sum = sum + vector[i];
    }
    return sum;
}

// =============================================================================
//  get index from ordered values
// =============================================================================
int get_index_from_ordered(
    const double value,
    const std::vector<double> &unique_values)
{
    int index = 0;
    int i = 0;
    int nV = unique_values.size();
    while((i < nV) && (value > unique_values[i])) {
        i++;
    }
    if (unique_values[i] == value) {
        index = i;
    }
    return index;
}

// =============================================================================
// define difference-functions
// =============================================================================
double get_deltasquared(
    const double c,
    const double k,
    const std::vector<double> &coding_values,
    const std::vector<double> &marginalSums,
    const int metric)
{
    int index_c = get_index_from_ordered(c, coding_values);
    int index_k = get_index_from_ordered(k, coding_values);

    double deltaSquared = 0;
    int U = 6;

    if (metric == 1) { // nominal
        if (c != k) {
            deltaSquared = 1;
        } else {
            deltaSquared = 0;
        }
    } else if (metric == 2) { // ordinal
        double n_c = marginalSums[index_c];
        double n_k = marginalSums[index_k];
        if (index_c <= index_k) {
            deltaSquared = std::pow(get_array_sum(marginalSums, index_c, index_k) - (n_c + n_k)/2, 2);
        } else {
            deltaSquared = std::pow(get_array_sum(marginalSums, index_k, index_c) - (n_c + n_k)/2, 2);
        }
    } else if (metric == 3) { // interval
        deltaSquared = std::pow(c - k, 2);
    } else if (metric == 4) { // ratio
        deltaSquared = std::pow((c - k) / (c + k), 2);
    } else if (metric == 5) { // circular
        deltaSquared = std::pow(sin(M_PI * (c - k) / U), 2);
    } else if (metric == 6) { // bipolar
        std::vector<double>::const_iterator c_min;
        std::vector<double>::const_iterator c_max;
        c_min = std::min_element(coding_values.begin(), coding_values.end());
        c_max = std::max_element(coding_values.begin(), coding_values.end());
        if (c != k) {
            deltaSquared = std::pow(c - k, 2) / ((c + k - 2 * (*c_min)) * (2 * (*c_max) - c - k));
        } else {
            deltaSquared = 0;
        }
    } else { // otherwise
        deltaSquared = NAN;
    }
    return deltaSquared;
}

// =============================================================================
// obtain delta2-matrix
// =============================================================================
std::vector<double> get_delta_matrix(
    const std::vector<double> &unique_values,
    const std::vector<double> &contributions,
    const int metric)
{
    int nV = unique_values.size();
    std::vector<double> delta_matrix (nV * nV, 0);
    double c = 0;
    double k = 0;
    for (int i = 0; i < nV; i++) {
        for (int j = 0; j < nV; j++) {
            c = unique_values[i];
            k = unique_values[j];
            delta_matrix[i * nV + j] = get_deltasquared(c, k, unique_values, contributions, metric);
        }
    }
    return delta_matrix;
}

// =============================================================================
// obtain expected disagreement
// =============================================================================
double get_expected_disagreement(
    const std::vector<double> contributions,
    const std::vector<double> delta_matrix,
    const int nV,
    const double pairable_units)
{
    double sum = 0;
    for (int i = 0; i < nV; i++) {
        for (int j = 0; j < nV; j++) {
            sum = sum + contributions[i] * contributions[j] * delta_matrix[i * nV + j];
        }
    }
    double exp_disagreement = sum * 1 / (pairable_units * (pairable_units - 1));
    return exp_disagreement;
}

// =============================================================================
// obtain observed disagreement
// =============================================================================
double get_observed_disagreement(
    const std::vector<double> coincidence_matrix,
    const std::vector<double> delta_matrix,
    const int nV,
    const double pairable_units)
{
    double sum = 0;
    for (int i = 0; i < nV; i++) {
        for (int j = 0; j < nV; j++) {
            sum = sum + coincidence_matrix[i * nV + j] * delta_matrix[i * nV + j];
        }
    }
    double obs_disagreement = sum / pairable_units;
    return obs_disagreement;
}

// =============================================================================
// compute alpha: get_alpha()
// =============================================================================
int get_alpha(
    const std::vector<double> &data,
    const int nC,
    const int nU,
    const int metric,
    resultsAlpha &results)
{
    //int U = 6;
    
    // get unique coding values and number of unique coding values
    std::vector<double> unique_values = get_unique_values(data);
    int nV = unique_values.size();

    // compute number of coded values in each unit
    std::vector<int> nonmissing_codings = check_nonmissing(data, nC, nU);

    // number of coders assigning non-missing values per unit
    std::vector<int> m_u = get_colsum_int(nonmissing_codings, nC, nU);

    // get coder-value pairs
    std::vector<double> coder_value_pairs = get_coder_value_pairs(data, nC, nU);

    // coincidence matrix
    std::vector<double> coincidence_matrix = get_coincidence_matrix(coder_value_pairs, unique_values, m_u, nC, nU);

    // contributions of each value within pairable observations
    std::vector<double> contributions = get_colsum_double(coincidence_matrix, nV, nV);

    // total number of pairable observations
    double pairable_units = std::accumulate(contributions.begin(), contributions.end(), 0);

    // compute delta^2-matrix
    std::vector<double> delta_matrix = get_delta_matrix(unique_values, contributions, metric);

    // compute expected disagreement (D_e)
    double D_e = get_expected_disagreement(contributions, delta_matrix, nV, pairable_units);

    // compute observed disagreement (D_o)
    double D_o = get_observed_disagreement(coincidence_matrix, delta_matrix, nV, pairable_units);

    // =========================================================================
    // compute alpha
    // =========================================================================
    if (D_e == 0) {
        results.alpha = 1;
    } else {
        results.alpha = 1 - D_o/D_e;
    }

    results.D_e = D_e;
    results.D_o = D_o;
    results.coincidence_matrix = coincidence_matrix;
    results.delta_matrix = delta_matrix;
    results.coding_values = unique_values;
    results.contributions = contributions;
    results.n_coders = nC;
    results.n_units = nU;
    results.metric = metric;

    return 0;
}
