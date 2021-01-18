/*
 * Part of the landsepi R package.
 * Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
 *                    Julien Papaix <julien.papaix@inrae.fr>
 *                    Jean-François Rey <jean-francois.rey@inrae.fr>
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
 * along with this program; if not, write to the Free Software Foundation, Inc.,i
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <Rcpp.h>
#include <gsl/gsl_rng.h> // gsl_rng_uniform, gsl_rng*...
#include <math.h> // pow...
#include <array>
#include <vector>

// Multi-dimensions vector
template <typename T>
using Vector2D = std::vector<std::vector<T> >;
template <typename T>
using Vector3D = std::vector<Vector2D<T> >;
template <typename T>
using Vector4D = std::vector<Vector3D<T> >;
template <typename T>
using Vector5D = std::vector<Vector4D<T> >;
template <typename T>
using Vector6D = std::vector<Vector5D<T> >;
template <typename T>
using Vector7D = std::vector<Vector6D<T> >;
template <typename T>
using Vector8D = std::vector<Vector7D<T> >;

// Matrix
template <typename T, int size>
using Matrix = std::array<std::array<T, size>, size>;

/****************************************************************/
/*                         functions.c                          */
/****************************************************************/
/* Useful functions to manipulate C objects */
/* Convert an integer to its binary equivalent using a recursive procedure */
int as_binary(const int& num);

/* Perform the sum of a table in 3 dimensions on its 1st dimension */
template <typename T>
const Vector2D<T> sum1_3(const int& z, const int& l, const int& c, const Vector3D<T>& t) {
    Vector2D<T> tsum1(l, std::vector<T>(c, 0));
    for(int i = 0; i < l; i++) {
        for(int j = 0; j < c; j++) {
            tsum1[i][j] = 0;
            for(int k = 0; k < z; k++)
                tsum1[i][j] += t[k][i][j];
        }
    }
    return tsum1;
}

/* Perform the sum of a table in 2 dimensions on its 2 dimensions */
template <typename T>
const T sum2_2(const Vector2D<T>& t) {
    T tsum = 0.0;
    for(unsigned int i = 0; i < t.size(); i++) {
        for(unsigned int j = 0; j < t[i].size(); j++)
            tsum += t[i][j];
    }
    return tsum;
}

/*-----------------------------------*/
/*    Gamma & beta distributions     */
/*-----------------------------------*/
/* Find the values of alpha1 and alpha2 from expectation and variance of a gamma distribution */
std::array<double, 2> find_paramGamma(const double& exp, const double& var);
/* Find the values of alpha1 and alpha2 from expectation and variance of a beta distribution */
std::array<double, 2> find_paramBeta(const double& exp, const double& var);

/*-------------------------*/
/*    Sigmoid function     */
/*-------------------------*/
/* Sigmoid function used for spore contamination and host dispersal */
double sigmoid(const double& plateau, const double& kappa, const double& sigma, const double& x);

/*--------------------------*/
/*    Trade-off function    */
/*--------------------------*/
/* Trade-off function on aggressiveness traits (cf Debarre et al. JEB 2010) */
std::vector<double> tradeoff(const std::vector<double>& x, const double& strength);

/*--------------------------*/
/*     Sample functions     */
/*--------------------------*/
/* Single random draw using a vector of probabilities   */
int sample_multinomial_once(const gsl_rng* gen, const std::vector<double>& cumProb);
/* Samples an integer array without replacement until entire array has been sampled once */
std::vector<int> sample(const gsl_rng* gen, const std::vector<int>& inArray);

#endif
