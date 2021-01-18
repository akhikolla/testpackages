// Copyright (C) 2012 - 2015  Philip Rinn
// Copyright (C) 2012 - 2015  Carl von Ossietzy Universität Oldenburg
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, see <https://www.gnu.org/licenses/gpl-2.0>.

// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include "linreg.h"
#ifdef _OPENMP
    #include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::export(".Langevin2D")]]
List Langevin2D(const arma::mat& data, const int& bins, const arma::vec& steps,
                const double& sf, const int& bin_min, int reqThreads) {
    // Set the number of threads
    #ifdef _OPENMP
        int haveCores = omp_get_num_procs();
        if(reqThreads <= 0 || reqThreads > haveCores)
            reqThreads = haveCores;
    #endif
    int nsteps = steps.n_elem;
    arma::mat U((bins+1), 2);
    U.col(0) = arma::linspace<arma::vec>(arma::min(data.col(0)),arma::max(data.col(0)),(bins+1));
    U.col(1) = arma::linspace<arma::vec>(arma::min(data.col(1)),arma::max(data.col(1)),(bins+1));
    arma::cube M1(bins, bins, 2*nsteps);
    arma::cube eM1(bins, bins, 2*nsteps);
    arma::cube M2(bins, bins, 3*nsteps);
    arma::cube D1(bins, bins, 2);
    arma::cube eD1(bins, bins, 2);
    arma::cube D2(bins, bins, 3);
    arma::mat dens(bins, bins);
    arma::cube mean_bin(bins, bins, 2);
    M1.fill(NA_REAL);
    eM1.fill(NA_REAL);
    M2.fill(NA_REAL);
    D1.fill(NA_REAL);
    eD1.fill(NA_REAL);
    D2.fill(NA_REAL);
    mean_bin.zeros();

#pragma omp parallel num_threads(reqThreads) default(shared)
{
    #pragma omp for collapse(2)
    for (int i = 0; i < bins; i++) {
        for (int j = 0; j < bins; j++) {
            arma::mat sum_m1(2, nsteps);
            arma::mat sum_m2(3, nsteps);
            arma::vec len_step(nsteps);
            sum_m1.zeros();
            sum_m2.zeros();
            len_step.zeros();
            double len_bin = 0;;
            for (int n = 0; n < data.n_rows - steps.max(); n++) {
                if(data(n,0) >= U(i,0) && data(n,0) < U(i+1,0) && data(n,1) >= U(j,1) && data(n,1) < U(j+1,1) && arma::is_finite(data(n,0)) && arma::is_finite(data(n,1))) {
                    for (int s = 0; s < nsteps; s++) {
                        if(arma::is_finite(data(n+steps(s),0)) && arma::is_finite(data(n+steps(s),1))) {
                            double inc0 = data(n+steps(s),0) - data(n,0);
                            double inc1 = data(n+steps(s),1) - data(n,1);
                            sum_m1(0,s) += inc0;
                            sum_m1(1,s) += inc1;
                            sum_m2(0,s) += inc0*inc0;
                            sum_m2(1,s) += inc0*inc1;
                            sum_m2(2,s) += inc1*inc1;
                            len_step(s)++;
                        }
                    }
                    mean_bin(i,j,0) += data(n,0);
                    mean_bin(i,j,1) += data(n,1);
                    len_bin++;
                }
            }
            mean_bin(i,j,0) /= len_bin;
            mean_bin(i,j,1) /= len_bin;
            dens(i,j) = arma::max(len_step);
            if(len_bin >= bin_min) {
                M1(arma::span(i),arma::span(j),arma::span(0,nsteps-1)) = sum_m1.row(0)/arma::trans(len_step);          // dim1
                M1(arma::span(i),arma::span(j),arma::span(nsteps,2*nsteps-1)) = sum_m1.row(1)/arma::trans(len_step);   // dim2
                M2(arma::span(i),arma::span(j),arma::span(0,nsteps-1)) = sum_m2.row(0)/arma::trans(len_step);          // M2_11
                M2(arma::span(i),arma::span(j),arma::span(nsteps,2*nsteps-1)) = sum_m2.row(1)/arma::trans(len_step);   // M2_12
                M2(arma::span(i),arma::span(j),arma::span(2*nsteps,3*nsteps-1)) = sum_m2.row(2)/arma::trans(len_step); // M2_22

                eM1(arma::span(i),arma::span(j),arma::span(0,nsteps-1)) = arma::sqrt((M2(arma::span(i),arma::span(j),arma::span(0,nsteps-1)) - arma::square(M1(arma::span(i),arma::span(j),arma::span(0,nsteps-1)))));
                eM1(arma::span(i),arma::span(j),arma::span(0,nsteps-1)) /= arma::sqrt(len_step);
                eM1(arma::span(i),arma::span(j),arma::span(nsteps,2*nsteps-1)) = arma::sqrt((M2(arma::span(i),arma::span(j),arma::span(2*nsteps,3*nsteps-1)) - arma::square(M1(arma::span(i),arma::span(j),arma::span(nsteps,2*nsteps-1)))));
                eM1(arma::span(i),arma::span(j),arma::span(nsteps,2*nsteps-1)) /= arma::sqrt(len_step);

                // linear regression with weights to get D1 and D2
                // D1
                arma::vec y = M1(arma::span(i),arma::span(j),arma::span(0,nsteps-1));
                arma::vec w = eM1(arma::span(i),arma::span(j),arma::span(0,nsteps-1));
                arma::vec coef = linreg(steps, y, 1/w);
                D1(i,j,0) = sf*coef(1);
                y = M1(arma::span(i),arma::span(j),arma::span(nsteps,2*nsteps-1));
                w = eM1(arma::span(i),arma::span(j),arma::span(nsteps,2*nsteps-1));
                coef = linreg(steps, y, 1/w);
                D1(i,j,1) = sf*coef(1);
                // D2
                // M2 - (tau * D1)^2 = tau * D2
                y = M2(arma::span(i),arma::span(j),arma::span(0,nsteps-1));
                y -= arma::square(steps*D1(i,j,0)/sf);
                coef = linreg(steps, y);
                D2(i,j,0) = sf*coef(1)/2;
                y = M2(arma::span(i),arma::span(j),arma::span(nsteps,2*nsteps-1));
                y -= arma::square(steps/sf)*D1(i,j,0)*D1(i,j,1);
                coef = linreg(steps, y);
                D2(i,j,1) = sf*coef(1)/2;
                y = M2(arma::span(i),arma::span(j),arma::span(2*nsteps,3*nsteps-1));
                y -= arma::square(steps*D1(i,j,1)/sf);
                coef = linreg(steps, y);
                D2(i,j,2) = sf*coef(1)/2;
                // calculate the error of D1
                eD1(i,j,0) = sqrt(sf*D2(i,j,0)/dens(i,j) - D1(i,j,0)*D1(i,j,0)/dens(i,j));
                eD1(i,j,1) = sqrt(sf*D2(i,j,2)/dens(i,j) - D1(i,j,1)*D1(i,j,1)/dens(i,j));
            }
        }
    }
}
    List ret;
    ret["D1"] = D1;
    ret["eD1"] = eD1;
    ret["D2"] = D2;
    ret["mean_bin"] = mean_bin;
    ret["density"] = dens;
    ret["M1"] = M1;
    ret["eM1"] = eM1;
    ret["M2"] = M2;
    ret["U"] = U;

    ret.attr("class") = CharacterVector::create("Langevin");
    return ret;
}
