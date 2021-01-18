// Copyright (C) 2015  Philip Rinn
// Copyright (C) 2015  Pedro G. Lind
// Copyright (C) 2015  Carl von Ossietzy Universität Oldenburg
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

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".timeseries2D")]]
NumericMatrix timeseries2D(const unsigned int& N, const double& startpointx,
                  const double& startpointy, const NumericMatrix& D1_1,
                  const NumericMatrix& D1_2, const NumericMatrix& g_11,
                  const NumericMatrix& g_12, const NumericMatrix& g_21,
                  const NumericMatrix& g_22, const double& sf, double dt) {
    NumericMatrix ts(N, 2);
    std::fill(ts.begin(), ts.end(), NA_REAL);
    ts(0,0) = startpointx;
    ts(0,1) = startpointy;

    // Calculate the integration time step and related values
    double stime = 1.0/sf;
    if(stime < dt || dt == 0) {
        dt = stime;
    }
    // Ration between sampling time and integration time step
    unsigned int m = ceil(stime/dt);
    dt = stime/m;

    // Integration
    double gamma1, gamma2;
    double sumD1_1, sumD1_2;
    double sumG_11, sumG_12, sumG_21, sumG_22;
    double x = ts(0,0);
    double y = ts(0,1);

    for (unsigned int i = 0; i < N; i++)  {
        // Integrate m steps and just save the last
        for (unsigned int j = 0; j < m; j++) {
            // Get a single Gaussian random number
            gamma1 = rnorm(1,0,std::sqrt(2))[0];
            gamma2 = rnorm(1,0,std::sqrt(2))[0];
            // Iterate with integration time step dt
            sumD1_1 = D1_1(0,0) + D1_1(1,0)*x + D1_1(0,1)*y + D1_1(2,0)*std::pow(x,2) + D1_1(0,2)*std::pow(y,2) + D1_1(3,0)*std::pow(x,3) + D1_1(0,3)*std::pow(y,3) + D1_1(1,1)*x*y + D1_1(2,1)*std::pow(x,2)*y + D1_1(1,2)*x*std::pow(y,2);
            sumD1_2 = D1_2(0,0) + D1_2(1,0)*x + D1_2(0,1)*y + D1_2(2,0)*std::pow(x,2) + D1_2(0,2)*std::pow(y,2) + D1_2(3,0)*std::pow(x,3) + D1_2(0,3)*std::pow(y,3) + D1_2(1,1)*x*y + D1_2(2,1)*std::pow(x,2)*y + D1_2(1,2)*x*std::pow(y,2);
            sumG_11 = g_11(0,0) + g_11(1,0)*x + g_11(0,1)*y + g_11(2,0)*std::pow(x,2) + g_11(0,2)*std::pow(y,2) + g_11(1,1)*x*y;
            sumG_12 = g_12(0,0) + g_12(1,0)*x + g_12(0,1)*y + g_12(2,0)*std::pow(x,2) + g_12(0,2)*std::pow(y,2) + g_12(1,1)*x*y;
            sumG_21 = g_21(0,0) + g_21(1,0)*x + g_21(0,1)*y + g_21(2,0)*std::pow(x,2) + g_21(0,2)*std::pow(y,2) + g_21(1,1)*x*y;
            sumG_22 = g_22(0,0) + g_22(1,0)*x + g_22(0,1)*y + g_22(2,0)*std::pow(x,2) + g_22(0,2)*std::pow(y,2) + g_22(1,1)*x*y;
            x += sumD1_1*dt + std::sqrt(sumG_11*dt)*gamma1 + std::sqrt(sumG_12*dt)*gamma2;
            y += sumD1_2*dt + std::sqrt(sumG_21*dt)*gamma1 + std::sqrt(sumG_22*dt)*gamma2;
        }
        // Save every mth step
        ts(i,0) = x;
        ts(i,1) = y;
    }
    ts.attr("class") = CharacterVector::create("mts", "ts", "matrix");
    ts.attr("tsp") = NumericVector::create(1, 1 + (N - 1)/sf, sf);
    return ts;
}
