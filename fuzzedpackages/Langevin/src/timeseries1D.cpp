// Copyright (C) 2012 - 2015  Philip Rinn
// Copyright (C) 2012  David Bastine
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

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".timeseries1D")]]
NumericVector timeseries1D(const unsigned int& N, const double& startpoint,
                           const double& d13, const double& d12,
                           const double& d11, const double& d10,
                           const double& d22, const double& d21,
                           const double& d20, const double& sf,
                           double dt) {
    NumericVector ts(N, NA_REAL);
    ts(0) = startpoint;

    // Calculate the integration time step and related values
    double stime = 1.0/sf;
    if(stime < dt || dt == 0) {
        dt = stime;
    }
    // Ration between sampling time and integration time step
    unsigned int m = ceil(stime/dt);
    dt = stime/m;

    // Integration
    double gamma;
    double x = ts(0);

    for (unsigned int i = 0; i < N; i++)  {
        // Integrate m steps and just save the last
        for (unsigned int j = 0; j < m; j++) {
            // Get a single Gaussian random number
            gamma = rnorm(1,0,std::sqrt(2))[0];
            // Iterate with integration time step dt
            x += (d13*std::pow(x,3) + d12*std::pow(x,2) + d11*x + d10)*dt + std::sqrt((d22*std::pow(x,2) + d21*x + d20)*dt)*gamma;
        }
        // Save every mth step
        ts(i) = x;
    }
    ts.attr("class") = CharacterVector::create("ts");
    ts.attr("tsp") = NumericVector::create(1, 1 + (N - 1)/sf, sf);
    return ts;
}
