/* -*- mode: C++; -*- */

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

#ifndef ALPHA_H
#define ALPHA_H

#include <vector>

class resultsAlpha {
    public:
        double alpha;
        std::vector<double> coincidence_matrix;
        std::vector<double> delta_matrix;
        double D_e;
        double D_o;
        std::vector<double> coding_values;
        std::vector<double> contributions;
        int n_coders;
        int n_units;
        int metric;
};

extern int get_alpha(
    const std::vector<double> &reliability_data,
    const int n_coders,
    const int n_units,
    const int metric,
    resultsAlpha &results);

extern double get_deltasquared(
    const double c,
    const double k,
    const std::vector<double> &coding_values,
    const std::vector<double> &contributions,
    const int metric);

#endif

/* M_PI not always defined */
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif
