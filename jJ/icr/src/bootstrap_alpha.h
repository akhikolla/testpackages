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

#ifndef BOOTSTRAP_ALPHA_H
#define BOOTSTRAP_ALPHA_H

#include <vector>

/* bootstrap functions */
extern int bootstrap_alpha(
    const double D_e,
    const std::vector<double> &reliability_data,
    const int nC,
    const int nU,
    const std::vector<double> &coding_values,
    const std::vector<double> &contributions,
    const int metric,
    const int bootstraps,
    const unsigned long seed[6],
    const int n_threads,
    std::vector<double> &alphas);

extern int bootstrap_alpha_nonparametric(
    const std::vector<double> &reliability_data,
    const int nC,
    const int nU,
    const int metric,
    const int bootstraps,
    const unsigned long seed[6],
    const int n_threads,
    std::vector<double> &alphas);

#endif
