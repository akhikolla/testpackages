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
#ifndef __GENE__
#define __GENE__

#include "functions.hpp"

class Gene {
    /* Initialisation of mutation matrix (only called by Gene constructor) */
    Vector2D<double> init_mutkernel(const double& mutation_prob);
    /* Initialisation of infectivity and aggressiveness matrices (only called by Gene constructor) */
    Vector2D<double> init_aggressiveness_matrix(const double& efficiency, const double& fitness_cost,
                                                const double& tradeoff_strength);

public:
    const double time_to_activ_exp;               /* Expected delay to resistance activation (for APRs) */
    const double time_to_activ_var;               /* Variance of the delay to resistance activation (for APRs) */
    const int Nlevels_aggressiveness;             /* Number of adaptation levels related to a given aggressiveness component */
    const std::string target_trait;               /* Aggressiveness component targeted by resistance (IR, LAT, IP, PR) */
    const Vector2D<double> mutkernel;             /* Mutation matrix (for pathogen evolution) */
    const Vector2D<double> aggressiveness_matrix; /* Aggressiveness matrix (plant-pathogen interaction) */

    Gene();
    Gene(const double& time_to_activ_exp, const double& time_to_activ_var, const int& Nlevels_aggressiveness,
         const std::string& target_trait, const double& mutation_prob, const double& efficiency,
         const double& fitness_cost, const double& tradeoff_strength);
    std::string to_string() const;
};

#endif
