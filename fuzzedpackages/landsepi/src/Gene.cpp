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

#include "Gene.hpp"

/* Initialisation of mutation matrix (only called by Gene constructor) */
Vector2D<double> Gene::init_mutkernel(const double& mutation_prob) {
    /* Mutation matrix for pathogen genes */
    Vector2D<double> mutkernel(Nlevels_aggressiveness, std::vector<double>(Nlevels_aggressiveness));
    for(int i = 0; i < Nlevels_aggressiveness; i++) {
        mutkernel[i][i] = 1 - mutation_prob; // Diagonal
    }
    if(Nlevels_aggressiveness > 1) {
        mutkernel[0][1] = mutation_prob; // Because only one possible direction for mutation
        mutkernel[Nlevels_aggressiveness - 1][Nlevels_aggressiveness - 2] = mutation_prob;
        for(int i = 1; i < Nlevels_aggressiveness - 1; i++) { // Two possible directions
            mutkernel[i][i + 1] = mutation_prob / 2;
            mutkernel[i][i - 1] = mutation_prob / 2;
        }
    }
    return mutkernel;
}

/* Initialisation of aggressiveness matrix (only called by Gene constructor) */
Vector2D<double> Gene::init_aggressiveness_matrix(const double& efficiency, const double& fitness_cost,
                                                  const double& tradeoff_strength) {
    /* Aggressiveness matrix */
    Vector2D<double> aggressiveness_matrix(Nlevels_aggressiveness, std::vector<double>(2));
    const double aggressiveness_0 = 1 - efficiency;

    /* (Nlevels_aggressiveness-1): intervals between Nlevels_aggressiveness values */
    const double step = (Nlevels_aggressiveness > 1) ? 1 / static_cast<double>(Nlevels_aggressiveness - 1) : 0;

    std::vector<double> gain(Nlevels_aggressiveness);
    for(int i = 0; i < Nlevels_aggressiveness; i++) {
        gain[i] = i * step;
    }

    const std::vector<double> cost = tradeoff(gain, tradeoff_strength);
    for(int i = 0; i < Nlevels_aggressiveness; i++) {
        aggressiveness_matrix[i][1] = aggressiveness_0 + gain[i] * efficiency;
        aggressiveness_matrix[i][0] = 1 - cost[i] * fitness_cost;
    }
    return aggressiveness_matrix;
}

// Default constructor
Gene::Gene()
    : time_to_activ_exp(0.0),
      time_to_activ_var(0.0),
      Nlevels_aggressiveness(0),
      target_trait(""),
      mutkernel(Vector2D<double>()),
      aggressiveness_matrix(Vector2D<double>()) {}

// Constructor
Gene::Gene(const double& time_to_activ_exp, const double& time_to_activ_var, const int& Nlevels_aggressiveness,
           const std::string& target_trait, const double& mutation_prob, const double& efficiency,
           const double& fitness_cost, const double& tradeoff_strength)
    : time_to_activ_exp(time_to_activ_exp),
      time_to_activ_var(time_to_activ_var),
      Nlevels_aggressiveness(Nlevels_aggressiveness),
      target_trait(target_trait),
      mutkernel(init_mutkernel(mutation_prob)),
      aggressiveness_matrix(init_aggressiveness_matrix(efficiency, fitness_cost, tradeoff_strength)) {}

// Transform Gene attributs to string
std::string Gene::to_string() const {
    std::string str("");
    str += "  time_to_activ_exp:      " + std::to_string(this->time_to_activ_exp) + "\n";
    str += "  time_to_activ_var:      " + std::to_string(this->time_to_activ_var) + "\n";
    str += "  Nlevels_aggressiveness: " + std::to_string(this->Nlevels_aggressiveness) + "\n";
    str += "  target_trait:           " + this->target_trait + "\n";
    str += "  mutkernel:\n";
    for(int i = 0; i < this->Nlevels_aggressiveness; i++) {
        str += "    " + std::to_string(i) + ": ";
        for(int j = 0; j < this->Nlevels_aggressiveness; j++) {
            str += std::to_string(mutkernel[i][j]) + " ";
        }
        str += "\n";
    }
    str += "  aggressiveness_matrix:\n";
    for(int i = 0; i < this->Nlevels_aggressiveness; i++) {
        str += "    " + std::to_string(i) + ": ";
        str += std::to_string(aggressiveness_matrix[i][0]) + " ";
        str += std::to_string(aggressiveness_matrix[i][1]) + "\n";
    }
    return str;
}
