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

#include "Cultivar.hpp"

// Default constructor
Cultivar::Cultivar()
    : initial_density(0.0),
      max_density(0.0),
      growth_rate(0.0),
      reproduction_rate(0.0),
      death_rate(0.0),
      genes_id(std::vector<int>()) {}

// Constructor
Cultivar::Cultivar(const double& initial_density, const double& max_density, const double& growth_rate,
                   const double& reproduction_rate, const double& death_rate, const std::vector<int>& genes_id)
    : initial_density(initial_density),
      max_density(max_density),
      growth_rate(growth_rate),
      reproduction_rate(reproduction_rate),
      death_rate(death_rate),
      genes_id(genes_id) {}

// Transform Cultivar attributs to string
std::string Cultivar::to_string() const {
    std::string str("");
    str += "    initial_density:   " + std::to_string(this->initial_density) + "\n";
    str += "    max_density:       " + std::to_string(this->max_density) + "\n";
    str += "    growth_rate:       " + std::to_string(this->growth_rate) + "\n";
    str += "    reproduction_rate: " + std::to_string(this->reproduction_rate) + "\n";
    str += "    death_rate:        " + std::to_string(this->death_rate) + "\n";
    str += "    genes_id:          ";
    for(const int& id : this->genes_id) {
        str += std::to_string(id) + " ";
    }
    str += "\n\n";

    return str;
}
