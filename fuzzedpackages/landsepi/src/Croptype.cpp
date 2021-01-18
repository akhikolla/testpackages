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

#include "Croptype.hpp"

Croptype::Croptype() : cultivar_proportion(std::vector<std::pair<int, double>>(0)) {}
Croptype::Croptype(std::vector<std::pair<int, double>>& cultivar_proportion)
    : cultivar_proportion(cultivar_proportion) {}

// Transform Croptype attributs to string
std::string Croptype::to_string() const {
    std::string str("");
    str += "cultivar_proportion {cultivar_id, proportion}:\n";
    for(const std::pair<int, double>& p : this->cultivar_proportion) {
        str += "    {" + std::to_string(p.first) + ", " + std::to_string(p.second) + "}\n";
    }

    return str;
}
