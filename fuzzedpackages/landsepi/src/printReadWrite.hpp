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

#ifndef __PRINT_READ_WRITE__
#define __PRINT_READ_WRITE__

#include <stdio.h> // fwrite, fprintf
#include <array>
#include <string>
#include <vector>

#include "Cultivar.hpp"
#include "Model.hpp"
#include "functions.hpp" // sum2_2, sum1_3...

/* ************************************************************************* */
/*                         printReadWrite.c                                  */
/* ************************************************************************* */

class Model; // Forward declaration of Model because printReadWrite compile before Model

/* Print a vector */
template <typename T>
void print_1d(FILE* f, const std::vector<T>& t, const std::string& title) {
    if(title != "") {
        fprintf(f, "%s : ", title.c_str());
    }
    for(unsigned int i = 0; i < t.size(); i++) {
        if(typeid(T) == typeid(int)) {
            fprintf(f, "%7d", t[i]);
        } else if(typeid(T) == typeid(double)) {
            fprintf(f, "%.3f ", t[i]);
        } else {
            fprintf(f, "NaN");
        }
    }
    fprintf(f, "\n");
}

/* Print an array */
template <typename T, size_t n>
void print_1d(FILE* f, const std::array<T, n>& t, const std::string& title) {
    if(title != "") {
        fprintf(f, "%s : \n", title.c_str());
    }
    for(unsigned int i = 0; i < n; i++) {
        if(typeid(T) == typeid(int)) {
            fprintf(f, "%7d ", t[i]);
        } else if(typeid(T) == typeid(double)) {
            fprintf(f, "%.3f ", t[i]);
        } else {
            fprintf(f, "NaN ");
        }
    }
    fprintf(f, "\n");
}

/* Print a table */
template <typename T>
void print_2d(FILE* f, const Vector2D<T>& t, const std::string& title) {
    if(title != "") {
        fprintf(f, "%s : \n", title.c_str());
    }
    for(unsigned int i = 0; i < t.size(); i++) {
        print_1d(f, t[i], "");
    }
    fprintf(f, "\n");
}

/* Print a 3-dimension table */
template <typename T>
void print_3d(FILE* f, const Vector3D<T>& t, const std::string& title) {
    if(title != "") {
        fprintf(f, "%s : \n", title.c_str());
    }
    for(unsigned int i = 0; i < t.size(); i++) {
        print_2d(f, t[i], "");
    }
    fprintf(f, "\n");
}
/* Print the sum of the 1st dimension a table of integer of dimension 3 */
void print_i2sum2(FILE* f, const Vector2D<int>& t, const std::string& title);
/* Print the sum of the 1st dimension a table of float of dimension 3 */
void print_d2sum2(FILE* f, const Vector2D<double>& t, const std::string& title);
/* Print the sum of the 1st dimension a table of integer of dimension 3 */
void print_i3sum1(FILE* f, const int& z, const int& l, const int& c, const Vector3D<int>& t, const std::string& title);
/* Print the sum of the 1st dimension a table of float of dimension 3 */
void print_d3sum1(FILE* f, const int& z, const int& l, const int& c, const Vector3D<double>& t,
                  const std::string& title);
#endif
