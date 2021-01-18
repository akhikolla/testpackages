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

/* ************************************************************************* */
/*                         printReadWrite.c                                  */
/* ************************************************************************* */
#include "printReadWrite.hpp"

/* Print the sum of the 1st dimension a table of integer of dimension 3 */
void print_i2sum2(FILE* f, const Vector2D<int>& t, const std::string& title) {
    fprintf(f, "%s = %5d\n", title.c_str(), sum2_2(t));
}

/* Print the sum of the 1st dimension a table of float of dimension 3 */
void print_d2sum2(FILE* f, const Vector2D<double>& t, const std::string& title) {
    fprintf(f, "%s = %.3f\n", title.c_str(), sum2_2(t));
}

/* Print the sum of the 1st dimension a table of integer of dimension 3 */
void print_i3sum1(FILE* f, const int& z, const int& l, const int& c, const Vector3D<int>& t, const std::string& title) {
    print_2d(f, sum1_3(z, l, c, t), title);
}

/* Print the sum of the 1st dimension a table of float of dimension 3 */
void print_d3sum1(FILE* f, const int& z, const int& l, const int& c, const Vector3D<double>& t,
                  const std::string& title) {
    print_2d(f, sum1_3(z, l, c, t), title);
}

/* Print the parameters in an output .txt file */
void Model::print_param(const int& seed, const std::vector<double>& mutation_prob,
                        const std::vector<double>& efficiency, const std::vector<double>& fitness_cost,
                        const std::vector<double>& tradeoff_strength) {
    std::ofstream param_file("parameters.txt");
    std::ofstream landscape_file("param_landscape.txt");
//    std::ofstream genes_file("param_genes.txt");
    std::ofstream dispP_file("param_disp_patho.txt");
//    std::ofstream dispH_file("param_disp_host.txt");

    param_file << "###     MODEL PARAMETERS     ###\n";
    param_file << "seed:                " << seed << "\n";

    param_file << "\n*****             Seasonality                 *****\n";
    param_file << "Nyears:              " << this->Nyears << "\n";
    param_file << "time_steps_per_year: " << this->time_steps_per_year << "\n";
    param_file << "pI0:                 " << this->pI0 << "\n";

    param_file << "\n*****     Landscape & deployment strategy     *****\n";
    param_file << "Npoly:               " << this->Npoly << "\n";

    landscape_file << "area, year 1 rotation, year 2 rotation... : \n";
    dispP_file << "Pathogen dispersal:\n";
//    dispH_file << "Host dispersal:\n";
    for(int i = 0; i < this->Npoly; i++) {
        landscape_file << "  poly " << i + 1 << ": ";
        landscape_file << area[i] << " / ";
        for(const int& croptype_id : rotation[i]) {
            landscape_file << croptype_id << " ";
        }
        landscape_file << "\n";

        dispP_file << i + 1 << ": ";
//        dispH_file << i + 1 << ": ";
        for(int j = 0; j < this->Npoly; j++) {
            dispP_file << std::to_string(this->disp_patho[i][j]) << " ";
//            dispH_file << std::to_string(this->disp_host[i][j]) << " ";
        }
        dispP_file << "\n";
//        dispH_file << "\n";
    }
    
    param_file << "Croptypes:\n";
    for(unsigned int i = 0; i < this->croptypes.size(); i++) {
        param_file << "  croptype " << i + 1 << ": " << this->croptypes[i].to_string() << "\n";
    }

    param_file << "\n*****                  Hosts                  *****\n";
    param_file << "Nhost:               " << this->Nhost << "\n";
    param_file << "cultivars: \n";
    for(int i = 0; i < this->Nhost; i++) {
        param_file << "  cultivar " << i + 1 << ":\n" << this->cultivars[i].to_string();
    }
    param_file << "sigmoid_kappa_host:   " << this->sigmoid_kappa_host << "\n";
    param_file << "sigmoid_sigma_host:   " << this->sigmoid_sigma_host << "\n";
    param_file << "sigmoid_plateau_host: " << this->sigmoid_plateau_host << "\n";

    param_file << "\n*****                  Genes                  *****\n";
    param_file << "Ngene: " << this->Ngene << "\n";
    for(int i = 0; i < Ngene; i++) {
        param_file << "Gene " << i + 1 << ":\n" << this->genes[i].to_string();
        param_file << "  mutation_prob:          " << mutation_prob[i] << "\n";
        param_file << "  efficiency:             " << efficiency[i] << "\n";
        param_file << "  fitness_cost:           " << fitness_cost[i] << "\n";
        param_file << "  tradeoff_strength:      " << tradeoff_strength[i] << "\n\n";
    }
    
    param_file << "\n*****                 Pathogen                *****\n";
    param_file << "Npatho: " << this->Npatho << "\n";
    param_file << "basic_patho:\n" << this->basic_patho.to_string();
}

/* Write model output in .txt files and print output on screen */
void Model::write_HHjuvPLIR(const Vector2D<int>& H, const Vector2D<int>& Hjuv, const Vector2D<int>& P,
                            const Vector3D<int>& L, const Vector3D<int>& I, const Vector3D<int>& R, FILE* fH,
                            FILE* fHjuv, FILE* fP, FILE* fL, FILE* fI, FILE* fR) {
    for(int poly = 0; poly < Npoly; poly++) {
        for(int patho = 0; patho < Npatho; patho++) {
            fwrite(&P[poly][patho], sizeof(int), 1, fP);
            for(int host = 0; host < Nhost; host++) {
                fwrite(&L[poly][patho][host], sizeof(int), 1, fL);
                fwrite(&I[poly][patho][host], sizeof(int), 1, fI);
                fwrite(&R[poly][patho][host], sizeof(int), 1, fR);
            }
        }
        for(int host = 0; host < Nhost; host++) {
            fwrite(&Hjuv[poly][host], sizeof(int), 1, fHjuv);
            fwrite(&H[poly][host], sizeof(int), 1, fH);
        }
    }
}
