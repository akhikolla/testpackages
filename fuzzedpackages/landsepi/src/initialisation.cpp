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

/****************************************************************/
/*                    initialisation.c                          */
/****************************************************************/
#include "Model.hpp"

/* Switch from indices of aggressiveness to index pathotype */
int Model::switch_aggr_to_patho(const std::vector<int>& aggr) {
    int index_patho = 0;

    for(int i = 0; i < this->Ngene; i++) {
        int prod = 1;
        for(int j = i + 1; j < this->Ngene; j++) {
            prod *= this->genes[j].Nlevels_aggressiveness;
        }
        index_patho += aggr[i] * prod;
    }
    return index_patho;
}

/* Switch from index pathotype to indices of aggressiveness (vector of dimension Ngene) */
std::vector<int> Model::switch_patho_to_aggr(const int& index_patho) {
    std::vector<int> aggr(this->Ngene);

    int remainder = index_patho;
    for(int i = 0; i < this->Ngene; i++) {
        int prod = 1;
        for(int j = i + 1; j < this->Ngene; j++) {
            prod *= this->genes[j].Nlevels_aggressiveness;
        }
        aggr[i] = remainder / prod; /* Quotient */
        remainder = remainder % prod;
    }
    return aggr;
}

/* Initialisation of activeR: matrix indicating for each poly and each gene, at which time step resistance starts to be
 * active */
Vector2D<int> Model::init_activeR() {
    Vector2D<int> activeR(this->Npoly, std::vector<int>(this->Ngene));
    // note: if Ngene==0, activeR is of size 0 and will never be used
    for(int poly = 0; poly < this->Npoly; poly++) {
        for(int gene = 0; gene < this->Ngene; gene++) {
            const double exp = this->genes[gene].time_to_activ_exp;
            const double var = this->genes[gene].time_to_activ_var;
            if(exp > 0 && var > 0) {
                const std::array<double, 2> alpha = find_paramGamma(exp, var);
                activeR[poly][gene] = static_cast<int>(ran_gamma(alpha[0], alpha[1]));
            } else if (var == 0){
                activeR[poly][gene] = (int) exp;  /* Deterministic delay if variance==0 */
            } else {
                activeR[poly][gene] = 0;
            }
        }
    }
    return activeR;
}

/* Indicates if a cultivar carries a given gene and this one is active at a given timestep t */
bool Model::get_resistance(const int& index_gene, const int& host, const int& t, const int& activeR) {
    std::vector<int> genes_id = this->cultivars[host].genes_id;
    if(std::find(genes_id.begin(), genes_id.end(), index_gene) != genes_id.end()) {
        if(t >= activeR) {
            return 1;
        }
    }
    return 0;
}

/* Initialisation of H, L, I, R and P at 0 */
void Model::init_HjuvPLIR(Vector2D<int>& Hjuv, Vector2D<int>& P, Vector3D<int>& L, Vector3D<int>& I, Vector3D<int>& R) {
    Hjuv = Vector2D<int>(this->Npoly, std::vector<int>(this->Nhost, 0));
    P = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0));
    L = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost, 0)));
    I = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost, 0)));
    R = Vector3D<int>(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost, 0)));
}

/* Initialise L2I and I2R with 0 */
void Model::init_L2I2R(Vector4D<int>& L2I, Vector4D<int>& I2R) {
    L2I = Vector4D<int>(
        this->Npoly,
        Vector3D<int>(this->Npatho, Vector2D<int>(this->Nhost, std::vector<int>(this->time_steps_per_year, 0))));
    I2R = Vector4D<int>(
        this->Npoly,
        Vector3D<int>(this->Npatho, Vector2D<int>(this->Nhost, std::vector<int>(this->time_steps_per_year, 0))));
}

/* Plantation of H in the beginning of a season */
Vector2D<int> Model::intro_H(const int& year) {
    Vector2D<int> H(this->Npoly, std::vector<int>(this->Nhost));
    for(int poly = 0; poly < this->Npoly; poly++) {
        // If there is no rotation (same croptype each year)
        int id_croptype = (this->rotation[poly].size() == 1) ? this->rotation[poly][0] : this->rotation[poly][year];

        for(std::pair<int, double> cultivar_prop : this->croptypes[id_croptype].cultivar_proportion) {
            int id_host = cultivar_prop.first;
            double prop = cultivar_prop.second;
            H[poly][id_host] = static_cast<int>(this->area[poly] * this->cultivars[id_host].initial_density * prop);
        }
    }
    return H;
}

/* Pathogen introduction : infectious sites from pathotype 0 in cultivar 0 */
void Model::intro_I(Vector2D<int>& H, Vector3D<int>& I, Vector4D<int>& I2R, const Vector2D<int>& activeR) {
    const int patho = 0; // Introduced pathotype is 0
    const std::vector<int> aggr = switch_patho_to_aggr(patho);
    const int host = 0; // Introduction on cultivar 0
    const int t = 0; // Introduction at t=0

    for(int poly = 0; poly < this->Npoly; poly++) {
        /* Infection */
        double infection_rate_exp = this->pI0;
        for(int g = 0; g < this->Ngene; g++) {
            if(this->genes[g].target_trait == "IR") {
                // Indicate if the cultivar has an active resistance gene
                bool Rgene = get_resistance(g, host, t, activeR[poly][g]);
                infection_rate_exp *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
            }
        }
        I[poly][patho][host] = ran_binomial(infection_rate_exp, H[poly][host]);
        H[poly][host] -= I[poly][patho][host];

        /* Calculation of infectious period for infected sites at t=0 */
        for(int i = 0; i < I[poly][patho][host]; i++) {
            double infectious_period_exp_new = this->basic_patho.infectious_period_exp;
            for(int g = 0; g < this->Ngene; g++) {
                if(this->genes[g].target_trait == "IP") {
                    // Indicate if the cultivar has an active resistance gene
                    bool Rgene = get_resistance(g, host, t, activeR[poly][g]);
                    infectious_period_exp_new *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
                }
            }
            // Security to avoid problem in alpha calculation
            infectious_period_exp_new += 0.001 * (infectious_period_exp_new == 0);
            std::array<double, 2> infectious_period_alpha =
                find_paramGamma(infectious_period_exp_new, this->basic_patho.infectious_period_var);
            int lag = static_cast<int>(ran_gamma(infectious_period_alpha[0], infectious_period_alpha[1]));
            lag += 1 * (lag == 0); // Security to avoid infectious duration of 0 day
            if(lag < this->time_steps_per_year) {
                I2R[poly][patho][host][lag]++;
            }
        }
    }
}
