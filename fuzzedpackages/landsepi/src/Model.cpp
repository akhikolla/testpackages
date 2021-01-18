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

/* REMINDER                              */
/*  | is bitwise OR,  || is logical OR   */
/*  & is bitwise AND, && is logical AND  */
/*  (int) a/b is the euclidian quotient  */
/*  (int) a%b is the euclidian remainder */
/*  ++variable == variable++ == variable=variable+1 */
/*  ++*pointer == (*pointer)++ == *pointer=*pointer+1    CAREFULL   *pointer++   DOES NOT WORK */

/* CULTIVAR GENOTYPE : mg are major resistant genes, qr are traits of quantitative resistance  */
/* PATHOGEN GENOTYPE : ig are infectivity genes,     ag are aggressiveness genes               */
/* Cultivar 0 is susceptible (no mg, no qr), Cultivars > 0 have one or more resistance sources */
/* Pathogen 0 is avirulent   (no ig, no ag), Pathogen > 0 have some ig or ag                   */

#include "Model.hpp"

// Model constructor
Model::Model(const int& Nyears, const int& time_steps_per_year, const int& Npoly, const int& Nhost, const int& Npatho,
             const int& Ngene, const std::vector<double>& area, const Vector2D<int>& rotation,
             const gsl_rng* random_generator, const std::vector<Cultivar>& cultivars, const std::vector<Gene>& genes,
             const Basic_patho& basic_patho, const std::vector<Croptype>& croptypes, const double& sigmoid_kappa_host,
             const double& sigmoid_sigma_host, const double& sigmoid_plateau_host, const double& pI0,
             const Vector2D<double>& disp_patho, const Vector2D<double>& disp_host, const int& seed)
    : Nyears(Nyears),
      time_steps_per_year(time_steps_per_year),
      Npoly(Npoly),
      Nhost(Nhost),
      Npatho(Npatho),
      Ngene(Ngene),
      area(area),
      rotation(rotation),
      random_generator(random_generator),
      cultivars(cultivars),
      genes(genes),
      basic_patho(basic_patho),
      croptypes(croptypes),
      sigmoid_kappa_host(sigmoid_kappa_host),
      sigmoid_sigma_host(sigmoid_sigma_host),
      sigmoid_plateau_host(sigmoid_plateau_host),
      pI0(pI0),
      disp_patho(disp_patho),
      disp_host(disp_host) {
    gsl_rng_set(random_generator, seed); // This function initializes (or `seeds') the random number generator.
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/* The following functions are nested in a { for poly } */

/* -------------------------------------- */
/*        MODE OF REPRODUCTION            */
/* -------------------------------------- */

/* Split reproducing infections (I) between clonal and sexual modes */
std::array<Vector2D<int>, 2> Model::split_IclonalIsex(const Vector2D<int>& I) {
    /* Random multinomial draw WITHOUT replacement */

    // Number of I exhibiting clonal reproduction in a given poly
    Vector2D<int> Iclonal_poly(this->Npatho, std::vector<int>(this->Nhost));
    // Number of I exhibiting sexual reproduction in a given poly
    Vector2D<int> Isex_poly(this->Npatho, std::vector<int>(this->Nhost));
    for(int patho = 0; patho < this->Npatho; patho++) {
        for(int host = 0; host < this->Nhost; host++) {
            Isex_poly[patho][host] = ran_binomial(this->basic_patho.repro_sex_prob, I[patho][host]);
            Iclonal_poly[patho][host] = I[patho][host] - Isex_poly[patho][host];
        }
    }
    return {Iclonal_poly, Isex_poly};
}

/* -------------------------------------- */
/*        CLONAL REPRODUCTION             */
/* -------------------------------------- */
/* Production of propagules through clonal reproduction */
/*  (i.e. same genotype as parental individual)         */
/*       Update P in a given poly                       */
void Model::reproClonal(const int& t, std::vector<int>& P, const Vector2D<int>& I, const std::vector<int>& activeR) {
    /* P and I are the numbers of individuals in a given poly */
    for(int patho = 0; patho < this->Npatho; patho++) {
        // Indicate if the pathogen has a aggressiveness gene on propagule production rate
        const std::vector<int> aggr = this->switch_patho_to_aggr(patho);
        // Poisson draw of the number of produced propagules, depending on the effect of aggressiveness on
        // propagule_prod_rate
        double Pprod_exp = 0.0;
        for(int host = 0; host < this->Nhost; host++) {
            double Pprod_tmp = this->basic_patho.propagule_prod_rate * I[patho][host];
            /* Interaction with resistance genes */
            for(int g = 0; g < this->Ngene; g++) {
                if(this->genes[g].target_trait == "PR") {
                    // Indicate if the cultivar has an active resistance gene
                    const bool Rgene = this->get_resistance(g, host, t, activeR[g]);
                    Pprod_tmp *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
                }
            }
            Pprod_exp += Pprod_tmp; // Expectation of the number of produced propagules
        }
        const int Pprod = ran_poisson(Pprod_exp);
        P[patho] += Pprod;
    }
}

/* -------------------------------------- */
/*        SEXUAL REPRODUCTION             */
/* -------------------------------------- */
/* Production of propagules through sexual reproduction (or recombination) */
/* Sex/recombination occurs between parents located in the same field and on the same host */
/* Update P in a given poly */
void Model::reproSex(const int& t, std::vector<int>& P, const Vector2D<int>& I, const std::vector<int>& activeR) {
    /* P and I are the numbers of individuals in a given poly */
    for(int host = 0; host < this->Nhost; host++) {
        int Itot = 0; // Total number of parents and couples
        for(int patho = 0; patho < this->Npatho; patho++) {
            Itot += I[patho][host];
        }
        const int coupleTot = Itot / 2;

        /* Sexual reproduction needs at least two parents */
        if(coupleTot > 0) {
            std::vector<int> pathoItot_tmp(Itot); // Vectors containing parent genotypes
            for(int patho = 0; patho < this->Npatho; patho++) {
                for(int i = 0; i < I[patho][host]; i++) {
                    pathoItot_tmp[i] = patho;
                }
            }
            const std::vector<int> pathoItot(sample(this->random_generator, pathoItot_tmp));

            for(int j = 0; j < coupleTot; j++) {
                const std::array<int, 2> pathoParents{pathoItot[j], pathoItot[j + coupleTot]};

                /* Propagule production */
                /* Reproduction rate per parent is the mean of both parents        */
                /* i.e. reproduction rate of the couple is the sum of both parents */
                double Pprod_exp = this->basic_patho.propagule_prod_rate;

                Vector2D<int> aggr_parents(2, std::vector<int>(this->Ngene));
                for(int k = 0; k < 2; k++) {
                    aggr_parents[k] = this->switch_patho_to_aggr(pathoParents[k]);
                }

                /* Interaction with resistance genes */
                for(int g = 0; g < this->Ngene; g++) {
                    if(this->genes[g].target_trait == "PR") {
                        // Indicate if the cultivar has an active resistance gene
                        const bool Rgene = this->get_resistance(g, host, t, activeR[g]);
                        Pprod_exp *= (this->genes[g].aggressiveness_matrix[aggr_parents[0][g]][Rgene] +
                                      this->genes[g].aggressiveness_matrix[aggr_parents[1][g]][Rgene]);
                    }
                }

                int Pprod = ran_poisson(Pprod_exp);
                if(pathoParents[0] == pathoParents[1]) {
                    P[pathoParents[0]] += Pprod;
                } else {
                    /* Random loci segregation for each propagule */
                    for(int s = 0; s < Pprod; s++) {
                        const int randNum = static_cast<int>(pow(2, static_cast<double>(this->Ngene)) * rng_uniform());
                        int randBin = as_binary(randNum);
                        std::vector<int> aggr_new(this->Ngene);
                        for(int locus = 0; locus < this->Ngene; locus++) { /* i.e. for each locus */
                            // Algorithm to decompose randBin in a vector of size Ngene
                            const int binFormula = randBin % 10;
                            randBin /= 10;
                            aggr_new[locus] = aggr_parents[binFormula][locus];
                        }
                        /* Increment P with the corresponding pathogen genotype */
                        const int patho_new = this->switch_aggr_to_patho(aggr_new);
                        P[patho_new]++;
                    }
                }
            }
        }
    }
}

/* -------------------------------------------------- */
/*      MUTATION OF THE PROPAGULES AT A GIVEN LOCUS   */
/* -------------------------------------------------- */

/* Update PpathoMut with mutation on "trait_mut" through multinomial draw */
void Model::mutation_locus(const int& patho, const int& trait_mut, Vector2D<int>& PpathoMut) {
    const int Nlevels = this->genes[trait_mut].Nlevels_aggressiveness;
    Vector2D<int> PaggrMut(this->Npatho, std::vector<int>(Nlevels));

    /* Mutation of trait_mut */
    for(int patho_old = 0; patho_old < this->Npatho; patho_old++) {
        const std::vector<int> aggr_old = this->switch_patho_to_aggr(patho_old);
        const int aggr_to_mutate = aggr_old[trait_mut]; // Aggressiveness index before mutation
        PaggrMut[patho_old] =
            ran_multinomial(PpathoMut[patho][patho_old], this->genes[trait_mut].mutkernel[aggr_to_mutate]);
        PpathoMut[patho][patho_old] = 0; // Re-initialisation of PpathoMut
    }

    /* Update PpathoMut with mutants */
    for(int patho_old = 0; patho_old < this->Npatho; patho_old++) {
        // Aggressiveness index relative to the different traits
        const std::vector<int> aggr_old = this->switch_patho_to_aggr(patho_old);
        for(int aggr_mut = 0; aggr_mut < Nlevels; aggr_mut++) {
            std::vector<int> aggr_new = aggr_old;
            aggr_new[trait_mut] = aggr_mut;
            const int id_patho_mut = this->switch_aggr_to_patho(aggr_new); // Pathotype index after mutation
            PpathoMut[patho][id_patho_mut] += PaggrMut[patho_old][aggr_mut]; // Add in index of new pathotype
        }
    }
}

/* -------------------------------------- */
/*               MUTATION                 */
/* -------------------------------------- */

/* Update P with mutation of propagules in a given poly */
void Model::mutation(std::vector<int>& P) {
    /* P is the numbers of individual propagules in a given poly */
    Vector2D<int> PpathoMut(this->Npatho, std::vector<int>(this->Npatho, 0));

    for(int patho = 0; patho < this->Npatho; patho++) {
        PpathoMut[patho][patho] = P[patho];
        /* Mutation of the different genes */
        for(int g = 0; g < this->Ngene; g++) {
            this->mutation_locus(patho, g, PpathoMut);
        }
    }

    /* Update P after mutation */
    for(int patho_mut = 0; patho_mut < this->Npatho; patho_mut++) {
        P[patho_mut] = 0;
        for(int patho = 0; patho < this->Npatho; patho++) {
            P[patho_mut] += PpathoMut[patho][patho_mut];
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*    functions affecting all poly    */

/* -------------------------------------- */
/*               DISPERSAL                */
/* -------------------------------------- */
/* Dispersal of new host and new pathogen propagules in the landscape */

/* Update Hjuv and P */
void Model::dispersal(const Vector2D<int>& H, Vector2D<int>& Hjuv, Vector2D<int>& P) {
    /* H, Hjuv, P are the numbers of individuals in a given poly */
    Vector3D<int> Pdisp(this->Npatho, Vector2D<int>(this->Npoly, std::vector<int>(this->Npoly)));
    Vector3D<int> Hjuvtmp(this->Nhost, Vector2D<int>(this->Npoly, std::vector<int>(this->Npoly)));

    /* Production and dispersal of the host & dispersal of pathogen propagules */
    for(int poly = 0; poly < this->Npoly; poly++) {
        /* Pathogen dispersal */
        for(int patho = 0; patho < this->Npatho; patho++) {
            Pdisp[patho][poly] = ran_multinomial(P[poly][patho], this->disp_patho[poly]);
        }
        /* Host reproduction: production and dispersal of Hjuv */
        for(int host = 0; host < this->Nhost; host++) {
            Hjuvtmp[host][poly] = ran_multinomial(
                static_cast<int>(this->cultivars[host].reproduction_rate * H[poly][host]), this->disp_host[poly]);
        }
    }

    for(int poly = 0; poly < this->Npoly; poly++) {
        /* Update the number of propagules (P) and Hjuv landing in each field */
        for(int patho = 0; patho < this->Npatho; patho++) {
            P[poly][patho] = 0;
            for(int polyE = 0; polyE < this->Npoly; polyE++) {
                P[poly][patho] += Pdisp[patho][polyE][poly];
            }
        }

        for(int host = 0; host < this->Nhost; host++) {
            Hjuv[poly][host] = 0;
            for(int polyE = 0; polyE < this->Npoly; polyE++) {
                Hjuv[poly][host] += Hjuvtmp[host][polyE][poly];
            }
        }
    }
}

/* --------------------------------------------- */
/*       BOTTLENECK AT THE END OF THE SEASON     */
/* --------------------------------------------- */
Vector3D<int> Model::bottleneck(const int& t, const Vector3D<int>& L, const Vector3D<int>& I,
                                const Vector2D<int>& activeR) {
    Vector3D<int> eqIsurv(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost)));

    for(int patho = 0; patho < this->Npatho; patho++) {
        const std::vector<int> aggr = this->switch_patho_to_aggr(patho);
        for(int host = 0; host < this->Nhost; host++) {
            for(int poly = 0; poly < this->Npoly; poly++) {
                /* Reduce the number of infected hosts (bottleneck) */
                eqIsurv[poly][patho][host] =
                    ran_binomial(this->basic_patho.survival_prob, L[poly][patho][host] + I[poly][patho][host]);
                /* Calculate the mean infectious period */
                double infectious_period_exp_new = this->basic_patho.infectious_period_exp;
                for(int g = 0; g < this->Ngene; g++) {
                    if(this->genes[g].target_trait == "IP") {
                        // Indicate if the cultivar has an active resistance gene
                        const bool Rgene = this->get_resistance(g, host, t, activeR[patho][g]);
                        infectious_period_exp_new *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
                    }
                }
                /* Calculate the equivalent number of infectious hosts */
                eqIsurv[poly][patho][host] *= static_cast<int>(infectious_period_exp_new);
            }
        }
    }
    return eqIsurv;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/* The following functions are nested in a { for poly } */

/* -------------------------- */
/*         HOST DYNAMIC       */
/* -------------------------- */

/* Compute host reproduction, death and growth and updtate the number of H in the concerned poly */
void Model::host_dynamic(const int& poly, const int& year, std::vector<int>& H, std::vector<int>& Hjuv,
                         const Vector2D<int>& L, const Vector2D<int>& I, const Vector2D<int>& R, std::vector<int>& N) {
    /* H, Hjuv, L, I and R are the number of host in a given poly */

    // If there is no rotation (same croptype each year), only take the first year rotation
    const int id_croptype = (this->rotation[poly].size() == 1) ? this->rotation[poly][0] : this->rotation[poly][year];

    for(std::pair<int, double> cultivar_prop : this->croptypes[id_croptype].cultivar_proportion) {
        const int id_host = cultivar_prop.first;
        const double prop = cultivar_prop.second;

        /* Calculation of totals for L, I, R and K */
        // Carrying capacity of the cultivar in the concerned paddock
        const int K = static_cast<int>(this->area[poly] * this->cultivars[id_host].max_density * prop);
        int L_host = 0, I_host = 0, R_host = 0;
        for(int patho = 0; patho < this->Npatho; patho++) {
            L_host += L[patho][id_host];
            I_host += I[patho][id_host];
            R_host += R[patho][id_host];
        }
        N[id_host] = H[id_host] + L_host + I_host + R_host; // Number of occupied sites

        /* HOST MORTALITY: H2M */
        /* ------------------- */
        const int H2M = ran_binomial(this->cultivars[id_host].death_rate, H[id_host]);
        const int H2Mjuv = ran_binomial(this->cultivars[id_host].death_rate, Hjuv[id_host]);
        Hjuv[id_host] -= H2Mjuv; // Update Hjuv

        /* HOST REPRODUCTION: Hnewind */
        /* -------------------------- */
        /* Hjuv settlement in the field */
        int availSites = K - N[id_host] - H2M; // Number of available sites
        if(availSites < 0) { /* Security */
            availSites = 0;
        }
        const double f1host = sigmoid(this->sigmoid_plateau_host, this->sigmoid_kappa_host, this->sigmoid_sigma_host,
                                      availSites / static_cast<double>(K));
        int siteaccess = ran_binomial(f1host, availSites);
        const int Hnewind = (siteaccess < Hjuv[id_host]) ? siteaccess : Hjuv[id_host];
  
        /* HOST GROWTH: H2H */
        /* ---------------- */
        int H2H = static_cast<int>(this->cultivars[id_host].growth_rate * (H[id_host] - H2M) *
                                   (1 - ((N[id_host] - H2M + Hnewind) / static_cast<double>(K))));
        if(H2H < 0) { /* Security */
        Rcpp::Rcerr << "hostID" << id_host << " growrate " << this->cultivars[id_host].growth_rate << " H " << H[id_host] 
                    << " H2M " << H2M << " N " << N[id_host] << " HnewId " << Hnewind << " K " << static_cast<double>(K) 
                    << std::endl;
            Rprintf("CAREFUL ! H2H < 0, one of the areas may be 0: check if Npoly, NpolyTot and idLAN are correct\n");
            H2H = 0;
        } else if((N[id_host] - H2M + Hnewind + H2H) > K) {
            Rprintf("CAREFUL ! H2H too big\n");
            H2H = K - (N[id_host] - H2M + Hnewind);
        }

        /* UPDATE NUMBER OF HOSTS */
        /* ---------------------- */
        H[id_host] = H[id_host] - H2M + Hnewind + H2H;
        N[id_host] = N[id_host] - H2M + Hnewind + H2H;
    }
}

/* ------------------------------------------------------------------ */
/*         CONTAMINATION : propagule deposition on healthy sites      */
/* ------------------------------------------------------------------ */

/* Calculation of the number of contaminated sites */
Vector2D<int> Model::contamination(const std::vector<int>& H, const std::vector<int>& P, const std::vector<int>& N) {
    /* H, P are the number of individuals in a given poly */
    Vector2D<int> Hcontaminated(this->Npatho, std::vector<int>(this->Nhost));

    /* Calculation of total for H and P */
    const int Htot = std::accumulate(H.begin(), H.end(), 0);
    const int Ptot = std::accumulate(P.begin(), P.end(), 0);

    std::vector<double> probaH(this->Nhost + 1); // Probability for H to belong to each host
    double probaHtot = 0;
    for(int host = 0; host < this->Nhost; host++) {
        probaH[host] = (Htot == 0) ? 0 : static_cast<double>(H[host]) / static_cast<double>(Htot);
        // (cond) ? res1 : res2;
        probaHtot += probaH[host];
    }
    probaH[this->Nhost] = 1 - probaHtot;

    std::vector<double> probaP(this->Npatho + 1); // Probability to belong to the different pathotypes
    double probaPtot = 0;
    for(int patho = 0; patho < this->Npatho; patho++) {
        // If there is no propagule, set everything to 0, else put each proportion
        probaP[patho] = (Ptot == 0) ? 0 : static_cast<double>(P[patho]) / static_cast<double>(Ptot);
        probaPtot += probaP[patho];
    }
    probaP[this->Npatho] = 1 - probaPtot;

    // Distribution of the propagules among the different cultivars
    const std::vector<int> P_host = ran_multinomial(Ptot, probaH);

    for(int host = 0; host < this->Nhost; host++) {
        // Distribution of the propagules among the different pathotypes
        const std::vector<int> P_host_patho = ran_multinomial(P_host[host], probaP);
        /* Calculation of the number of contaminable sites */
        const double f1patho = (N[host] > 0)
                                   ? sigmoid(this->basic_patho.sigmoid_plateau, this->basic_patho.sigmoid_kappa,
                                             this->basic_patho.sigmoid_sigma, (H[host] / static_cast<double>(N[host])))
                                   : 0.0;

        // Contaminable site, where a propagule may deposit
        const int Hcontaminable = ran_binomial(f1patho, H[host]);
        // Distribution of the contaminable sites among the different pathotypes
        const std::vector<int> Hcontaminated_tmp = ran_multinomial(Hcontaminable, probaP);
        /* The true number of contaminated sites is the minimum between sites and propagules */
        for(int patho = 0; patho < this->Npatho; patho++) {
            Hcontaminated[patho][host] = std::min(Hcontaminated_tmp[patho], P_host_patho[patho]);
        }
    }
    return Hcontaminated;
}

/* ----------------------------------------------------------- */
/*         INFECTIOUS CYCLE : transitions H -> L -> I -> R     */
/* ----------------------------------------------------------- */

/* Calculate the number of contaminated sites that become infected, infectious or removed and update H, L, I, R */
void Model::infection(const int& t, std::vector<int>& H, const Vector2D<int>& Hcontaminated, Vector2D<int>& L,
                      Vector2D<int>& I, Vector2D<int>& R, Vector3D<int>& L2I, Vector3D<int>& I2R,
                      const std::vector<int>& activeR) {
    /* H, Hcontaminated, L, I, R, L2I and I2R are the number of individuals in a given poly */
    for(int patho = 0; patho < this->Npatho; patho++) {
        const std::vector<int> aggr = this->switch_patho_to_aggr(patho);
        for(int host = 0; host < this->Nhost; host++) {
            /* Infection of healthy sites: H2L */
            /* ------------------------------- */
            double infection_rate_exp = this->basic_patho.infection_rate;
            /* Interaction with resistance genes */
            for(int g = 0; g < this->Ngene; g++) {
                if(this->genes[g].target_trait == "IR") {
                    // Indicate if the cultivar has an active resistance gene
                    bool Rgene = this->get_resistance(g, host, t, activeR[g]);
                    infection_rate_exp *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
                }
            }
            const int H2L = ran_binomial(infection_rate_exp, Hcontaminated[patho][host]);

            /* Latent and infectious periods */
            /* ----------------------------- */
            /* Find parameters of gamma distributions from mean and variance */
            double latent_period_exp_new = this->basic_patho.latent_period_exp;
            double infectious_period_exp_new = this->basic_patho.infectious_period_exp;
            std::array<double, 2> latent_period_alpha;
            std::array<double, 2> infectious_period_alpha;
            for(int g = 0; g < this->Ngene; g++) {
                /* Latent period */
                if(this->genes[g].target_trait == "LAT") {
                    // Indicate if the cultivar has an active resistance gene
                    bool Rgene = this->get_resistance(g, host, t, activeR[g]);
                    latent_period_exp_new /= (this->genes[g].aggressiveness_matrix[aggr[g]][Rgene] +
                                              0.001 * (this->genes[g].aggressiveness_matrix[aggr[g]][Rgene] == 0));
                }
                /* Security to avoid problem in alpha calculation */
                latent_period_exp_new += 0.001 * (latent_period_exp_new == 0);

                /* Infectious period */
                if(this->genes[g].target_trait == "IP") {
                    // Indicate if the cultivar has an active resistance gene
                    bool Rgene = this->get_resistance(g, host, t, activeR[g]);
                    infectious_period_exp_new *= this->genes[g].aggressiveness_matrix[aggr[g]][Rgene];
                }
                // Security to avoid problem in alpha calculation
                infectious_period_exp_new += 0.001 * (infectious_period_exp_new == 0);
            }
            latent_period_alpha = find_paramGamma(latent_period_exp_new, this->basic_patho.latent_period_var);
            infectious_period_alpha = find_paramGamma(infectious_period_exp_new, this->basic_patho.infectious_period_var);

            /* Recently infected hosts */
            for(int h2l = 0; h2l < H2L; h2l++) {
                /* Latent period */
                const int lag1 = static_cast<int>(ran_gamma(latent_period_alpha[0], latent_period_alpha[1]));
                if((t + lag1) < this->time_steps_per_year) {
                    L2I[patho][host][t + lag1]++;
                }
                /* Infectious period */
                const int lag2 = lag1 + static_cast<int>(ran_gamma(infectious_period_alpha[0], infectious_period_alpha[1]));
                if((t + lag2) < this->time_steps_per_year) {
                    I2R[patho][host][t + lag2]++;
                }
            }

            /* Update H, L, I, R */
            /* ----------------- */
            H[host] -= H2L;
            L[patho][host] += H2L;
            L[patho][host] -= L2I[patho][host][t];
            I[patho][host] += L2I[patho][host][t];
            I[patho][host] -= I2R[patho][host][t];
            R[patho][host] += I2R[patho][host][t];
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

/* --------------------------- */
/*   DYNAMIC OF THE EPIDEMIC   */
/* --------------------------- */

void Model::dynepi() {
    char name_fH[20], name_fHjuv[20], name_fP[20], name_fL[20], name_fI[20], name_fR[20];
    Vector2D<int> Hcontaminated; // Contaminated sites (where a propagule is deposited)
    Vector2D<int> Hjuv, P;
    Vector3D<int> L(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost)));
    Vector3D<int> I(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost)));
    Vector3D<int> R(this->Npoly, Vector2D<int>(this->Npatho, std::vector<int>(this->Nhost)));
    std::vector<int> N(this->Nhost);
    Vector4D<int> L2I;
    Vector4D<int> I2R;

    /* Initialisation (t=0) */
    this->init_HjuvPLIR(Hjuv, P, L, I, R);
    this->init_L2I2R(L2I, I2R);

    Vector2D<int> H = this->intro_H(0);
    Vector2D<int> activeR = this->init_activeR();
    intro_I(H, I, I2R, activeR);

    for(int year = 1; year < (this->Nyears + 1); year++) {
        Rprintf("----------------------------- YEAR %d -----------------------------\n", year);
        /* Create the files to write the output */
        sprintf(name_fH, "H-%02d.bin", year);
        sprintf(name_fHjuv, "Hjuv-%02d.bin", year);
        sprintf(name_fP, "P-%02d.bin", year);
        sprintf(name_fL, "L-%02d.bin", year);
        sprintf(name_fI, "I-%02d.bin", year);
        sprintf(name_fR, "R-%02d.bin", year);

        FILE* fH = fopen(name_fH, "wb");
        FILE* fHjuv = fopen(name_fHjuv, "wb");
        FILE* fL = fopen(name_fL, "wb");
        FILE* fI = fopen(name_fI, "wb");
        FILE* fP = fopen(name_fP, "wb");
        FILE* fR = fopen(name_fR, "wb");

        /* Loop for all the timesteps of the cropping season */
        for(int t = 1; t < this->time_steps_per_year; t++) {
            // Writing model output for timestep t
            this->write_HHjuvPLIR(H, Hjuv, P, L, I, R, fH, fHjuv, fP, fL, fI, fR);
            P = Vector2D<int>(this->Npoly, std::vector<int>(this->Npatho, 0)); // Re-initialisation at 0
            for(int poly = 0; poly < this->Npoly; poly++) {
                const std::array<Vector2D<int>, 2> splited_I(this->split_IclonalIsex(I[poly]));
                this->reproClonal(t, P[poly], splited_I[0], activeR[poly]);
                this->reproSex(t, P[poly], splited_I[1], activeR[poly]);
                this->mutation(P[poly]);
            }
            this->dispersal(H, Hjuv, P);

            for(int poly = 0; poly < this->Npoly; poly++) {
                this->host_dynamic(poly, year - 1, H[poly], Hjuv[poly], L[poly], I[poly], R[poly], N);
                Hcontaminated = this->contamination(H[poly], P[poly], N);
                this->infection(t, H[poly], Hcontaminated, L[poly], I[poly], R[poly], L2I[poly], I2R[poly],
                                activeR[poly]);
            }
        }

        /* Last time-step of the season: bottleneck before starting a new season */

        // Writing model output for last timestep
        this->write_HHjuvPLIR(H, Hjuv, P, L, I, R, fH, fHjuv, fP, fL, fI, fR);
        // Calculation of the equivalent number of I that survive and produce propagules for the next season
        const Vector3D<int> eqIsurv = this->bottleneck(this->time_steps_per_year, L, I, activeR);

        /* Re-initialisation at 0 */
        this->init_HjuvPLIR(Hjuv, P, L, I, R);
        this->init_L2I2R(L2I, I2R);

        /* Generate P issued from eqIsurv = (remaining L+I) * Tspo */
        for(int poly = 0; poly < this->Npoly; poly++) {
            const std::array<Vector2D<int>, 2> splited_I(this->split_IclonalIsex(eqIsurv[poly]));
            this->reproClonal(this->time_steps_per_year, P[poly], splited_I[0], activeR[poly]);
            this->reproSex(this->time_steps_per_year, P[poly], splited_I[1], activeR[poly]);
            this->mutation(P[poly]);
        }
        this->dispersal(H, Hjuv, P);

        /* Re-plantation --> regenerate H */
        H = this->intro_H(year);
        activeR = this->init_activeR();

        /* Infection of newly planted hosts to generate the inoculum of the next season */
        for(int poly = 0; poly < this->Npoly; poly++) {
            /* N = H[poly] in beginning of next season */
            Hcontaminated = this->contamination(H[poly], P[poly], H[poly]);
            this->infection(0, H[poly], Hcontaminated, L[poly], I[poly], R[poly], L2I[poly], I2R[poly], activeR[poly]);
        }

        fclose(fH);
        fclose(fHjuv);
        fclose(fL);
        fclose(fI);
        fclose(fP);
        fclose(fR);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void model_landsepi(Rcpp::List time_param, Rcpp::NumericVector area_vector, Rcpp::NumericMatrix rotation_matrix, Rcpp::NumericMatrix croptypes_cultivars_prop, Rcpp::List dispersal, Rcpp::List inits, int seed,
                    Rcpp::List cultivars_param, Rcpp::List basic_patho_param, Rcpp::List genes_param) {

    /*------------------*/
    /* Seed generation  */
    /*------------------*/
    /* Makoto Matsumoto and Takuji Nishimura, “Mersenne Twister: A 623-dimensionally equidistributed uniform
     * pseudorandom number generator”. ACM Transactions on Modeling and Computer Simulation, Vol. 8, No. 1 (Jan. 1998),
     * Pages 3–30 */
    const gsl_rng_type* gen_type = gsl_rng_mt19937;
    // This function returns a pointer to a newly-created instance of a random number generator of type gen_type
    const gsl_rng* gen = gsl_rng_alloc(gen_type);

    /*-------------------*/
    /*  Time parameters  */
    /*-------------------*/
    const int Nyears = Rcpp::as<int>(time_param["Nyears"]);
    const int nTSpY = Rcpp::as<int>(time_param["nTSpY"]);

    /*-------------------------------------*/
    /*  Landscape & deployment parameters  */
    /*-------------------------------------*/

    std::vector<double> area = Rcpp::as<std::vector<double>>(area_vector);
    Vector2D<int> rotation(0);
    for( int i = 0; i < rotation_matrix.nrow() ; i++) {
      Rcpp::NumericVector poly_rotation = rotation_matrix.row(i);
      //poly_rotation.push_back(poFeature_rotation->GetFieldAsInteger(year));
      rotation.push_back(Rcpp::as<std::vector<int>>(poly_rotation));
    }
    
    const int Npoly = area.size();
    
    std::vector<Croptype> croptypes(0);
    
    for(int i=0 ; i < croptypes_cultivars_prop.nrow() ; i++) {
      if( croptypes_cultivars_prop(i,0) >= croptypes.size()) croptypes.push_back(Croptype());
      croptypes[croptypes_cultivars_prop(i,0)].cultivar_proportion.push_back(
          std::pair<int, double>(
              croptypes_cultivars_prop(i,1),
              croptypes_cultivars_prop(i,2) ));                        
    }

    /*------------------------*/
    /*  Dispersal parameters  */
    /*------------------------*/
    const std::vector<double> disp_patho_tmp = Rcpp::as<std::vector<double>>(dispersal["disp_patho"]);
    const std::vector<double> disp_host_tmp = Rcpp::as<std::vector<double>>(dispersal["disp_host"]);
    Vector2D<double> disp_patho(Npoly, std::vector<double>(Npoly));
    Vector2D<double> disp_host(Npoly, std::vector<double>(Npoly));

    for(int i = 0; i < Npoly; i++) {
        for(int j = 0; j < Npoly; j++) {
            disp_patho[j][i] = disp_patho_tmp[j + i * Npoly];
            disp_host[j][i] = disp_host_tmp[j + i * Npoly];
        }
    }
    
    /*----------------------*/
    /*  Initial conditions  */
    /*----------------------*/
    const double pI0 = Rcpp::as<double>(inits["pI0"]);

    /*-------------------*/
    /*  Host parameters  */
    /*-------------------*/
    //const int Nhost = Rcpp::as<int>(cultivars_param["Nhost"]);
    const std::vector<double> initial_density = Rcpp::as<std::vector<double>>(cultivars_param["initial_density"]);
    const std::vector<double> max_density = Rcpp::as<std::vector<double>>(cultivars_param["max_density"]);
    const std::vector<double> growth_rate = Rcpp::as<std::vector<double>>(cultivars_param["growth_rate"]);
    const std::vector<double> reproduction_rate = Rcpp::as<std::vector<double>>(cultivars_param["reproduction_rate"]);
    const std::vector<double> death_rate = Rcpp::as<std::vector<double>>(cultivars_param["death_rate"]);
    const int Nhost = initial_density.size();
    Rcpp::List cultivars_genes_list = Rcpp::as<Rcpp::List>(cultivars_param["cultivars_genes_list"]);
    std::vector<Cultivar> cultivars;
    std::vector<int> total_genes_id(0); // Contains all the genes_id used by the cultivars
    for(int i = 0; i < Nhost; i++) {
        const std::vector<int> genes_id = Rcpp::as<std::vector<int>>(cultivars_genes_list[i]);
        cultivars.push_back(Cultivar(initial_density[i], max_density[i], growth_rate[i], reproduction_rate[i],
                                     death_rate[i], genes_id));
        total_genes_id.insert(total_genes_id.end(), genes_id.begin(), genes_id.end()); // Appends a vector to another
    }

    // Remove duplicate value in "total_genes_id"
    std::sort(total_genes_id.begin(), total_genes_id.end());
    auto last = std::unique(total_genes_id.begin(), total_genes_id.end());
    total_genes_id.erase(last, total_genes_id.end());

    double sigmoid_kappa_host = Rcpp::as<double>(cultivars_param["sigmoid_kappa_host"]);
    // security to avoid kappa = 0
    sigmoid_kappa_host += 1E-6 * (sigmoid_kappa_host == 0);
    const double sigmoid_sigma_host = Rcpp::as<double>(cultivars_param["sigmoid_sigma_host"]);
    const double sigmoid_plateau_host = Rcpp::as<double>(cultivars_param["sigmoid_plateau_host"]);

    /*---------------------*/
    /* Pathogen parameters */
    /*---------------------*/
    const double infection_rate = Rcpp::as<double>(basic_patho_param["infection_rate"]);
    const double propagule_prod_rate = Rcpp::as<double>(basic_patho_param["propagule_prod_rate"]);
    const double latent_period_exp = Rcpp::as<double>(basic_patho_param["latent_period_exp"]);
    double latent_period_var = Rcpp::as<double>(basic_patho_param["latent_period_var"]);
    // security to avoid variance = 0
    latent_period_var += 1E-6 * (latent_period_var == 0);
    const double infectious_period_exp = Rcpp::as<double>(basic_patho_param["infectious_period_exp"]);
    double infectious_period_var = Rcpp::as<double>(basic_patho_param["infectious_period_var"]);
    // security to avoid variance = 0
    infectious_period_var += 1E-6 * (infectious_period_var == 0);
    const double survival_prob = Rcpp::as<double>(basic_patho_param["survival_prob"]);
    const double repro_sex_prob = Rcpp::as<double>(basic_patho_param["repro_sex_prob"]);
    double sigmoid_kappa = Rcpp::as<double>(basic_patho_param["sigmoid_kappa"]);
    // security to avoid kappa = 0
    sigmoid_kappa += 1E-6 * (sigmoid_kappa == 0);
    const double sigmoid_sigma = Rcpp::as<double>(basic_patho_param["sigmoid_sigma"]);
    const double sigmoid_plateau = Rcpp::as<double>(basic_patho_param["sigmoid_plateau"]);

    Basic_patho basic_patho(infection_rate, propagule_prod_rate, latent_period_exp, latent_period_var,
                            infectious_period_exp, infectious_period_var, survival_prob, repro_sex_prob, sigmoid_kappa,
                            sigmoid_sigma, sigmoid_plateau);

    /*------------------------*/
    /*  Evolution parameters  */
    /*------------------------*/
    const std::vector<double> time_to_activ_exp = Rcpp::as<std::vector<double>>(genes_param["time_to_activ_exp"]);
    const std::vector<double> time_to_activ_var = Rcpp::as<std::vector<double>>(genes_param["time_to_activ_var"]);
    const std::vector<int> Nlevels_aggressiveness = Rcpp::as<std::vector<int>>(genes_param["Nlevels_aggressiveness"]);
    const std::vector<std::string> target_trait = Rcpp::as<std::vector<std::string>>(genes_param["target_trait"]);
    const std::vector<double> mutation_prob = Rcpp::as<std::vector<double>>(genes_param["mutation_prob"]);
    const std::vector<double> efficiency = Rcpp::as<std::vector<double>>(genes_param["efficiency"]);
    const std::vector<double> fitness_cost = Rcpp::as<std::vector<double>>(genes_param["fitness_cost"]);
    const std::vector<double> tradeoff_strength = Rcpp::as<std::vector<double>>(genes_param["tradeoff_strength"]);
    
    // ! BE CAREFUL ! 
    // total_genes_id -> GENE ID in DB
    // genes_param -> liste sorted by GENE ID order ?
    std::vector<Gene> genes(0);
    //for( int g=0; g < total_genes_id.size(); g++) {
    for(int g : total_genes_id) { // Only keep the gene used by the cultivars (Reduce Npatho)
        genes.push_back(Gene(time_to_activ_exp[g], time_to_activ_var[g], Nlevels_aggressiveness[g], target_trait[g],
                             mutation_prob[g], efficiency[g], fitness_cost[g], tradeoff_strength[g]));
    }
    const int Ngene = genes.size();
    
    int Npatho = 1;
    for(int g = 0; g < Ngene; g++) {
        Npatho *= genes[g].Nlevels_aggressiveness;
    }

    // Create the model
    Model model(Nyears, nTSpY, Npoly, Nhost, Npatho, Ngene, area, rotation, gen, cultivars, genes, basic_patho,
                croptypes, sigmoid_kappa_host, sigmoid_sigma_host, sigmoid_plateau_host, pI0, disp_patho, disp_host,
                seed);
    
    /*--------------------------------------*/
    /* Write and Print the model parameters */
    /*--------------------------------------*/
    model.print_param(seed, mutation_prob, efficiency, fitness_cost, tradeoff_strength);

    /* -------------- */
    /* Epidemic model */
    /* -------------- */
    Rprintf("\n*** SPATIOTEMPORAL MODEL SIMULATING THE SPREAD AND EVOLUTION OF A PATHOGEN IN A LANDSCAPE ***\n\n");
    model.dynepi();
}
