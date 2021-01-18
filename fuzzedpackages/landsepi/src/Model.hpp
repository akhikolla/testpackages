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

#ifndef __MODEL_SIMUL__
#define __MODEL_SIMUL__

#include <Rcpp.h>
#include <gsl/gsl_randist.h> // gsl_ran_binomial, gsl_ran_poisson...
#include <gsl/gsl_rng.h> // gsl_rng_mt19937, gsl_rng*...
#include <math.h> // pow...
#include <stdio.h> // FILE*...
#include <array>
#include <fstream> // ofstream
#include <string>
#include <vector>

#include "Basic_patho.hpp"
#include "Croptype.hpp"
#include "Cultivar.hpp"
#include "Gene.hpp"
#include "functions.hpp" // find_paramGamma, sigmoid, sample...
#include "printReadWrite.hpp" // write_HHjuvPLIR

struct Model {
    const int Nyears;               // Duration of the simulation (number of cropping seasons)
    const int time_steps_per_year;  // Number of timesteps per cropping season
    const int Npoly;                // Number of polygons (i.e. pieces of land where crops may be cultivated)
    const int Nhost;                // Number of host genotypes
    const int Npatho;               // Number of pathogen genotypes
    const int Ngene;                // Number of genes
    const std::vector<double> area;        // Vector of areas of polygons
    const Vector2D<int> rotation;          // Array of croptypes for each polygon and each year
    const gsl_rng* random_generator;       // Randomness generator
    const std::vector<Cultivar> cultivars; // Array of available cultivars
    const std::vector<Gene> genes;         // Array of available genes
    const Basic_patho basic_patho;         // Aggressiveness of a wild-type pathogen on a S cultivar
    const std::vector<Croptype> croptypes; // List of pairs {cultivar, relative proportion} for each croptype
    const double sigmoid_kappa_host;       // Kappa parameter for the sigmoid invasion function (for host dispersal)
    const double sigmoid_sigma_host;       // Sigma parameter for the sigmoid invasion function (for host dispersal)
    const double sigmoid_plateau_host;     // Plateau parameter for the sigmoid invasion function (for host dispersal)
    const double pI0;                      // Initial inoculum: probability to be I at t=0
    const Vector2D<double> disp_patho;     // Pathogen dispersal matrix
    const Vector2D<double> disp_host;      // Host dispersal matrix

    Model(const int& Nyears, const int& time_steps_per_year, const int& Npoly, const int& Nhost, const int& Npatho,
          const int& Ngene, const std::vector<double>& area, const Vector2D<int>& rotation,
          const gsl_rng* random_generator, const std::vector<Cultivar>& cultivars, const std::vector<Gene>& genes,
          const Basic_patho& basic_patho, const std::vector<Croptype>& croptypes, const double& sigmoid_kappa_host,
          const double& sigmoid_sigma_host, const double& sigmoid_plateau_host, const double& pI0,
          const Vector2D<double>& disp_patho, const Vector2D<double>& disp_host, const int& seed);

    void dynepi();
    void infection(const int& t, std::vector<int>& H, const Vector2D<int>& Hcontaminated, Vector2D<int>& L,
                   Vector2D<int>& I, Vector2D<int>& R, Vector3D<int>& L2I, Vector3D<int>& I2R,
                   const std::vector<int>& activeR);
    Vector2D<int> contamination(const std::vector<int>& H, const std::vector<int>& P, const std::vector<int>& N);
    void host_dynamic(const int& poly, const int& year, std::vector<int>& H, std::vector<int>& Hjuv,
                      const Vector2D<int>& L, const Vector2D<int>& I, const Vector2D<int>& R, std::vector<int>& N);
    Vector3D<int> bottleneck(const int& t, const Vector3D<int>& L, const Vector3D<int>& I,
                             const Vector2D<int>& activeQR);
    void dispersal(const Vector2D<int>& H, Vector2D<int>& Hjuv, Vector2D<int>& P);
    void mutation(std::vector<int>& P);
    void mutation_locus(const int& patho, const int& trait_mut, Vector2D<int>& PpathoMut);
    void reproSex(const int& t, std::vector<int>& P, const Vector2D<int>& I, const std::vector<int>& activeQR);
    void reproClonal(const int& t, std::vector<int>& P, const Vector2D<int>& I, const std::vector<int>& activeR);
    std::array<Vector2D<int>, 2> split_IclonalIsex(const Vector2D<int>& I);
    bool get_resistance(const int& index_gene, const int& host, const int& t, const int& activeR);

    /* Init functions */
    int switch_aggr_to_patho(const std::vector<int>& aggr);
    std::vector<int> switch_patho_to_aggr(const int& index_patho);
    Vector2D<int> init_activeR();
    void init_HjuvPLIR(Vector2D<int>& Hjuv, Vector2D<int>& P, Vector3D<int>& L, Vector3D<int>& I, Vector3D<int>& R);
    void init_L2I2R(Vector4D<int>& L2I, Vector4D<int>& I2R);
    Vector2D<int> intro_H(const int& year);
    void intro_I(Vector2D<int>& H, Vector3D<int>& I, Vector4D<int>& I2R, const Vector2D<int>& activeR);

    /* Random functions */
    double rng_uniform();
    int ran_poisson(const double& mu);
    double ran_gamma(const double& a, const double& b);
    int ran_binomial(const double& p, const int& n);
    std::vector<int> ran_multinomial(const int& N, const std::vector<double>& p);

    /* Print parameters in an output .txt file */
    void print_param(const int& seed, const std::vector<double>& mutation_prob, const std::vector<double>& efficiency,
                     const std::vector<double>& fitness_cost, const std::vector<double>& tradeoff_strength);

    /* Write model output in .txt files and print output on screen */
    void write_HHjuvPLIR(const Vector2D<int>& H, const Vector2D<int>& Hjuv, const Vector2D<int>& P,
                         const Vector3D<int>& L, const Vector3D<int>& I, const Vector3D<int>& R, FILE* fH, FILE* fHjuv,
                         FILE* fP, FILE* fL, FILE* fI, FILE* fR);
};

//' @title Model for Landscape Epidemiology & Evolution
//' @name model_landsepi
//' @description Stochastic, spatially-explicit, demo-genetic model simulating the spread and evolution of a 
//' plant pathogen in a heterogeneous landscape.
//' 
//' @param time_param list of simulation parameters:\itemize{ 
//' \item Nyears = number cropping seasons, 
//' \item nTSpY = number of time-steps per cropping season.
//' }
//' @param area_vector a vector containing areas of polygons (i.e. fields), in surface units.
//' @param rotation_matrix a matrix containing for each field (rows) and year (columns, named "year_1", "year_2", etc.), 
//' the index of the cultivated croptype. Importantly, the matrix must contain 1 more column than the real number 
//' of simulated years. 
//' @param croptypes_cultivars_prop a matrix with three columns named 'croptypeID' for croptype index, 
//' 'cultivarID' for cultivar index and 'proportion' for the proportion of the cultivar within the croptype. 
//' @param dispersal list of dispersal parameters:\itemize{ 
//' \item disp_patho = vectorised dispersal matrix of the pathogen, 
//' \item disp_host = vectorised dispersal matrix of the host.
//' }
//' @param inits list of initial conditions:\itemize{
//' \item pI0 = initial probability for the first host (whose index is 0) to be infectious (i.e. state I) at t=0.
//' }
//' @param seed seed (for random number generation).
//' @param cultivars_param list of parameters associated with each host genotype (i.e. cultivars) 
//' when cultivated in pure crops:\itemize{   
//' \item initial_density = vector of host densities (per surface unit) at the beginning of the cropping season,  
//' \item max_density = vector of maximum host densities (per surface unit) at the end of the cropping season, 
//' \item growth rate = vector of host growth rates, 
//' \item reproduction rate = vector of host reproduction rates, 
//' \item death rate = vector of host death rates, 
//' \item sigmoid_kappa_host = kappa parameter for the sigmoid invasion function (for host dispersal),
//' \item sigmoid_sigma_host = sigma parameter for the sigmoid invasion function (for host dispersal),
//' \item sigmoid_plateau_host = plateau parameter for the sigmoid invasion function (for host dispersal),
//' \item cultivars_genes_list = a list containing, for each host genotype, the indices of carried resistance genes.
//' } 
//' @param basic_patho_param list of pathogen aggressiveness parameters on a susceptible host 
//' for a pathogen genotype not adapted to resistance: \itemize{
//' \item infection_rate = maximal expected infection rate of a propagule on a healthy host, 
//' \item propagule_prod_rate = maximal expected reproduction_rate of an infectious host per timestep, 
//' \item latent_period_exp = minimal expected duration of the latent period, 
//' \item latent_period_var = variance of the latent period duration, 
//' \item infectious_period_exp = maximal expected duration of the infectious period, 
//' \item infectious_period_var = variance of the infectious period duration,
//' \item survival_prob = probability for a propagule to survive the off-season, 
//' \item repro_sex_prob = probability for an infectious host to reproduce via sex rather than via cloning, 
//' \item sigmoid_kappa = kappa parameter of the sigmoid contamination function, 
//' \item sigmoid_sigma = sigma parameter of the sigmoid contamination function, 
//' \item sigmoid_plateau = plateau parameter of the sigmoid contamination function.
//' }
//' @param genes_param list of parameters associated with each resistance gene and with the evolution of 
//' each corresponding pathogenicity gene:\itemize{ 
//' \item target_trait = vector of aggressiveness components (IR, LAT, IP, or PR) targeted by resistance genes, 
//' \item efficiency = vector of resistance gene efficiencies (percentage of reduction of the targeted 
//' aggressiveness component: IR, 1/LAT, IP and PR), 
//' \item time_to_activ_exp = vector of expected delays to resistance activation (for APRs), 
//' \item time_to_activ_var = vector of variances of the delay to resistance activation (for APRs),  
//' \item mutation_prob = vector of mutation probabilities for pathogenicity genes (each of them corresponding to a resistance gene), 
//' \item Nlevels_aggressiveness = vector of number of adaptation levels related to each resistance gene (i.e. 1 + number 
//' of required mutations for a pathogenicity gene to fully adapt to the corresponding resistance gene), 
//' \item fitness_cost = vector of fitness penalties paid by pathogen genotypes fully adapted 
//' to the considered resistance genes on hosts that do not carry this gene, 
//' \item tradeoff_strength = vector of strengths of the trade-off relationships between the level of aggressiveness 
//' on hosts that do and do not carry the resistance genes.
//' }
//' 
//' @details See ?landsepi for details on the model and assumptions. 
//' Briefly, the model is stochastic, spatially explicit (the basic spatial unit is an individual field), based on a SEIR
//' (‘susceptible-exposed-infectious-removed’, renamed HLIR for 'healthy-latent-infectious-removed' to avoid confusions 
//' with 'susceptible host') structure with a discrete time step. It simulates the spread and 
//' evolution of a pathogen in a heterogeneous cropping landscape, across cropping seasons split by host harvests which impose
//' potential bottlenecks to the pathogen. A wide array of resistance deployment strategies can be simulated.
//'  
//' @return A set of binary files is generated for every year of simulation and every compartment: 
//' \itemize{
//'  \item H: healthy hosts,
//'  \item Hjuv: juvenile healthy hosts,
//'  \item L: latently infected hosts,
//'  \item I: infectious hosts,
//'  \item R: removed hosts,
//'  \item P: propagules.}
//' Each file indicates for every time-step the number of individuals in each field, and when 
//' appropriate for each host and pathogen genotypes).
//' 
//' @examples
//' \dontrun{
//' #### Spatially-implicit simulation with 2 patches (S + R) during 3 years ####
//' 
//' ## Simulation parameters
//' time_param <- list(Nyears=3, nTSpY=120)
//' Npoly=2
//' Npatho=2
//' area <- c(100000, 100000)
//' cultivars <- as.list(rbind(loadCultivar(name="Susceptible", type="growingHost")
//' , loadCultivar(name="Resistant", type="growingHost")))
//' names(cultivars)[names(cultivars)=="cultivarName"] <- "name"
//' cultivars <- c(cultivars, list(sigmoid_kappa_host=0.002, sigmoid_sigma_host=1.001,
//'  sigmoid_plateau_host=1, cultivars_genes_list=list(numeric(0),0)))
//' rotation <- data.frame(year_1=c(0,1), year_2=c(0,1), year_3=c(0,1), year_4=c(0,1))
//' croptypes_cultivars_prop <- data.frame(croptypeID=c(0,1), cultivarID=c(0,1), proportion=c(1,1))
//' genes <- as.list(loadGene(name="MG", type="majorGene"))
//'     
//' ## run simulation
//' model_landsepi(seed=1,
//'                time_param = time_param,
//'                basic_patho_param = loadPathogen(disease = "rust"),
//'                inits = list(pI0=0.01), area_vector = area,
//'                dispersal = list(disp_patho=c(0.99,0.01,0.01,0.99),
//'                disp_host=c(1,0,0,1)),
//'                rotation_matrix = as.matrix(rotation),
//'                croptypes_cultivars_prop = as.matrix(croptypes_cultivars_prop),
//'                cultivars_param = cultivars, genes_param = genes)  
//' 
//' ## Compute outputs
//' eco_param <- list(yield_perHa = cbind(H = as.numeric(cultivars$yield_H),
//'              L = as.numeric(cultivars$yield_L),
//'              I = as.numeric(cultivars$yield_I),
//'              R = as.numeric(cultivars$yield_R)),
//'              production_cost_perHa = as.numeric(cultivars$production_cost),
//'              market_value = as.numeric(cultivars$market_value))
//'             
//' evol_res <- evol_output(, time_param, Npoly, cultivars, genes)
//' epid_output(, time_param, Npatho, area, rotation
//' , croptypes_cultivars_prop, cultivars, eco_param)
//' 
//' 
//' #### 1-year simulation of a rust epidemic in pure susceptible crop in a single 1 ha patch ####
//' time_param <- list(Nyears=1, nTSpY=120)
//' Npoly=1
//' Npatho=1
//' area <- c(100000)
//' cultivars <- as.list(rbind(loadCultivar(name="Susceptible", type="growingHost")))
//' names(cultivars)[names(cultivars)=="cultivarName"] <- "name"
//' cultivars <- c(cultivars, list(sigmoid_kappa_host=0.002, sigmoid_sigma_host=1.001,
//'                                sigmoid_plateau_host=1, cultivars_genes_list=list(numeric(0))))
//' rotation <- data.frame(year_1=c(0), year_2=c(0))
//' croptypes_cultivars_prop <- data.frame(croptypeID=c(0), cultivarID=c(0), proportion=c(1))
//' genes <-   list(geneName = character(0) , fitness_cost = numeric(0)
//' , mutation_prob = numeric(0)
//' , efficiency = numeric(0) , tradeoff_strength = numeric(0)
//' , Nlevels_aggressiveness = numeric(0)
//' , time_to_activ_exp = numeric(0) , time_to_activ_var = numeric(0)
//' , target_trait = character(0))
//'     
//' ## run simulation
//' model_landsepi(seed=1, time_param = time_param
//' , basic_patho_param = loadPathogen(disease = "rust"),
//' inits = list(pI0=0.01), area_vector = area, dispersal = list(disp_patho=c(1), disp_host=c(1)),
//' rotation_matrix = as.matrix(rotation),
//' croptypes_cultivars_prop = as.matrix(croptypes_cultivars_prop),
//' cultivars_param = cultivars,  genes_param = genes) 
//' }
//' @references Rimbaud L., Papaïx J., Rey J.-F., Barrett L. G. and Thrall P. H. (2018).
//' Assessing the durability andefficiency of landscape-based strategies to deploy plant 
//' resistance to pathogens. \emph{PLoS Computational Biology} 14(4):e1006067.
//' 
//' @export
// [[Rcpp::export]]
void model_landsepi(Rcpp::List time_param
                        , Rcpp::NumericVector area_vector
                        , Rcpp::NumericMatrix rotation_matrix
                        , Rcpp::NumericMatrix croptypes_cultivars_prop
                        , Rcpp::List dispersal
                        , Rcpp::List inits
                        , int seed
                        , Rcpp::List cultivars_param
                        , Rcpp::List basic_patho_param
                        , Rcpp::List genes_param);
#endif
