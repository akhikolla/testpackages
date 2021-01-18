#include "Simulator.h"
#include <string.h>

//' Simulates species tree using constant rate birth-death process
//'
//' @description Forward simulates to a number of tips. This function does so using
//'     the general algorithm of Hartmann et al. 2010.
//' @param sbr species birth rate (i.e. speciation rate)
//' @param sdr species death rate (i.e. extinction rate)
//' @param numbsim number of species trees to simulate
//' @param n_tips number of tips to simulate to
//' @param gsa_stop_mult number of tips to simulate the GSA tip to
//' @return List of objects of the tree class (as implemented in APE)
//' @references
//' K. Hartmann, D. Wong, T. Stadler. Sampling trees from evolutionary models.
//'     Syst. Biol., 59(4): 465-476, 2010.
//'
//' T. Stadler. Simulating trees on a fixed number of extant species.
//'     Syst. Biol., 60: 676-684, 2011.
//' @examples
//' mu <- 0.5 # death rate
//' lambda <- 2.0 # birth rate
//' numb_replicates <- 10
//' numb_extant_tips <- 4
//' # simulate trees under the GSA so first simulates a tree with
//' # numb_extant_tips * 100 tips counting each time we have a tree with 10 tips
//' # then randomly picks one of those trees
//'
//' tree_list <- sim_sptree_bdp(sbr = lambda,
//'                 sdr = mu,
//'                 numbsim = numb_replicates,
//'                 n_tips = numb_extant_tips)
// [[Rcpp::export]]
Rcpp::List sim_sptree_bdp(SEXP sbr,
                          SEXP sdr,
                          SEXP numbsim,
                          Rcpp::NumericVector n_tips,
                          Rcpp::NumericVector gsa_stop_mult = 10){
    double sbr_ = as<double>(sbr);
    double sdr_ = as<double>(sdr);
    unsigned numbsim_ = as<int>(numbsim);
    unsigned n_tips_ = as<int>(n_tips);
    unsigned gsa_stop_ = as<int>(gsa_stop_mult);
    unsigned gsa_stop = gsa_stop_ * n_tips_;
    RNGScope scope;
    if(sbr_ <= 0.0)
        stop("'sbr' must be bigger than 0.0.");
    if(sbr_ < sdr_)
        stop("'sbr' must be greater than 'sdr'");
    if(numbsim_ < 1)
        stop("'numbsim' must be more than 1");
    if(sdr_ < 0.0)
        stop("'sdr' must be 0.0 or greater.");
    if(n_tips_ < 1)
        stop("'n_tips' must be greater than 1.");
    if(gsa_stop_ < 1)
        stop("'gsa_stop_mult' must be greater than 1.");
    return bdsim_species_tree(sbr_, sdr_, numbsim_, n_tips_, gsa_stop);
}
//' Simulates species tree using constant rate birth-death process to a time
//'
//' @description Forward simulates a tree until a provided time is reached.
//' @param sbr species birth rate (i.e. speciation rate)
//' @param sdr species death rate (i.e. extinction rate)
//' @param numbsim number of species trees to simulate
//' @param t time to simulate to
//' @return List of objects of the tree class (as implemented in APE)
//' @references
//' K. Hartmann, D. Wong, T. Stadler. Sampling trees from evolutionary models.
//'     Syst. Biol., 59(4): 465-476, 2010.
//'
//' T. Stadler. Simulating trees on a fixed number of extant species.
//'     Syst. Biol., 60: 676-684, 2011.
//' @examples
//' mu <- 0.5 # death rate
//' lambda <- 2.0 # birth rate
//' numb_replicates <- 10
//' time <- 1
//'
//' tree_list <- sim_sptree_bdp_time(sbr = lambda,
//'                 sdr = mu,
//'                 numbsim = numb_replicates,
//'                 t = time)
// [[Rcpp::export]]
Rcpp::List sim_sptree_bdp_time(SEXP sbr, SEXP sdr, SEXP numbsim, SEXP t){
    double sbr_ = as<double>(sbr);
    double sdr_ = as<double>(sdr);
    unsigned numbsim_ = as<int>(numbsim);
    double t_ = as<double>(t);
    if(sbr_ <= 0.0)
        stop("'sbr' must be bigger than 0.0.");
    if(sbr_ < sdr_)
        stop("'sbr' must be greater than 'sdr'");
    if(numbsim_ < 1)
        stop("'numbsim' must be more than 1");
    if(sdr_ < 0.0)
        stop("'sdr' must be 0.0 or greater.");
    if(t_ <= 0.0)
        stop("'t' must be greater than 0.");
    return sim_bdsimple_species_tree(sbr_, sdr_, numbsim_, t_);
}
//' Simulates locus tree using constant rate birth-death-transfer process
//'
//' @description Given a species tree simulates a locus or gene family tree along
//'     the species tree.
//' @param species_tree species tree to simulate along
//' @param gbr gene birth rate
//' @param gdr gene death rate
//' @param lgtr gene transfer rate
//' @param num_loci number of locus trees to simulate
//' @param transfer_type The type of transfer input. Acceptable options: "cladewise" or "random"
//' @return List of objects of the tree class (as implemented in APE)
//' @details Given a species tree will perform a birth-death process coupled with transfer.
//' The simulation runs along the species tree speciating and going extinct in addition to locus tree birth and deaths.
//' Thus with parameters set to 0.0 a tree identical to the species tree is returned (it is relabel however).
//'
//' Transfers are implemented as a birth-death process.
//' One descendant lineage retains its species identity the other gains a new identity.
//' At present, two types of transfers are implemented: "random" an "cladewise".
//' The random transfer mode transfers one randomly chooses a contemporaneous lineage.
//' Cladewise transfers choose lineages based on relatedness with more closely related lineages being more likely.
//' @references
//' Rasmussen MD, Kellis M. Unified modeling of gene duplication, loss, and
//'     coalescence using a locus tree. Genome Res. 2012;22(4):755–765.
//'     doi:10.1101/gr.123901.111
//' @examples
//' # first simulate a species tree
//' mu <- 0.5 # death rate
//' lambda <- 2.0 # birth rate
//' numb_replicates <- 10
//' numb_extant_tips <- 4
//' # simulate trees under the GSA so first simulates a tree with
//' # numb_extant_tips * 100 tips counting each time we have a tree with 10 tips
//' # then randomly picks one of those trees
//'
//' sp_tree <- sim_sptree_bdp(sbr = lambda,
//'                 sdr = mu,
//'                 numbsim = numb_replicates,
//'                 n_tips = numb_extant_tips)
//'
//' gene_br <- 1.0
//' gene_dr <- 0.2
//' transfer_rate <- 0.2
//' sim_locustree_bdp(species_tree = sp_tree[[1]],
//'                   gbr = gene_br,
//'                   gdr = gene_dr,
//'                   lgtr = transfer_rate,
//'                   num_loci = 10)
// [[Rcpp::export]]
Rcpp::List sim_locustree_bdp(Rcpp::List species_tree,
                             SEXP gbr,
                             SEXP gdr,
                             SEXP lgtr,
                             SEXP num_loci,
                             Rcpp::String transfer_type = "random"){
    RNGScope scope;
    double gbr_ = as<double>(gbr);
    double gdr_ = as<double>(gdr);
    double lgtr_ = as<double>(lgtr);
    unsigned numLoci = as<int>(num_loci);
    std::string trans_type = transfer_type;

    if(gbr_ < 0.0)
        stop("'gbr' must be a positive number or 0.0");
    if(gbr_ < gdr_){
        if(gbr_ == gdr_)
            stop("'gbr' can only be equal to 'gdr' if both are 0.0");
        stop("'gbr' must be greater than 'gdr'");
    }
    if(numLoci < 1)
        stop("'num_loci' must be greater than 1");
    if(lgtr_ < 0.0)
        stop("'lgtr' must be a positive number or 0.0");
    if(gdr_ < 0.0)
        stop("'gdr' must be greater than or equal to 0.0");
    if(strcmp(species_tree.attr("class"), "phylo") != 0)
        stop("species_tree must be an object of class phylo'.");

    if(trans_type !=  "cladewise" && trans_type != "random")
        stop("the transfer_type must be set to 'cladewise' or 'random'");

    std::shared_ptr<SpeciesTree> specTree = std::shared_ptr<SpeciesTree>(new SpeciesTree(species_tree));
    return sim_locus_tree(specTree, gbr_, gdr_, lgtr_, numLoci, trans_type);
}
//' Simulates a cophylogenetic system using a paired birth-death process
//'
//' @details Simulates a cophylogenetic system using birth-death processes. The
//'     host tree is simulated following a constant rate birth-death process
//'     with an additional parameter - the cospeciation rate. This rate works as
//'     the speciation rate with the additional effect that if cospeciation
//'     occurs the symbiont tree also speciates. The symbiont tree is related to
//'     the host tree via an association matrix that describes which lineages
//'     are associated with which. The symbiont tree has an independent
//'     birth-death process with the addition of a host shift speciation rate
//'     that allows for the addition of more associated hosts upon symbiont
//'     speciation.
//' @param hbr host tree birth rate
//' @param hdr host tree death rate
//' @param sbr symbiont tree birth rate
//' @param sdr symbiont tree death rate
//' @param host_exp_rate host shift speciation rate
//' @param cosp_rate cospeciation rate
//' @param time_to_sim time units to simulate until
//' @param numbsim number of replicates
//' @return A list containing the `host_tree`, the `symbiont_tree`, the
//'     association matrix at present, and the history of events that have
//'     occurred.
//' @examples
//'
//' host_mu <- 0.5 # death rate
//' host_lambda <- 2.0 # birth rate
//' numb_replicates <- 10
//' time <- 1.0
//' symb_mu <- 0.2
//' symb_lambda <- 0.4
//' host_shift_rate <- 0.0
//' cosp_rate <- 2.0
//'
//' cophylo_pair <- sim_cophylo_bdp(hbr = host_lambda,
//'                            hdr = host_mu,
//'                            cosp_rate = cosp_rate,
//'                            host_exp_rate = host_shift_rate,
//'                            sdr = symb_mu,
//'                            sbr = symb_lambda,
//'                            numbsim = numb_replicates,
//'                            time_to_sim = time)
//'
// [[Rcpp::export]]
Rcpp::List sim_cophylo_bdp(SEXP hbr,
                           SEXP hdr,
                           SEXP sbr,
                           SEXP sdr,
                           SEXP host_exp_rate,
                           SEXP cosp_rate,
                           SEXP time_to_sim,
                           SEXP numbsim){
    double hbr_ = as<double>(hbr);
    double hdr_ = as<double>(hdr);
    double sbr_ = as<double>(sbr);
    double sdr_ = as<double>(sdr);
    double cosp_rate_ = as<double>(cosp_rate);
    double host_exp_rate_ = as<double>(host_exp_rate);
    double timeToSimTo_ = as<double>(time_to_sim);
    int numbsim_ = as<int>(numbsim);
    RNGScope scope;
    if(hbr_ < 0.0){
         stop("'hbr' must be positive or 0.0.");
    }
    if((hbr_ + cosp_rate_) < hdr_){
        stop("'hbr + cosp_rate' must be greater than 'hdr'.");
    }
    if(hdr_ < 0.0){
        stop("'hdr' must be a positive value or 0.0.");
    }
    if(host_exp_rate_ < 0.0){
        stop("'host_exp_rate' must be a positive value or 0.0.");
    }
    if(numbsim_ < 1){
        stop("'numbsim' must be larger than 1.");
    }
    if(cosp_rate_ < 0.0)
        stop("'cosp_rate' must be a positive value or 0.0.");
    if(timeToSimTo_ < 0.0)
        stop("'timeToSimTo' must be a positive value or 0.0.");

    return sim_host_symb_treepair(hbr_,
                                  hdr_,
                                  sbr_,
                                  sdr_,
                                  host_exp_rate_,
                                  cosp_rate_,
                                  timeToSimTo_,
                                  numbsim_);
}
//' Simulate multispecies coalescent on a species tree
//'
//' @description Simulates the multispecies coalescent on a species tree.
//' @param species_tree input species tree of class "phylo"
//' @param ne Effective population size
//' @param generation_time The number of time units per generation
//' @param num_sampled_individuals number of individuals sampled within each lineage
//' @param num_genes number of genes to simulate within each locus
//' @param mutation_rate The rate of mutation per generation
//' @param rescale Rescale the tree into coalescent units (otherwise assumes it is in those units)
//' @details
//' This a multispecies coalescent simulator with two usage options.
//' The function can rescale the given tree into coalescent units given the `mutation_rate`, `ne`, and the `generation_time`.
//' These result in a tree with coalescent times in units of expected number of mutations per site.
//' The generation_time parameter default is 1 time unit per generation if the units of the tree are in millions of years
//' The mutation_rate parameter is by default set to 1 mutations per site per generation (which is nonsensical).
//' Rescale is set to true by default.
//'
//' If rescale is set to false the tree is assumed to be in coalescent units and `ne` is used as the population
//' genetic parameter theta.
//' @return A list of coalescent trees
//' @seealso sim_locustree_bdp, sim_sptree_bdp, sim_sptree_bdp_time
//'
//' @examples
//' # first simulate a species tree
//' mu <- 0.5
//' lambda <- 1.0
//' nt <- 6
//' tr <- sim_sptree_bdp(sbr = lambda, sdr = mu, numbsim = 1, n_tips = nt)
//' # for a locus tree with 100 genes sampled per locus tree
//' gentrees <- sim_multispecies_coal(tr[[1]],
//'                                     ne = 10000,
//'                                     mutation_rate = 1e-9,
//'                                     generation_time = 1e-6,
//'                                     num_sampled_individuals = 1,
//'                                     num_genes = 100)
//'
//' @references
//' Bruce Rannala and Ziheng Yang (2003) Bayes Estimation of Species Divergence Times and Ancestral Population Sizes Using DNA Sequences From Multiple Loci Genetics August 1, 2003 vol. 164 no. 4 1645-1656
//' Mallo D, de Oliveira Martins L, Posada D (2015) SimPhy: Phylogenomic Simulation of Gene, Locus and Species Trees. Syst. Biol. doi: http://dx.doi.org/10.1093/sysbio/syv082
// [[Rcpp::export]]
Rcpp::List sim_multispecies_coal(SEXP species_tree,
                                 SEXP ne,
                                 SEXP num_sampled_individuals,
                                 SEXP num_genes,
                                 Rcpp::LogicalVector rescale = true,
                                 Rcpp::NumericVector mutation_rate = 1,
                                 Rcpp::NumericVector generation_time = 1){
    Rcpp::List species_tree_ = as<Rcpp::List>(species_tree);
    if(strcmp(species_tree_.attr("class"), "phylo") != 0)
        stop("species_tree must be an object of class phylo'.");
    auto specTree = std::shared_ptr<SpeciesTree>(new SpeciesTree(species_tree_));

    RNGScope scope;
    int num_sampled_individuals_ = as<int>(num_sampled_individuals);
    double ne_ = as<double>(ne);
    int num_genes_ = as<int>(num_genes);
    double mutation_rate_ = as<double>(mutation_rate);
    double generation_time_ = as<double>(generation_time);
    bool rescale_ = as<bool>(rescale);
    double u = std::exp(std::log(1) - std::log(generation_time_) + std::log(mutation_rate_)); //mut per site per gen x unit time per gen
    double theta = ne_;
    if(rescale_){
        theta = 4 * ne_ * u;
        specTree->scaleTree(theta);
    }
    if(mutation_rate_ <= 0.0)
        stop("'mutation_rate' must be greater than 0.0.");
    if(generation_time_ <= 0.0)
        stop("'generation_time' must be greater than 0.0.");
    if(ne_ <= 0.0)
        stop("'ne' must be greater than 0.0.");
    if(num_genes_ < 1)
        stop("'num_genes' must be greater than or equal to 1");
    if(num_sampled_individuals_ < 1)
        stop("'num_sampled_individuals' must be greater than or equal to 1");

    return sim_genetree_msc(specTree,
                            theta,
                            num_sampled_individuals_,
                            num_genes_);
}

