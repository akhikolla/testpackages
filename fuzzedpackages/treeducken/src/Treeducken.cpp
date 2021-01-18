#include <iostream>
#include <string>
#include "SpeciesTree.h"
#include "Simulator.h"
#include <sstream>
#include <RcppArmadillo.h>

using namespace Rcpp;


Rcpp::List bdsim_species_tree(double sbr,
                        double sdr,
                        int numbsim,
                        int n_tips,
                        int gsa_stop){

    List multiphy(numbsim);
    for(int i = 0; i < numbsim; i++){
        std::shared_ptr<Simulator> phySimulator = std::shared_ptr<Simulator>(new Simulator(n_tips,
                                                                                            sbr,
                                                                                            sdr,
                                                                                            1));
        phySimulator->setGSAStop(gsa_stop);
        phySimulator->simSpeciesTree();

        List phy = List::create(Named("edge") = phySimulator->getSpeciesEdges(),
                                Named("edge.length") = phySimulator->getSpeciesEdgeLengths(),
                                Named("Nnode") = phySimulator->getSpeciesNnodes(),
                                Named("tip.label") = phySimulator->getSpeciesTipNames(),
                                Named("root.edge") = phySimulator->getSpeciesTreeRootEdge());;
        phy.attr("class") = "phylo";
        multiphy[i] = phy;
    }


    multiphy.attr("class") = "multiPhylo";


    return multiphy;
}

Rcpp::List sim_bdsimple_species_tree(double sbr,
                                     double sdr,
                                     int numbsim,
                                     double timeToSimTo){
    List multiphy(numbsim);
    for(int i = 0; i < numbsim; i++){

        auto phySimulator = std::shared_ptr<Simulator>(new Simulator(1,
                                                                      sbr,
                                                                      sdr,
                                                                      1));
        phySimulator->setTimeToSim(timeToSimTo);
        phySimulator->simSpeciesTreeTime();

        List phy = List::create(Named("edge") = phySimulator->getSpeciesEdges(),
                                Named("edge.length") = phySimulator->getSpeciesEdgeLengths(),
                                Named("Nnode") = phySimulator->getSpeciesNnodes(),
                                Named("tip.label") = phySimulator->getSpeciesTipNames(),
                                Named("root.edge") = phySimulator->getSpeciesTreeRootEdge());
        phy.attr("class") = "phylo";
        multiphy[i] = phy;
    }

    multiphy.attr("class") = "multiPhylo";
    return multiphy;
}

Rcpp::List sim_locus_tree(std::shared_ptr<SpeciesTree> species_tree,
                          double gbr,
                          double gdr,
                          double lgtr,
                          int numbsim,
                          std::string trans_type){
    Rcpp::List multiphy;
    int ntax = species_tree->getNumExtant();
    double lambda = 0.0;
    double mu = 0.0;
    double rho = 0.0;
    unsigned numLociToSim = numbsim;
    for(int i = 0; i < numbsim; i++){

        auto phySimulator = std::shared_ptr<Simulator>(new Simulator( ntax,
                                                        lambda,
                                                        mu,
                                                        rho,
                                                        numLociToSim,
                                                        gbr,
                                                        gdr,
                                                        lgtr,
                                                        trans_type));
        phySimulator->setSpeciesTree(species_tree);

        phySimulator->simLocusTree();
        List phy = List::create(Named("edge") = phySimulator->getLocusEdges(),
                                Named("edge.length") = phySimulator->getLocusEdgeLengths(),
                                Named("Nnode") = phySimulator->getLocusNnodes(),
                                Named("tip.label") = phySimulator->getLocusTipNames(),
                                Named("root.edge") = phySimulator->getLocusTreeRootEdge(),
                                Named("node.label") = phySimulator->getLocusTreeNodeLabels());

        phy.attr("class") = "phylo";

        multiphy.push_back(phy);
    }
    multiphy.attr("class") = "multiPhylo";


    return multiphy;
}

Rcpp::List sim_locus_tree_gene_tree(std::shared_ptr<SpeciesTree> species_tree,
                                    double gbr,
                                    double gdr,
                                    double lgtr,
                                    int numLoci,
                                    double popsize,
                                    int samples_per_lineage,
                                    int numGenesPerLocus){
    Rcpp::List multiphy;
    int ntax = species_tree->getNumExtant();
    double lambda = 0.0;
    double mu = 0.0;
    double rho = 1.0;
    unsigned numLociToSim = numLoci;
    double genTime = 1.0;
    double ts = 1.0;
    bool sout = false;
    double og = 0.0;
    for(int i = 0; i < numLoci; i++){

        auto phySimulator = std::shared_ptr<Simulator>(new Simulator(ntax,
                                                        lambda,
                                                        mu,
                                                        rho,
                                                        numLociToSim,
                                                        gbr,
                                                        gdr,
                                                        lgtr,
                                                        samples_per_lineage,
                                                        popsize,
                                                        genTime,
                                                        numGenesPerLocus,
                                                        og,
                                                        ts,
                                                        sout));


        phySimulator->setSpeciesTree(species_tree);
        if(gbr + gdr + lgtr > 0.0){
            phySimulator->simLocusTree();
        }
        else{
            phySimulator->setLocusTree(std::shared_ptr<LocusTree>(new LocusTree(*species_tree, 
                                                                                 ntax, 
                                                                                 0.0, 
                                                                                 0.0, 
                                                                                 0.0)));
        }
        List phyGenesPerLoc(numGenesPerLocus);

        for(int j=0; j<numGenesPerLocus; j++){
            phySimulator->simGeneTree(j);

            List phyGene = List::create(Named("edge") = phySimulator->getGeneEdges(j),
                         _("edge.length") = phySimulator->getGeneEdgeLengths(j),
                         _("Nnode") = phySimulator->getGeneNnodes(j),
                         _("tip.label") = phySimulator->getGeneTipNames(j),
                         _("root.edge") = phySimulator->getGeneTreeRootEdge(j));
            phyGene.attr("class") = "phylo";
            phyGenesPerLoc[j] = phyGene;
        }

        // THE PROBLEM IS HERE I DON'T EVEN THINK I NEED THIS...
        List phyLoc = List::create(Named("edge") = phySimulator->getLocusEdges(),
                                   Named("edge.length") = phySimulator->getLocusEdgeLengths(),
                                   Named("Nnode") = phySimulator->getLocusNnodes(),
                                   Named("tip.label") = phySimulator->getLocusTipNames(),
                                   Named("root.edge") = phySimulator->getLocusTreeRootEdge());

        phyLoc.attr("class") = "phylo";
        List locusGeneSet = List::create(Named("container.tree") = phyLoc,
                                         Named("gene.trees") = phyGenesPerLoc);
        multiphy.push_back(locusGeneSet);
        
    }



    return multiphy;
}


Rcpp::List sim_genetree_msc(std::shared_ptr<SpeciesTree> species_tree,
                            double popsize,
                            int samples_per_lineage,
                            int numbsim){
    return sim_locus_tree_gene_tree(species_tree,
                             0.0,
                             0.0,
                             0.0,
                             1,
                             popsize,
                             samples_per_lineage,
                             numbsim);
    // this one is a wrapper for above function with locus tree parameters set to 0
}
