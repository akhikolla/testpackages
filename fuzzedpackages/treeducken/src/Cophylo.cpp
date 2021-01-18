#include <RcppArmadillo.h>
#include "SymbiontTree.h"
#include "Simulator.h"
using namespace Rcpp;


Rcpp::List sim_host_symb_treepair(double hostbr,
                                  double hostdr,
                                  double symbbr,
                                  double symbdr,
                                  double switchRate,
                                  double cospeciationRate,
                                  double timeToSimTo,
                                  int numbsim){

    double rho = 1.0;
    Rcpp::List multiphy;
    Rcpp::List hostSymbPair;
    for(int i = 0; i < numbsim; i++){
        auto phySimulator = std::shared_ptr<Simulator>(new Simulator( timeToSimTo,
                                                 hostbr,
                                                 hostdr,
                                                 symbbr,
                                                 symbdr,
                                                 switchRate,
                                                 cospeciationRate,
                                                 rho,
                                                 1));

        phySimulator->simHostSymbSpeciesTreePair();





        List phyHost = List::create(Named("edge") = phySimulator->getSpeciesEdges(),
                                    Named("edge.length") = phySimulator->getSpeciesEdgeLengths(),
                                    Named("Nnode") = phySimulator->getSpeciesNnodes(),
                                    Named("tip.label") = phySimulator->getSpeciesTipNames(),
                                    Named("root.edge") = phySimulator->getSpeciesTreeRootEdge());
        phyHost.attr("class") = "phylo";


        List phySymb = List::create(Named("edge") = phySimulator->getSymbiontEdges(),
                                Named("edge.length") = phySimulator->getSymbiontEdgeLengths(),
                                Named("Nnode") = phySimulator->getSymbiontNnodes(),
                                Named("tip.label") = phySimulator->getSymbiontTipNames(),
                                Named("root.edge") = phySimulator->getSymbiontTreeRootEdge());
        phySymb.attr("class") = "phylo";

        hostSymbPair = List::create(Named("host_tree") = phyHost,
                                    Named("symb_tree") = phySymb,
                                    Named("association_mat") = phySimulator->getAssociationMatrix(),
                                    Named("event_history") = phySimulator->createEventDF());
        hostSymbPair.attr("class") = "cophy";
        multiphy.push_back(hostSymbPair);
    }


    multiphy.attr("class") = "multiCophy";

    return multiphy;
}
