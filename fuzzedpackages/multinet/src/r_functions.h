/*
 * r_functions.h
 *
 * Created on: Jun 19, 2014
 * Author: matteomagnani
 * Version: 0.0.1
 */

#ifndef _R_FUNCTIONS_H_
#define _R_FUNCTIONS_H_

#include <Rcpp.h>
#include "networks/MultilayerNetwork.hpp"
#include "generation/EvolutionModel.hpp"
#include <unordered_set>
#include <vector>
#include <memory>

using namespace Rcpp;

class RMLNetwork
{
  private:
    std::shared_ptr<uu::net::MultilayerNetwork> ptr;

  public:

    std::string
    name(
    ) const
    {
        return ptr->name;
    }

    RMLNetwork(std::shared_ptr<uu::net::MultilayerNetwork> ptr) : ptr(ptr)
    {
        // @todo check not null?
    }

    uu::net::MultilayerNetwork*
    get_mlnet() const
    {
        return ptr.get();
    }

};

class REvolutionModel
{
  private:
    std::shared_ptr<uu::net::EvolutionModel<uu::net::MultilayerNetwork>> ptr;
    std::string desc;

  public:
    std::string
    description(
    ) const
    {
        return desc;
    }

    REvolutionModel(
        std::shared_ptr<uu::net::EvolutionModel<uu::net::MultilayerNetwork>> ptr,
        const std::string description
    ) : ptr(ptr), desc(description) {}

    uu::net::EvolutionModel<uu::net::MultilayerNetwork>*
    get_model() const
    {
        return ptr.get();
    }
};

// CREATION AND STORAGE


RMLNetwork
emptyMultilayer(
    const std::string& name
);


RMLNetwork
readMultilayer(
               const std::string& input_file,
               const std::string& name, char sep, bool vertex_aligned);

void
writeMultilayer(
                const RMLNetwork& mnet,
                const std::string& output_file,
                const std::string& format,
                const CharacterVector& layer_names, char sep, bool merge_actors, bool all_actors);

REvolutionModel
ba_evolution_model(
    size_t m0,
    size_t m
);

REvolutionModel
er_evolution_model(
    size_t n
);

RMLNetwork
growMultiplex(
    size_t num_actors,
    long num_of_steps,
    const GenericVector& evolution_model,
    const NumericVector& pr_internal_event,
    const NumericVector& pr_external_event,
    const NumericMatrix& dependency
);


// INFORMATION ON NETWORKS

CharacterVector
layers(
    const RMLNetwork& mnet
);

CharacterVector
actors(
    const RMLNetwork& mnet,
    const CharacterVector& layer_names
);

DataFrame
vertices(
    const RMLNetwork& mnet,
    const CharacterVector& layer_names
);

DataFrame
edges(
    const RMLNetwork& mnet,
    const CharacterVector& layer_names1,
    const CharacterVector& layer_names2
);

DataFrame
edges_idx(
    const RMLNetwork& rmnet
);

size_t
numLayers(
    const RMLNetwork& mnet
);

size_t
numActors(
    const RMLNetwork& mnet,
    const CharacterVector& layers
);

size_t
numNodes(
    const RMLNetwork& mnet,
    const CharacterVector& layers
);
size_t
numEdges(
    const RMLNetwork& mnet,
    const CharacterVector& layer_names1,
    const CharacterVector& layer_names2
);

DataFrame
isDirected(
    const RMLNetwork& mnet,
    const CharacterVector& layer_names1,
    const CharacterVector& layer_names2
);

std::unordered_set<std::string>
actor_neighbors(
    const RMLNetwork& rmnet,
    const std::string& actor_name,
    const CharacterVector& layer_names,
    const std::string& mode_name
);

std::unordered_set<std::string>
actor_xneighbors(
    const RMLNetwork& rmnet,
    const std::string& actor_name,
    const CharacterVector& layer_names,
    const std::string& mode_name
);


// NETWORK MANIPULATION

void
addLayers(
    RMLNetwork& rmnet,
    const CharacterVector& layer_names,
    const LogicalVector& directed
);

void
addActors(
    RMLNetwork& rmnet,
    const CharacterVector& actor_names
);

void
addNodes(
    RMLNetwork& rmnet,
    const DataFrame& vertices
);

void
addEdges(
    RMLNetwork& rmnet,
    const DataFrame& edges);

void
setDirected(
    const RMLNetwork&,
    const DataFrame& directionalities
);

void
deleteLayers(
    RMLNetwork& rmnet,
    const CharacterVector& layer_names
);

void
deleteActors(
    RMLNetwork& rmnet,
    const CharacterVector& actor_names
);

void
deleteNodes(
    RMLNetwork& rmnet,
    const DataFrame& vertices
);

void
deleteEdges(
    RMLNetwork& rmnet,
    const DataFrame& edges
);



void
newAttributes(
    RMLNetwork& rmnet,
    const CharacterVector& attribute_names,
    const std::string& type,
    const std::string& target,
    const std::string& layer_name,
    const std::string& layer_name1,
    const std::string& layer_name2
);

DataFrame
getAttributes(
    const RMLNetwork&,
    const std::string& target
);

DataFrame
getValues(
    RMLNetwork& rmnet,
    const std::string& attribute_name,
    const CharacterVector& actor_names,
    const DataFrame& vertex_matrix,
    const DataFrame& edge_matrix
);


void
setValues(
    RMLNetwork& rmnet,
    const std::string& attribute_name,
    const CharacterVector& actor_names,
    const DataFrame& vertex_matrix,
    const DataFrame& edge_matrix,
    const GenericVector& values
);


// TRANSFORMATION

void
flatten(
    RMLNetwork& rmnet,
    const std::string& new_layer,
    const CharacterVector& layer_names,
    const std::string& method,
    bool force_directed,
    bool all_actors
);


void
project(
    RMLNetwork& rmnet,
    const std::string& new_layer,
    const std::string& layer1,
    const std::string& layer2,
    const std::string& method);

// MEASURES


NumericVector
degree_ml(
    const RMLNetwork&,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type
);


NumericVector
degree_deviation_ml(
    const RMLNetwork&,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type
);


/*
NumericVector
occupation_ml(
    const RMLNetwork&, const NumericMatrix& transitions, double teleportation, long steps);

 */

NumericVector
neighborhood_ml(
    const RMLNetwork& mnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type
);

NumericVector
xneighborhood_ml(
    const RMLNetwork& mnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type
);


NumericVector
connective_redundancy_ml(
    const RMLNetwork& mnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type
);

NumericVector
relevance_ml(
    const RMLNetwork& mnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type
);

NumericVector
xrelevance_ml(
    const RMLNetwork& mnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type
);



double
summary_ml(
    const RMLNetwork&,
    const std::string& layer,
    const std::string& method,
    const std::string& type
);

DataFrame
comparison_ml(
    const RMLNetwork&,
    const CharacterVector& layer_names,
    const std::string& method,
    const std::string& type,
    int K
);



DataFrame
distance_ml(const RMLNetwork& mnet,
            const std::string& from,
            const CharacterVector& to,
            const std::string& method);

// CLUSTERING

DataFrame
cliquepercolation_ml(
    const RMLNetwork& rmnet,
    int k,
    int m
);


DataFrame
infomap_ml(const RMLNetwork& mnet,
           bool overlapping,
           bool directed,
           bool include_self_links
          );


DataFrame
glouvain_ml(
    const RMLNetwork&,
    double gamma,
    double omega,
    int limit
);

DataFrame
glouvain2_ml(
    const RMLNetwork&,
    double omega
);

DataFrame
abacus_ml(
    const RMLNetwork&,
    int min_actors,
    int min_layers
);

double
modularity_ml(
              const RMLNetwork& rmnet,
              const DataFrame& com, double gamma,
              double omega
              );


List
to_list(
        const DataFrame& cs,
        const RMLNetwork& mnet
        );

// Layout

DataFrame
multiforce_ml(
    const RMLNetwork& mnet,
    const NumericVector& w_in,
    const NumericVector& w_out,
    const NumericVector& gravity,
    int iterations
);

DataFrame
circular_ml(
    const RMLNetwork& mnet
);


/*

// SPREADING
// NumericMatrix sir_ml(const RMLNetwork& mnet, double beta, int tau, long num_iterations);

 */

#endif
