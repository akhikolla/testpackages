#include <Rcpp.h>
#include "r_functions.h"

using namespace Rcpp;

//RCPP_EXPOSED_CLASS_NODECL(MLNetwork)
RCPP_EXPOSED_CLASS(RMLNetwork)
RCPP_EXPOSED_CLASS(REvolutionModel)


void
show_multilayer_network(
    RMLNetwork *mnet
)
{
    Rcpp::Rcout << mnet->get_mlnet()->summary() << std::endl;
}

void
show_evolution_model(
    REvolutionModel *evm
)
{
    Rcpp::Rcout << evm->description() << std::endl;
}

RCPP_MODULE(multinet)
{

    class_<RMLNetwork>("RMLNetwork")
    .method("name", &RMLNetwork::name, "name of the multilayer network" )
    .method("show", &show_multilayer_network);

    class_<REvolutionModel>("REvolutionModel")
    .method("show", &show_evolution_model);

    /******************************/
    /* CREATION AND STORAGE       */
    /******************************/

    function("ml_empty", &emptyMultilayer, List::create( _["name"]=""), "Creates an empty multilayer network");

    // OTHER FUNCTIONS TO CREATE MLNETWORKS FROM THE LITERATURE (E.G., AUCS) ARE DEFINED IN functions.R

    // Generation
    function("evolution_pa_ml",
             &ba_evolution_model,
             List::create(
                 _["m0"],
                 _["m"]),
             "Creates a layer evolutionary model based on preferential attachment");

    function("evolution_er_ml",
             &er_evolution_model,
             List::create( _["n"]),
             "Creates a layer evolutionary model based on random edge creation, as in the ER model"
            );

    function("grow_ml",
             &growMultiplex,
             List::create(
                 _["num.actors"],
                 _["num.steps"],
                 _["models"],
                 _["pr.internal"],
                 _["pr.external"],
                 _["dependency"]),
             "Grows a multiplex network"
            );

    // IO

    function("read_ml", &readMultilayer, List::create( _["file"], _["name"]="unnamed", _["sep"]=',', _["aligned"]=false), "Reads a multilayer network from a file");

    function("write_ml", &writeMultilayer, List::create( _["n"], _["file"], _["format"]="multilayer", _["layers"]=CharacterVector(), _["sep"]=',', _["merge.actors"]=true, _["all.actors"]=false), "Writes a multilayer network to a file");



    /**************************************/
    /* INFORMATION ON MULTILAYER NETWORKS */
    /**************************************/

    function("layers_ml",
             &layers,
             List::create( _["n"]),
             "Returns the list of layers in the input multilayer network");

    function("actors_ml",
             &actors,
             List::create( _["n"], _["layers"]=CharacterVector()),
             "Returns the list of actors present in the input layers, or in the whole multilayer network if no layers are specified");

    function("vertices_ml",
             &vertices,
             List::create( _["n"], _["layers"]=CharacterVector()),
             "Returns the list of vertices in the input layers, or in the whole multilayer network if no layers are specified");

    function("edges_ml",
             &edges,
             List::create( _["n"], _["layers1"]=CharacterVector(), _["layers2"]=CharacterVector()),
             "Returns the list of edges among vertices in the input layers (if only one set of layers is specified), or from the first set of input layers to the second set of input layers, or in the whole multilayer network if no layers are specified");

    function("edges_idx_ml",
             &edges_idx,
             List::create( _["n"]),
             "Returns the list of edges, where vertex ids are used instead of vertex names");

    function("num_layers_ml",
             &numLayers,
             List::create( _["n"]),
             "Returns the number of layers in the input mlnetwork");
    function("num_actors_ml",
             &numActors,
             List::create( _["n"], _["layers"]=CharacterVector()),
             "Returns the number of actors in the set of input layers, or in the whole mlnetwork if no layers are specified");
    function("num_vertices_ml",
             &numNodes,
             List::create( _["n"], _["layers"]=CharacterVector()),
             "Returns the number of vertices in the set of input layers, or in the whole mlnetwork if no layers are specified");
    function("num_edges_ml",
             &numEdges,
             List::create( _["n"], _["layers1"]=CharacterVector(), _["layers2"]=CharacterVector()),
             "Returns the number of edges in the set of input layers, or in the whole mlnetwork if no layers are specified");

    function("is_directed_ml",
             &isDirected,
             List::create( _["n"], _["layers1"]=CharacterVector(), _["layers2"]=CharacterVector()),
             "Returns a logical vector indicating for each pair of layers if it is directed or not");

    // NAVIGATION

    function("neighbors_ml", &actor_neighbors, List::create( _["n"], _["actor"], _["layers"]=CharacterVector(), _["mode"] = "all"), "Returns the neighbors of a global identity on the set of input layers");
    function("xneighbors_ml", &actor_xneighbors, List::create( _["n"], _["actor"], _["layers"]=CharacterVector(), _["mode"] = "all"), "Returns the exclusive neighbors of a global identity on the set of input layers");


    // NETWORK MANIPULATION

    function("add_layers_ml", &addLayers, List::create( _["n"], _["layers"], _["directed"]=false), "Adds one or more layers to a multilayer network");
    function("add_actors_ml_v3", &addActors, List::create( _["n"], _["actors"]), "Adds one or more actors to a multilayer network");
    function("add_vertices_ml", &addNodes, List::create( _["n"], _["vertices"]), "Adds one or more vertices to a layer of a multilayer network");
    function("add_edges_ml", &addEdges, List::create( _["n"], _["edges"]), "Adds one or more edges to a multilayer network - each edge is a quadruple [actor,layer,actor,layer]");

    function("set_directed_ml", &setDirected, List::create( _["n"], _["directionalities"]), "Set the directionality of one or more pairs of layers");

    function("delete_layers_ml", &deleteLayers, List::create( _["n"], _["layers"]), "Deletes one or more layers from a multilayer network");
    function("delete_actors_ml", &deleteActors, List::create( _["n"], _["actors"]), "Deletes one or more actors from a multilayer network");


    function("delete_vertices_ml", &deleteNodes, List::create( _["n"], _["vertices"]), "Deletes one or more vertices from a layer of a multilayer network");


    function("delete_edges_ml", &deleteEdges, List::create( _["n"], _["edges"]), "Deletes one or more edges from a multilayer network - each edge is a quadruple [actor,layer,actor,layer]");


    // ATTRIBUTE HANDLING

    function("add_attributes_ml", &newAttributes, List::create( _["n"], _["attributes"], _["type"]="string", _["target"]="actor", _["layer"]="", _["layer1"]="", _["layer2"]=""), "Creates a new attribute so that values can be associated to actors, layers, vertices or edges");
    function("attributes_ml", &getAttributes, List::create( _["n"], _["target"]="actor"), "Returns the list of attributes defined for the input multilayer network");

    function("get_values_ml", &getValues, List::create( _["n"], _["attribute"], _["actors"]=CharacterVector(), _["vertices"]=CharacterMatrix(0,0), _["edges"]=CharacterMatrix(0,0)), "Returns the value of an attribute on the specified actors, layers, vertices or edges");
    function("set_values_ml", &setValues, List::create( _["n"], _["attribute"], _["actors"]=CharacterVector(), _["vertices"]=CharacterMatrix(0,0), _["edges"]=CharacterMatrix(0,0), _["values"]), "Sets the value of an attribute for the specified actors/vertexes/edges");


    // TRANSFORMATION

    function("flatten_ml", &flatten, List::create( _["n"], _["new.layer"]="flattening", _["layers"]=CharacterVector(), _["method"] = "weighted", _["force.directed"] = false, _["all.actors"] = false), "Adds a new layer with the actors in the input layers and an edge between A and B if they are connected in any of the merged layers");
    
     function("project_ml", &project, List::create( _["n"], _["new.layer"]="projection", _["layer1"], _["layer2"], _["method"] = "clique"), "Adds a new layer with the actors in layer 1, and edges between actors A and B if they are connected to a common object in layer 2");

    // MEASURES
    
    function("degree_ml", &degree_ml, List::create( _["n"], _["actors"]=CharacterVector(), _["layers"]=CharacterVector(), _["mode"] = "all"), "Returns the degree of each actor");

    function("degree_deviation_ml", &degree_deviation_ml, List::create( _["n"], _["actors"]=CharacterVector(), _["layers"]=CharacterVector(), _["mode"] = "all"), "Returns the standard deviation of the degree of each actor on the specified layers");
    /*
    function("occupation_ml", &occupation_ml, List::create( _["n"], _["transitions"], _["teleportation"]=.2, _["steps"]=0), "Returns the occupation centrality value of each actor");
     */
    function("neighborhood_ml", &neighborhood_ml, List::create( _["n"], _["actors"]=CharacterVector(), _["layers"]=CharacterVector(), _["mode"] = "all"), "Returns the neighborhood of each actor");
    function("xneighborhood_ml", &xneighborhood_ml, List::create( _["n"], _["actors"]=CharacterVector(), _["layers"]=CharacterVector(), _["mode"] = "all"), "Returns the exclusive neighborhood of each actor");

    function("connective_redundancy_ml", &connective_redundancy_ml, List::create( _["n"], _["actors"]=CharacterVector(), _["layers"]=CharacterVector(), _["mode"] = "all"), "Returns the connective redundancy of each actor");
    function("relevance_ml", &relevance_ml, List::create( _["n"], _["actors"]=CharacterVector(), _["layers"]=CharacterVector(), _["mode"] = "all"), "Returns the layer relevance of each actor");
    function("xrelevance_ml", &xrelevance_ml, List::create( _["n"], _["actors"]=CharacterVector(), _["layers"]=CharacterVector(), _["mode"] = "all"), "Returns the exclusive layer relevance of each actor");

    function("layer_summary_ml", &summary_ml, List::create( _["n"], _["layer"], _["method"] = "entropy.degree", _["mode"] = "all"), "Computes a summary of the input layer");

    function("layer_comparison_ml", &comparison_ml, List::create( _["n"], _["layers"]=CharacterVector(), _["method"] = "jaccard.edges", _["mode"] = "all", _["K"] = 0), "Computes the similarity between the input layers");


    function("distance_ml", &distance_ml, List::create( _["n"], _["from"], _["to"]=CharacterVector(), _["method"] = "multiplex"), "Computes the distance between two actors");


    // CLUSTERING
    function("clique_percolation_ml",
             &cliquepercolation_ml,
             List::create(
                          _["n"],
                          _["k"]=3,
                          _["m"]=1
            ), "Extension of the clique percolation method");
    

    function("glouvain_ml",
             &glouvain_ml,
             List::create(_["n"],_["gamma"]=1,_["omega"]=1,_["limit"]=0),
             "Extension of the louvain method");
    
    function("glouvain2_ml",
             &glouvain2_ml,
             List::create(_["n"], _["omega"]=1),
             "Extension of the louvain method");
    
    function("abacus_ml", &abacus_ml,List::create(_["n"],_["min.actors"]=3,_["min.layers"]=1),
            "Community extraction based on frequent itemset mining");
    
    function("infomap_ml",
             &infomap_ml, List::create(_["n"],_["overlapping"]=false,_["directed"]=false,_["self.links"]=true),
             "Community extraction based on the flow equation");

    
    function("modularity_ml",
             &modularity_ml,
             List::create(_["n"], _["comm.struct"],_["gamma"]=1,_["omega"]=1),
             "Generalized modularity");
    
    /*
    function("lart_ml", &lart_ml, List::create( _["n"], _["t"]=-1, _["eps"]=1, _["gamma"]=1), "Community extraction based on locally adaptive random walks");


    // FOR VISUALIZATION
     */
    function("layout_multiforce_ml", &multiforce_ml, List::create( _["n"], _["w_in"]=1, _["w_inter"]=1, _["gravity"]=0, _["iterations"]=100), "Multiforce method: computes vertex coordinates");

    function("layout_circular_ml", &circular_ml, List::create( _["n"] ), "Circular method: computes vertex coordinates arranging actors on a circle");

    // plotting function defined in functions.R


    function("get_community_list_ml", &to_list, List::create( _["comm.struct"], _["n"]), "Converts a community structure (data frame) into a list of communities, layer by layer");

    /*

     //function("sir_ml", &sir_ml, List::create( _["n"], _["beta"], _["tau"], _["num_iterations"] = 1000), "Executes a SIR spreading process, returning the number of vertices in each status at each iteration");
     */

}


