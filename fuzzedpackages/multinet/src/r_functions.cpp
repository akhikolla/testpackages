/*
 * r_functions.cpp
 *
 * Created on: Jun 19, 2014
 * Author: matteomagnani
 * Version: 0.0.1
 */

#include <sstream>
#include "r_functions.h"
#include "rcpp_utils.h"

#include "operations/union.hpp"
#include "operations/project.hpp"
#include "community/glouvain.hpp"
#include "community/glouvain2.hpp"
#include "community/abacus.hpp"
#include "community/infomap.hpp"
#include "community/ml-cpm.hpp"
#include "community/modularity.hpp"
#include "io/read_multilayer_network.hpp"
#include "io/write_multilayer_network.hpp"
#include "measures/degree_ml.hpp"
#include "measures/neighborhood.hpp"
#include "measures/relevance.hpp"
#include "measures/redundancy.hpp"
#include "measures/layer.hpp"
#include "measures/distance.hpp"
#include "generation/evolve.hpp"
#include "generation/PAEvolutionModel.hpp"
#include "generation/EREvolutionModel.hpp"
#include "core/datastructures/propertymatrix/summarization.hpp"
#include "layout/multiforce.hpp"
#include "layout/circular.hpp"

using namespace Rcpp;

RCPP_EXPOSED_CLASS(RMLNetwork)
RCPP_EXPOSED_CLASS(REvolutionModel)

using M = uu::net::MultilayerNetwork;
using G = uu::net::Network;
// CREATION AND STORAGE

RMLNetwork
emptyMultilayer(
    const std::string& name
)
{
    return RMLNetwork(std::make_shared<uu::net::MultilayerNetwork>(name));
}

/*

void renameMultilayer(
    RMLNetwork& rmnet,
    const std::string& new_name) {
    rmnet.get_mlnet()->name = new_name;
}

 */

RMLNetwork
readMultilayer(const std::string& input_file,
               const std::string& name, char sep, bool vertex_aligned)
{
    return RMLNetwork(uu::net::read_attributed_homogeneous_multilayer_network(input_file,name,sep,vertex_aligned));
}


void
writeMultilayer(
    const RMLNetwork& rmnet,
    const std::string& output_file,
    const std::string& format,
    const CharacterVector& layer_names,
    char sep,
    bool merge_actors,
    bool all_actors
)
{
    auto mnet = rmnet.get_mlnet();
    auto layers = resolve_layers_unordered(mnet,layer_names);

    if (format=="multilayer")
    {
        write_attributed_homogeneous_multilayer_network(mnet,layers.begin(),layers.end(),output_file,sep);
    }

    else if (format=="graphml")
    {
        if (!merge_actors && all_actors)
        {
            Rcout << "option all.actors not used when merge.actors=FALSE" << std::endl;
        }

        write_graphml(mnet,layers.begin(),layers.end(),output_file,merge_actors,all_actors);
    }

    else
    {
        stop("unexpected value: format " + format);
    }
}


REvolutionModel
ba_evolution_model(
    size_t m0,
    size_t m
)
{
    auto pa = std::make_shared<uu::net::PAEvolutionModel<uu::net::MultilayerNetwork>>(m0, m);

    return REvolutionModel(pa,"Preferential attachment evolution model (" + std::to_string(m0) + "," + std::to_string(m) + ")");
}


REvolutionModel
er_evolution_model(
    size_t n
)
{
    auto er = std::make_shared<uu::net::EREvolutionModel<uu::net::MultilayerNetwork>>(n);

    return REvolutionModel(er, "Uniform evolution model (" + std::to_string(n) + ")");
}


RMLNetwork
growMultiplex(
    size_t num_actors,
    long num_of_steps,
    const GenericVector& evolution_model,
    const NumericVector& pr_internal_event,
    const NumericVector& pr_external_event,
    const NumericMatrix& dependency
)
{

    if (num_actors<=0)
    {
        stop("The number of actors must be positive");
    }

    if (num_of_steps<=0)
    {
        stop("The number of steps must be positive");
    }

    size_t num_layers = evolution_model.size();

    if (dependency.nrow()!=num_layers || dependency.ncol()!=num_layers ||
            pr_internal_event.size()!=num_layers || pr_external_event.size()!=num_layers)
    {
        stop("The number of evolution models, evolution probabilities and the number of rows/columns of the dependency matrix must be the same");
    }

    std::vector<double> pr_int(pr_internal_event.size());

    for (size_t i=0; i<pr_internal_event.size(); i++)
    {
        pr_int[i] = pr_internal_event.at(i);
    }

    std::vector<double> pr_ext(pr_external_event.size());

    for (size_t i=0; i<pr_external_event.size(); i++)
    {
        pr_ext[i] = pr_external_event.at(i);
    }

    std::vector<std::vector<double> > dep(dependency.nrow());

    for (size_t i=0; i<dependency.nrow(); i++)
    {
        std::vector<double> row(dependency.ncol());

        for (size_t j=0; j<dependency.ncol(); j++)
        {
            row[j] = dependency(i,j);
        }

        dep[i] = row;
    }

    std::vector<uu::net::EvolutionModel<uu::net::MultilayerNetwork>*> models(evolution_model.size());

    for (size_t i=0; i<models.size(); i++)
    {
        models[i] = (as<REvolutionModel>(evolution_model[i])).get_model();

    }

    auto res = std::make_shared<uu::net::MultilayerNetwork>("synth");

    for (size_t a=0; a<num_actors; a++)
    {
        res->actors()->add("a"+std::to_string(a));
    }

    std::vector<std::string> layer_names;

    for (size_t l=0; l<num_layers; l++)
    {
        std::string layer_name = "l"+std::to_string(l);
        auto layer = std::make_unique<uu::net::Network>(layer_name, uu::net::EdgeDir::UNDIRECTED, true);
        res->layers()->add(std::move(layer));
        layer_names.push_back(layer_name);
    }

    uu::net::evolve(res.get(), layer_names, pr_int, pr_ext, dep, models, num_of_steps);
    return RMLNetwork(res);
}


// INFORMATION ON NETWORKS

CharacterVector
layers(
    const RMLNetwork& rmnet
)
{
    auto mnet = rmnet.get_mlnet();
    CharacterVector res(0);

    for (auto layer: *mnet->layers())
    {
        res.push_back(layer->name);
    }

    return res;
}

CharacterVector
actors(
    const RMLNetwork& rmnet,
    const CharacterVector& layer_names
)
{
    CharacterVector actors(0);
    auto mnet = rmnet.get_mlnet();

    if (layer_names.size()==0)
    {
        for (auto actor: *mnet->actors())
        {
            actors.push_back(actor->name);
        }
    }

    else
    {
        auto layers = resolve_layers(mnet,layer_names);

        for (auto layer: layers)
        {
            for (auto vertex: *layer->vertices())
            {
                actors.push_back(vertex->name);
            }
        }
    }

    //for (auto actor: actors)
    //  res.push_back(actor->name);
    return actors;
}

DataFrame
vertices(
    const RMLNetwork& rmnet,
    const CharacterVector& layer_names
)
{
    auto mnet = rmnet.get_mlnet();
    auto layers = resolve_layers_unordered(mnet,layer_names);
    CharacterVector actor, layer;

    for (auto l: *mnet->layers())
    {

        if (layers.count(l)==0)
        {
            continue;
        }

        for (auto vertex: *l->vertices())
        {
            actor.push_back(vertex->name);
            layer.push_back(l->name);
        }
    }

    return DataFrame::create(_["actor"] = actor, _["layer"] = layer);
}

DataFrame
edges_idx(
    const RMLNetwork& rmnet
)
{
    auto mnet = rmnet.get_mlnet();
    NumericVector from, to, directed;

    // stores at which index vertices start in a layer
    std::unordered_map<const G*, size_t> offset;
    size_t num_vertices = 0;
    for (auto layer: *mnet->layers())
    {
        offset[layer] = num_vertices;
        num_vertices += layer->vertices()->size();
    }
    
    // intralayer
    
    for (auto l: *mnet->layers())
    {
        auto vertices = l->vertices();

        for (auto edge: *l->edges())
        {
            from.push_back(vertices->index_of(edge->v1)+offset[l]+1);
            to.push_back(vertices->index_of(edge->v2)+offset[l]+1);
            directed.push_back((edge->dir==uu::net::EdgeDir::DIRECTED)?1:0);
        }
    }

    // interlayer
    for (auto edge: *mnet->interlayer_edges())
    {
        from.push_back(edge->l1->vertices()->index_of(edge->v1)+offset[edge->l1]+1);
        to.push_back(edge->l2->vertices()->index_of(edge->v2)+offset[edge->l2]+1);
        directed.push_back((edge->dir==uu::net::EdgeDir::DIRECTED)?1:0);
    }


    return DataFrame::create(_["from"] = from, _["to"] = to, _["dir"] = directed );
}

DataFrame
edges(
    const RMLNetwork& rmnet,
    const CharacterVector& layer_names1,
    const CharacterVector& layer_names2
)
{
    auto mnet = rmnet.get_mlnet();
    std::vector<uu::net::Network*> layers1 = resolve_layers(mnet,layer_names1);
    std::vector<uu::net::Network*> layers2;

    if (layer_names2.size()==0)
    {
        layers2 = layers1;
    }

    else
    {
        layers2 = resolve_layers(mnet,layer_names2);
    }

    CharacterVector from_a, from_l, to_a, to_l;
    NumericVector directed;

    for (auto layer1: layers1)
    {
        for (auto layer2: layers2)
        {
            if (layer2<layer1)
            {
                continue;
            }

            else if (layer1==layer2)
            {

                for (auto edge: *layer1->edges())
                {
                    from_a.push_back(edge->v1->name);
                    from_l.push_back(layer1->name);
                    to_a.push_back(edge->v2->name);
                    to_l.push_back(layer1->name);
                    directed.push_back((edge->dir==uu::net::EdgeDir::DIRECTED)?1:0);
                }
            }

            else
            {
                for (auto edge: *mnet->interlayer_edges()->get(layer1,layer2))
                {
                    from_a.push_back(edge->v1->name);
                    from_l.push_back(edge->l1->name);
                    to_a.push_back(edge->v2->name);
                    to_l.push_back(edge->l2->name);
                    directed.push_back((edge->dir==uu::net::EdgeDir::DIRECTED)?1:0);
                }
            }

        }
    }

    return DataFrame::create(_["from_actor"] = from_a, _["from_layer"] = from_l, _["to_actor"] = to_a, _["to_layer"] = to_l, _["dir"] = directed );
}

size_t
numLayers(
    const RMLNetwork& rmnet
)
{
    auto mnet = rmnet.get_mlnet();
    return mnet->layers()->size();
}

size_t
numActors(
    const RMLNetwork& rmnet,
    const CharacterVector& layer_names
)
{
    auto mnet = rmnet.get_mlnet();

    if (layer_names.size()==0)
    {
        return mnet->actors()->size();
    }

    std::vector<uu::net::Network*> layers = resolve_layers(mnet,layer_names);
    std::unordered_set<const uu::net::Vertex*> actors;

    for (auto layer: layers)
    {
        for (auto vertex: *layer->vertices())
        {
            actors.insert(vertex);
        }
    }

    return actors.size();
}

size_t
numNodes(
    const RMLNetwork& rmnet,
    const CharacterVector& layer_names
)
{
    auto mnet = rmnet.get_mlnet();
    std::vector<uu::net::Network*> layers = resolve_layers(mnet,layer_names);
    size_t num_vertices = 0;

    for (auto layer: layers)
    {
        num_vertices += layer->vertices()->size();
    }

    return num_vertices;
}

size_t
numEdges(
    const RMLNetwork& rmnet,
    const CharacterVector& layer_names1,
    const CharacterVector& layer_names2
)
{
    auto mnet = rmnet.get_mlnet();
    std::unordered_set<const uu::net::Network*> layers1 = resolve_const_layers_unordered(mnet,layer_names1);
    std::unordered_set<const uu::net::Network*> layers2;

    if (layer_names2.size()==0)
    {
        layers2 = layers1;
    }

    else
    {
        layers2 = resolve_const_layers_unordered(mnet,layer_names2);
    }

    size_t num_edges = 0;



    for (auto layer1: layers1)
    {
        for (auto layer2: layers2)
        {
            if (layer2<layer1)
            {
                continue;
            }

            else if (layer1==layer2)
            {
                num_edges += layer1->edges()->size();
            }

            else
            {
                num_edges += mnet->interlayer_edges()->get(layer1,layer2)->size();
            }

        }
    }

    return num_edges;
}

DataFrame
isDirected(
    const RMLNetwork& rmnet,
    const CharacterVector& layer_names1,
    const CharacterVector& layer_names2)
{
    auto mnet = rmnet.get_mlnet();
    std::vector<uu::net::Network*> layers1 = resolve_layers(mnet,layer_names1);
    std::vector<uu::net::Network*> layers2;

    if (layer_names2.size()==0)
    {
        layers2 = layers1;
    }

    else
    {
        layers2 = resolve_layers(mnet,layer_names2);
    }

    CharacterVector l1, l2;
    NumericVector directed;

    for (auto layer1: layers1)
    {
        for (auto layer2: layers2)
        {
            l1.push_back(layer1->name);
            l2.push_back(layer2->name);
            
            if (layer1==layer2)
            {
                directed.push_back(layer1->is_directed()?1:0);
            }
            else
            {
                directed.push_back((mnet->interlayer_edges()->is_directed(layer1,layer2))?1:0);
            }
        }
    }

    return DataFrame::create(_["layer1"] = l1, _["layer2"] = l2, _["dir"] = directed );
}

std::unordered_set<std::string>
actor_neighbors(
    const RMLNetwork& rmnet,
    const std::string& actor_name,
    const CharacterVector& layer_names,
    const std::string& mode_name
)
{
    std::unordered_set<std::string> res_neighbors;
    auto mnet = rmnet.get_mlnet();
    auto actor = mnet->actors()->get(actor_name);

    if (!actor)
    {
        stop("actor " + actor_name + " not found");
    }

    auto layers = resolve_layers_unordered(mnet, layer_names);
    auto mode = resolve_mode(mode_name);
    auto actors = uu::net::neighbors(layers.begin(), layers.end(), actor, mode);

    for (auto neigh: actors)
    {
        res_neighbors.insert(neigh->name);
    }

    return res_neighbors;
}

std::unordered_set<std::string>
actor_xneighbors(
    const RMLNetwork& rmnet,
    const std::string& actor_name,
    const CharacterVector& layer_names,
    const std::string& mode_name
)
{
    std::unordered_set<std::string> res_xneighbors;
    auto mnet = rmnet.get_mlnet();
    auto actor = mnet->actors()->get(actor_name);

    if (!actor)
    {
        stop("actor " + actor_name + " not found");
    }

    auto layers = resolve_layers_unordered(mnet,layer_names);
    auto mode = resolve_mode(mode_name);
    auto actors = uu::net::xneighbors(mnet, layers.begin(), layers.end(), actor, mode);

    for (auto neigh: actors)
    {
        res_xneighbors.insert(neigh->name);
    }

    return res_xneighbors;
}


// NETWORK MANIPULATION

void
addLayers(
    RMLNetwork& rmnet,
    const CharacterVector& layer_names,
    const LogicalVector& directed
)
{
    auto mnet = rmnet.get_mlnet();

    if (directed.size()==1)
    {
        for (size_t i=0; i<layer_names.size(); i++)
        {
            auto layer_name = std::string(layer_names[i]);
            auto dir = directed[0]?uu::net::EdgeDir::DIRECTED:uu::net::EdgeDir::UNDIRECTED;
            auto layer = std::make_unique<G>(layer_name, dir, true);
            mnet->layers()->add(std::move(layer));
        }
    }

    else if (layer_names.size()!=directed.size())
    {
        stop("Same number of layer names and layer directionalities expected");
    }

    else
    {
        for (size_t i=0; i<layer_names.size(); i++)
        {
            auto layer_name = std::string(layer_names[i]);
            auto dir = directed[i]?uu::net::EdgeDir::DIRECTED:uu::net::EdgeDir::UNDIRECTED;
            auto layer = std::make_unique<G>(layer_name, dir, true);
            mnet->layers()->add(std::move(layer));
        }
    }
}

void
addActors(
    RMLNetwork& rmnet,
    const CharacterVector& actor_names
)
{
    auto mnet = rmnet.get_mlnet();

    for (size_t i=0; i<actor_names.size(); i++)
    {
        auto actor_name = std::string(actor_names[i]);
        mnet->actors()->add(actor_name);
    }
}

void
addNodes(
    RMLNetwork& rmnet,
    const DataFrame& vertices)
{
    auto mnet = rmnet.get_mlnet();

    CharacterVector a = vertices(0);
    CharacterVector l = vertices(1);

    for (size_t i=0; i<a.size(); i++)
    {
        auto actor_name = std::string(a[i]);
        mnet->actors()->add(actor_name);
    }
    
    for (size_t i=0; i<vertices.nrow(); i++)
    {
        auto actor = mnet->actors()->get(std::string(a(i)));

        if (!actor)
        {
            actor = mnet->actors()->add(std::string(a(i)));
        }
        
        auto layer = mnet->layers()->get(std::string(l(i)));

        if (!layer)
        {
            auto dir = uu::net::EdgeDir::UNDIRECTED;
            auto layer_ = std::make_unique<G>(std::string(l(i)), dir, true);
            layer = mnet->layers()->add(std::move(layer_));
        }

        layer->vertices()->add(actor);
    }
}

void
addEdges(
    RMLNetwork& rmnet,
    const DataFrame& edges)
{
    auto mnet = rmnet.get_mlnet();

    CharacterVector a_from = edges(0);
    CharacterVector l_from = edges(1);
    CharacterVector a_to = edges(2);
    CharacterVector l_to = edges(3);

    
    for (size_t i=0; i<edges.nrow(); i++)
    {
        auto actor1 = mnet->actors()->get(std::string(a_from(i)));

        if (!actor1)
        {
            actor1 = mnet->actors()->add(std::string(a_from(i)));
        }
        
        auto actor2 = mnet->actors()->get(std::string(a_to(i)));

        if (!actor2)
        {
            actor2 = mnet->actors()->add(std::string(a_to(i)));
        }
        
        auto layer1 = mnet->layers()->get(std::string(l_from(i)));

        if (!layer1)
        {
            auto dir = uu::net::EdgeDir::UNDIRECTED;
            auto layer_ = std::make_unique<G>(std::string(l_from(i)), dir, true);
            layer1 = mnet->layers()->add(std::move(layer_));
        }
        
        layer1->vertices()->add(actor1);

        auto layer2 = mnet->layers()->get(std::string(l_to(i)));

        if (!layer2)
        {
            auto dir = uu::net::EdgeDir::UNDIRECTED;
            auto layer_ = std::make_unique<G>(std::string(l_to(i)), dir, true);
            layer2 = mnet->layers()->add(std::move(layer_));
        }

        layer2->vertices()->add(actor2);
        

        if (layer1==layer2)
        {
            layer1->edges()->add(actor1, actor2);
        }

        else
        {
            mnet->interlayer_edges()->add(actor1, layer1, actor2, layer2);
        }
    }
}

void
setDirected(
    const RMLNetwork& rmnet,
    const DataFrame& layers_dir)
{
    auto mnet = rmnet.get_mlnet();
    CharacterVector l1 = layers_dir(0);
    CharacterVector l2 = layers_dir(1);
    NumericVector dir = layers_dir(2);

    for (size_t i=0; i<layers_dir.nrow(); i++)
    {
        auto layer1 = mnet->layers()->get(std::string(l1(i)));

        if (!layer1)
        {
            stop("cannot find layer " + std::string(l1(i)));
        }

        auto layer2 = mnet->layers()->get(std::string(l2(i)));

        if (!layer2)
        {
            stop("cannot find layer " + std::string(l2(i)));
        }

        int directed = (int)dir(i);

        if (directed!=0 && directed!=1)
        {
            stop("directionality can only be 0 or 1");
        }

        if (layer1==layer2)
        {
            // @todo do nothing?
        }
        else
        {
            mnet->interlayer_edges()->set_directed(layer1, layer2, directed);
        }
    }
}

void
deleteLayers(
    RMLNetwork& rmnet,
    const CharacterVector& layer_names)
{
    auto mnet = rmnet.get_mlnet();

    for (size_t i=0; i<layer_names.size(); i++)
    {
        auto layer = mnet->layers()->get(std::string(layer_names(i)));
        mnet->layers()->erase(layer);
    }
}

void
deleteActors(
    RMLNetwork& rmnet,
    const CharacterVector& actor_names)
{
    auto mnet = rmnet.get_mlnet();

    for (size_t i=0; i<actor_names.size(); i++)
    {
        auto actor = mnet->actors()->get(std::string(actor_names(i)));
        mnet->actors()->erase(actor);
    }
}

void
deleteNodes(
    RMLNetwork& rmnet,
    const DataFrame& vertex_matrix
)
{
    auto mnet = rmnet.get_mlnet();
    auto vertices = resolve_vertices(mnet, vertex_matrix);

    for (auto vertex: vertices)
    {
        auto actor = vertex.first;
        auto layer = vertex.second;
        layer->vertices()->erase(actor);
    }
}

void
deleteEdges(
    RMLNetwork& rmnet,
    const DataFrame& edge_matrix
)
{
    auto mnet = rmnet.get_mlnet();
    auto edges = resolve_edges(mnet, edge_matrix);

    for (auto edge: edges)
    {
        auto actor1 = std::get<0>(edge);
        auto layer1 = std::get<1>(edge);
        auto actor2 = std::get<2>(edge);
        auto layer2 = std::get<3>(edge);
        if (layer1 == layer2)
        {
            auto e = layer1->edges()->get(actor1, actor2);
            layer1->edges()->erase(e);
        }
        else
        {
            auto e = mnet->interlayer_edges()->get(actor1, layer1, actor2, layer2);
            mnet->interlayer_edges()->erase(e);
        }
    }
}


void
newAttributes(
    RMLNetwork& rmnet,
    const CharacterVector& attribute_names,
    const std::string& type,
    const std::string& target,
    const std::string& layer_name,
    const std::string& layer_name1,
    const std::string& layer_name2
)
{
    auto mnet = rmnet.get_mlnet();

    uu::core::AttributeType a_type;

    if (type=="string")
    {
        a_type = uu::core::AttributeType::STRING;
    }

    else if (type=="numeric")
    {
        a_type = uu::core::AttributeType::DOUBLE;
    }

    else
    {
        stop("Wrong type");
    }

    if (target=="actor")
    {
        if (layer_name!="" || layer_name1!="" || layer_name2!="")
        {
            stop("No layers should be specified for target 'actor'");
        }

        for (size_t i=0; i<attribute_names.size(); i++)
        {
            mnet->actors()->attr()->add(std::string(attribute_names[i]),a_type);
        }
    }

    else if (target=="layer")
    {
        stop("layer attributes are not available in this version of the library");

        /*
        if (layer_name!="" || layer_name1!="" || layer_name2!="")
            stop("No layers should be specified for target 'layer'");
        for (size_t i=0; i<attribute_names.size(); i++) {
            mnet->layer_features()->add(std::string(attribute_names[i]),a_type);
        }*/
    }

    else if (target=="node" || target=="vertex")
    {
        if (target=="node")
        {
            Rf_warning("target 'node' deprecated: use 'vertex' instead");
        }

        if (layer_name1!="" || layer_name2!="")
        {
            stop("layer1 and layer2 should not be specified for target '" + target + "'");
        }

        auto layer = mnet->layers()->get(layer_name);

        if (!layer)
        {
            stop("layer " + layer_name + " not found");
        }

        for (size_t i=0; i<attribute_names.size(); i++)
        {
            layer->vertices()->attr()->add(std::string(attribute_names[i]),a_type);
        }
    }

    else if (target=="edge")
    {
        if (layer_name!="" && (layer_name1!="" || layer_name2!=""))
        {
            stop("either layers (for intra-layer edges) or layers1 and layers2 (for inter-layer edges) must be specified for target 'edge'");
        }

        uu::net::Network* layer1;
        uu::net::Network* layer2;

        if (layer_name1=="")
        {
            layer1 = mnet->layers()->get(layer_name);
            layer2 = layer1;

            if (!layer1)
            {
                stop("layer " + layer_name + " not found");
            }
        }

        else if (layer_name2!="")
        {
            layer1 = mnet->layers()->get(layer_name1);
            layer2 = mnet->layers()->get(layer_name2);
        }

        else
        {
            stop("if layer1 is specified, also layer2 is required");
        }

        if (layer1 == layer2)
        {
            for (size_t i=0; i<attribute_names.size(); i++)
            {
                layer1->edges()->attr()->add(std::string(attribute_names[i]),a_type);
            }
        }

        else
        {
            stop("attributes on inter-layer edges are not available in this version of the library");
        }

    }

    else
    {
        stop("wrong target: " + target);
    }
}


DataFrame
getAttributes(
    const RMLNetwork& rmnet,
    const std::string& target
)
{
    auto mnet = rmnet.get_mlnet();

    if (target=="actor")
    {
        auto attributes = mnet->actors()->attr();
        CharacterVector a_name, a_type;

        for (auto att: *attributes)
        {
            a_name.push_back(att->name);
            a_type.push_back(uu::core::to_string(att->type));
        }

        return DataFrame::create(_["name"] = a_name, _["type"] = a_type);
    }

    else if (target=="layer")
    {
        stop("layer attributes are not available in this version of the library");
        /*
        auto store = mnet->layer_features();
        CharacterVector a_name, a_type;
        for (auto att: store->attributes()) {
            a_name.push_back(att->name);
            a_type.push_back(uu::core::to_string(att->type));
        }
        return DataFrame::create(_["name"] = a_name, _["type"] = a_type);
         */
    }

    else if (target=="node" || target=="vertex")
    {
        if (target=="node")
        {
            Rf_warning("target 'node' deprecated: use 'vertex' instead");
        }

        CharacterVector a_layer, a_name, a_type;

        for (auto layer: *mnet->layers())
        {
            auto attributes = layer->vertices()->attr();

            for (auto att: *attributes)
            {
                a_layer.push_back(layer->name);
                a_name.push_back(att->name);
                a_type.push_back(uu::core::to_string(att->type));
            }
        }

        return DataFrame::create(_["layer"] = a_layer, _["name"] = a_name, _["type"] = a_type);
    }

    else if (target=="edge")
    {
        /*
        CharacterVector a_layer1, a_layer2, a_name, a_type;
        for (auto layer1: *mnet->get_layers()) {
            for (auto layer2: *mnet->get_layers()) {
                if (!mnet->is_directed(layer1,layer2) && layer1->name>layer2->name) continue;
                auto store=mnet->edge_features(layer1,layer2);
                for (auto att: store->attributes()) {
                    a_layer1.push_back(layer1->name);
                    a_layer2.push_back(layer2->name);
                    a_name.push_back(att->name);
                    a_type.push_back(uu::core::to_string(att->type));
                }
            }
        }
        return DataFrame::create(_["layer1"] = a_layer1, _["layer2"] = a_layer2, _["name"] = a_name, _["type"] = a_type);
         */

        CharacterVector a_layer, a_name, a_type;

        for (auto layer: *mnet->layers())
        {
            auto attributes = layer->edges()->attr();

            for (auto att: *attributes)
            {
                a_layer.push_back(layer->name);
                a_name.push_back(att->name);
                a_type.push_back(uu::core::to_string(att->type));
            }
        }

        auto attributes = mnet->interlayer_edges()->attr();
        
        for (auto att: *attributes)
        {
            a_layer.push_back("--");
            a_name.push_back(att->name);
            a_type.push_back(uu::core::to_string(att->type));
        }
        
        return DataFrame::create(_["layer"] = a_layer, _["name"] = a_name, _["type"] = a_type);
    }

    else
    {
        stop("wrong target: " + target);
    }

    return 0; // never gets here
}


DataFrame
getValues(
    RMLNetwork& rmnet,
    const std::string& attribute_name,
    const CharacterVector& actor_names,
    const DataFrame& vertex_matrix,
    const DataFrame& edge_matrix
)
{
    auto mnet = rmnet.get_mlnet();

    if (actor_names.size() != 0)
    {
        if (vertex_matrix.size() > 0)
        {
            Rcout << "Warning: unused parameter: \"vertices\"" << std::endl;
        }

        if (edge_matrix.size() > 0)
        {
            Rcout << "Warning: unused parameter: \"edges\"" << std::endl;
        }

        auto actors = resolve_actors(mnet,actor_names);
        auto attributes = mnet->actors()->attr();
        auto att = attributes->get(attribute_name);

        if (!att)
        {
            stop("cannot find attribute: " + attribute_name + " for actors");
        }

        if (att->type==uu::core::AttributeType::DOUBLE)
        {
            NumericVector value;

            for (auto actor: actors)
            {
                value.push_back(attributes->get_double(actor,att->name).value);
            }

            return DataFrame::create(_["value"] = value);
        }

        else if (att->type==uu::core::AttributeType::STRING)
        {
            CharacterVector value;

            for (auto actor: actors)
            {
                value.push_back(attributes->get_string(actor,att->name).value);
            }

            return DataFrame::create(_["value"] = value);
        }

        else
        {
            stop("attribute type not supported: " + uu::core::to_string(att->type));
        }
    }

    // local attributes: vertices
    // must be from a single layer
    else if (vertex_matrix.size() > 0)
    {
        if (edge_matrix.size() > 0)
        {
            Rcout << "Warning: unused parameter: \"edges\"" << std::endl;
        }

        auto vertices = resolve_vertices(mnet,vertex_matrix);

        // Check if the attribute is available for all layers
        const uu::core::Attribute* att;
        std::set<const G*> layers;
        std::set<uu::core::AttributeType> types;
        for (auto vertex: vertices)
        {
            auto layer = vertex.second;
            layers.insert(layer);
        }
        for (auto layer: layers)
        {
            auto attributes = layer->vertices()->attr();
            att = attributes->get(attribute_name);
            if (!att)
            {
                stop("cannot find attribute: " + attribute_name + " for vertices on layer " + layer->name);
            }
            types.insert(att->type);
        }
        if (types.size() > 1)
        {
            stop("different attribute types on different layers");
        }
        
        auto attribute_type = *types.begin();
        
        if (attribute_type==uu::core::AttributeType::DOUBLE)
        {
            NumericVector value;
            
            for (auto vertex: vertices)
            {
                
                auto layer = vertex.second;
                auto attributes = layer->vertices()->attr();

            value.push_back(attributes->get_double(vertex.first,attribute_name).value);
            
        }
            
            return DataFrame::create(_["value"] = value);
        }

        else if (attribute_type == uu::core::AttributeType::STRING)
        {
            CharacterVector value;
            
            for (auto vertex: vertices)
            {
                
                auto layer = vertex.second;
                auto attributes = layer->vertices()->attr();
                
                value.push_back(attributes->get_string(vertex.first,attribute_name).value);
                
            }
            
            return DataFrame::create(_["value"] = value);
        }
        
        else
        {
            stop("attribute type not supported: " + uu::core::to_string(attribute_type));
        }
        

    }

    else if (edge_matrix.size() > 0)
    {
        auto edges = resolve_edges(mnet,edge_matrix);
        
        // Check if the attribute is available for all combinations of layers
        const uu::core::Attribute* att;
        std::set<std::pair<const G*,const G*>> layers;
        std::set<uu::core::AttributeType> types;
        for (auto edge: edges)
        {
            auto layer1 = std::get<1>(edge);
            auto layer2 = std::get<3>(edge);
            layers.insert(std::pair<const G*,const G*>(layer1, layer2));
        }
        for (auto p: layers)
        {
            auto layer1 = p.first;
            auto layer2 = p.second;
            if (layer1 == layer2)
            {
                auto attributes = layer1->edges()->attr();
                att = attributes->get(attribute_name);
                if (!att)
                {
                    stop("cannot find attribute: " + attribute_name + " for edges on layer " + layer1->name);
                }
                types.insert(att->type);
                
            }
            else
            {
                auto attributes = mnet->interlayer_edges()->attr();
                att = attributes->get(attribute_name);
                if (!att)
                {
                    stop("cannot find attribute: " + attribute_name + " for edges on layers " + layer1->name + ", " + layer2->name);
                }
                types.insert(att->type);
            }
        }
        if (types.size() > 1)
        {
            stop("different attribute types on different combinations of layers");
        }
        
        auto attribute_type = *types.begin();
        
        if (attribute_type == uu::core::AttributeType::DOUBLE)
        {
            NumericVector value;

            for (auto edge: edges)
            {
                auto actor1 = std::get<0>(edge);
                auto layer1 = std::get<1>(edge);
                auto actor2 = std::get<2>(edge);
                auto layer2 = std::get<3>(edge);
                if (layer1 == layer2)
                {
                    auto attributes = layer1->edges()->attr();
                    auto e = layer1->edges()->get(actor1, actor2);
                    value.push_back(attributes->get_double(e, attribute_name).value);
                }
                else
                {
                    auto attributes = mnet->interlayer_edges()->attr();
                    auto e = mnet->interlayer_edges()->get(actor1, layer1, actor2, layer2);
                    value.push_back(attributes->get_double(e, attribute_name).value);
                }
            }
            return DataFrame::create(_["value"] = value);
        }

        else if (attribute_type == uu::core::AttributeType::STRING)
        {
            CharacterVector value;

            for (auto edge: edges)
            {
                auto actor1 = std::get<0>(edge);
                auto layer1 = std::get<1>(edge);
                auto actor2 = std::get<2>(edge);
                auto layer2 = std::get<3>(edge);
                if (layer1 == layer2)
                {
                    auto attributes = layer1->edges()->attr();
                    auto e = layer1->edges()->get(actor1, actor2);
                    value.push_back(attributes->get_string(e, attribute_name).value);
                }
                else
                {
                    auto attributes = mnet->interlayer_edges()->attr();
                    auto e = mnet->interlayer_edges()->get(actor1, layer1, actor2, layer2);
                    value.push_back(attributes->get_string(e, attribute_name).value);
                }
            }
            return DataFrame::create(_["value"] = value);
        }

        else
        {
            stop("attribute type not supported: " + uu::core::to_string(attribute_type));
        }
    }

    else
    {
        stop("Required at least one parameter: \"actors\", \"vertices\" or \"edges\"");
    }

    // Never gets here
    return 0;
}

void
setValues(
    RMLNetwork& rmnet,
    const std::string& attribute_name,
    const CharacterVector& actor_names,
    const DataFrame& vertex_matrix,
    const DataFrame& edge_matrix,
    const GenericVector& values
)
{
    auto mnet = rmnet.get_mlnet();

    if (actor_names.size() != 0)
    {
        if (actor_names.size() != values.size() && values.size()!=1)
        {
            stop("wrong number of values");
        }

        if (vertex_matrix.size() > 0)
        {
            Rcout << "Warning: unused parameter: \"vertices\"" << std::endl;
        }

        if (edge_matrix.size() > 0)
        {
            Rcout << "Warning: unused parameter: \"edges\"" << std::endl;
        }

        auto actors = resolve_actors(mnet,actor_names);
        auto attributes = mnet->actors()->attr();
        auto att = attributes->get(attribute_name);

        if (!att)
        {
            stop("cannot find attribute: " + attribute_name + " for actors");
        }

        size_t i=0;

        for (auto actor: actors)
        {
            switch (att->type)
            {
            case uu::core::AttributeType::NUMERIC:
            case uu::core::AttributeType::DOUBLE:
                if (values.size()==1)
                {
                    attributes->set_double(actor,att->name,as<double>(values[0]));
                }

                else
                {
                    attributes->set_double(actor,att->name,as<double>(values[i]));
                }

                break;

            case uu::core::AttributeType::STRING:
                if (values.size()==1)
                {
                    attributes->set_string(actor,att->name,as<std::string>(values[0]));
                }

                else
                {
                    attributes->set_string(actor,att->name,as<std::string>(values[i]));
                }

                break;

            case uu::core::AttributeType::TEXT:
            case uu::core::AttributeType::TIME:
            case uu::core::AttributeType::INTEGER:
                stop("attribute type not supported: " + uu::core::to_string(att->type));

            }

            i++;
        }
    }

    // local attributes: vertices
    else if (vertex_matrix.size() > 0)
    {
        if (edge_matrix.size() > 0)
        {
            Rcout << "Warning: unused parameter: \"edges\"" << std::endl;
        }

        auto vertices = resolve_vertices(mnet,vertex_matrix);

        if (vertices.size() != values.size() && values.size()!=1)
        {
            stop("wrong number of values");
        }

        size_t i=0;

        for (auto vertex: vertices)
        {
            auto attributes = vertex.second->vertices()->attr();
            auto att = attributes->get(attribute_name);

            if (!att)
            {
                stop("cannot find attribute: " + attribute_name + " for vertices on layer " + vertex.second->name);
            }

            switch (att->type)
            {
            case uu::core::AttributeType::NUMERIC:
            case uu::core::AttributeType::DOUBLE:
                if (values.size()==1)
                {
                    attributes->set_double(vertex.first,att->name,as<double>(values[0]));
                }

                else
                {
                    attributes->set_double(vertex.first,att->name,as<double>(values[i]));
                }

                    i++;
                    
                break;

            case uu::core::AttributeType::STRING:
                if (values.size()==1)
                {
                    attributes->set_string(vertex.first,att->name,as<std::string>(values[0]));
                }

                else
                {
                    attributes->set_string(vertex.first,att->name,as<std::string>(values[i]));
                }

                    i++;
                    
                break;

            case uu::core::AttributeType::TEXT:
            case uu::core::AttributeType::TIME:
            case uu::core::AttributeType::INTEGER:
                stop("attribute type not supported: " + uu::core::to_string(att->type));

            }

        }
    }

    else if (edge_matrix.size() > 0)
    {
        auto edges = resolve_edges(mnet,edge_matrix);

        if (edges.size() != values.size() && values.size()!=1)
        {
            stop("wrong number of values");
        }

        size_t i=0;

        for (auto edge: edges)
        {
         
         auto actor1 = std::get<0>(edge);
         auto layer1 = std::get<1>(edge);
         auto actor2 = std::get<2>(edge);
         auto layer2 = std::get<3>(edge);
         
            if (layer1 == layer2)
            {
                
            
            auto attributes = layer1->edges()->attr();
            auto att = attributes->get(attribute_name);
                auto e = layer1->edges()->get(actor1, actor2);

            if (!att)
            {
                stop("cannot find attribute: " + attribute_name + " for edges on layer " + layer1->name);
            }

            switch (att->type)
            {
            case uu::core::AttributeType::NUMERIC:
            case uu::core::AttributeType::DOUBLE:
                if (values.size()==1)
                {
                    attributes->set_double(e,att->name,as<double>(values[0]));
                }

                else
                {
                    attributes->set_double(e,att->name,as<double>(values[i]));
                }

                    i++;
                    
                break;

            case uu::core::AttributeType::STRING:
                if (values.size()==1)
                {
                    attributes->set_string(e,att->name,as<std::string>(values[0]));
                }

                else
                {
                    attributes->set_string(e,att->name,as<std::string>(values[i]));
                }

                    i++;
                    
                break;

            case uu::core::AttributeType::TEXT:
            case uu::core::AttributeType::TIME:
            case uu::core::AttributeType::INTEGER:
                stop("attribute type not supported: " + uu::core::to_string(att->type));

            }
            }

            else
            {
                
                auto attributes = mnet->interlayer_edges()->attr();
                auto att = attributes->get(attribute_name);
                auto e = mnet->interlayer_edges()->get(actor1, layer1, actor2, layer2);
                
                if (!att)
                {
                    stop("cannot find attribute: " + attribute_name + " for edges on layers " + layer1->name +
                         ", " + layer2->name);
                }
                
                switch (att->type)
                {
                    case uu::core::AttributeType::NUMERIC:
                    case uu::core::AttributeType::DOUBLE:
                    if (values.size()==1)
                    {
                        attributes->set_double(e,att->name,as<double>(values[0]));
                    }
                    
                    else
                    {
                        attributes->set_double(e,att->name,as<double>(values[i]));
                    }
                    
                        i++;
                        
                    break;
                    
                    case uu::core::AttributeType::STRING:
                    if (values.size()==1)
                    {
                        attributes->set_string(e,att->name,as<std::string>(values[0]));
                    }
                    
                    else
                    {
                        attributes->set_string(e,att->name,as<std::string>(values[i]));
                    }
                    
                        i++;
                        
                    break;
                    
                    case uu::core::AttributeType::TEXT:
                    case uu::core::AttributeType::TIME:
                    case uu::core::AttributeType::INTEGER:
                    stop("attribute type not supported: " + uu::core::to_string(att->type));
                    
                }
            }
        }
    }

    else
    {
        stop("Required at least one parameter: \"actors\", \"vertices\" or \"edges\"");
    }
}

// TRANSFORMATION

void
flatten(
    RMLNetwork& rmnet,
    const std::string& new_layer_name,
    const CharacterVector& layer_names,
    const std::string& method,
    bool force_directed,
    bool all_actors
)
{

    // @todo
    if (all_actors)
    {
        stop("option to include all actors not currently implemented");
    }

    auto mnet = rmnet.get_mlnet();

    auto layers = resolve_layers_unordered(mnet,layer_names);

    bool directed = force_directed;

    if (!force_directed)
    {
        for (auto layer: layers)
        {
            if (layer->is_directed())
            {
                directed = true;
                break;
            }
        }
    }

    auto edge_directionality = directed?uu::net::EdgeDir::DIRECTED:uu::net::EdgeDir::UNDIRECTED;


    auto new_layer = std::make_unique<G>(new_layer_name, edge_directionality, true);
    new_layer->edges()->attr()->add("weight", uu::core::AttributeType::DOUBLE);

    auto target = mnet->layers()->add(std::move(new_layer));

    if (method=="weighted")
    {
        uu::net::weighted_graph_union(layers.begin(),layers.end(),target,"weight");
    }

    else if (method=="or")
    {
        // todo replace with new union
        for (auto g=layers.begin(); g!=layers.end(); ++g)
        {
            uu::net::graph_add(*g, target);
        }
    }

    else
    {
        stop("Unexpected value: method");
    }
}



void project(
    RMLNetwork& rmnet,
    const std::string& new_layer,
    const std::string& layer_name1,
    const std::string& layer_name2,
    const std::string& method) {
    auto mnet = rmnet.get_mlnet();
    auto layer1 = mnet->layers()->get(layer_name1);
    auto layer2 = mnet->layers()->get(layer_name2);
    if (!layer1 || !layer2)
        stop("Layer not found");
    if (method=="clique")
    {
        auto target = std::make_unique<G>(new_layer, uu::net::EdgeDir::UNDIRECTED, true);
        auto target_ptr = mnet->layers()->add(std::move(target));
        uu::net::project_unweighted(mnet, layer1, layer2, target_ptr);
    }
    else stop("Unexpected value: algorithm");
}

// MEASURES

NumericVector
degree_ml(
    const RMLNetwork& rmnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type
)
{
    auto mnet = rmnet.get_mlnet();

    auto actors = resolve_actors(mnet,actor_names);
    auto layers = resolve_layers_unordered(mnet,layer_names);
    NumericVector res(actors.size());

    size_t i = 0;
    for (auto actor: actors)
    {
        long deg = 0;
        auto mode = resolve_mode(type);
        deg = degree(layers.begin(), layers.end(), actor, mode);

        if (deg==0)
        {
            // check if the actor is missing from all layer_names
            bool is_missing = true;

            for (auto layer: layers)
            {
                if (layer->vertices()->contains(actor))
                {
                    is_missing = false;
                }
            }

            if (is_missing)
            {
                res[i] = NA_REAL;
            }

            else
            {
                res[i] = 0;
            }
        }

        else
        {
            res[i] = deg;
        }
        i++;
    }

    return res;
}


NumericVector
degree_deviation_ml(
    const RMLNetwork& rmnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type)
{
    auto mnet = rmnet.get_mlnet();

    auto actors = resolve_actors(mnet,actor_names);
    auto layers = resolve_layers_unordered(mnet,layer_names);
    NumericVector res(actors.size());

    size_t i = 0;
    for (auto actor: actors)
    {
        double deg = 0;
        auto mode = resolve_mode(type);
        deg = degree_deviation(layers.begin(), layers.end(), actor, mode);

        if (deg==0)
        {
            // check if the actor is missing from all layer_names
            bool is_missing = true;

            for (auto layer: layers)
            {
                if (layer->vertices()->contains(actor))
                {
                    is_missing = false;
                }
            }

            if (is_missing)
            {
                res[i] = NA_REAL;
            }

            else
            {
                res[i] = 0;
            }
        }

        else
        {
            res[i] = deg;
        }
        i++;
    }

    return res;
}

/*
NumericVector occupation_ml(
    const RMLNetwork& rmnet, const NumericMatrix& transitions, double teleportation, long steps) {
auto mnet = rmnet.get_mlnet();

if (steps==0) {
    // completely arbitrary value :)
    steps = 100*mnet->get_edges()->size();
}
if (transitions.nrow()!=transitions.ncol()) {
    stop("expected NxN matrix");
}
if (transitions.nrow()!=mnet->get_layers()->size()) {
    stop("dimensions of transition probability matrix do not match the number of layers in the network");
}
matrix<double> m(transitions.nrow());
for (size_t i=0; i<transitions.nrow(); i++) {
    std::vector<double> row(transitions.ncol());
    for (size_t j=0; j<transitions.ncol(); j++) {
        row[j] = transitions(i,j);
    }
    m[i] = row;
}

NumericVector res(actors.size());
std::unordered_map<const uu::net::Vertex*, int > occ = occupation(mnet,teleportation,m,steps);

for (const auto &p : occ) {
    res[p.first->name] = p.second;
}
return res;
}

*/

NumericVector
neighborhood_ml(
    const RMLNetwork& rmnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type
)
{
    auto mnet = rmnet.get_mlnet();

    auto actors = resolve_actors(mnet,actor_names);
    auto layers = resolve_layers_unordered(mnet,layer_names);
    NumericVector res(actors.size());

    size_t i = 0;
    for (auto actor: actors)
    {
        long neigh = 0;
        auto mode = resolve_mode(type);
        neigh = neighbors(layers.begin(), layers.end(), actor, mode).size();

        if (neigh==0)
        {
            // check if the actor is missing from all layer_names
            bool is_missing = true;

            for (auto layer: layers)
            {
                if (layer->vertices()->contains(actor))
                {
                    is_missing = false;
                }
            }

            if (is_missing)
            {
                res[i] = NA_REAL;
            }

            else
            {
                res[i] = 0;
            }
        }

        else
        {
            res[i] = neigh;
        }
        i++;
    }

    return res;
}



NumericVector
xneighborhood_ml(
    const RMLNetwork& rmnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type)
{
    auto mnet = rmnet.get_mlnet();

    auto actors = resolve_actors(mnet,actor_names);
    auto layers = resolve_layers_unordered(mnet,layer_names);
    NumericVector res(actors.size());

    size_t i = 0;
    for (auto actor: actors)
    {
        long neigh = 0;
        auto mode = resolve_mode(type);
        neigh = xneighbors(mnet, layers.begin(), layers.end(), actor, mode).size();

        if (neigh==0)
        {
            // check if the actor is missing from all layer_names
            bool is_missing = true;

            for (auto layer: layers)
            {
                if (layer->vertices()->contains(actor))
                {
                    is_missing = false;
                }
            }

            if (is_missing)
            {
                res[i] = NA_REAL;
            }

            else
            {
                res[i] = 0;
            }
        }

        else
        {
            res[i] = neigh;
        }
        i++;
    }

    return res;
}

NumericVector
connective_redundancy_ml(
    const RMLNetwork& rmnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type)
{
    auto mnet = rmnet.get_mlnet();

    auto actors = resolve_actors(mnet,actor_names);
    auto layers = resolve_layers_unordered(mnet,layer_names);
    NumericVector res(actors.size());
    double cr = 0;

    size_t i = 0;
    for (auto actor: actors)
    {
        auto mode = resolve_mode(type);

        cr = uu::net::connective_redundancy(mnet, layers.begin(), layers.end(), actor, mode);

        if (cr==0)
        {
            // check if the actor is missing from all layer_names
            bool is_missing = true;

            for (auto layer: layers)
            {
                if (layer->vertices()->contains(actor))
                {
                    is_missing = false;
                }
            }

            if (is_missing)
            {
                res[i] = NA_REAL;
            }

            else
            {
                res[i] = 0;
            }
        }

        else
        {
            res[i] = cr;
        }
        i++;
    }

    return res;
}

NumericVector
relevance_ml(
    const RMLNetwork& rmnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type)
{
    auto mnet = rmnet.get_mlnet();

    auto actors = resolve_actors(mnet,actor_names);
    auto layers = resolve_layers_unordered(mnet,layer_names);
    NumericVector res(actors.size());

    size_t i = 0;
    for (auto actor: actors)
    {
        double rel = 0;
        auto mode = resolve_mode(type);
        rel = uu::net::relevance(mnet, layers.begin(), layers.end(), actor, mode);

        if (rel==0)
        {
            // check if the actor is missing from all layer_names
            bool is_missing = true;

            for (auto layer: layers)
            {
                if (layer->vertices()->contains(actor))
                {
                    is_missing = false;
                }
            }

            if (is_missing)
            {
                res[i] = NA_REAL;
            }

            else
            {
                res[i] = 0;
            }
        }

        else
        {
            res[i] = rel;
        }
        i++;
    }

    return res;
}


NumericVector
xrelevance_ml(
    const RMLNetwork& rmnet,
    const CharacterVector& actor_names,
    const CharacterVector& layer_names,
    const std::string& type)
{
    auto mnet = rmnet.get_mlnet();

    auto actors = resolve_actors(mnet,actor_names);
    auto layers = resolve_layers_unordered(mnet,layer_names);

    NumericVector res(actors.size());

    size_t i = 0;
    for (auto actor: actors)
    {
        double rel = 0;
        auto mode = resolve_mode(type);
        rel = uu::net::xrelevance(mnet, layers.begin(), layers.end(), actor, mode);

        if (rel==0)
        {
            // check if the actor is missing from all layer_names
            bool is_missing = true;

            for (auto layer: layers)
            {
                if (layer->vertices()->contains(actor))
                {
                    is_missing = false;
                }
            }

            if (is_missing)
            {
                res[i] = NA_REAL;
            }

            else
            {
                res[i] = 0;
            }
        }

        else
        {
            res[i] = rel;
        }
        i++;
    }

    return res;
}

DataFrame
comparison_ml(
    const RMLNetwork& rmnet,
    const CharacterVector& layer_names,
    const std::string& method,
    const std::string& type,
    int K
)
{

    auto mnet = rmnet.get_mlnet();
    std::vector<uu::net::Network*> layers = resolve_layers(mnet,layer_names);
    std::vector<NumericVector> values;

    for (size_t i=0; i<layers.size(); i++)
    {
        NumericVector v;
        values.push_back(v);
    }

    DataFrame res = DataFrame::create();

    if (method=="jaccard.actors")
    {
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,bool> P = uu::net::actor_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::jaccard<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="coverage.actors")
    {
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,bool> P = uu::net::actor_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::coverage<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="kulczynski2.actors")
    {
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,bool> P = uu::net::actor_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::kulczynski2<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="sm.actors")
    {
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,bool> P = uu::net::actor_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::simple_matching<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="rr.actors")
    {
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,bool> P = uu::net::actor_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::russell_rao<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="hamann.actors")
    {
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,bool> P = uu::net::actor_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::hamann<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="jaccard.edges")
    {
        auto P = uu::net::edge_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::jaccard<std::pair<const typename M::vertex_type*,const typename M::vertex_type*>, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="coverage.edges")
    {
        auto P = uu::net::edge_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::coverage<std::pair<const typename M::vertex_type*,const typename M::vertex_type*>, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="kulczynski2.edges")
    {
        auto P = uu::net::edge_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::kulczynski2<std::pair<const typename M::vertex_type*,const typename M::vertex_type*>, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="sm.edges")
    {
        auto P = uu::net::edge_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::simple_matching<std::pair<const typename M::vertex_type*,const typename M::vertex_type*>, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="rr.edges")
    {
        auto P = uu::net::edge_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::russell_rao<std::pair<const typename M::vertex_type*,const typename M::vertex_type*>, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="hamann.edges")
    {
        auto P = uu::net::edge_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::hamann<std::pair<const typename M::vertex_type*,const typename M::vertex_type*>, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="jaccard.triangles")
    {
        uu::core::PropertyMatrix<uu::net::Triad, const uu::net::Network*,bool> P = uu::net::triangle_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::jaccard<uu::net::Triad, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="coverage.triangles")
    {
        uu::core::PropertyMatrix<uu::net::Triad, const uu::net::Network*,bool> P = uu::net::triangle_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::coverage<uu::net::Triad, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="kulczynski2.triangles")
    {
        uu::core::PropertyMatrix<uu::net::Triad, const uu::net::Network*,bool> P = uu::net::triangle_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::kulczynski2<uu::net::Triad, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="sm.triangles")
    {
        uu::core::PropertyMatrix<uu::net::Triad, const uu::net::Network*,bool> P = uu::net::triangle_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::simple_matching<uu::net::Triad, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="rr.triangles")
    {
        uu::core::PropertyMatrix<uu::net::Triad, const uu::net::Network*,bool> P = uu::net::triangle_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::russell_rao<uu::net::Triad, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="hamann.triangles")
    {
        uu::core::PropertyMatrix<uu::net::Triad, const uu::net::Network*,bool> P = uu::net::triangle_existence_property_matrix(mnet);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::hamann<uu::net::Triad, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="dissimilarity.degree")
    {
        auto mode = resolve_mode(type);
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,double> P = uu::net::actor_degree_property_matrix(mnet,mode);

        if (K<=0)
        {
            K=std::ceil(std::log2(P.num_structures) + 1);
        }

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::dissimilarity_index<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j],K));
            }
        }
    }

    else if (method=="KL.degree")
    {
        auto mode = resolve_mode(type);
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,double> P = uu::net::actor_degree_property_matrix(mnet,mode);

        if (K<=0)
        {
            K=std::ceil(std::log2(P.num_structures) + 1);
        }

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::KL_divergence<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j],K));
            }
        }
    }

    else if (method=="jeffrey.degree")
    {
        auto mode = resolve_mode(type);
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,double> P = uu::net::actor_degree_property_matrix(mnet,mode);

        if (K<=0)
        {
            K=std::ceil(std::log2(P.num_structures) + 1);
        }

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::jeffrey_divergence<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j],K));
            }
        }
    }

    else if (method=="pearson.degree")
    {
        auto mode = resolve_mode(type);
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,double> P = uu::net::actor_degree_property_matrix(mnet,mode);

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::pearson<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else if (method=="rho.degree")
    {
        auto mode = resolve_mode(type);
        uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,double> P = uu::net::actor_degree_property_matrix(mnet,mode);
        P.rankify();

        for (size_t j=0; j<layers.size(); j++)
        {
            for (size_t i=0; i<layers.size(); i++)
            {
                values[j].push_back(uu::core::pearson<const uu::net::Vertex*, const uu::net::Network*>(P,layers[i],layers[j]));
            }
        }
    }

    else
    {
        stop("Unexpected value: method parameter");
    }

    if (layer_names.size()==0)
    {
        CharacterVector names;

        for (auto l: layers)
        {
            names.push_back(l->name);
        }

        for (size_t i=0; i<layers.size(); i++)
        {
            res.push_back(values[i],std::string(names[i]));
        }

        res.attr("class") = "data.frame";
        res.attr("row.names") = names;
    }

    else
    {
        for (size_t i=0; i<layers.size(); i++)
        {
            res.push_back(values[i],std::string(layer_names[i]));
        }

        res.attr("class") = "data.frame";
        res.attr("row.names") = layer_names;
    }

    return res;
}

double
summary_ml(
    const RMLNetwork& rmnet,
    const std::string& layer_name,
    const std::string& method,
    const std::string& type
)
{

    auto mnet = rmnet.get_mlnet();
    auto layer = mnet->layers()->get(layer_name);

    if (!layer)
    {
        stop("no layer named " + layer_name);
    }

    auto mode = resolve_mode(type);
    uu::core::PropertyMatrix<const uu::net::Vertex*, const uu::net::Network*,double> P = uu::net::actor_degree_property_matrix(mnet,mode);

    if (method=="min.degree")
    {
        return uu::core::min<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else if (method=="max.degree")
    {
        return uu::core::max<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else if (method=="sum.degree")
    {
        return uu::core::sum<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else if (method=="mean.degree")
    {
        return uu::core::mean<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else if (method=="sd.degree")
    {
        return uu::core::sd<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else if (method=="skewness.degree")
    {
        return uu::core::skew<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else if (method=="kurtosis.degree")
    {
        return uu::core::kurt<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else if (method=="entropy.degree")
    {
        return uu::core::entropy<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else if (method=="CV.degree")
    {
        return uu::core::CV<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else if (method=="jarque.bera.degree")
    {
        return uu::core::jarque_bera<const uu::net::Vertex*, const uu::net::Network*>(P,layer);
    }

    else
    {
        stop("Unexpected value: method parameter");
    }

    return 0;
}



DataFrame
distance_ml(
    const RMLNetwork& rmnet,
    const std::string& from_actor,
    const CharacterVector& to_actors,
    const std::string& method)
{
    auto mnet = rmnet.get_mlnet();
    std::vector<const uu::net::Vertex*> actors_to = resolve_actors(mnet,to_actors);
    auto actor_from = mnet->actors()->get(from_actor);

    if (!actor_from)
    {
        stop("no actor named " + from_actor);
    }

    if (method=="multiplex")
    {
        auto dists = uu::net::pareto_distance(mnet, actor_from);

        CharacterVector from, to;
        std::vector<NumericVector> lengths;

        for (size_t i=0; i<mnet->layers()->size(); i++)
        {
            NumericVector v;
            lengths.push_back(v);
        }

        for (auto actor: actors_to)
        {
            for (auto d: dists[actor])
            {
                from.push_back(from_actor);
                to.push_back(actor->name);

                for (size_t i=0; i<mnet->layers()->size(); i++)
                {
                    lengths[i].push_back(d.length(mnet->layers()->at(i)));
                }
            }
        }

        DataFrame res = DataFrame::create(_["from"] = from, _["to"] = to);

        for (size_t i=0; i<mnet->layers()->size(); i++)
        {
            res.push_back(lengths[i],mnet->layers()->at(i)->name);
        }

        return DataFrame(res);
    }

    else
    {
        stop("Unexpected value: method");
    }

    return 0;
}

/*
NumericMatrix sir_ml(
    const RMLNetwork& rmnet, double beta, int tau, long num_iterations) {
auto mnet = rmnet.get_mlnet();
matrix<long> stats = sir(mnet, beta, tau, num_iterations);

NumericMatrix res(3,num_iterations+1);

CharacterVector colnames(0);
CharacterVector rownames(3);
rownames(0) = "S";
rownames(1) = "I";
rownames(2) = "R";
res.attr("dimnames") = List::create(rownames, colnames);

for (size_t i=0; i<3; i++) {
for (long j=0; j<num_iterations+1; j++) {
res(i,j) = stats[i][j];
}
}
return res;
} ///////
*/

// COMMUNITY DETECTION

DataFrame
cliquepercolation_ml(
    const RMLNetwork& rmnet,
    int k,
    int m
)
{
    auto mnet = rmnet.get_mlnet();

    auto com_struct = mlcpm(mnet, k, m);
    return to_dataframe(com_struct.get());
}


DataFrame
glouvain_ml(
    const RMLNetwork& rmnet,
    double gamma,
    double omega,
    int limit
)
{
    auto mnet = rmnet.get_mlnet();

    auto com_struct = uu::net::generalized_louvain<M,G>(mnet, gamma, omega, limit);

    return to_dataframe(com_struct.get());
}

DataFrame
glouvain2_ml(
    const RMLNetwork& rmnet,
    double omega
)
{
    auto mnet = rmnet.get_mlnet();
    
    auto com_struct = uu::net::glouvain2<M>(mnet, omega);
    
    return to_dataframe(com_struct.get());
}

DataFrame
infomap_ml(const RMLNetwork& rmnet,
           bool overlapping,
           bool directed,
           bool include_self_links
          )
{
    auto mnet = rmnet.get_mlnet();

    try
    {
        auto com_struct = uu::net::infomap(mnet, overlapping, directed, include_self_links);
        return to_dataframe(com_struct.get());
    }

    catch (std::exception& e)
    {
        Rcout << "Warning: could not run external library: " << e.what() << std::endl;
        Rcout << "Returning empty community set." << std::endl;
    }

    auto com_struct = std::make_unique<uu::net::CommunityStructure<uu::net::VertexLayerCommunity<const G>>>();
    return to_dataframe(com_struct.get());
}

DataFrame
abacus_ml(
    const RMLNetwork& rmnet,
    int min_actors,
    int min_layers
)
{
    auto mnet = rmnet.get_mlnet();

    try
    {
        auto pillar_com_struct = uu::net::abacus<M,G>(mnet, min_actors, min_layers);
        auto com_struct = to_vertex_layer_community_structure(pillar_com_struct.get());
        return to_dataframe(com_struct.get());
    }

    catch (std::exception& e)
    {
        Rcout << "Warning: could not run external library: " << e.what() << std::endl;
        Rcout << "Returning empty community set." << std::endl;
    }

    auto pillar_com_struct = std::make_unique<uu::net::CommunityStructure<uu::net::PillarCommunity<const G>>>();
    auto com_struct = to_vertex_layer_community_structure(pillar_com_struct.get());
    return to_dataframe(com_struct.get());

}
/*
DataFrame lart_ml(
   const RMLNetwork& rmnet, int t, double eps, double gamma) {
auto mnet = rmnet.get_mlnet();

lart k;
if (t<0) t = mnet->get_layers()->size()*3;
CommunityStructureSharedPtr community_structure = k.fit(mnet, t, eps, gamma);
return to_dataframe(community_structure);
}


*/

double
modularity_ml(
              const RMLNetwork& rmnet,
              const DataFrame& com,
              double gamma,
              double omega
              )
{
    auto mnet = rmnet.get_mlnet();
    auto communities = to_communities(com, mnet);
    return uu::net::modularity(mnet, communities.get(), omega);
}


List
to_list(
    const DataFrame& cs,
    const RMLNetwork& rmnet
)
{
    auto mnet = rmnet.get_mlnet();

    // stores at which index vertices start in a layer
    std::unordered_map<const G*, size_t> offset;
    size_t num_vertices = 0;

    for (auto layer: *mnet->layers())
    {
        offset[layer] = num_vertices;
        num_vertices += layer->vertices()->size();
    }


    std::map<int, std::map<int, std::vector<int> > > list;
    CharacterVector cs_actor = cs["actor"];
    CharacterVector cs_layer = cs["layer"];
    NumericVector cs_cid = cs["cid"];

    for (size_t i=0; i<cs.nrow(); i++)
    {
        int comm_id = cs_cid[i];
        auto layer = mnet->layers()->get(std::string(cs_layer[i]));

        if (!layer)
        {
            stop("cannot find layer " + std::string(cs_layer[i]) + " (community structure not compatible with this network?)");
        }

        int l = mnet->layers()->index_of(layer);
        auto actor = mnet->actors()->get(std::string(cs_actor[i]));

        if (!actor)
        {
            stop("cannot find actor " + std::string(cs_actor[i]) + " (community structure not compatible with this network?)");
        }

        int vertex_idx = layer->vertices()->index_of(actor);

        if (vertex_idx==-1)
        {
            stop("cannot find vertex " + std::string(cs_actor[i]) + "::" + std::string(cs_layer[i]) + " (community structure not compatible with this network?)");
        }

        int n = vertex_idx+offset[layer]+1;
        list[comm_id][l].push_back(n);
    }

    List res = List::create();

    for (auto clist: list)
    {
        for (auto llist: clist.second)
        {
            res.push_back(List::create(_["cid"]=clist.first, _["lid"]=llist.first, _["aid"]=llist.second));
        }
    }

    return res;
}


// LAYOUT

DataFrame
multiforce_ml(
    const RMLNetwork& rmnet,
    const NumericVector& w_in,
    const NumericVector& w_inter,
    const NumericVector& gravity,
    int iterations
)
{
    auto mnet = rmnet.get_mlnet();
    std::unordered_map<const G*,double> weight_in, weight_inter, weight_gr;
    auto layers = mnet->layers();

    if (w_in.size()==1)
    {
        for (size_t i=0; i<layers->size(); i++)
        {
            weight_in[layers->at(i)] = w_in[0];
        }
    }

    else if (w_in.size()==layers->size())
    {
        for (size_t i=0; i<layers->size(); i++)
        {
            weight_in[layers->at(i)] = w_in[i];
        }
    }

    else
    {
        stop("wrong dimension: internal weights (should contain 1 or num.layers.ml weights)");
    }

    if (w_inter.size()==1)
    {
        for (size_t i=0; i<layers->size(); i++)
        {
            weight_inter[layers->at(i)] = w_inter[0];
        }
    }

    else if (w_inter.size()==layers->size())
    {
        for (size_t i=0; i<layers->size(); i++)
        {
            weight_inter[layers->at(i)] = w_inter[i];
        }
    }

    else
    {
        stop("wrong dimension: external weights (should contain 1 or num.layers.ml weights)");
    }

    if (gravity.size()==1)
    {
        for (size_t i=0; i<layers->size(); i++)
        {
            weight_gr[layers->at(i)] = gravity[0];
        }
    }

    else if (gravity.size()==layers->size())
    {
        for (size_t i=0; i<layers->size(); i++)
        {
            weight_gr[layers->at(i)] = gravity[i];
        }
    }

    else
    {
        stop("wrong dimension: gravity (should contain 1 or num.layers.ml weights)");
    }

    auto coord = uu::net::multiforce(mnet, 10, 10, weight_in, weight_inter, weight_gr, iterations);

    std::unordered_map<const G*, size_t> offset;
    size_t num_rows = 0;

    for (auto layer: *mnet->layers())
    {
        num_rows += layer->vertices()->size();
    }

    CharacterVector actor_n(num_rows);
    CharacterVector layer_n(num_rows);
    NumericVector x_n(num_rows);
    NumericVector y_n(num_rows);
    NumericVector z_n(num_rows);
    int current_row=0;

    for (auto l: *mnet->layers())
    {
        for (auto a: *l->vertices())
        {

            auto n = std::make_pair(a, l);
            actor_n(current_row) = a->name;
            layer_n(current_row) = l->name;
            auto c = coord.at(n);
            x_n(current_row) = c.x;
            y_n(current_row) = c.y;
            z_n(current_row) = c.z;
            current_row++;
        }
    }

    DataFrame vertices = DataFrame::create(
                             Named("actor")=actor_n,
                             Named("layer")=layer_n,
                             Named("x")=x_n,
                             Named("y")=y_n,
                             Named("z")=z_n);

    return vertices;
}


DataFrame
circular_ml(
    const RMLNetwork& rmnet)
{
    auto mnet = rmnet.get_mlnet();

    auto coord = uu::net::circular(mnet, 10.0);

    std::unordered_map<const G*, size_t> offset;
    size_t num_rows = 0;

    for (auto layer: *mnet->layers())
    {
        num_rows += layer->vertices()->size();
    }

    CharacterVector actor_n(num_rows);
    CharacterVector layer_n(num_rows);
    NumericVector x_n(num_rows);
    NumericVector y_n(num_rows);
    NumericVector z_n(num_rows);
    int current_row=0;

    for (auto l: *mnet->layers())
    {
        for (auto a: *l->vertices())
        {

            auto n = std::make_pair(a, l);
            actor_n(current_row) = a->name;
            layer_n(current_row) = l->name;
            auto c = coord.at(n);
            x_n(current_row) = c.x;
            y_n(current_row) = c.y;
            z_n(current_row) = c.z;
            current_row++;
        }
    }

    DataFrame vertices = DataFrame::create(
                             Named("actor")=actor_n,
                             Named("layer")=layer_n,
                             Named("x")=x_n,
                             Named("y")=y_n,
                             Named("z")=z_n);

    return vertices;
}

