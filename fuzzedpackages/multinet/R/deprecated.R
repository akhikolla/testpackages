loadModule("multinet",TRUE)

# Deprecated from version 3.1.0

add_actors_ml <- function(...)
{
    .Deprecated(NULL, package="multinet",
    msg="'add_actors_ml' is deprecated. From version 3.1 add_vertices_ml automatically adds the actors")
    add_actors_ml_v3(...)
}

# Deprecated from version 3.0.0

add.igraph.layer.ml <- function(...)
{
    .Deprecated("add_igraph_layer_ml",package="multinet")
    add_igraph_layer_ml(...)
}

ml.aucs <- function(...)
{
    .Deprecated("ml_aucs",package="multinet")
    ml_aucs(...)
}

ml.bankwiring <- function(...)
{
    .Deprecated("ml_bankwiring",package="multinet")
    ml_bankwiring(...)
}

ml.florentine <- function(...)
{
    .Deprecated("ml_florentine",package="multinet")
    ml_florentine(...)
}

ml.monastery <- function(...)
{
    .Deprecated("ml_monastery",package="multinet")
    ml_monastery(...)
}

ml.tailorshop <- function(...)
{
    .Deprecated("ml_tailorshop",package="multinet")
    ml_tailorshop(...)
}

ml.toy <- function(...)
{
    .Deprecated("ml_toy",package="multinet")
    ml_toy(...)
}

ml.empty <- function(...)
{
    .Deprecated("ml_empty",package="multinet")
    ml_empty(...)
}

evolution.pa.ml <- function(...)
{
    .Deprecated("evolution_pa_ml",package="multinet")
    evolution_pa_ml(...)
}

evolution.er.ml <- function(...)
{
    .Deprecated("evolution_er_ml",package="multinet")
    evolution_er_ml(...)
}

grow.ml <- function(...)
{
    .Deprecated("grow_ml",package="multinet")
    grow_ml(...)
}

read.ml <- function(...)
{
    .Deprecated("read_ml",package="multinet")
    read_ml(...)
}

write.ml <- function(...)
{
    .Deprecated("write_ml",package="multinet")
    write_ml(...)
}

layers.ml <- function(...)
{
    .Deprecated("layers_ml",package="multinet")
    layers_ml(...)
}

actors.ml <- function(...)
{
    .Deprecated("actors_ml",package="multinet")
    actors_ml(...)
}

vertices.ml <- function(...)
{
    .Deprecated("vertices_ml",package="multinet")
    vertices_ml(...)
}

edges.ml <- function(...)
{
    .Deprecated("edges_ml",package="multinet")
    edges_ml(...)
}

edges.idx.ml <- function(...)
{
    .Deprecated("edges_idx_ml",package="multinet")
    edges_idx_ml(...)
}

num.layers.ml <- function(...)
{
    .Deprecated("num_layers_ml",package="multinet")
    num_layers_ml(...)
}

num.actors.ml <- function(...)
{
    .Deprecated("num_actors_ml",package="multinet")
    num_actors_ml(...)
}

num.vertices.ml <- function(...)
{
    .Deprecated("num_vertices_ml",package="multinet")
    num_vertices_ml(...)
}

num.edges.ml <- function(...)
{
    .Deprecated("num_edges_ml",package="multinet")
    num_edges_ml(...)
}

is.directed.ml <- function(...)
{
    .Deprecated("is_directed_ml",package="multinet")
    is_directed_ml(...)
}

neighbors.ml <- function(...)
{
    .Deprecated("neighbors_ml",package="multinet")
    neighbors_ml(...)
}

xneighbors.ml <- function(...)
{
    .Deprecated("xneighbors_ml",package="multinet")
    xneighbors_ml(...)
}

add.layers.ml <- function(...)
{
    .Deprecated("add_layers_ml",package="multinet")
    add_layers_ml(...)
}

add.actors.ml <- function(...)
{
    .Deprecated("add_actors_ml",package="multinet")
    add_actors_ml(...)
}

add.vertices.ml <- function(...)
{
    .Deprecated("add_vertices_ml", package="multinet")
    add_vertices_ml(...)
}

add.edges.ml <- function(...)
{
    .Deprecated("add_edges_ml",package="multinet")
    add_edges_ml(...)
}

set.directed.ml <- function(...)
{
    .Deprecated("set_directed_ml",package="multinet")
    set_directed_ml(...)
}

delete.layers.ml <- function(...)
{
    .Deprecated("delete_layers_ml",package="multinet")
    delete_layers_ml(...)
}

delete.actors.ml <- function(...)
{
    .Deprecated("delete_actors_ml",package="multinet")
    delete_actors_ml(...)
}

delete.vertices.ml <- function(...)
{
    .Deprecated("delete_vertices_ml",package="multinet")
    delete_vertices_ml(...)
}

delete.edges.ml <- function(...)
{
    .Deprecated("delete_edges_ml",package="multinet")
    delete_edges_ml(...)
}

add.attributes.ml <- function(...)
{
    .Deprecated("add_attributes_ml",package="multinet")
    add_attributes_ml(...)
}

attributes.ml <- function(...)
{
    .Deprecated("attributes_ml",package="multinet")
    attributes_ml(...)
}

get.values.ml <- function(...)
{
    .Deprecated("get_values_ml",package="multinet")
    get_values_ml(...)
}

set.values.ml <- function(...)
{
    .Deprecated("set_values_ml",package="multinet")
    set_values_ml(...)
}

flatten.ml <- function(...)
{
    .Deprecated("flatten_ml",package="multinet")
    flatten_ml(...)
}

degree.ml <- function(...)
{
    .Deprecated("degree_ml",package="multinet")
    degree_ml(...)
}

degree.deviation.ml <- function(...)
{
    .Deprecated("degree_deviation_ml",package="multinet")
    degree_deviation_ml(...)
}

neighborhood.ml <- function(...)
{
    .Deprecated("neighborhood_ml",package="multinet")
    neighborhood_ml(...)
}

xneighborhood.ml <- function(...)
{
    .Deprecated("xneighborhood_ml",package="multinet")
    xneighborhood_ml(...)
}

connective.redundancy.ml <- function(...)
{
    .Deprecated("connective_redundancy_ml",package="multinet")
    connective_redundancy_ml(...)
}

relevance.ml <- function(...)
{
    .Deprecated("relevance_ml",package="multinet")
    relevance_ml(...)
}

xrelevance.ml <- function(...)
{
    .Deprecated("xrelevance_ml",package="multinet")
    xrelevance_ml(...)
}

layer.summary.ml <- function(...)
{
    .Deprecated("layer_summary_ml",package="multinet")
    layer_summary_ml(...)
}

layer.comparison.ml <- function(...)
{
    .Deprecated("layer_comparison_ml",package="multinet")
    layer_comparison_ml(...)
}

distance.ml <- function(...)
{
    .Deprecated("distance_ml",package="multinet")
    distance_ml(...)
}

clique.percolation.ml <- function(...)
{
    .Deprecated("clique_percolation_ml",package="multinet")
    clique_percolation_ml(...)
}

glouvain.ml <- function(...)
{
    .Deprecated("glouvain_ml",package="multinet")
    glouvain_ml(...)
}

abacus.ml <- function(...)
{
    .Deprecated("abacus_ml",package="multinet")
    abacus_ml(...)
}

infomap.ml <- function(...)
{
    .Deprecated("infomap_ml",package="multinet")
    infomap_ml(...)
}

layout.multiforce.ml <- function(...)
{
    .Deprecated("layout_multiforce_ml",package="multinet")
    layout_multiforce_ml(...)
}

layout.circular.ml <- function(...)
{
    .Deprecated("layout_circular_ml",package="multinet")
    layout_circular_ml(...)
}

get.community.list.ml <- function(...)
{
    .Deprecated("get_community_list_ml",package="multinet")
    get_community_list_ml(...)
}

nodes.ml <- function(...)
{
    .Deprecated("vertices_ml",package="multinet")
    vertices_ml(...)
}

num.nodes.ml <- function(...)
{
    .Deprecated("num_vertices_ml",package="multinet")
    num_vertices_ml(...)
}

add.nodes.ml <- function(...)
{
    .Deprecated("add_vertices_ml",package="multinet")
    add_vertices_ml(...)
}

delete.nodes.ml <- function(...)
{
    .Deprecated("delete_vertices_ml",package="multinet")
    delete_vertices_ml(...)
}

new.attributes.ml <- function(...)
{
    .Deprecated("add_attributes_ml",package="multinet")
    add_attributes_ml(...)
}

list.attributes.ml <- function(...)
{
    .Deprecated("attributes_ml",package="multinet")
    attributes_ml(...)
}


