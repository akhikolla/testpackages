loadModule("multinet",TRUE)

# Casting of (a portion of) a multilayer network into an igraph (multi)graph. This is done by creating an intermediate graphml file and loading it as an igraph file
as.igraph.Rcpp_RMLNetwork <- function (x, layers=NULL, merge.actors=TRUE, all.actors=FALSE, ...) {
    if (is.null(layers)) {
        layers <- layers_ml(x)
    }
    temp_f <- tempfile()
    write_ml(x, temp_f, format = "graphml", layers, ",", merge.actors, all.actors)
    g <- read.graph(temp_f, format = "graphml")
    unlink(temp_f)
    g
}


# A function to convert the network into a list of igraph objects
as.list.Rcpp_RMLNetwork <- function(x, ...) {
    layer.names = sort(layers_ml(x))
    layers <- vector("list",num_layers_ml(x)+1)
    layers[[1]] <- as.igraph(x)
    names(layers)[1] <- "_flat_"
    if (num_layers_ml(x)>0) {
        for (i in 1 : num_layers_ml(x)) {
            layers[[i+1]] <- as.igraph(x,layer.names[i])
            bad.vs<-V(layers[[i+1]])[degree(layers[[i+1]]) == 0]
            # remove isolated nodes
            layers[[i+1]] <-delete.vertices(layers[[i+1]], bad.vs)
            names(layers)[i+1] <- layer.names[i]
        }
    }
    layers
}

#
add_igraph_layer_ml <- function(n, g, name)
{
    if (is.null(vertex_attr(g)$name))
    {
        stop("the igraph object must have a vertex attribute 'name' with the names of the actors")
    }
    
    add_layers_ml(n, name, is.directed(g))
    
    # no longer necessary: add_actors_ml(n, vertex_attr(g)$name) # remove from next version
    
    vertices = data.frame(actor=vertex_attr(g)$name, layer=name)
    add_vertices_ml(n, vertices)
    
    for (attr in names(vertex_attr(g)))
    {
        if (is.numeric(vertex_attr(g)[[attr]]))
        {
            add_attributes_ml(n, attributes=attr, type="numeric", target="vertex", layer=name)
            set_values_ml(n, attr, vertices=vertices, values=vertex_attr(g)[[attr]])
        }
        if (is.character(vertex_attr(g)[[attr]]))
        {
            add_attributes_ml(n, attributes=attr, type="string", target="vertex", layer=name)
            set_values_ml(n, attr, vertices=vertices, values=vertex_attr(g)[[attr]])
        }
        if (is.logical(vertex_attr(g)[[attr]]))
        {
            add_attributes_ml(n, attributes=attr, type="numeric", target="vertex", layer=name)
            set_values_ml(n, attr, vertices=vertices, values=as.numeric(vertex_attr(g)[[attr]]))
        }
    }
    
    
    edges = data.frame(
    actor1=get.edgelist(g)[,1], layer1=name,
    actor2=get.edgelist(g)[,2], layer2=name)
    
    add_edges_ml(n, edges)
    
    for (attr in names(edge_attr(g)))
    {
        if (is.numeric(edge_attr(g)[[attr]]))
        {
            add_attributes_ml(n, attributes=attr, type="numeric", target="edge", layer=name)
            set_values_ml(n, attr, edges=edges, values=edge_attr(g)[[attr]])
        }
        if (is.character(edge_attr(g)[[attr]]))
        {
            add_attributes_ml(n, attributes=attr, type="string", target="edge", layer=name)
            set_values_ml(n, attr, edges=edges, values=edge_attr(g)[[attr]])
        }
        if (is.logical(vertex_attr(g)[[attr]]))
        {
            add_attributes_ml(n, attributes=attr, type="numeric", target="vertex", layer=name)
            set_values_ml(n, attr, edges=edges, values=as.numeric(edge_attr(g)[[attr]]))
        }
    }
}
