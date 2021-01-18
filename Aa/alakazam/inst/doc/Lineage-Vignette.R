## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Load required packages
library(alakazam)
library(igraph)
library(dplyr)

# Select a clone from the example database
data(ExampleDb)
sub_db <- subset(ExampleDb, clone_id == 3138)

## ---- eval=TRUE---------------------------------------------------------------
# This example data set does not have ragged ends
# Preprocess clone without ragged end masking (default)
clone <- makeChangeoClone(sub_db, text_fields=c("sample_id", "c_call"), 
                          num_fields="duplicate_count")

# Show combined annotations
clone@data[, c("sample_id", "c_call", "duplicate_count")]

## ---- eval=FALSE--------------------------------------------------------------
#  # Run PHYLIP and parse output
#  phylip_exec <- "~/apps/phylip-3.69/dnapars"
#  graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)

## ---- echo=FALSE, warning=FALSE, message=FALSE--------------------------------
# Load data insted of running phylip
# Clone 3138 is at index 23
graph <- ExampleTrees[[23]]

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# The graph has shared annotations for the clone
data.frame(clone_id=graph$clone,
           junction_length=graph$junc_len,
           v_gene=graph$v_gene,
           j_gene=graph$j_gene)

# The vertices have sequence specific annotations
data.frame(sequence_id=V(graph)$name, 
           c_call=V(graph)$c_call,
           duplicate_count=V(graph)$duplicate_count)

## ---- eval=TRUE---------------------------------------------------------------
# Plot graph with defaults
plot(graph)

## ---- eval=TRUE---------------------------------------------------------------
# Modify graph and plot attributes
V(graph)$color <- "steelblue"
V(graph)$color[V(graph)$name == "Germline"] <- "black"
V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
V(graph)$label <- V(graph)$c_call
E(graph)$label <- ""

# Remove large default margins
par(mar=c(0, 0, 0, 0) + 0.1)
# Plot graph
plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
     vertex.label.color="black", vertex.size=40)
# Add legend
legend("topleft", c("Germline", "Inferred", "Sample"), 
       fill=c("black", "white", "steelblue"), cex=0.75)

## ---- eval=TRUE, warning=FALSE, results="hide"--------------------------------
# Preprocess clones
clones <- ExampleDb %>%
    group_by(clone_id) %>%
    do(CHANGEO=makeChangeoClone(., text_fields=c("sample_id", "c_call"), 
                                num_fields="duplicate_count"))

## ---- eval=FALSE--------------------------------------------------------------
#  # Build lineages
#  phylip_exec <- "~/apps/phylip-3.69/dnapars"
#  graphs <- lapply(clones$CHANGEO, buildPhylipLineage,
#                   phylip_exec=phylip_exec, rm_temp=TRUE)

## ---- echo=FALSE, warning=FALSE, message=FALSE--------------------------------
# Load data insted of running phylip
graphs <- ExampleTrees

## ---- eval=TRUE---------------------------------------------------------------
# Note, clones with only a single sequence will not be processed.
# A warning will be generated and NULL will be returned by buildPhylipLineage
# These entries may be removed for clarity
graphs[sapply(graphs, is.null)] <- NULL

# The set of tree may then be subset by node count for further 
# analysis, if desired.
graphs <- graphs[sapply(graphs, vcount) >= 5]

## ---- eval=TRUE, show=FALSE---------------------------------------------------
# Modify graph and plot attributes
V(graph)$color <- categorical_pal(8)[1]
V(graph)$label <- V(graph)$name
E(graph)$label <- E(graph)$weight

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
##plot lineage tree using igraph
plot(graph, layout=layout_as_tree)

# convert to phylo
phylo <- graphToPhylo(graph)

#plot using ape
plot(phylo, show.node.label=TRUE)

#write tree file in Newick format
ape::write.tree(phylo, file="example.tree")

## ---- eval=TRUE---------------------------------------------------------------
#read in tree as phylo object
phylo_r <- ape::read.tree("example.tree")

#convert to graph object
graph_r <- phyloToGraph(phylo_r, germline="Germline")

#plot converted form using igraph - it's the same as before
plot(graph_r,layout=layout_as_tree)

