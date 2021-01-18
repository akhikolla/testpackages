## ----pkg-load, echo = FALSE, message = FALSE----------------------------------
library (dodgr)

## -----------------------------------------------------------------------------
graph <- weight_streetnet (hampi, wt_profile = "foot")
from <- sample (graph$from_id, size = 10)
to <- sample (graph$to_id, size = 10)
flows <- matrix (10 * runif (length (from) * length (to)),
                 nrow = length (from))

## -----------------------------------------------------------------------------
graph_f <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
head (graph_f)

## -----------------------------------------------------------------------------
summary (graph_f$flow)

## -----------------------------------------------------------------------------
dens <- rep (1, length (from)) # uniform densities
graph_f <- dodgr_flows_disperse (graph, from = from, dens = dens)
summary (graph_f$flow)

## -----------------------------------------------------------------------------
graph_undir <- merge_directed_graph (graph_f)

## -----------------------------------------------------------------------------
nrow (graph_f); nrow (graph_undir) # the latter is much smaller

## -----------------------------------------------------------------------------
graph <- graph [graph_undir$edge_id, ]
graph$flow <- graph_undir$flow

## ----flowmap, eval = FALSE----------------------------------------------------
#  graph_f <- graph_f [graph_f$flow > 0, ]
#  dodgr_flowmap (graph_f, linescale = 5)

