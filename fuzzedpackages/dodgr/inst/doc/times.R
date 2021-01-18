## ----pkg-load, echo = FALSE, message = FALSE----------------------------------
library (dodgr)

## ----streetnet_times, eval = FALSE--------------------------------------------
#  dat_sc <- dodgr_streetnet_sc ("ogbomosho nigeria")
#  graph <- weight_streetnet (dat_sc, wt_profile = "bicycle")
#  graph_t <- weight_streetnet (dat_sc, wt_profile = "bicycle", turn_penalty = TRUE)
#  nrow (graph); nrow (graph_t)

## ----streetnet_times-out1, echo = FALSE---------------------------------------
c (164168, 173160)

## ----streetnet-vertices, eval = FALSE-----------------------------------------
#  graphc <- dodgr_contract_graph (graph) # not graph_t!
#  v <- dodgr_vertices (graphc)
#  n <- 100 # number of desired vertices
#  from <- sample (v$id, size = n)
#  to <- sample (v$id, size = n)

## ----routing, eval = FALSE----------------------------------------------------
#  graph_tc <- dodgr_contract_graph (graph_t)
#  nrow (graph_tc); nrow (graph_t)

## ----routing-out, echo = FALSE------------------------------------------------
c (35808, 176160)

## ----weighting_profiles, eval = TRUE------------------------------------------
lapply (dodgr::weighting_profiles, function (i) i [i$name == "bicycle", ])

## ----shortest-vs-fastest------------------------------------------------------
graph <- weight_streetnet (hampi, wt_profile = "foot")
n <- 100 # number of sample routing vertices
from <- sample (graph$from_id, size = n)
to <- sample (graph$from_id, size = n)
d_dist <- dodgr_dists (graph, from = from, to = to, shortest = TRUE) # default
d_time <- dodgr_dists (graph, from = from, to = to, shortest = FALSE) # fastest paths
plot (d_dist / 1000, d_time / 1000, col = "orange",
      xlab = "distances along shortest paths (km)",
      ylab = "distances along fastest paths (km)")
lines (0:100, 0:100, col = "red", lty = 2)

## ----shortest-vs-fastest2-----------------------------------------------------
mean (abs (d_time - d_dist), na.rm = TRUE)

## ----fastest-is-shorter-------------------------------------------------------
index <- which (!is.na (d_time) & !is.na (d_dist))
length (which (d_time [index] < d_dist [index])) / length (index)

## ----shortest-vs-fastest-times------------------------------------------------
t_dist <- dodgr_times (graph, from = from, to = to, shortest = TRUE) # default
t_time <- dodgr_times (graph, from = from, to = to, shortest = FALSE) # fastest paths
plot (t_dist / 3600, t_time / 3600, col = "orange",
      xlab = "times along shortest paths (hours)",
      ylab = "times along fastest paths (hours)")
lines (0:100, 0:100, col = "red", lty = 2)
mean (abs (t_time - t_dist), na.rm = TRUE)

## ----fastest-is-slower--------------------------------------------------------
index <- which (!is.na (t_time) & !is.na (t_dist))
length (which (t_dist [index] < t_time [index])) / length (index)

## ----colswap, eval = FALSE----------------------------------------------------
#  graph$d_weighted <- graph$time_weighted

