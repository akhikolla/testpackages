## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
library(stplanr)
library(sf)

## ---- out.width="40%", fig.show='hold', fig.width=5, message=FALSE------------
library(stplanr)
library(sf)
sample_routes <- routes_fast_sf[2:6, 1]
sample_routes$value <- rep(1:3, length.out = 5)
rnet <- overline(sample_routes, attrib = "value")
plot(sample_routes["value"], lwd = sample_routes$value, main = "Routes")
plot(rnet["value"], lwd = rnet$value, main = "Route network")

## ----rnets1, message=FALSE, warning=FALSE, out.width="100%", fig.width=6, fig.height=6, echo=FALSE----
# knitr::include_graphics("route-networks.png")

## -----------------------------------------------------------------------------
touching_list = st_intersects(sample_routes)
g = igraph::graph.adjlist(touching_list)
igraph::is_connected(g)

## ---- eval=FALSE--------------------------------------------------------------
#  # piggyback::pb_download_url("r_key_roads_test.Rds")
#  u = "https://github.com/ropensci/stplanr/releases/download/0.6.0/r_key_roads_test.Rds"
#  rnet_disconnected = readRDS(url(u))
#  touching_list = sf::st_intersects(rnet_disconnected)
#  g = igraph::graph.adjlist(touching_list)
#  igraph::is_connected(g)
#  #> [1] FALSE
#  sf:::plot.sfc_LINESTRING(rnet_disconnected$geometry)

## ---- eval=FALSE--------------------------------------------------------------
#  rnet_disconnected$group = rnet_igroup(rnet_disconnected)

## ----rnet-routing1------------------------------------------------------------
sln <- SpatialLinesNetwork(rnet)
class(sln)

## -----------------------------------------------------------------------------
class(sln@sl)
nrow(sln@sl)
class(sln@g)
length(igraph::edge.attributes(sln@g)[["weight"]])
class(sln@nb)
length(unique(unlist(sln@nb)))
identical(sln@sl$geometry, rnet$geometry)

## -----------------------------------------------------------------------------
sln_nodes <- sln2points(sln)
nrow(sln_nodes)
length(sln@nb)

## -----------------------------------------------------------------------------
rnet_coordinates <- sf::st_coordinates(rnet)
set.seed(85)
x <- runif(n = 2, min = min(rnet_coordinates[, 1]), max = max(rnet_coordinates[, 1]))
y <- runif(n = 2, min = min(rnet_coordinates[, 2]), max = max(rnet_coordinates[, 2]))
crs <- sf::st_crs(rnet)
xy_sf <- sf::st_as_sf(data.frame(n = 1:2, x, y), coords = c("x", "y"), crs = crs)
xy_nodes <- stplanr::find_network_nodes(sln = sln, x = x, y = y)

## ---- out.width="49%", fig.show='hide'----------------------------------------
# plot(rnet$geometry)
# plot(sln_nodes, add = TRUE)
# xy_path <- sum_network_routes(sln = sln, start = xy_nodes[1], end = xy_nodes[2], sumvars = "length")
# # xy_path = sum_network_links(sln = sln, start = xy_nodes[1], end = xy_nodes[2])
# plot(rnet$geometry)
# plot(xy_sf$geometry, add = TRUE)
# plot(xy_path$geometry, add = TRUE, lwd = 5)

## ----netpoint-----------------------------------------------------------------
new_point_coordinates <- c(-1.540, 53.826)
p <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(new_point_coordinates)), crs = crs)

## ---- fig.show='hold', out.width="49%"----------------------------------------
sln_new <- sln_add_node(sln = sln, p = p)
route_new <- route_local(sln = sln_new, from = p, to = xy_sf[1, ])
plot(sln_new)
plot(p, add = TRUE)
plot(route_new, lwd = 5, add = TRUE)

## ---- message=FALSE, warning=FALSE, out.width="100%", fig.width=6, fig.height=6, echo=FALSE, eval=FALSE----
#  # Show solutions to https://github.com/ropensci/stplanr/issues/237
#  
#  library(tmap)
#  # sample_routes = routes_fast_sf
#  # sample_routes$id <- 1:nrow(sample_routes)
#  # sample_routes <- sample_routes[!is.na(sample_routes$plan),]
#  sample_routes <- routes_fast_sf[c(22, 38, 39, 46, 47), NULL]
#  sample_routes$type <- " Routes"
#  v <- 1:5
#  n <- c(2, 3, 5)
#  # route_list = purrr::map(n, ~sample_routes[1:., ])
#  route_list <- lapply(n, function(x) {
#    l <- sample_routes[1:x, ]
#    l$n <- x
#    l$value <- rep(v, length.out = x)
#    l
#  })
#  routes_all <- do.call(rbind, route_list)
#  # rnet_list = purrr::map(route_list, function(x) {
#  rnet_list <- lapply(route_list, function(x) {
#    l <- overline2(x, "value")
#    l$n <- mean(x$n)
#    l
#  })
#  rnet_all <- do.call(rbind, rnet_list)
#  rnet_all$type <- "Route network"
#  all_routes <- rbind(routes_all, rnet_all)
#  p <- sf::st_centroid(all_routes)
#  tmap_mode("plot")
#  m <- tm_shape(all_routes, bbox = tmaptools::bb(all_routes)) +
#    tm_lines(
#      col = "value", lwd = 1, palette = "OrRd", scale = 8, alpha = 0.8, breaks = 0:6,
#      legend.lwd.show = FALSE, labels = as.character(1:6),
#      legend.col.show = FALSE
#    ) +
#    tm_text("value") +
#    tm_facets(by = c("n", "type")) +
#    tm_layout(scale = 1.5)
#  # m
#  tmap_save(m, "vignettes/route-networks.png")

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  # test code:
#  
#  sample_routes2 <- sample_routes5[2:3, ]
#  sample_routes3 <- sample_routes5[2:4, ]
#  rnet2 <- overline2(sample_routes2, attrib = "value")
#  rnet3 <- overline2(sample_routes3, attrib = "value")
#  rnet5 <- overline2(sample_routes5, attrib = "value")
#  
#  b <- 0:6
#  bb <- tmaptools::bb(rnet, ext = 1.1)
#  
#  rnet5$n <- 5
#  rnet5$type <- "Route network"
#  sample_routes5$n <- 5
#  sample_routes5$type <- " Routes"
#  
#  all_routes <- rbind(rnet5, sample_routes5)
#  
#  m2 <- tm_shape(sample_routes[1:2, ], bbox = bb) +
#    tm_lines(col = "value", lwd = "value", palette = "magma", scale = 8, alpha = 0.5, breaks = b) +
#    tm_layout(title = "2 Routes", legend.show = FALSE)
#  r2 <- tm_shape(rnet2, bbox = bb) +
#    tm_lines(
#      col = "value", palette = "viridis", scale = 10, alpha = 0.5, breaks = b,
#      legend.lwd.show = FALSE, labels = as.character(1:6)
#    )
#  m3 <- tm_shape(sample_routes[1:3, ], bbox = bb) +
#    tm_lines(col = "value", lwd = "value", palette = "viridis", scale = 15, alpha = 0.5, breaks = b) +
#    tm_layout(title = "3 Routes", legend.show = FALSE)
#  r3 <- tm_shape(rnet3, bbox = bb) +
#    tm_lines(col = "value", palette = "viridis", scale = 8, alpha = 0.5, breaks = b) +
#    tm_layout(legend.show = FALSE)
#  m6 <- tm_shape(sample_routes, bbox = bb) +
#    tm_lines(col = "value", lwd = "value", palette = "viridis", scale = 10, alpha = 0.5, breaks = b) +
#    tm_layout(title = "6 Routes", legend.show = FALSE)
#  r6 <- tm_shape(rnet, bbox = bb) +
#    tm_lines(col = "value", palette = "viridis", scale = 8, alpha = 0.5, breaks = b) +
#    tm_layout(legend.show = FALSE)
#  
#  
#  
#  tmap_arrange(m2, r2, m3, r3, m6, r6, nrow = 3)

