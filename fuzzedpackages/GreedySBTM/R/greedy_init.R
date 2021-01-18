
GreedyInit <- function(adj_cube, Kup, list_of_inactive_nodes = NULL)
{
  tframes <- dim(adj_cube)[3]
  N <- dim(adj_cube)[2]
  aggregated <- adj_cube[,,1]
  for (t in 2:tframes) aggregated <- cbind(aggregated,adj_cube[,,t])
  kmeans_out <- kmeans(x = aggregated, centers = Kup)$cluster
  allocations_init <- matrix(kmeans_out,tframes,N,T)
  if (!is.null(list_of_inactive_nodes)) if (nrow(list_of_inactive_nodes) > 0) for (index in 1:nrow(list_of_inactive_nodes)) allocations_init[list_of_inactive_nodes[index,1],list_of_inactive_nodes[index,2]] = 0
  CollapseLabels(allocations_init)
}
