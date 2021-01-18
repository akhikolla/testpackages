
CollapseLabels <- function(allocations) 
{
  tframes <- nrow(allocations)
  N <- ncol(allocations)
  index <- 1
  vec <- rep(NA,tframes*N)
  for (t in 1:tframes) for (i in 1:N) if (allocations[t,i] != 0)
  {
    vec[index] = allocations[t,i]
    index = index + 1
  }
  vec_out <- as.numeric(cpp_CollapseLabels(vec)$vec) + 1
  allocations_out <- allocations
  index <- 1
  for (t in 1:tframes) for (i in 1:N) if (allocations[t,i] != 0)
  {
    allocations_out[t,i] = vec_out[index]
    index = index + 1
  }
  allocations_out
}

