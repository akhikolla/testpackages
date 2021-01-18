# Initialization ----------------------------------------------------------

initializePopulation  = function(N, method="uniform", ...) {
  out = switch (method,
                "uniform" = .uniformSquare(N=N, ...),
                .uniformSquare(N=N, ...)
                )
  return(out)
}


# Initialization methods --------------------------------------------------

.uniformSquare = function(N, n, min=0, max=1, maxpop=1e6, ...) {
  out = matrix(NA, nrow = min(5*N, maxpop), ncol = n)
  out[seq_len(N), ] = runif(n*N, min=min, max=max) 
  return(out)
}
