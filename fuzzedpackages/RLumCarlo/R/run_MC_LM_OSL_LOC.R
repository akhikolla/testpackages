#' @title Run Monte-Carlo Simulation for LM-OSL (localized transitions)
#'
#' @description Runs a Monte-Carlo (MC) simulation of linearly modulated optically stimulated
#' luminescence (LM-OSL) using the generalized one trap (GOT) model. Localized transitions refer to
#' transitions which do not involve the conduction or valence band. These transitions take place
#' between the ground state and an excited state of the trap, and also involve
#' a an energy state of the recombination centre.
#'
#' @details
#'
#' **The model**
#'
#' \deqn{
#' I_{LOC}(t) = -dn/dt = (A * t/P) * (n^2 / (r + n))
#' }
#'
#' Where in the function: \cr
#'  A := optical excitation rate from the ground state into the excited state of the trap (1/s)\cr
#'  P := total excitation time (s) \cr
#'  t := time (s) \cr
#'  n := `n_filled`, the instantaneous number of electrons \cr
#'  r := the retrapping ratio for localized transitions
#'
#' @param A [numeric] (**required**): The optical excitation rate from the ground state into the excited
#' state of the trap (s^-1)
#'
#' @param times [numeric] (**required**): The sequence of time steps within the simulation (s)
#'
#' @param clusters [numeric] (*with default*): The number of created clusters for the MC runs. The input can be the output of [create_ClusterSystem]. In that case `n_filled` indicate absolute numbers of a system.
#'
#' @param n_filled [integer] (*with default*): The number of filled electron traps at the
#' beginning of the simulation (dimensionless). Can be a vector of `length(clusters)`, shorter values are recycled.
#'
#' @param r [numeric] (**required**): The retrapping ratio for localized transitions
#'
#' @param method [character] (*with default*): Sequential `'seq'` or parallel `'par'`processing. In
#' the parallel mode the function tries to run the simulation on multiple CPU cores (if available) with
#' a positive effect on the computation time.
#'
#' @param output [character] (*with default*): output is either the `'signal'` (the default)
#' or `'remaining_e'` (the remaining charges, electrons, in the trap)
#'
#' @param \dots further arguments, such as `cores` to control the number of used CPU cores or `verbose` to silence the terminal
#'
#' @return This function returns an object of class `RLumCarlo_Model_Output` which
#' is a [list] consisting of an [array] with dimension length(times) x clusters
#' and a [numeric] time vector.
#'
#' @section Function version: 0.1.0
#'
#' @author Sebastian Kreutzer, Geography & Earth Sciences, Aberystwyth University (United Kingdom)
#'
#' @references
#' Pagonis, V., Friedrich, J., Discher, M., Müller-Kirschbaum, A., Schlosser, V.,
#' Kreutzer, S., Chen, R. and Schmidt, C., 2019. Excited state luminescence signals
#' from a random distribution of defects: A new Monte Carlo simulation approach for feldspar.
#' Journal of Luminescence 207, 266–272. \doi{10.1016/j.jlumin.2018.11.024}
#'
#' @examples
#' ## short example
#' run_MC_LM_OSL_LOC(
#'  A = 1,
#'  times = 0:40,
#'  clusters = 10,
#'  n_filled = 100,
#'  r = 1e-7,
#'  method = "seq",
#'  output = "signal") %>%
#' plot_RLumCarlo(legend = TRUE)
#'
#' \dontrun{
#' ## the long (meaningful) example
#' results <- run_MC_LM_OSL_LOC(
#'  A = 1,
#'  times = 0:100,
#'  clusters = 100,
#'  n_filled = 100,
#'  r = 1e-7,
#'  method = "par",
#'  output = "signal")
#'
#' ## plot
#' plot_RLumCarlo(results, legend = TRUE)
#' }
#'
#' @keywords models data
#' @encoding UTF-8
#' @md
#' @export
run_MC_LM_OSL_LOC <- function(
  A,
  times,
  clusters = 10,
  n_filled = 100,
  r,
  method = "par",
  output = "signal",
  ...){

# Integrity checks ----------------------------------------------------------------------------
  if(!output %in% c("signal", "remaining_e"))
    stop("[run_MC_LM_OSL_LOC()] Allowed keywords for 'output' are either 'signal' or 'remaining_e'!", call. = FALSE)

# Register multi-core back end ----------------------------------------------------------------
cl <- .registerClusters(method, ...)
on.exit(parallel::stopCluster(cl))

# Enable dosimetric cluster system -----------------------------------------
if(class(clusters)[1] == "RLumCarlo_ClusterSystem"){
  n_filled <- .distribute_electrons(
    clusters = clusters,
    N_system = n_filled[1])[["e_in_cluster"]]
  clusters <- clusters$cl_groups

}

# Expand parameters -------------------------------------------------------
n_filled <- rep(n_filled, length.out = max(clusters))

# Run model -----------------------------------------------------------------------------------
  temp <- foreach(c = 1:max(clusters),
                  .packages = 'RLumCarlo',
                  .combine = 'comb_array',
                  .multicombine = TRUE) %dopar% {

    results <- MC_C_LM_OSL_LOC(
      times = times,
      n_filled = n_filled[c],
      r = r[1],
      A = A[1]
      )

    return(results[[output]])

  }  # end c-loop

# Return --------------------------------------------------------------------------------------
  .return_ModelOutput(signal = temp, time = times)
}
