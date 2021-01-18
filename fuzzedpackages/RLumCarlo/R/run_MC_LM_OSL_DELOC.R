#' @title Run Monte-Carlo Simulation for LM-OSL (delocalized transitions)
#'
#' @description Runs a Monte-Carlo (MC) simulation of linearly modulated optically stimulated
#' luminescence (LM-OSL) using the one trap one recombination centre (OTOR) model.
#' Delocalised refers to involvement of the conduction band.
#'
#' @details
#'
#' **The model**
#'
#' \deqn{
#' I_{DELOC}(t) = -dn/dt = A * t/P * (n^2 / (N*R + n(1-R)))
#' }
#'
#' Where in the function: \cr
#'  t := time (s) \cr
#'  A := the optical excitation rate from trap to conduction band (1/s)\cr
#'  n := `n_filled`, the instantaneous number of electrons \cr
#'  R :=  the retrapping ratio for delocalized transitions \cr
#'  N := `N_e`, the total number of electron traps available (dimensionless) \cr
#'  P := total stimulation time (s)
#'
#' @param A [numeric] (**required**): The optical excitation rate from trap to conduction band (s^-1)
#'
#' @param times [numeric] (**required**): The sequence of time steps within the simulation (s)
#'
#' @param clusters [numeric] (*with default*): The number of created clusters for the MC runs. The input can be the output of [create_ClusterSystem]. In that case `n_filled` indicate absolute numbers of a system.
#'
#' @param N_e [integer] (*with default*): The total number of electron traps available (dimensionless). Can be a vector of `length(clusters)`, shorter values are recycled.
#'
#' @param n_filled [integer] (*with default*): The number of filled electron traps at the beginning
#' of the simulation (dimensionless). Can be a vector of `length(clusters)`, shorter values are recycled.
#'
#' @param R [numeric] (**required**): The retrapping ratio for delocalized transitions
#'
#' @param method [character] (*with default*): Sequential `'seq'` or parallel `'par'`processing. In
#' the parallel mode the function tries to run the simulation on multiple CPU cores (if available) with
#' a positive effect on the computation time.
#'
#' @param output [character] (*with default*): output is either the `'signal'` (the default)
#' or `'remaining_e'` (the remaining charges/electrons in the trap)
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
#' Pagonis, V., Friedrich, J., Discher, M., Müller-Kirschbaum, A., Schlosser, V., Kreutzer, S.,
#' Chen, R. and Schmidt, C., 2019. Excited state luminescence signals from a random distribution of
#' defects: A new Monte Carlo simulation approach for feldspar.
#' Journal of Luminescence 207, 266–272. \doi{10.1016/j.jlumin.2018.11.024}
#'
#' **Further reading**
#'
#' Chen, R., McKeever, S.W.S., 1997. Theory of Thermoluminescence and Related Phenomena.
#' WORLD SCIENTIFIC. \doi{10.1142/2781}
#'
#' @examples
#' run_MC_LM_OSL_DELOC(
#'  A = 0.12,
#'  R = 0.1,
#'  times = 0:50,
#'  method = "seq",
#'  clusters = 10) %>%
#' plot_RLumCarlo(legend = TRUE)
#'
#' @keywords models data
#' @encoding UTF-8
#' @md
#' @export
run_MC_LM_OSL_DELOC <- function(
  A,
  times,
  clusters = 10,
  N_e = 200,
  n_filled = N_e,
  R,
  method = "par",
  output = "signal",
  ...){

# Integrity checks ----------------------------------------------------------------------------
  if(!output %in% c("signal", "remaining_e"))
    stop("[run_MC_LM_OSL_DELOC()] Allowed keywords for 'output' are either 'signal' or 'remaining_e'!", call. = FALSE)

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
N_e <- rep(N_e, length.out = max(clusters))


# Run model -----------------------------------------------------------------------------------
  temp <- foreach(c = 1:max(clusters),
                  .packages = 'RLumCarlo',
                  .combine = 'comb_array',
                  .multicombine = TRUE) %dopar% {

    results <- MC_C_LM_OSL_DELOC(
      times = times,
      N_e = N_e[c],
      n_filled = n_filled[c],
      R = R[1],
      A = A[1])

    return(results[[output]])

  }  # end c-loop

# Return --------------------------------------------------------------------------------------
.return_ModelOutput(signal = temp, time = times)
}

