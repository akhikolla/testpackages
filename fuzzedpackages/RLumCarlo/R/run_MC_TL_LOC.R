#' @title Run Monte-Carlo Simulation for TL (localized transitions)
#'
#' @description Runs a Monte-Carlo (MC) simulation of thermoluminescence (TL) using
#' the generalized one trap (GOT) model. Localized transitions refer to transitions
#' which do not involve the conduction or valence band. These transitions take place between the
#' ground state and an excited state of the trapped charge, and also involve an energy
#' state of the recombination centre. The heating rate in this function is assumed to be 1 K/s.
#'
#' @details
#'
#' **The model**
#'
#' \deqn{
#' I_{LOC}(t) = -dn/dt = (s * exp(-E/(k_{B} * T))) * (n^2 / (r + n))
#' }
#'
#'Where in the function: \cr
#' E := the thermal activation energy (eV) \cr
#' s := the frequency factor for the trap (s^-1) \cr
#' t := time (s) \cr
#' \eqn{k_{B}} := Boltzmann constant (8.617 x 10^-5 eV K^-1)\cr
#' T := temperature (°C) \cr
#' n := the instantaneous number of electrons \cr
#' r := the retrapping ratio for localized transitions
#'
#' @param E [numeric] (**required**): Thermal activation energy of the trap (eV)
#'
#' @param s [numeric] (**required**): The frequency factor of the trap (s^-1)
#'
#' @param times [numeric] (**required**): The sequence of temperature steps within the simulation (s).
#' The default heating rate is set to 1 K/s. The final temperature is `max(times) * b`
#'
#' @param b [numeric] (*with default*): the heating rate in K/s
#'
#' @param clusters [numeric] (*with default*): The number of created clusters for the MC runs. The input can be the output of [create_ClusterSystem]. In that case `n_filled` indicate absolute numbers of a system.
#'
#' @param n_filled [integer] (*with default*): The number of filled electron traps at
#' the beginning of the simulation (dimensionless). Can be a vector of `length(clusters)`, shorter values are recycled.
#'
#' @param r [numeric] (**required**): The localized retrapping ratio (dimensionless)
#'
#' @param method [character] (*with default*): Sequential `'seq'` or parallel `'par'`processing. In
#' the parallel mode the function tries to run the simulation on multiple CPU cores (if available) with
#' a positive effect on the computation time.
#'
#' @param output [character] (*with default*): output is either the `'signal'` (the default) or
#' `'remaining_e'` (the remaining charges/electrons in the trap)
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
#' ## the short example
#' run_MC_TL_LOC(
#'  s = 1e14,
#'  E = 0.9,
#'  times = 50:100,
#'  b = 1,
#'  method = "seq",
#'  clusters = 30,
#'  r = 1) %>%
#' plot_RLumCarlo()
#'
#' \dontrun{
#' ## the long (meaningful) example
#' results <- run_MC_TL_LOC(
#'  s = 1e14,
#'  E = 0.9,
#'  times = 50:100,
#'  method = "par",
#'  clusters = 100,
#'  r = 1)
#'
#' ## plot
#' plot_RLumCarlo(results)
#'
#' }
#'
#' @keywords models data
#' @encoding UTF-8
#' @md
#' @export
run_MC_TL_LOC <- function(
  s,
  E,
  times,
  b = 1,
  clusters = 10,
  n_filled = 100,
  r,
  method = "par",
  output = "signal",
  ...){

# Integrity checks ----------------------------------------------------------------------------
  if(!output %in% c("signal", "remaining_e"))
    stop("[run_MC_TL_LOC()] Allowed keywords for 'output' are either 'signal' or 'remaining_e'!", call. = FALSE)

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

    results <- MC_C_TL_LOC(
      times = times,
      b = b[1],
      n_filled = n_filled[c],
      r = r[1],
      E = E[1],
      s = s[1])

    return(results[[output]])

  }  # end c-loop

# Return --------------------------------------------------------------------------------------
if (output == "signal") temp <- temp / b
.return_ModelOutput(time = times * b, signal = temp)
}
