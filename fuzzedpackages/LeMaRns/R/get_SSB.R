#' @include run_LeMans.R
#' @include param_setup.R
NULL

#' Calculate Spawning Stock Biomass (SSB) or total biomass
#'
#' @description Calculates the Spawning Stock Biomass (SSB) or the total biomass of each species in the model.
#' @param inputs A \linkS4class{LeMans_param} object containing the parameter values of the current LeMans model.
#' @param outputs A \linkS4class{LeMans_outputs} object containing the outputs of the model run.
#' @param wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class.
#' @param mature A matrix with dimensions \code{ncs} and \code{nfish} with elements in the range 0-1 representing the proportion of mature individuals of each species in each length class.
#' @param N An array with dimensions \code{nsc}, \code{nfish} and \code{tot_time} representing the number of individuals in each length class for each time step.
#' @param species A numeric value or vector or a character string or vector of the species that you wish to calculate the mean maximum length. The default is \code{1:dim(N)[2]}.
#' @param time_steps A numeric vector of the time steps that you wish to calculate the mean maximum length. The default is \code{1:dim(N)[3]}.
#' @param species_names A character vector of length \code{nfish} that denotes the names of the species in the model.
#' @param ... Additional arguments.
#' @details The SSB for each species in \code{species} is equal to:
#' @details \code{colSums(N[, species]*wgt[, species]*mature[, species])}.
#' @details The biomass for each species in \code{species} is equal to:
#' @details \code{colSums(N[, species]*wgt[, species])}
#' @return If \code{length(species)>1}, \code{get_SSB} returns a matrix with dimensions \code{length(time_steps)} by \code{length(species)} where the \code{i,j}th element is the SSB (g) of the \code{j}th \code{species} in the \code{i}th \code{time_steps}. If \code{length(species)==1}, the function will return a numeric vector of length \code{length(time_steps)}.
#' @return If \code{length(species)>1}, \code{get_biomass} returns a matrix with dimensions \code{length(time_steps)} by \code{length(species)} where the \code{i,j}th element is the biomass (g) of the \code{j}th \code{species} in the \code{i}th \code{time_steps}. If \code{length(species)==1}, the function will return a numeric vector of length \code{length(time_steps)}.
#' @examples
#' # Set up and run the model
#' NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=1e12)
#' effort <- matrix(0.5, 10, dim(NS_params@Qs)[3])
#' model_run <- run_LeMans(NS_params, years=10, effort=effort)
#'
#' # Calculate SSB
#' get_SSB(inputs=NS_params, outputs=model_run)
#'
#' # Calculate biomass
#' get_biomass(inputs=NS_params, outputs=model_run)
#' @export
setGeneric('get_SSB', function(inputs, outputs, ...)
  standardGeneric('get_SSB'))

#' @rdname get_SSB
setMethod('get_SSB', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3]) {
            return(get_SSB(N=outputs@N, wgt=inputs@wgt, mature=inputs@mature, species=species, time_steps=time_steps, species_names=inputs@species_names))
          })

#' @rdname get_SSB
setMethod('get_SSB', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs, N, species=1:dim(N)[2], time_steps=1:dim(N)[3]) {
            return(get_SSB(N=N, wgt=inputs@wgt, mature=inputs@mature, species=species, time_steps=time_steps, species_names=inputs@species_names))
          })

#' @rdname get_SSB
setMethod('get_SSB', signature(inputs="missing", outputs='LeMans_outputs'),
          function(wgt, mature, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], species_names=NULL) {
            return(get_SSB(N=outputs@N, wgt=wgt, mature=mature, species=species, time_steps=time_steps, species_names=species_names))
          })

#' @rdname get_SSB
setMethod('get_SSB', signature(inputs="missing", outputs='missing'),
          function(wgt, mature, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], species_names=NULL) {
            if (is.character(species)) {
              species <- sapply(species, function(species){which(species==species_names)})}
            if (length(time_steps)>1) {
              new_n <- N[, species, time_steps]
              if (length(species)==1) {
                return(apply(new_n, length(dim(new_n)), calc_SSB, wgt=wgt[, species], mature=mature[, species]))
              }
              return(t(apply(new_n, length(dim(new_n)), calc_SSB, wgt=wgt[, species], mature=mature[, species])))
            }
            new_n <- N[, species, time_steps]
            return(calc_SSB(new_n, wgt=wgt[, species], mature=mature[, species]))
          })

#' @rdname get_SSB
#' @export
setGeneric('get_biomass', function(inputs, outputs, ...)
  standardGeneric('get_biomass'))

#' @rdname get_SSB
setMethod('get_biomass', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3]) {
            return(get_biomass(wgt=inputs@wgt, N=outputs@N, species=species, time_steps=time_steps, species_names=inputs@species_names))
          })

#' @rdname get_SSB
setMethod('get_biomass', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs, N, species=1:dim(N)[2], time_steps=1:dim(N)[3]) {
            return(get_biomass(wgt=inputs@wgt, N=N, species=species, time_steps=time_steps, species_names=inputs@species_names))
          })

#' @rdname get_SSB
setMethod('get_biomass', signature(inputs="missing", outputs='LeMans_outputs'),
          function(wgt, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], species_names=NULL) {
            return(get_biomass(wgt=wgt, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names))
          })

#' @rdname get_SSB
setMethod('get_biomass', signature(inputs="missing", outputs='missing'),
          function(wgt, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], species_names=NULL) {
            if (is.character(species)) {
              species <- sapply(species, function(species){which(species==species_names)})
            }
            if (length(time_steps)>1) {
              new_n <- N[, species, time_steps]
              if (length(species)==1) {
                return(apply(new_n, length(dim(new_n)), calc_biomass, wgt=wgt[, species]))
              }
              return(t(apply(new_n, length(dim(new_n)), calc_biomass, wgt=wgt[, species])))
            }
            new_n <- N[, species, time_steps]
            return(calc_biomass(new_n, wgt=wgt[, species]))
          })
