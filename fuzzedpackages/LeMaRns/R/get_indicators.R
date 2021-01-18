#' @include run_LeMans.R
#' @include param_setup.R
NULL

#' Calculate indicators
#'
#' @description Calculates Mean Maximum Length (MML), the Large Fish Indicator (LFI), Typical Length (TyL) and Length Quantiles (LQ) of the whole community or a subset of the species.
#' @param inputs A \code{\linkS4class{LeMans_param}} object containing the parameter values of the current LeMans model.
#' @param outputs A \code{\linkS4class{LeMans_outputs}} object containing the outputs of the model run.
#' @param wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class.
#' @param mid A numeric vector of length \code{nfish} representing the mid-point of the length classes in the model.
#' @param l_bound A numeric vector of length \code{nsc} representing the lower bounds of the length classes.
#' @param u_bound A numeric vector of length \code{nsc} representing the upper bounds of the length classes.
#' @param Linf A numeric vector of length \code{nfish} representing the asymptotic length of each species.
#' @param N An array with dimensions \code{nsc}, \code{nfish} and \code{tot_time} representing the number of individuals in each length class for each time step.
#' @param species A numeric value or vector or a character string representing the species that you wish to use to calculate the indicators. The default is \code{1:dim(N)[2]}.
#' @param time_steps A numeric vector of the time steps that you wish to use to calculate the indicators. The default is \code{1:dim(N)[3]}.
#' @param species_names A character vector of length \code{nfish} that denotes the names of the species in the model.
#' @param prob A numeric value or vector between 0 and 1 denoting the length quantiles to be calculated. The default is \code{0.5}.
#' @param length_LFI A numeric vector representing the thresholds to be used to calculate the LFI. The default value is \code{40}.
#' @param ... Additional arguments.
#' @details The LFI represents the proportion of biomass with a length larger than \code{length_LFI}. The MML is the biomass weighted mean of \code{Linf}:
#' @details \code{sum(biomass[species]*Linf[species])/sum(biomass[species])}
#' @details where \code{biomass} is a numeric vector of length \code{nfish} representing the biomass of each species. TyL is the biomass-weighted geometric mean length of the community:
#' @details \code{exp(sum(biomass_*log(mid))/sum(Bio_l))}
#' @details where \code{biomass_} is a numeric vector of length \code{nsc} representing the biomass of all the species in each length class. The LQ is the length at which the biomass exceeds a given proportion \code{prob} of the total biomass.
#' @return \code{get_indicators} returns a list object with names `LFI`, `MML`, `TYL` and `LQ`. If \code{length(length_LFI)>1}, `LFI` is a matrix with dimensions \code{length(time_steps)} by \code{length(length_LFI)} where the \code{i,j}th element represents the LFI using the \code{j}th \code{length_LFI} in the \code{i}th \code{time_steps}. If \code{length(length_LFI)==1}, the function will return a numeric vector of length \code{length(time_steps)}. `MML` is a numeric vector of length \code{time_steps} where each element is the MML for the species in \code{species}. `TYL` is a numeric vector of length \code{time_steps} where each element is the TyL for the species in \code{species}. If \code{length(prob)==1}, `LQ` is a matrix with dimensions \code{length(time_steps)} by \code{length(prob)} where the \code{i,j}th element is the LQ using the\code{j}th \code{prob} in the \code{i}th \code{time_steps}. If \code{length(prob)==1}, the function will return a numeric vector of length \code{length(time_steps)}.
#' @return If \code{length(length_LFI)==1}, \code{get_LFI} returns a matrix with dimensions \code{length(time_steps)} by \code{length(length_LFI)} where the \code{i,j}th element is the LFI using the \code{j}th \code{length_LFI} in the \code{i}th \code{time_steps}. If \code{length(length_LFI)==1}, the function will return a numeric vector of length \code{length(time_steps)}.
#' @return \code{get_MML} returns a numeric vector of length \code{time_steps} where each element is the MML for the species in \code{species}.
#' @return If \code{length(prob)>1}, \code{get_LQ} returns a matrix with dimensions \code{length(time_steps)} and \code{length(prob)} where the \code{i,j}th element is the LQ using the the \code{j}th \code{prob} in the \code{i}th \code{time_steps}. If \code{length(prob)==1}, the function will return a numeric vector of length \code{length(time_steps)}.
#' @seealso \code{\link{plot_indicators}}
#' @examples
#' # Set up and run the model
#' NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=1e12)
#' effort <- matrix(0.5, 10, dim(NS_params@Qs)[3])
#' model_run <- run_LeMans(NS_params, years=10, effort=effort)
#'
#' # Calculate the indicators
#' get_indicators(inputs=NS_params, outputs=model_run)
#'
#' # Calculate the LFI
#' get_LFI(inputs=NS_params, outputs=model_run)
#'
#' # Calculate MML
#' get_MML(inputs=NS_params, outputs=model_run)
#'
#' # Calculate TyL
#' get_TyL(inputs=NS_params, outputs=model_run)
#'
#' # Calculate LQs
#' get_LQ(inputs=NS_params, outputs=model_run)
#' @export
setGeneric('get_indicators', function(inputs, outputs, ...)
  standardGeneric('get_indicators'))

#' @rdname get_indicators
setMethod('get_indicators', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], prob=0.5, length_LFI=40) {
            return(get_indicators(wgt=inputs@wgt, mid=inputs@mid, l_bound=inputs@l_bound, u_bound=inputs@u_bound, Linf=inputs@Linf, N=outputs@N, species=species, time_steps=time_steps, species_names=inputs@species_names, prob=prob, length_LFI=length_LFI))
          })

#' @rdname get_indicators
setMethod('get_indicators', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], prob=0.5, length_LFI=40) {
            return(get_indicators(wgt=inputs@wgt, mid=inputs@mid, l_bound=inputs@l_bound, u_bound=inputs@u_bound, Linf=inputs@Linf, N=N, species=species, time_steps=time_steps, species_names=inputs@species_names, prob=prob, length_LFI=length_LFI))
          })

#' @rdname get_indicators
setMethod('get_indicators', signature(inputs="missing", outputs='LeMans_outputs'),
          function(wgt, mid, l_bound, u_bound, Linf, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], species_names=NULL, prob=0.5, length_LFI=40) {
            return(get_indicators(wgt=wgt, mid=mid, l_bound=l_bound, u_bound=u_bound, Linf=Linf, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, prob=prob, length_LFI=length_LFI))
          })

#' @rdname get_indicators
setMethod('get_indicators', signature(inputs="missing", outputs='missing'),
          function(wgt, mid, l_bound, u_bound, Linf, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], species_names=NULL, prob=0.5, length_LFI=40) {
            ret <- list()
            ret[["LFI"]] <- get_LFI(wgt=wgt, l_bound=l_bound, u_bound=u_bound, N=N, species=species, time_steps=time_steps, species_names=species_names, length_LFI=length_LFI)
            ret[["MML"]] <- get_MML(wgt=wgt, Linf=Linf, N=N, species=species, time_steps=time_steps, species_names=species_names)
            ret[["TYL"]] <- get_TyL(wgt=wgt, mid=mid, N=N, species=species, time_steps=time_steps, species_names=species_names)
            ret[["LQ"]] <- get_LQ(wgt=wgt, u_bound=u_bound, N=N, species=species, time_steps=time_steps, species_names=species_names, prob=prob)
            return(ret)
          })

#' @export
#' @rdname get_indicators
setGeneric('get_LFI', function(inputs, outputs, ...)
  standardGeneric('get_LFI'))

#' @rdname get_indicators
setMethod('get_LFI', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], length_LFI=40) {
            return(get_LFI(wgt=inputs@wgt, l_bound=inputs@l_bound, u_bound=inputs@u_bound, N=outputs@N,species=species, time_steps=time_steps, species_names=inputs@species_names, length_LFI=length_LFI))
          })

#' @rdname get_indicators
setMethod('get_LFI', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], length_LFI=40) {
            return(get_LFI(wgt=inputs@wgt, l_bound=inputs@l_bound, u_bound=inputs@u_bound, N=N, species=species, time_steps=time_steps, species_names=inputs@species_names, length_LFI=length_LFI))
          })

#' @rdname get_indicators
setMethod('get_LFI', signature(inputs="missing", outputs='LeMans_outputs'),
          function(wgt, l_bound, u_bound, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], species_names=NULL, length_LFI=40) {
            return(get_LFI(wgt=wgt, l_bound=l_bound, u_bound=u_bound, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, length_LFI=length_LFI))
          })

#' @rdname get_indicators
setMethod('get_LFI', signature(inputs="missing", outputs='missing'),
          function(wgt, l_bound, u_bound, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], species_names=NULL, length_LFI=40) {
            if (length(length_LFI)>1){
              return(sapply(length_LFI, function(wgt, l_bound, u_bound, N, species, time_steps, species_names, length_LFI){
                get_LFI(wgt=wgt, l_bound=l_bound, u_bound=u_bound, N=N, species=species, time_steps=time_steps, species_names=species_names, length_LFI=length_LFI)},
                wgt=wgt, l_bound=l_bound, u_bound=u_bound, N=N, species=species, time_steps=time_steps, species_names=species_names))
            }
            if (is.character(species)){
              species <- sapply(species, function(species){which(species==species_names)})
            }
            bry <- which(l_bound<=length_LFI & u_bound>length_LFI)
            prop <- (length_LFI-l_bound[bry])/(u_bound[bry]-l_bound[bry])
            if (length(time_steps)>1){
              new_n <- N[, species, time_steps]
              return(apply(new_n, length(dim(new_n)), calc_LFI, wgt=wgt[,species], bry=bry, prop=prop))
            }
            new_n <- N[, species, time_steps]
            return(calc_LFI(new_n, wgt=wgt[, species], bry=bry, prop=prop))
          })

#' @export
#' @rdname get_indicators
setGeneric('get_MML', function(inputs, outputs, ...)
  standardGeneric('get_MML'))

#' @rdname get_indicators
setMethod('get_MML', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs,outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3]) {
            return(get_MML(wgt=inputs@wgt, Linf=inputs@Linf, N=outputs@N, species=species, time_steps=time_steps, species_names=inputs@species_names))
          })

#' @rdname get_indicators
setMethod('get_MML', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs,N, species=1:dim(N)[2], time_steps=1:dim(N)[3]) {
            return(get_MML(wgt=inputs@wgt, Linf=inputs@Linf, N=N, species=species, time_steps=time_steps, species_names=inputs@species_names))
          })

#' @rdname get_indicators
setMethod('get_MML', signature(inputs="missing", outputs='LeMans_outputs'),
          function(wgt, Linf, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], species_names=NULL) {
            return(get_MML(wgt=wgt, Linf=Linf, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names))
          })

#' @rdname get_indicators
setMethod('get_MML', signature(inputs="missing", outputs='missing'),
          function(wgt, Linf, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], species_names=NULL) {
            if (is.character(species)){
              species <- sapply(species, function(species){which(species==species_names)})
            }
            if (length(species)==1){
              return(rep(Linf[species], length(time_steps)))
            }
            if (length(time_steps)>1){
              new_n <- N[, species, time_steps]
              biomass <- apply(new_n, length(dim(new_n)), calc_biomass, wgt=wgt[, species])
              return(apply(biomass, 2, weighted.mean, x=Linf[species]))
            }
            new_n <- N[, species, time_steps]
            return(weighted.mean(Linf[species], w=calc_biomass(new_n, wgt=wgt[, species])))
          })

#' @export
#' @rdname get_indicators
setGeneric('get_TyL', function(inputs, outputs, ...)
  standardGeneric('get_TyL'))

#' @rdname get_indicators
setMethod('get_TyL', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3]) {
            return(get_TyL(wgt=inputs@wgt, mid=inputs@mid, N=outputs@N, species=species, time_steps=time_steps, species_names=inputs@species_names))
          })

#' @rdname get_indicators
setMethod('get_TyL', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs, N, species=1:dim(N)[2], time_steps=1:dim(N)[3]) {
            return(get_TyL(wgt=inputs@wgt, mid=inputs@mid, N=N, species=species, time_steps=time_steps, species_names=inputs@species_names))
          })

#' @rdname get_indicators
setMethod('get_TyL', signature(inputs="missing", outputs='LeMans_outputs'),
          function(wgt, mid, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], species_names=NULL) {
            return(get_TyL(wgt=wgt, mid=mid, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names))
          })

#' @rdname get_indicators
setMethod('get_TyL', signature(inputs="missing", outputs='missing'),
          function(wgt, mid, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], species_names=NULL) {
            if(is.character(species)){
              species <- sapply(species, function(species){which(species==species_names)})
            }
            if (length(time_steps)>1){
              new_n <- N[, species, time_steps]
              return(apply(new_n, length(dim(new_n)), calc_TyL, wgt=wgt[, species], mid=mid))
            }
            new_n <- N[, species, time_steps]
            return(calc_TyL(new_n, wgt=wgt[, species], mid=mid))
          })

#' @export
#' @rdname get_indicators
setGeneric('get_LQ', function(inputs, outputs, ...)
  standardGeneric('get_LQ'))

#' @rdname get_indicators
setMethod('get_LQ', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], prob=0.5) {
            return(get_LQ(wgt=inputs@wgt, u_bound=inputs@u_bound, N=outputs@N, species=species, time_steps=time_steps, species_names=inputs@species_names, prob=prob))
          })

#' @rdname get_indicators
setMethod('get_LQ', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], prob=0.5) {
            return(get_LQ(wgt=inputs@wgt, u_bound=inputs@u_bound, N=N, species=species, time_steps=time_steps, species_names=inputs@species_names, prob=prob))
          })

#' @rdname get_indicators
setMethod('get_LQ', signature(inputs="missing", outputs='LeMans_outputs'),
          function(wgt, u_bound, outputs, species=1:dim(outputs@N)[2], time_steps=1:dim(outputs@N)[3], species_names=NULL, prob=0.5) {
            return(get_LQ(wgt=wgt, u_bound=u_bound, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, prob=prob))
          })

#' @rdname get_indicators
setMethod('get_LQ', signature(inputs="missing", outputs='missing'),
          function(wgt, u_bound, N, species=1:dim(N)[2], time_steps=1:dim(N)[3], species_names=NULL, prob=0.5) {
            if (length(prob)>1){
              return(sapply(prob, function(wgt, u_bound, N, species, time_steps, species_names, prob)
                get_LQ(wgt=wgt, u_bound=u_bound, N=N, species=species, time_steps=time_steps, species_names=species_names, prob=prob),
                wgt=wgt, u_bound=u_bound, N=N, species=species, time_steps=time_steps, species_names=species_names))
            }
            if (is.character(species)){
              species <- sapply(species, function(species){which(species==species_names)})
            }
            if (length(time_steps)>1){
              new_n <- N[,species,time_steps]
              return(apply(new_n, length(dim(new_n)), calc_LQ, wgt=wgt[, species], u_bound=u_bound, prob=prob))
            }
            new_n <- N[,species,time_steps]
            return(calc_LQ(new_n, wgt=wgt[,species], u_bound=u_bound, prob=prob))
          })
