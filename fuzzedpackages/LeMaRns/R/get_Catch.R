#' @include run_LeMans.R
#' @include param_setup.R
NULL

#' Get annual catch for each species, catch per unit effort or catch per gear
#'
#' @description Get annual catch for each species, the Catch Per Unit Effort (CPUE) or the Catch Per Gear (CPG)
#' @param inputs A \code{\linkS4class{LeMans_param}} object containing the parameter values of the current LeMans model.
#' @param outputs A \code{\linkS4class{LeMans_outputs}} object containing the outputs of the model run.
#' @param Catch An array with dimensions \code{nsc}, \code{nfish} and \code{tot_time} representing the number of individuals in each length class for each time step.
#' @param years A numeric value representing the number of years included in \code{Catch}.
#' @param phi_min A numeric value representing the time step of the model.
#' @param inc_first A logical statement indicating whether the first time step of \code{Catch} should be included. The default is \code{FALSE}.
#' @param Qs An array of dimensions \code{nsc}, \code{nfish} and \code{gear} representing the catchability of each species by each of the fishing gears.
#' @param effort A matrix with dimensions \code{years} and the number of fishing gears, representing fishing effort in each year for each gear.
#' @param ... Additional arguments.
#' @return \code{get_annual_catch} returns a matrix with dimensions \code{years} and \code{length(species)} where the \code{i,j}th element represents the total catch (g) of the \code{j}th \code{species} in the \code{i}th year.
#' @return \code{get_CPUE} returns a matrix with dimensions \code{tot_time} and \code{nfish} where the \code{i,j}th element represents the CPUE in the \code{i}th time step for the \code{j}th species.
#' @return \code{get_CPG} returns an array with dimensions \code{nfish}, \code{dim(Qs[3])} and the number of time steps, where the \code{i,j,k}th element denotes the total catch of the \code{i}th species by the \code{j}th gear in the \code{k}th time step.
#' @examples
#' # Set up and run the model
#' NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=1e12)
#' effort <- matrix(0.5, 10, dim(NS_params@Qs)[3])
#' model_run <- run_LeMans(NS_params, years=10, effort=effort)
#'
#' # Get annual catch
#' get_annual_catch(inputs=NS_params, outputs=model_run)
#'
#' # Calculate the CPUE
#' get_CPUE(inputs=NS_params, outputs=model_run, effort=effort)
#'
#' # Calculate Catch Per Gear (CPG)
#' get_CPG(inputs=NS_params, outputs=model_run, effort=effort)
#' @export
setGeneric("get_annual_catch", function(inputs, outputs, ...)
  standardGeneric('get_annual_catch'))

#' @rdname get_annual_catch
setMethod("get_annual_catch", signature(inputs="missing", outputs="missing"),
          function(Catch, years=(dim(Catch)[3]+(inc_first-1))*phi_min, phi_min=0.1, inc_first=FALSE) {
            if (years<=0) {
              stop("years must be positive")
            }
            if ((1/phi_min)%%1!=0) {
              warning("The time step is not a multiple of a year")
            }
            year <- ceiling(((1:(years/phi_min))*phi_min))
            if (inc_first==FALSE) {
              Catch <- Catch[,,-1]
            }
            tmp <- apply(Catch, 3, colSums)
            return(t(sapply(sort(unique(year)), function(x, tmp, year){
              tele <- which(year==x)
              if (length(tele)>1) {
                return(rowSums(tmp[, tele]))
              }
              return(tmp[, tele])
            }, tmp=tmp, year=year)))
          })

#' @rdname get_annual_catch
setMethod("get_annual_catch", signature(inputs="LeMans_param", outputs="missing"),
          function(inputs, Catch, years=(dim(Catch)[3]+(inc_first-1))*inputs@phi_min, inc_first=FALSE)
          {
            return(get_annual_catch(Catch=Catch, years=years, phi_min=inputs@phi_min, inc_first=inc_first))
          })

#' @rdname get_annual_catch
setMethod("get_annual_catch", signature(inputs="LeMans_param", outputs="LeMans_outputs"),
          function(inputs, outputs, years=(dim(outputs@Catch)[3]+(inc_first-1))*inputs@phi_min, inc_first=FALSE)
          {
            return(get_annual_catch(Catch=outputs@Catch, years=years, phi_min=inputs@phi_min, inc_first=inc_first))
          })

#' @rdname get_annual_catch
setMethod("get_annual_catch", signature(inputs="missing", outputs="LeMans_outputs"),
          function(outputs, years=(dim(outputs@Catch)[3]+(inc_first-1))*phi_min, phi_min=0.1, inc_first=FALSE)
          {
            return(get_annual_catch(Catch=outputs@Catch, years=years, phi_min=phi_min, inc_first=inc_first))
          })

#' @export
#' @rdname get_annual_catch
setGeneric('get_CPUE', function(inputs, outputs, ...)
  standardGeneric('get_CPUE'))

#' @rdname get_annual_catch
setMethod('get_CPUE', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs, outputs, effort, years=nrow(effort), inc_first=FALSE){
            get_CPUE(Catch=outputs@Catch, Qs=inputs@Qs, effort=effort, years=years, phi_min=inputs@phi_min, inc_first=inc_first)
          })

#' @rdname get_annual_catch
setMethod('get_CPUE', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs, Catch, effort, years=nrow(effort), inc_first=FALSE){
            get_CPUE(Catch=Catch, Qs=inputs@Qs, effort=effort, years=years, phi_min=inputs@phi_min, inc_first=inc_first)
          })

#' @rdname get_annual_catch
setMethod('get_CPUE', signature(inputs="missing", outputs='LeMans_outputs'),
          function(outputs, Qs, effort, years=nrow(effort), phi_min=0.1, inc_first=FALSE){
            get_CPUE(Catch=outputs@Catch, Qs=Qs, effort=effort, years=years, phi_min=phi_min, inc_first=inc_first)
          })

#' @rdname get_annual_catch
setMethod('get_CPUE', signature(inputs="missing", outputs='missing'),
          function(Catch, Qs, effort, years=nrow(effort), phi_min=0.1, inc_first=FALSE) {
            if (inc_first==FALSE){
              Catch <- Catch[,,-1]
            }
            if (years<=0) {
              stop("years must be positive")
            }
            tot_time <- dim(Catch)[3]
            if ((1/phi_min)%%1!=0) {
              warning("The time step is not a multiple of a year")
            }
            Fs <- array(0, dim=c(dim(Catch)[1], dim(Catch)[2], tot_time, dim(Qs)[3]))
            if (ncol(effort)!=dim(Qs)[3]) {
              stop("Incorrect number of gears in effort")
            }
            year <- ceiling(1:tot_time*phi_min)
            if (is.null(colnames(effort))) {
              for (j in 1:ncol(effort)) {
                for (ts in 1:tot_time) {
                  Fs[,,ts,j] <- effort[year[ts], j]*Qs[,,j]
                }
              }
            } else {
              if (identical(colnames(effort), dimnames(Qs)[[3]])) {
                dimnames(Fs)[[4]] <- colnames(effort)
                for (j in colnames(effort)){
                  for (ts in 1:tot_time){
                    Fs[,,ts,j] <- effort[year[ts], j]*Qs[,,j]
                  }
                }
              } else {
                stop("Incorrect gear names in effort")
              }
            }
            CPUE <- matrix(0, dim(Catch)[3], dim(Catch)[2])
            for (ts in 1:tot_time) {
              tmp <- apply(Fs[,,ts,], c(1,2), sum)
              for(j in 1:dim(Qs)[[3]]) {
                tele <- Fs[,,ts,j]/tmp
                if (any(is.na(tele))){
                  tele[is.na(tele)] <- 0
                }
                telly <- colSums(Catch[,,ts]*tele)/effort[year[ts], j]
                if (any(is.na(telly))){telly[any(is.na(telly))] <- 0}
                CPUE[ts, ] <- CPUE[ts, ]+telly
              }
            }
            return(CPUE)
          })

#' @export
#' @rdname get_annual_catch
setGeneric('get_CPG', function(inputs, outputs, ...)
  standardGeneric('get_CPG'))

#' @rdname get_annual_catch
setMethod('get_CPG', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs, outputs, effort, years=nrow(effort), inc_first=FALSE){
            get_CPG(Catch=outputs@Catch, Qs=inputs@Qs, effort=effort, years=years, phi_min=inputs@phi_min, inc_first=inc_first)
          })

#' @rdname get_annual_catch
setMethod('get_CPG', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs, Catch, effort, years=nrow(effort), inc_first=FALSE){
            get_CPG(Catch=Catch, Qs=inputs@Qs, effort=effort, years=years, phi_min=inputs@phi_min, inc_first=inc_first)
          })

#' @rdname get_annual_catch
setMethod('get_CPG', signature(inputs="missing", outputs='LeMans_outputs'),
          function(outputs, Qs, effort, years=nrow(effort), phi_min=0.1, inc_first=FALSE){
            get_CPG(Catch=outputs@Catch, Qs=Qs, effort=effort, years=years, phi_min=phi_min, inc_first=inc_first)
          })

#' @rdname get_annual_catch
setMethod('get_CPG', signature(inputs="missing", outputs='missing'),
          function(Catch, Qs, effort, years=nrow(effort), phi_min=0.1, inc_first=FALSE) {
            if (inc_first==FALSE) {
              Catch <- Catch[,,-1]
            }
            if (years<=0) {
              stop("years must be positive")
            }
            tot_time <- dim(Catch)[3]
            if ((1/phi_min)%%1!=0) {
              warning("The time step is not a multiple of a year")
            }
            Fs <- array(0, dim=c(dim(Catch)[1], dim(Catch)[2], tot_time,dim(Qs)[3]))
            if (ncol(effort)!=dim(Qs)[3]) {
              stop("Incorrect number of gears in effort")
            }
            year <- ceiling(1:tot_time*phi_min)
            if (is.null(colnames(effort))) {
              for (j in 1:ncol(effort)) {
                for (ts in 1:tot_time) {
                  Fs[,,ts,j] <- effort[year[ts], j]*Qs[,,j]
                }
              }
            } else {
              if (identical(colnames(effort), dimnames(Qs)[[3]])) {
                dimnames(Fs)[[4]] <- colnames(effort)
                for (j in colnames(effort)) {
                  for (ts in 1:tot_time) {
                    Fs[,,ts,j] <- effort[year[ts], j]*Qs[,,j]
                  }
                }
              } else {
                stop("Incorrect gear names in effort")
              }
            }
            CPG <- array(0, dim=c(dim(Catch)[2], dim(Qs)[3], dim(Catch)[3]))
            for (ts in 1:tot_time) {
              tmp <- apply(Fs[,,ts,], c(1,2), sum)
              for(j in 1:dim(CPG)[2]) {
                tele <- Fs[,,ts,j]/tmp
                if (any(is.na(tele))) {
                  tele[is.na(tele)] <- 0
                }
                CPG[,j,ts] <- colSums(Catch[,,ts]*tele)
              }
            }
            return(CPG)
          })
