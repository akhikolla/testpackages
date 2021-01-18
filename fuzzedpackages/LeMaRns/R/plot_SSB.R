#' @include run_LeMans.R
#' @include param_setup.R
NULL

#' Plot Spawning Stock Biomass (SSB)
#'
#' @description Plots community and/or species-specific Spawning Stock Biomass (SSB) or total biomass.
#' @param inputs A \code{\linkS4class{LeMans_param}} object containing the parameter values of the current LeMans model. This option is only required if \code{SSB} is not provided.
#' @param outputs A \code{\linkS4class{LeMans_outputs}} object containing the outputs of the model run. This option is only required if \code{SSB} is not provided.
#' @param wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class. This option is only required if \code{inputs} and \code{SSB} are not provided.
#' @param mature A matrix with dimensions \code{nsc} and \code{nfish} and elements in the range 0-1 representing the proportion of individuals that are mature for each species and length class, where \code{nsc} represents the number of length classes in the model and \code{nfish} represents the number of species in the model. This option is only required if \code{outputs} and \code{SSB} are not provided.
#' @param N An array with dimensions \code{nsc}, \code{nfish} and \code{tot_time} representing the number of individuals in each length class for each time step, where \code{nsc} represents the number of length classes in the model, \code{nfish} represents the number of species in the model and \code{tot_time} represents the number of time steps that the model was run for. This option is only required if \code{outputs} and \code{SSB} are not provided.
#' @param species A numeric value or vector or a character string or vector denoting the species to be used in the plot(s). The default is \code{1:dim(N)[2]}.
#' @param time_steps A numeric vector denoting the time steps to be used to calculate and/or plot SSB. The default is \code{1:dim(N)[3]} or \code{1:length(SSB)}.
#' @param species_names A character vector of length \code{1:dim(N)[2]} that denotes the names of the species in the model. This option is only required if \code{inputs} is not provided.
#' @param SSB A numeric vector or a matrix representing the outputs of the function \code{get_SSB()}. This option is only required if \code{inputs} and \code{outputs} are not provided.
#' @param biomass A numeric vector or a matrix representing the outputs of the function \code{get_biomass()}. This option is only required if \code{inputs} and \code{outputs} are not provided.
#' @param full_plot_only A logical statement indicating whether a single plot depicting the SSB of all of the selected species should be produced (\code{full_plot_only=TRUE}) or multiple plots depicting the SSB of individual species should be produced (\code{full_plot_only=FALSE}). The default is \code{TRUE}.
#' @param units A character string denoting the units of weight used in the model. The default is \code{"g"}.
#' @param ... Additional arguments.
#' @return \code{plot_SSB} returns line plots of the in SSB of the selected species through time.
#' @return \code{plot_biomass} returns line plots of the changes in biomass of the selected species through time.
#' @seealso \code{\link{get_SSB}}, \code{\link{get_biomass}}
#' @examples
#' # Set up and run the model
#' NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=1e12)
#' effort <- matrix(0.5, 10, dim(NS_params@Qs)[3])
#' model_run <- run_LeMans(NS_params, years=10, effort=effort)
#'
#' # Calculate SSB
#' SSB <- get_SSB(NS_params, model_run)
#'
#' # Plot SSB
#' plot_SSB(SSB=SSB)
#'
#' # Calculate biomass
#' biomass <- get_biomass(NS_params, model_run)
#'
#' # Plot biomass
#' plot_biomass(biomass=biomass)
#' @export
setGeneric('plot_SSB', function(inputs, outputs, ...)
  standardGeneric('plot_SSB'))

#' @rdname plot_SSB
setMethod('plot_SSB', signature(inputs="LeMans_param", outputs="LeMans_outputs"),
          function(inputs, outputs, species, time_steps, species_names, SSB, ...) {
            if (missing(SSB)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              return(plot_SSB(wgt=inputs@wgt, mature=inputs@mature, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (is.matrix(SSB)) {
              if (missing(species)) {
                species <- 1:ncol(SSB)
              }
              if (missing(time_steps)) {
                time_steps <- 1:nrow(SSB)
              }
            } else if (is.vector(SSB)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(SSB)
              }
            }
            if(missing(species_names)){
              species_names <- inputs@species_names
            }
            return(plot_SSB(species=species, time_steps=time_steps, species_names=species_names, SSB=SSB, ...))
          })

#' @rdname plot_SSB
setMethod('plot_SSB', signature(inputs="LeMans_param", outputs="missing"),
          function(inputs, N, species, time_steps, species_names, SSB, ...) {
            if (missing(SSB)) {
              if (missing(species)) {
                species <- 1:dim(N)[2]
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(N)[3]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              return(plot_SSB(wgt=inputs@wgt, mature=inputs@mature, N=N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (is.matrix(SSB)) {
              if (missing(species)) {
                species <- 1:ncol(SSB)
              }
              if (missing(time_steps)) {
                time_steps <- 1:nrow(SSB)
              }
            } else if (is.vector(SSB)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(SSB)
              }
            }
            if(missing(species_names)){
              species_names <- inputs@species_names
            }
            return(plot_SSB(species=species, time_steps=time_steps, species_names=species_names, SSB=SSB, ...))
          })

#' @rdname plot_SSB
setMethod('plot_SSB', signature(inputs="missing", outputs="LeMans_outputs"),
          function(wgt,mature, outputs, species, time_steps, species_names, SSB, ...) {
            if (missing(SSB)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              if(missing(species_names)){
                species_names <- paste("Species", 1:length(species), sep=" ")
              }
              return(plot_SSB(wgt=wgt, mature=mature, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (is.matrix(SSB)) {
              if (missing(species)) {
                species <- 1:ncol(SSB)
              }
              if (missing(time_steps)) {
                time_steps <- 1:nrow(SSB)
              }
            } else if (is.vector(SSB)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(SSB)
              }
            }
            if(missing(species_names)){
              species_names <- paste("Species", 1:length(species), sep=" ")
            }
            return(plot_SSB(wgt=wgt, mature=mature, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
          })

#' @rdname plot_SSB
setMethod('plot_SSB', signature(inputs="missing", outputs="missing"),
          function(wgt, mature, N, species, species_names, time_steps, SSB, full_plot_only=TRUE, units="g", ...) {
            # If SSB is present but species and time_steps are missing, give default values
            if (!missing(SSB)) {
              if (is.matrix(SSB)) {
                if (missing(species)) {
                  species <- 1:ncol(SSB)
                }
                if (missing(species_names)) {
                  species_names <- paste("Species", species, sep=" ")
                }
                else{
                  if(is.character(species)){
                    species <- sapply(species, function(species){which(species==species_names)})
                  }
                  species_names <- species_names[species]
                }
                if (missing(time_steps)) {
                  time_steps <- 1:nrow(SSB)
                }
                SSB <- SSB[time_steps, species]
              } else if (is.vector(SSB)) {
                if (missing(time_steps)) {
                  time_steps <- 1:length(SSB)
                }
                SSB <- SSB[time_steps]
              }
            }
            # If SSB is missing, calculate it using N and wgt
            if (missing(SSB)) {
              if (missing(species)) {
                species <- 1:dim(N)[2]
              }
              if (missing(species_names)) {
                species_names <- paste("Species", 1:length(species), sep=" ")
              }
              if (is.character(species)) {
                species <- sapply(species, function(species) {which(species==species_names)})
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(N)[3]
              }
              SSB <- get_SSB(wgt=wgt, mature=mature, N=N, species=species, time_steps=time_steps, species_names=species_names)
              species_names <- species_names[species]
            }

            # Save users plot settings
            def.par <- par(no.readonly=TRUE)

            # Extract required data and plot
            if (is.vector(SSB)) {
              par(mfrow=c(1,1), mar=c(5,5,5,0))
              layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), nrow=3, byrow=TRUE))
              plot(time_steps, SSB, xlab="Time steps", ylab=paste("SSB (", units, ")", sep = ""), type="l", font.lab=2, cex.lab=1.5, cex.axis=1.5, log="y")
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=species_names, col="black", lty=1, cex=1.2, bty="n")
            } else if (is.matrix(SSB)) {
              if (full_plot_only==FALSE) {
                # Plot species individually in a 2x2 matrix
                par(mfrow=c(2,2), mar=c(5,5,5,5))
                for (i in 1:ncol(SSB)) {
                  plot(c(min(time_steps), max(time_steps)), c(min(SSB), max(SSB)), type="n", xlab="Time step", ylab=paste("SSB (", units, ")", sep = ""), font.lab=2, log="y", main=species_names[i])
                  lines(time_steps, SSB[, i])
                }
              }

              # Plot all species together
              par(mfrow=c(1,1), mar=c(5,5,5,0))
              layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), nrow=3, byrow=TRUE))
              rainbowcols <- rainbow(ncol(SSB), s=0.75)
              plot(c(min(time_steps), max(time_steps)), c(min(SSB), max(SSB)), type="n", xlab="Time step", ylab=paste("SSB (", units, ")", sep = ""), font.lab=2, cex.lab=1.5, cex.axis=1.5, log="y")
              for (i in 1:ncol(SSB)) {
                lines(time_steps, SSB[, i], col=rainbowcols[i])
              }
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=species_names, col=rainbowcols, lty=1, cex=1.2, bty="n")
            }

            # Reset plot settings
            par(def.par)
          })

#' @export
#' @rdname plot_SSB
setGeneric('plot_biomass', function(inputs, outputs, ...)
  standardGeneric('plot_biomass'))

#' @rdname plot_SSB
setMethod('plot_biomass', signature(inputs="LeMans_param", outputs="LeMans_outputs"),
          function(inputs, outputs, species, time_steps, species_names, biomass, ...) {
            if (missing(biomass)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              return(plot_biomass(wgt=inputs@wgt, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (is.matrix(biomass)) {
              if (missing(species)) {
                species <- 1:ncol(biomass)
              }
              if (missing(time_steps)) {
                time_steps <- 1:nrow(biomass)
              }
            } else if (is.vector(biomass)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(biomass)
              }
            }
            if(missing(species_names)){
              species_names <- inputs@species_names
            }
            return(plot_biomass(species=species, time_steps=time_steps, species_names=species_names, biomass=biomass, ...))
          })

#' @rdname plot_SSB
setMethod('plot_biomass', signature(inputs="LeMans_param", outputs="missing"),
          function(inputs, N, species, time_steps, biomass, species_names, ...) {
            if (missing(biomass)) {
              if (missing(species)) {
                species <- 1:dim(N)[2]
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(N)[3]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              return(plot_biomass(wgt=inputs@wgt, N=N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (is.matrix(biomass)) {
              if (missing(species)) {
                species <- 1:ncol(biomass)
              }
              if (missing(time_steps)) {
                time_steps <- 1:nrow(biomass)
              }
            } else if (is.vector(biomass)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(biomass)
              }
            }
            if(missing(species_names)){
              species_names <- inputs@species_names
            }
            return(plot_biomass(species=species, time_steps=time_steps, species_names=species_names, biomass=biomass, ...))
          })

#' @rdname plot_SSB
setMethod('plot_biomass', signature(inputs="missing", outputs="LeMans_outputs"),
          function(wgt, outputs, species, time_steps, biomass, species_names, ...) {
            if (missing(biomass)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              if(missing(species_names)){
                species_names <- paste("Species", 1:length(species), sep=" ")
              }
              return(plot_biomass(wgt=wgt, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (is.matrix(biomass)) {
              if (missing(species)) {
                species <- 1:ncol(biomass)
              }
              if (missing(time_steps)) {
                time_steps <- 1:nrow(biomass)
              }
            } else if (is.vector(biomass)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(biomass)
              }
            }
            if(missing(species_names)){
              species_names <- paste("Species", 1:length(species), sep=" ")
            }
            return(plot_biomass(wgt=wgt, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
          })

#' @rdname plot_SSB
setMethod('plot_biomass', signature(inputs="missing", outputs="missing"),
          function(wgt, N, species, time_steps, species_names, biomass, full_plot_only=TRUE, units="g", ...) {
            # If biomass is present but species and time_steps are missing, give default values
            if (!missing(biomass)) {
              if (is.matrix(biomass)) {
                if (missing(species)) {
                  species <- 1:ncol(biomass)
                }
                if (missing(species_names)) {
                  species_names <- paste("Species", species, sep=" ")
                }
                else{
                  if(is.character(species)){
                    species <- sapply(species, function(species){which(species==species_names)})
                  }
                  species_names <- species_names[species]
                }
                if (missing(time_steps)) {
                  time_steps <- 1:nrow(biomass)
                }
                biomass <- biomass[time_steps, species]
              } else if (is.vector(biomass)) {
                if (missing(time_steps)) {
                  time_steps <- 1:length(biomass)
                }
                biomass <- biomass[time_steps]
              }
            }
            # If biomass is missing, calculate it using N and wgt
            if (missing(biomass)) {
              if (missing(species)) {
                species <- 1:dim(N)[2]
              }
              if (missing(species_names)) {
                species_names <- paste("Species", 1:length(species), sep=" ")
              }
              if (is.character(species)) {
                species <- sapply(species, function(species) {which(species==species_names)})
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(N)[3]
              }
              biomass <- get_biomass(wgt=wgt, N=N, species=species, time_steps=time_steps, species_names=species_names)
              species_names <- species_names[species]
            }

            # Save users plot settings
            def.par <- par(no.readonly=TRUE)

            # Extract required data and plot
            if (is.vector(biomass)) {
              par(mfrow=c(1,1), mar=c(5,5,5,0))
              layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), nrow=3, byrow=TRUE))
              plot(time_steps, biomass, xlab="Time steps", ylab=paste("Biomass (", units, ")", sep = ""), type="l", font.lab=2, cex.lab=1.5, cex.axis=1.5, log="y")
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=species_names, col="black", lty=1, cex=1.2, bty="n")
            } else if (is.matrix(biomass)) {
              if (full_plot_only==FALSE) {
                # Plot species individually in a 2x2 matrix
                par(mfrow=c(2,2), mar=c(5,5,5,5))
                for (i in 1:ncol(biomass)) {
                  plot(c(min(time_steps), max(time_steps)), c(min(biomass), max(biomass)), type="n", xlab="Time step", ylab=paste("Biomass (", units, ")", sep = ""), font.lab=2, log="y", main=species_names[i])
                  lines(time_steps, biomass[, i])
                }
              }

              # Plot all species together
              par(mfrow=c(1,1), mar=c(5,5,5,0))
              layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), nrow=3, byrow=TRUE))
              rainbowcols <- rainbow(ncol(biomass), s=0.75)
              plot(c(min(time_steps), max(time_steps)), c(min(biomass), max(biomass)), type="n", xlab="Time step", ylab=paste("Biomass (", units, ")", sep = ""), font.lab=2, cex.lab=1.5, cex.axis=1.5, log="y")
              for (i in 1:ncol(biomass)) {
                lines(time_steps, biomass[, i], col=rainbowcols[i])
              }
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              if (!missing(species_names)) {
                legend("center", legend=species_names, col=rainbowcols, lty=1, cex=1.2, bty="n")
              }
            }

            # Reset plot settings
            par(def.par)
          })
