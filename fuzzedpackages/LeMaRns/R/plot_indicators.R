#' @include run_LeMans.R
#' @include param_setup.R
NULL

#' Plot indicators
#'
#' @description Plots Mean Maximum Length (MML), the Large Fish Indicator (LFI), Typical Length (TyL) and the Length Quantiles (LQ) of the whole community or a subset of the species.
#' @param inputs A \code{\linkS4class{LeMans_param}} object containing the parameter values of the current LeMans model.
#' @param outputs A \code{\linkS4class{LeMans_outputs}} object containing the outputs of the model run.
#' @param wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class. This option is only required if \code{inputs} or \code{MML}, \code{LFI}, \code{TyL} and \code{LQ} are not provided.
#' @param mid A numeric vector of length \code{nfish} representing the mid-point of the length classes in the model, where \code{nfish} represents the number of species in the model. This option is only required if \code{inputs} or \code{TyL} are not provided.
#' @param l_bound A numeric vector of length \code{nsc} representing the lower bounds of the length classes, where \code{nsc} represents the number of length classes in the model. This option is only required if \code{inputs} or \code{LFI} are not provided.
#' @param u_bound A numeric vector of length \code{nsc} representing the upper bounds of the length classes, where \code{nsc} represents the number of length classes in the model. This option is only required if \code{inputs}, \code{LFI} and/or \code{LQ} are not provided.
#' @param Linf A numeric vector of length \code{nfish} representing the asymptotic length of each species, where \code{nfish} represents the number of species in the model. This option is only required if \code{inputs} or \code{MML} are not provided.
#' @param N An array with dimensions \code{nsc}, \code{nfish} and \code{tot_time} representing the number of individuals in each length class for each time step, where \code{nsc} represents the number of length classes in the model, \code{nfish} represents the number of species in the model and \code{tot_time} represents the number of time steps that the model was run for. This option is only required if \code{outputs} or \code{MML}, \code{LFI}, \code{TyL} and \code{LQ} are not provided.
#' @param species A numeric value or vector or a character string or vector denoting the species to be used to calculate the indicators. This option is only required if \code{MML}, \code{LFI}, \code{TyL} and \code{LQ} are not provided. The default is \code{1:dim(N)[2]}.
#' @param species_names A character vector of length \code{1:dim(N)[2]} that denotes the names of the species in the model. This option is only required if \code{inputs} is not provided.
#' @param time_steps A numeric vector denoting the time steps to be used to calculate and/or plot the indicators. The default is \code{1:dim(N)[3]} or \code{1:length(indicator)}.
#' @param length_LFI A numeric value or vector representing the threshold(s) to be used to calculate the LFI. This option is only required if \code{LFI} is not provided. The default is \code{40}.
#' @param prob A numeric value or vector between 0 and 1 denoting the LQ(s) to be calculated. This option is only required if \code{LQ} is not provided. The default is \code{0.5}.
#' @param MML A numeric vector representing the outputs of the function \code{get_MML()}. This option is only required if \code{inputs} and \code{outputs} are not provided.
#' @param LFI A numeric vector or matrix representing the outputs of the function \code{get_LFI()}. This option is only required if \code{inputs} and \code{outputs} are not provided.
#' @param TyL A numeric vector representing the outputs of the function \code{get_TyL()}. This option is only required if \code{inputs} and \code{outputs} are not provided.
#' @param LQ A numeric vector or matrix representing the outputs of the function \code{get_LQ()}. This option is only required if \code{inputs} and \code{outputs} are not provided.
#' @param units A character string denoting the units of length used in the model. The default is \code{"cm"}.
#' @param ... Additional arguments.
#' @return \code{plot_indicators} returns a set of four line plots depicting changes in the MML, LFI, TyL and LQ(s) of the community (including only the selected species) over time.
#' @return \code{plot_LFI} returns a line plot depicting the LFI of the community (including only the selected species) over time.
#' @return \code{plot_MML} returns a line plot depicting the MML of the community (including only the selected species) over time.
#' @return \code{plot_TYL} returns a line plot depicting the TyL of the community (including only the selected species) over time.
#' @return \code{plot_LQ} returns a line plot depicting the LQ(s) of the community (including only the selected species) over time.
#' @seealso \code{\link{get_indicators}}
#' @examples
#' # Set up and run the model
#' NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=1e12)
#' effort <- matrix(0.5, 10, dim(NS_params@Qs)[3])
#' model_run <- run_LeMans(NS_params, years=10, effort=effort)
#'
#' # Calculate the indicators
#' tmp <- get_indicators(NS_params, model_run)
#' MML <- tmp$MML
#' LFI <- tmp$LFI
#' TyL <- tmp$TYL
#' LQ <- tmp$LQ
#'
#' # Plot the indicators
#' plot_indicators(MML=MML, LFI=LFI, TyL=TyL, LQ=LQ)
#'
#' # Plot the LFI
#' plot_LFI(LFI=LFI)
#'
#' # Plot MML
#' plot_MML(MML=MML)
#'
#' # Plot the TyL
#' plot_TyL(TyL=TyL)
#'
#' # Plot the LQs
#' plot_LQ(LQ=LQ)
#' @export
setGeneric('plot_indicators', function(inputs, outputs, ...)
  standardGeneric('plot_indicators'))

#' @rdname plot_indicators
setMethod('plot_indicators', signature(inputs="LeMans_param", outputs='LeMans_outputs'),
          function(inputs, outputs, species, time_steps, prob, length_LFI, ...) {
            return(plot_indicators(N=outputs@N, wgt=inputs@wgt, Linf=inputs@Linf, species=species, time_steps=time_steps, l_bound=inputs@l_bound, u_bound=inputs@u_bound, mid=inputs@mid, species_names=inputs@species_names, prob=prob, length_LFI=length_LFI, ...))
          })

#' @rdname plot_indicators
setMethod('plot_indicators', signature(inputs="LeMans_param", outputs='missing'),
          function(inputs, N, species, time_steps, prob, length_LFI, ...) {
            return(plot_indicators(N=N, wgt=inputs@wgt, Linf=inputs@Linf, species=species, time_steps=time_steps, l_bound=inputs@l_bound, u_bound=inputs@u_bound, mid=inputs@mid, species_names=inputs@species_names, prob=prob, length_LFI=length_LFI, ...))
          })

#' @rdname plot_indicators
setMethod('plot_indicators', signature(inputs="missing", outputs='LeMans_outputs'),
          function(wgt, mid, l_bound, u_bound, Linf, species, outputs, time_steps, species_names=NULL, prob, length_LFI, ...) {
            return(plot_indicators(wgt=wgt, mid=mid,l_bound=l_bound, u_bound=u_bound, Linf=Linf, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, prob=prob, length_LFI=length_LFI, ...))
          })

#' @rdname plot_indicators
setMethod('plot_indicators', signature(inputs="missing", outputs='missing'),
          function(wgt, mid, l_bound, u_bound, Linf, N, species, time_steps, species_names=NULL, prob, length_LFI, MML, TyL, LFI, LQ, units="cm", ...) {
            # Save users plot settings
            def.par <- par(no.readonly=TRUE)

            # Set up plot layout
            layout(matrix(c(1,1,1,2,2,2,3,1,1,1,2,2,2,3,4,4,4,5,5,5,6,4,4,4,5,5,5,6), nrow=4, byrow=TRUE))

            ############################ MML ############################
            if (missing(MML)) {
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
              MML <- get_MML(N=N, wgt=wgt, species=species, time_steps=time_steps, Linf=Linf)
            } else {
              # If MML is present but time_steps are missing, give default values
              if (missing(time_steps)) {
                time_steps <- 1:length(MML)
              }
              MML <- MML[time_steps]
            }

            # Plot
            par(mar=c(5,5,5,0))
            plot(x=time_steps, y=MML, ylab=paste("MML (", units, ")", sep = ""), xlab="Time steps", type="l", font.lab=2, cex.lab=1.5, cex.axis=1.5)

            ############################ LFI ############################
            if (missing(LFI)) {
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
              if (missing(length_LFI)) {
                length_LFI <- 40
              }
              LFI <- get_LFI(N=N, wgt=wgt, species=species, time_steps=time_steps, species_names=species_names, length_LFI=length_LFI, l_bound=l_bound, u_bound=u_bound)
            } else {
              if (is.matrix(LFI)) {
                if (missing(time_steps)) {
                  time_steps <- 1:nrow(LFI)
                }
                if (missing(length_LFI)) {
                  length_LFI <- paste("LFI", 1:ncol(LFI), sep=" ")
                }
                LFI <- LFI[time_steps, ]
              } else if (is.vector(LFI)) {
                if (missing(time_steps)) {
                  time_steps <- 1:length(LFI)
                }
                if (missing(length_LFI)) {
                  length_LFI <- "LFI 1"
                }
                LFI <- LFI[time_steps]
              }
            }

            # Extract required data and plot
            if (is.vector(LFI)) {
              par(mar=c(5,5,5,0))
              plot(x=time_steps, y=LFI, xlab="Time steps", ylab="LFI", type="l", font.lab=2, cex.lab=1.5, cex.axis=1.5)
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=paste(length_LFI, " ", units, sep=""), col="black", lty=1, cex=1.5, bty="n")
            } else if (is.matrix(LFI)) {
              # Plot the quantiles together
              rainbowcols <- rainbow(ncol(LFI), s=0.75)
              par(mar=c(5,5,5,0))
              plot(c(min(time_steps), max(time_steps)), c(min(LFI), max(LFI)), type="n", xlab="Time steps", ylab="LFI", font.lab=2, cex.lab=1.5, cex.axis=1.5)
              for (i in 1:ncol(LFI)) {
                lines(time_steps, LFI[, i], col=rainbowcols[i])
              }
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=paste(length_LFI, " ", units, sep=""), col=rainbowcols, lty=1, cex=1.5, bty="n")
            }

            ############################ TyL ############################
            if (missing(TyL)) {
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
              TyL <- get_TyL(N=N, wgt=wgt, species=species, time_steps=time_steps, mid=mid)
            } else {
              if (missing(time_steps)) {
                time_steps <- 1:length(TyL)
              }
              TyL <- TyL[time_steps]
            }

            # Plot
            par(mar=c(5,5,5,0))
            plot(x=time_steps, y=TyL, ylab=paste("TyL (", units, ")", sep = ""), xlab="Time steps", type="l", font.lab=2, cex.lab=1.5, cex.axis=1.5)

            ############################ LQ ############################
            if (missing(LQ)) {
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
              if (missing(prob)) {
                prob <- 0.5
              }
              LQ <- get_LQ(N=N, wgt=wgt, species=species, time_steps=time_steps, species_names=species_names, u_bound=u_bound, prob=prob)
            } else {
              if (is.matrix(LQ)) {
                if (missing(time_steps)) {
                  time_steps <- 1:nrow(LQ)
                }
                if (missing(prob)) {
                  prob <- paste("Quantile", 1:ncol(LQ), sep=" ")
                }
                LQ <- LQ[time_steps, ]
              } else if (is.vector(LQ)) {
                if (missing(time_steps)) {
                  time_steps <- 1:length(LQ)
                }
                if (missing(prob)) {
                  prob <- "Quantile 1"
                }
                LQ <- LQ[time_steps]
              }
            }

            # Save users plot settings
            def.par <- par(no.readonly=TRUE)

            # Extract required data and plot
            if (is.vector(LQ)) {
              par(mar=c(5,5,5,0))
              plot(x=time_steps, y=LQ, xlab="Time steps", ylab=paste("LQ (", units, ")", sep = ""), type="l", font.lab=2, cex.lab=1.5, cex.axis=1.5)
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=prob, col="black", lty=1, cex=1.5, bty="n")
            } else if (is.matrix(LQ)) {
              # Plot the quantiles together
              rainbowcols <- rainbow(ncol(LQ), s=0.75)
              par(mar=c(5,5,5,0))
              plot(c(min(time_steps), max(time_steps)), c(min(LQ), max(LQ)), type="n", xlab="Time steps", ylab=paste("LQ (", units, ")", sep = ""), font.lab=2,  cex.lab=1.5, cex.axis=1.5)
              for (i in 1:ncol(LQ)) {
                lines(time_steps, LQ[, i], col=rainbowcols[i])
              }
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=prob, col=rainbowcols, lty=1, cex=1.5, bty="n")
            }

            # Reset plot settings
            par(def.par)
          })

#' @export
#' @rdname plot_indicators
setGeneric('plot_LFI', function(inputs, outputs, ...)
  standardGeneric('plot_LFI'))

#' @rdname plot_indicators
setMethod('plot_LFI', signature(inputs="LeMans_param", outputs="LeMans_outputs"),
          function(inputs, outputs, species, time_steps, species_names, length_LFI, LFI, ...) {
            if (missing(LFI)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              if (missing(length_LFI)) {
                length_LFI <- 40
              }
              return(plot_LFI(wgt=inputs@wgt, l_bound=inputs@l_bound, u_bound=inputs@u_bound, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, length_LFI=length_LFI, ...))
            }
            if (is.matrix(LFI)) {
              if (missing(time_steps)) {
                time_steps <- 1:nrow(LFI)
              }
              if (missing(length_LFI)) {
                length_LFI <- paste("LFI", ncol(LFI), sep=" ")
              }
            } else if (is.vector(LFI)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(LFI)
              }
              if (missing(length_LFI)) {
                length_LFI <- "LFI 1"
              }
            }
            return(plot_LFI(time_steps=time_steps, LFI=LFI, length_LFI=length_LFI))
          })

#' @rdname plot_indicators
setMethod('plot_LFI', signature(inputs="LeMans_param", outputs="missing"),
          function(inputs, N, species, time_steps, species_names, length_LFI, LFI, ...) {
            if (missing(LFI)) {
              if (missing(species)) {
                species <- 1:dim(N)[2]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(N)[3]
              }
              if (missing(length_LFI)) {
                length_LFI <- 40
              }
              return(plot_LFI(wgt=inputs@wgt, l_bound=inputs@l_bound, u_bound=inputs@u_bound, N=N, species=species, time_steps=time_steps, species_names=species_names, length_LFI=length_LFI, ...))
            }
            if (is.matrix(LFI)) {
              if (missing(time_steps)) {
                time_steps <- 1:nrow(LFI)
              }
              if (missing(length_LFI)) {
                length_LFI <- paste("LFI", ncol(LFI), sep=" ")
              }
            } else if (is.vector(LFI)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(LFI)
              }
              if (missing(length_LFI)) {
                length_LFI <- "LFI 1"
              }
            }
            return(plot_LFI(time_steps=time_steps, LFI=LFI, length_LFI=length_LFI))
          })

#' @rdname plot_indicators
setMethod('plot_LFI', signature(inputs="missing", outputs="LeMans_outputs"),
          function(wgt, l_bound, u_bound, outputs, species, time_steps, species_names, length_LFI, LFI, ...) {
            if (missing(LFI)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if(missing(species_names)){
                species_names <- paste("Species", 1:length(species), sep=" ")
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              if (missing(length_LFI)) {
                length_LFI <- 40
              }
              return(plot_LFI(wgt=wgt, l_bound=l_bound, u_bound=u_bound, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, length_LFI=length_LFI, ...))
            }
            if (is.matrix(LFI)) {
              if (missing(time_steps)) {
                time_steps <- 1:nrow(LFI)
              }
              if (missing(length_LFI)) {
                length_LFI <- paste("LFI", ncol(LFI), sep=" ")
              }
            } else if (is.vector(LFI)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(LFI)
              }
              if (missing(length_LFI)) {
                length_LFI <- "LFI 1"
              }
            }
            return(plot_LFI(time_steps=time_steps, LFI=LFI, length_LFI=length_LFI))
          })

#' @rdname plot_indicators
setMethod('plot_LFI', signature(inputs="missing", outputs="missing"),
          function(wgt, l_bound, u_bound, N, species, time_steps, species_names, length_LFI, LFI, units="cm", ...) {
            # If LFI is present but species and time_steps are missing, give default values
            if (!missing(LFI)) {
              if (is.matrix(LFI)) {
                if (missing(time_steps)) {
                  time_steps <- 1:nrow(LFI)
                }
                if (missing(length_LFI)) {
                  length_LFI <- paste("LFI", ncol(LFI), sep=" ")
                }
                LFI <- LFI[time_steps, ]
              } else if (is.vector(LFI)) {
                if (missing(time_steps)) {
                  time_steps <- 1:length(LFI)
                }
                if (missing(length_LFI)) {
                  length_LFI <- "LFI 1"
                }
                LFI <- LFI[time_steps]
              }
            }

            # If LFI is missing, calculate it using N and wgt
            if (missing(LFI)) {
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
              if (missing(length_LFI)) {
                length_LFI <- 40
              }
              LFI <- get_LFI(wgt=wgt, l_bound=l_bound, u_bound=u_bound, N=N, species=species, time_steps=time_steps, species_names=species_names, length_LFI=length_LFI)
            }

            # Save users plot settings
            def.par <- par(no.readonly=TRUE)

            # Extract required data and plot
            if (is.vector(LFI)) {
              par(mfrow=c(1,1), mar=c(5,5,5,0))
              layout(matrix(c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2), nrow=3, byrow=TRUE))
              plot(time_steps, LFI, xlab="Time steps", ylab="LFI", type="l", font.lab=2, cex.lab=1.5, cex.axis=1.5)
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=paste(length_LFI, " ", units, sep=""), col="black", lty=1, cex=1.5, bty="n")
            } else if (is.matrix(LFI)) {
              # Plot the quantiles together
              par(mfrow=c(1,1), mar=c(5,5,5,0))
              layout(matrix(c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2), nrow=3, byrow=TRUE))
              rainbowcols <- rainbow(ncol(LFI), s=0.75)
              plot(c(min(time_steps), max(time_steps)), c(min(LFI), max(LFI)), type="n", xlab="Time steps", ylab="LFI", font.lab=2, cex.lab=1.5, cex.axis=1.5)
              for (i in 1:ncol(LFI)) {
                lines(time_steps, LFI[, i], col=rainbowcols[i])
              }
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=paste(length_LFI, " ", units, sep=""), col=rainbowcols, lty=1, cex=1.5, bty="n")
            }

            # Reset plot settings
            par(def.par)
          })

#' @export
#' @rdname plot_indicators
setGeneric('plot_MML', function(inputs, outputs, ...)
  standardGeneric('plot_MML'))

#' @rdname plot_indicators
setMethod('plot_MML', signature(inputs="LeMans_param", outputs="LeMans_outputs"),
          function(inputs, outputs, species, species_names, time_steps, MML, ...) {
            if (missing(MML)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              return(plot_MML(wgt=inputs@wgt, Linf=inputs@Linf, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (missing(time_steps)) {
              time_steps <- 1:length(MML)
            }
            return(plot_MML(time_steps=time_steps, MML=MML))
          })

#' @rdname plot_indicators
setMethod('plot_MML', signature(inputs="LeMans_param", outputs="missing"),
          function(inputs, N, species, time_steps, species_names, MML, ...) {
            if (missing(MML)) {
              if (missing(species)) {
                species <- 1:dim(N)[2]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(N)[3]
              }
              return(plot_MML(wgt=inputs@wgt, Linf=inputs@Linf, N=N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (missing(time_steps)) {
              time_steps <- 1:length(MML)
            }
            return(plot_MML(time_steps=time_steps, MML=MML))
          })

#' @rdname plot_indicators
setMethod('plot_MML', signature(inputs="missing", outputs="LeMans_outputs"),
          function(wgt, Linf, outputs, species, time_steps, species_names, MML, ...) {
            if (missing(MML)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if(missing(species_names)){
                species_names <- paste("Species", 1:length(species), sep=" ")
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              return(plot_MML(wgt=wgt, Linf=Linf, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (missing(time_steps)) {
              time_steps <- 1:length(MML)
            }
            return(plot_MML(time_steps=time_steps, MML=MML))
          })

#' @rdname plot_indicators
setMethod('plot_MML', signature(inputs="missing", outputs="missing"),
          function(wgt, Linf, N, species, time_steps, species_names, MML, units="cm", ...) {
            # Save users plot settings
            def.par <- par(no.readonly=TRUE)

            # If MML is missing, calculate it using N and wgt
            if (missing(MML)) {
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
              MML <- get_MML(wgt=wgt, Linf=Linf, N=N, species=species, time_steps=time_steps)
            } else {
              # If MML is present but time_steps are missing, give default values
              if (missing(time_steps)) {
                time_steps <- 1:length(MML)
              }
              MML <- MML[time_steps]
            }

            # Plot
            par(mfrow=c(1,1), mar=c(5,5,5,5))
            plot(time_steps, MML, ylab=paste("MML (", units, ")", sep = ""), xlab="Time steps", type="l", font.lab=2, cex.lab=1, cex.axis=1)

            # Reset plot settings
            par(def.par)
          })

#' @export
#' @rdname plot_indicators
setGeneric('plot_TyL', function(inputs, outputs, ...)
  standardGeneric('plot_TyL'))

#' @rdname plot_indicators
setMethod('plot_TyL', signature(inputs="LeMans_param", outputs="LeMans_outputs"),
          function(inputs, outputs, species, time_steps, species_names, TyL, ...) {
            if (missing(TyL)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              return(plot_TyL(wgt=inputs@wgt, mid=inputs@mid, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (missing(time_steps)) {
              time_steps <- 1:length(TyL)
            }
            return(plot_TyL(time_steps=time_steps, TyL=TyL))
          })

#' @rdname plot_indicators
setMethod('plot_TyL', signature(inputs="LeMans_param", outputs="missing"),
          function(inputs, N, species, time_steps, species_names, TyL, ...) {
            if (missing(TyL)) {
              if (missing(species)) {
                species <- 1:dim(N)[2]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(N)[3]
              }
              return(plot_TyL(wgt=inputs@wgt, mid=inputs@mid, N=N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (missing(time_steps)) {
              time_steps <- 1:length(TyL)
            }
            return(plot_TyL(time_steps=time_steps, TyL=TyL))
          })

#' @rdname plot_indicators
setMethod('plot_TyL', signature(inputs="missing", outputs="LeMans_outputs"),
          function(wgt, mid, outputs, species, time_steps, species_names, TyL, ...) {
            if (missing(TyL)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if(missing(species_names)){
                species_names <- paste("Species", 1:length(species), sep=" ")
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              return(plot_TyL(wgt=wgt, mid=mid, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, ...))
            }
            if (missing(time_steps)) {
              time_steps <- 1:length(TyL)
            }
            return(plot_TyL(time_steps=time_steps, TyL=TyL))
          })

#' @rdname plot_indicators
setMethod('plot_TyL', signature(inputs="missing", outputs="missing"),
          function(wgt, mid, N, species, time_steps, species_names, TyL, units="cm", ...) {
            # Save users plot settings
            def.par <- par(no.readonly=TRUE)

            # If TyL is missing, calculate it using N and wgt
            if (missing(TyL)) {
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
              TyL <- get_TyL(wgt=wgt, mid=mid, N=N, species=species, time_steps=time_steps)
            } else {
              # If TyL is present but time_steps are missing, give default values
              if (missing(time_steps)) {
                time_steps <- 1:length(TyL)
              }
              TyL <- TyL[time_steps]
            }

            # Plot
            par(mfrow=c(1,1), mar=c(5,5,5,5))
            plot(time_steps, TyL, ylab=paste("TyL (", units, ")", sep = ""), xlab="Time steps", type="l", font.lab=2, cex.lab=1, cex.axis=1)

            # Reset plot settings
            par(def.par)
          })

#' @export
#' @rdname plot_indicators
setGeneric('plot_LQ', function(inputs, outputs, ...)
  standardGeneric('plot_LQ'))

#' @rdname plot_indicators
setMethod('plot_LQ', signature(inputs="LeMans_param", outputs="LeMans_outputs"),
          function(inputs, outputs, species, time_steps, species_names, LQ, prob, ...) {
            if (missing(LQ)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              if (missing(prob)) {
                prob <- 0.5
              }
              return(plot_LQ(wgt=inputs@wgt, u_bound=inputs@u_bound, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, prob=prob, ...))
            }
            if (is.matrix(LQ)) {
              if (missing(time_steps)) {
                time_steps <- 1:nrow(LQ)
              }
              if (missing(prob)) {
                prob <- paste("Quantile", 1:ncol(LQ), sep=" ")
              }
            } else if (is.vector(LQ)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(LQ)
              }
              if (missing(prob)) {
                prob <- "Quantile 1"
              }
            }
            return(plot_LQ(time_steps=time_steps, LQ=LQ, prob=prob))
          })

#' @rdname plot_indicators
setMethod('plot_LQ', signature(inputs="LeMans_param", outputs="missing"),
          function(inputs, N, species, species_names, time_steps, LQ, prob, ...) {
            if (missing(LQ)) {
              if (missing(species)) {
                species <- 1:dim(N)[2]
              }
              if(missing(species_names)){
                species_names <- inputs@species_names
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(N)[3]
              }
              if (missing(prob)) {
                prob <- 0.5
              }
              return(plot_LQ(wgt=inputs@wgt, u_bound=inputs@u_bound, N=N, species=species, time_steps=time_steps, species_names=species_names, prob=prob, ...))
            }
            if (is.matrix(LQ)) {
              if (missing(time_steps)) {
                time_steps <- 1:nrow(LQ)
              }
              if (missing(prob)) {
                prob <- paste("Quantile", 1:ncol(LQ), sep=" ")
              }
            } else if (is.vector(LQ)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(LQ)
              }
              if (missing(prob)) {
                prob <- "Quantile 1"
              }
            }
            return(plot_LQ(time_steps=time_steps, LQ=LQ, prob=prob))
          })

#' @rdname plot_indicators
setMethod('plot_LQ', signature(inputs="missing", outputs="LeMans_outputs"),
          function(wgt, u_bound, outputs, species, time_steps, species_names, LQ, prob, ...) {
            if (missing(LQ)) {
              if (missing(species)) {
                species <- 1:dim(outputs@N)[2]
              }
              if(missing(species_names)){
                species_names <- paste("Species", 1:length(species), sep=" ")
              }
              if (missing(time_steps)) {
                time_steps <- 1:dim(outputs@N)[3]
              }
              if (missing(prob)) {
                prob <- 0.5
              }
              return(plot_LQ(wgt=wgt, u_bound=u_bound, N=outputs@N, species=species, time_steps=time_steps, species_names=species_names, prob=prob, ...))
            }
            if (is.matrix(LQ)) {
              if (missing(time_steps)) {
                time_steps <- 1:nrow(LQ)
              }
              if (missing(prob)) {
                prob <- paste("Quantile", 1:ncol(LQ), sep=" ")
              }
            } else if (is.vector(LQ)) {
              if (missing(time_steps)) {
                time_steps <- 1:length(LQ)
              }
              if (missing(prob)) {
                prob <- "Quantile 1"
              }
            }
            return(plot_LQ(time_steps=time_steps, LQ=LQ, prob=prob))
          })

#' @rdname plot_indicators
setMethod('plot_LQ', signature(inputs="missing", outputs="missing"),
          function(wgt, u_bound, N, species, time_steps, species_names, LQ, prob, units="cm", ...) {
            # If LQ is present but species and time_steps are missing, give default values
            if (!missing(LQ)) {
              if (is.matrix(LQ)) {
                if (missing(time_steps)) {
                  time_steps <- 1:nrow(LQ)
                }
                if (missing(prob)) {
                  prob <- paste("Quantile", 1:ncol(LQ), sep=" ")
                }
                LQ <- LQ[time_steps, ]
              } else if (is.vector(LQ)) {
                if (missing(time_steps)) {
                  time_steps <- 1:length(LQ)
                }
                if (missing(prob)) {
                  prob <- "Quantile 1"
                }
                LQ <- LQ[time_steps]
              }
            }

            # If LQ is missing, calculate it using N and wgt
            if (missing(LQ)) {
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
              if (missing(prob)) {
                prob <- 0.5
              }
              LQ <- get_LQ(wgt=wgt, u_bound=u_bound, N=N, species=species, species_names=species_names, time_steps=time_steps, prob=prob)
            }

            # Save users plot settings
            def.par <- par(no.readonly=TRUE)

            # Extract required data and plot
            if (is.vector(LQ)) {
              par(mfrow=c(1,1), mar=c(5,5,5,0))
              layout(matrix(c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2), nrow=3, byrow=TRUE))
              plot(time_steps, LQ, xlab="Time steps", ylab=paste("LQ (", units, ")", sep = ""), type="l", font.lab=2, cex.lab=1.5, cex.axis=1.5)
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=prob, col="black", lty=1, cex=1.5, bty="n")
            } else if (is.matrix(LQ)) {
              # Plot the quantiles together
              par(mfrow=c(1,1), mar=c(5,5,5,0))
              layout(matrix(c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2), nrow=3, byrow=TRUE))
              rainbowcols <- rainbow(ncol(LQ), s=0.75)
              plot(c(min(time_steps), max(time_steps)), c(min(LQ), max(LQ)), type="n", xlab="Time steps", ylab=paste("LQ (", units, ")", sep = ""), font.lab=2, cex.lab=1.5, cex.axis=1.5)
              for (i in 1:ncol(LQ)) {
                lines(time_steps, LQ[, i], col=rainbowcols[i])
              }
              par(mar=c(0,0,0,0))
              plot(0, 0, axes=F, type="n", xlab="", ylab="")
              legend("center", legend=prob, col=rainbowcols, lty=1, cex=1.5, bty="n")
            }

            # Reset plot settings
            par(def.par)
          })
