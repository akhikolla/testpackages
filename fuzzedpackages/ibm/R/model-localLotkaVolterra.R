# Main function -----------------------------------------------------------
#' @title Lotka-Volterra with local predation interactions 
#' @description This function simulates several trajectories for a 
#' Lotka-Volterra model with local predation interactions as decribed
#' in Brigatti et al. (2009).
#' 
#' @param par A list containing the parameters to run the model, currently
#' the growth rate of prey (r), the mortality rate of predator (l), predation
#' interaction parameters (alpha and beta), diffusion rates (D), diameters of 
#' local interaction (L) and initial population size (initial). For D, L and
#' initial population, a list with two values (named N and P) is required.
#' @param T Time horizon, number of time steps to be simulated. 
#' @param replicates Number of replicates (trajectories) to be simulated.
#' @param  dim Spatial dimension for the space. Can be 1, 2 or 3.
#' @param periodic Spatial boundary conditions. If \code{periodic} is set to 
#' \code{TRUE}, the space is a torus. If set to \code{FALSE}, the boundaries 
#' are reflective.
#' @param spatial Boolean, should spatial outputs (position of individuals) to 
#' be saved?
#' @param verbose Boolean, to print population sizes by step?
#' @param maxpop Maximum population size. If predator or prey population size
#' get bigger, the simulation ends.
#' @return A list with the following elements:
#' \item{N}{A matrix with prey population sizes by time (rows) and replicates (columns)}
#' \item{P}{A matrix with predator population sizes by time (rows) and replicates (columns)} 
#' \item{pop}{Prey and predator positions by time, if \code{spatial} is \code{TRUE}} 
#' @author Ricardo Oliveros--Ramos
#' @references Brigatti et al. 2009.
#' @keywords Lotka-volterra local interactions
#' @examples
#' \dontrun{
#' set.seed(880820)
#' par = list(alpha=5e-4, beta=5e-4, r=0.1, m=0.05, D=list(N=8e-5, P=8e-5), 
#' L=list(N=0.2, P=0.2))
#' N0 = with(par, m/(2*beta*L$P))
#' P0 = with(par, r/(2*alpha*L$N))
#' par$initial = list(N=round(N0), P=round(P0))
#' sim = localLotkaVolterra(par, T=240, replicates=100, maxpop = 1e4)
#' plot(sim)
#' } 
#' @export
localLotkaVolterra = function(par, T, replicates=1, dim=1, periodic=TRUE, 
                              spatial=FALSE, verbose=FALSE, maxpop=1e6) {
  
  if(isTRUE(spatial) & replicates!=1)
    stop("Spatial outputs only available for 1 replicate.")
  
  N = array(dim=c(T+1, replicates))
  P = array(dim=c(T+1, replicates))
  
  for(irep in seq_len(replicates)) {
    xtime = Sys.time()
    sim = .localLotkaVolterra(par=par, T=T, dim=dim, periodic=periodic, 
                              spatial=spatial, verbose=FALSE, maxpop=1e6) 
    xtime = c(Sys.time() - xtime)
    N[, irep] = sim$N  
    P[, irep] = sim$P  
    cat(sprintf("Replicate %d - Ellapsed %0.2fs\n", irep, xtime))
    
  }
  
  output = list(N=N, P=P)
  class(output) = c("ibm.LLV", class(output))
  return(output)
}

# One trajectory ----------------------------------------------------------
.localLotkaVolterra = function(par, T, dim=1, periodic=TRUE, spatial=FALSE, 
                               verbose=FALSE, maxpop=1e6) {
  

  par$sd = lapply(par$D, FUN = function(x) sqrt(2*x))
  
  # Initializing population vectors
  
  N = array(dim=T+1)
  P = array(dim=T+1)
  pop = NULL
  
  N[1] = par$initial$N
  P[1] = par$initial$P
  
  # Initializing individuals positions
  tpop = list()
  tpop$N = initializePopulation(N=par$initial$N, n=dim, maxpop=maxpop)
  tpop$P = initializePopulation(N=par$initial$P, n=dim, maxpop=maxpop)

  if(isTRUE(spatial)) {
    pop = array(dim=c(maxpop, dim, 2, T+1))
    pop[, , 1, 1] = tpop$N
    pop[, , 2, 1] = tpop$P
  }
  
  for(t in 1:T) {
   
    # diffusion
    tpop = diffusion(tpop, sd=par$sd, N=par$initial)
    tpop = boundaries(tpop, periodic=periodic)
    
    # (local) predation
    predation = localPredationInteractions(pop=tpop, R=par$L, N=N[t], P=P[t])
    
    newN  = reproduction(N[t], rates=par$r)    
    newP  = reproduction(P[t], rates=par$beta*predation$PN)
    
    survN = mortality(N[t], rates=par$alpha*predation$NP, survivors=TRUE)
    survP = mortality(P[t], rates=par$m, survivors=TRUE)
  
    # saving outputs
    N[t+1] = length(survN) + length(newN) 
    P[t+1] = length(survP) + length(newP)
    
    if(all(N[t+1]==0, P[t+1]==0)) break
    if(any(N[t+1]>maxpop, P[t+1]>maxpop)) break
    
    if(N[t+1]>nrow(tpop$N)) tpop$N = updateMatrixSize(x=tpop$N, n=N[t+1], max=maxpop)
    if(P[t+1]>nrow(tpop$P)) tpop$P = updateMatrixSize(x=tpop$P, n=P[t+1], max=maxpop)
    
    tpop$N[seq_len(N[t+1]), ] = tpop$N[c(survN, newN), ]
    tpop$P[seq_len(P[t+1]), ] = tpop$P[c(survP, newP), ]
   
    if(N[t+1]<N[t]) tpop$N[seq(N[t+1]+1, N[t]), ] = NA
    if(P[t+1]<P[t]) tpop$P[seq(P[t+1]+1, P[t]), ] = NA
      
    if(isTRUE(spatial)) {
      pop[, , 1, t+1] = tpop$N
      pop[, , 2, t+1] = tpop$P
    }
    
    if(isTRUE(verbose)) {
      cat("t=", t, ", N=", N[t+1], ", P=", P[t+1],"\n")
    }
    
  }
  
  xlim = max(N, P, na.rm=TRUE)
   
  output = list(N = as.numeric(N), P = as.numeric(P), pop=pop[xlim, , , ])
  return(output)
  
}


# methods -----------------------------------------------------------------

#' @export
plot.ibm.LLV = function(x, alpha=0.95, nmax=10, ...) {
  opar = par(no.readonly = TRUE)
  par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(4,4,1,4))
  prey = .summary(x$N, alpha=alpha, nmax=nmax)
  pred = .summary(x$P, alpha=alpha, nmax=nmax)
  # prey
  plot.new()
  plot.window(xlim=prey$xlim, ylim=prey$ylim)
  matplot(prey$rep, col="grey", type="l", lty=1, add=TRUE)
  lines(prey$median, col="red", lwd=2)
  lines(prey$ll, col="red", lty=3)
  lines(prey$ul, col="red", lty=3)
  axis(2)
  box()
  mtext("PREY", 2, line=3)
  # pred
  plot.new()
  plot.window(xlim=pred$xlim, ylim=pred$ylim)
  matlines(pred$rep, col="grey", lty=1)
  lines(pred$median, col="blue", lwd=2)
  lines(pred$ll, col="blue", lty=3)
  lines(pred$ul, col="blue", lty=3) 
  axis(4)
  axis(1)
  box()
  mtext("PREDATOR", 4, line=3)
  
  par(opar)
  
  return(invisible())
  
}



