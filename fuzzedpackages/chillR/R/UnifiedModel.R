#' ChuineCF
#'
#' chilling and forcing response function for
#' the unified model by Chuine
#'
#' @param x temperature
#' @param a numeric. paramter
#' @param b numeric. paramter
#' @param c numeric. paramter
#'
#' @references
#' Isabelle Chuine, A Unified Model for Budburst of Trees, J. theor. Biol. (2000) 207
#' 
#' @return
#' Returns a numeric vector.
ChuineCF <- function(x, a, b, c) {
  return(1./(1 + exp(a*(x - c)^2 + b * (x - c))))
}

#' ChuineFstar
#'
#' Critical forcing value
#'
#' @param Ctot numeric. total state of chilling
#' @param k numeric < 0. parameter.
#' @param w numeric > 0. parameter.
#' 
#' @references
#' Isabelle Chuine, A Unified Model for Budburst of Trees, J. theor. Biol. (2000) 207
#' 
#' @return
#' Returns a numeric vector.
ChuineFstar <- function(Ctot, w, k) {
  return(w*exp(k*Ctot))
}

#' UnifiedModel_Wrapper
#'
#' @param x data.frame with at least columns `Temp` and `JDay`
#' @param par numeric vector of length 9 with the parameters of the
#' unified model: 1. ac, 2. bc, 3. cc, 4. bf, 5. cf, 6. w, 7. k,
#' 8. Cstar and 9. tc.
#' 
#' @references
#' Isabelle Chuine, A Unified Model for Budburst of Trees, J. theor. Biol. (2000) 207
#'
#' @author Carsten Urbach <urbach@hiskp.uni-bonn.de>
#' @return
#' A single numeric value with the JDay prediction for the
#' temperaturs in `x$Temp` and the Unified Model parameters
#' in `par`.
#'
#' @export
UnifiedModel_Wrapper <- function(x, par) {
  tend <- length(x$Temp)
  if(tend < par[9] || par[6] <= 0 || par[7] >= 0) return(NA)
  chilling  <- cumsum(ChuineCF(x=x$Temp, a=par[1], b=par[2], c=par[3]))
  if(chilling[length(chilling)] < par[8]) return(NA)
  t1 <- min(which(chilling >= par[8]))
  Ctot <- sum(ChuineCF(x=x$Temp[c(1:round(par[9]))], a=par[1], b=par[2], c=par[3]))
  Fstar <- ChuineFstar(Ctot=Ctot, w=par[6], k=par[7])
  forcing <- ChuineCF(x=x$Temp[c(t1:tend)], a=0, b=par[4], c=par[5])
  cs <- cumsum(forcing)
  if(cs[length(cs)] < Fstar) return(NA)
  return(x$JDay[min(which(cs >= Fstar))])
}


#' UniChill_Wrapper
#'
#' @param x data.frame with at least columns `Temp` and `JDay`
#' @param par numeric vector of length 7 with the parameters of the
#' UniChill model: 1. ac, 2. bc, 3. cc, 4. bf, 5. cf, 6. Cstar and 7. Fstar.
#' 
#' @references
#' Isabelle Chuine, A Unified Model for Budburst of Trees, J. theor. Biol. (2000) 207
#'
#' @author Carsten Urbach <urbach@hiskp.uni-bonn.de>
#' @return
#' A single numeric value with the JDay prediction for the
#' temperaturs in `x$Temp` and the model parameters
#' in `par`.
#'
#' @export
UniChill_Wrapper <- function(x, par) {
  tend <- length(x$Temp)
  chilling  <- cumsum(ChuineCF(x=x$Temp, a=par[1], b=par[2], c=par[3]))
  if(chilling[length(chilling)] < par[6]) return(NA)
  t1 <- min(which(chilling >= par[6]))
  forcing <- ChuineCF(x=x$Temp[c(t1:tend)], a=0, b=par[4], c=par[5])
  cs <- cumsum(forcing)
  if(cs[length(cs)] < par[7]) return(NA)
  return(x$JDay[min(which(cs >= par[7]))])
}

#' StepChill_Wrapper
#'
#' Same as UniChill_Wrapper, but with a step function for chilling
#' 
#' @param x data.frame with at least columns `Temp` and `JDay`
#' @param par numeric vector of length 7 with the parameters of the
#' StepChill model: 1. Tc, 2. bf, 3. cf, 4. Cstar and 5. Fstar.
#' 
#' @references
#' Isabelle Chuine, A Unified Model for Budburst of Trees, J. theor. Biol. (2000) 207
#' 
#' Asse et al., Process-based models outcompete correlative models in projecting spring
#' phenology of trees in a future warmer climate,Agricultural and Forest Meteorology (2020) 107913
#'
#' @author Carsten Urbach <urbach@hiskp.uni-bonn.de>
#' @return
#' A single numeric value with the JDay prediction for the
#' temperaturs in `x$Temp` and the model parameters
#' in `par`.
#'
#' @export
StepChill_Wrapper <- function(x, par) {
  tend <- length(x$Temp)
  tmp <- rep(0, times=tend)
  tmp[which(x$Temp <= par[1])]  <- 1
  chilling  <- cumsum(tmp)
  if(chilling[length(chilling)] < par[4]) return(NA)
  t1 <- min(which(chilling >= par[4]))
  forcing <- ChuineCF(x=x$Temp[c(t1:tend)], a=0, b=par[2], c=par[3])
  cs <- cumsum(forcing)
  if(cs[length(cs)] < par[5]) return(NA)
  return(x$JDay[min(which(cs >= par[5]))])
}

#' UniForce_Wrapper
#'
#' @param x data.frame with at least columns `Temp` and `JDay`
#' @param par numeric vector of length 4 with the parameters of the
#' UniForce model: 1. bf, 2. cf, 3. Fstar, 4. t1.
#' 
#' @references
#' Isabelle Chuine, A Unified Model for Budburst of Trees, J. theor. Biol. (2000) 207
#'
#' @author Carsten Urbach <urbach@hiskp.uni-bonn.de>
#' @return
#' A single numeric value with the JDay prediction for the
#' temperaturs in `x$Temp` and the Unified Model parameters
#' in `par`.
#'
#' @export
UniForce_Wrapper <- function(x, par) {
  tend <- length(x$Temp)
  t1 <- round(par[4])
  if(t1 >= tend) return(NA)
  forcing <- ChuineCF(x=x$Temp[c(t1:tend)], a=0, b=par[1], c=par[2])
  cs <- cumsum(forcing)
  if(cs[length(cs)] < par[3]) return(NA)
  return(x$JDay[min(which(cs >= par[3]))])
}
