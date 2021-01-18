#' A function to extract specific model output parameters from result
#'
#' @param mod Output from cnp_mod_mcmc()
#' @param par Character vector specifying which output parameter that should be returned.
#' 
#' @return Main model output parameters:
#' \itemize{
#' \item{F0c:} C-specific minimal inorganic flux (g/day)
#' \item{F0n:} N-specific minimal inorganic flux (g/day)
#' \item{F0p:} P-specific minimal inorganic flux (g/day)
#' \item{Gc:} Carbon-specific growth rate (g/day)
#' \item{Gn:} Nitrogen-specific growth rate (g/day)
#' \item{Gp:} Phosphorus-specific growth rate (g/day)
#' \item{Sc:} C-specific minimal supply rate (g/day)
#' \item{Sn:} N-specific minimal supply rate (g/day)
#' \item{Sp:} P-specific minimal supply rate (g/day)
#' \item{Ic:} Ingestion rate of C (g/day)
#' \item{In:} Ingestion rate of N (g/day)
#' \item{Ip:} Ingestion rate of P (g/day)
#' \item{Wc:} Egestion rate of C (g/day)
#' \item{Wn:} Egestion rate of N (g/day)
#' \item{Wp:} Egestion rate of P (g/day)
#' \item{Fc:} Total inorganic flux of C (respiration) (g/day)
#' \item{Fn:} Total inorganic flux of N (excretion) (g/day)
#' \item{Fp:} Total inorganic flux of P (excretion) (g/day)
#' }
#' 
#' @details Returns a data.frame with a summary of the selected output parameters
#' 
#' @keywords fish stoichiometry excretion mcmc
#'
#' @importFrom dplyr filter
#' 
#' @examples
#' model <- cnp_model_mcmc(TL = 5:10, param = list(Qc_m = 40, Qn_m = 10, Qp_m = 4))
#' extract(model, c("Fn","Fp"))
#'
#' @export
extract <- function(mod, par) {
  summary <- mod$summary
  TL <- unique(summary$TL)
  list <- lapply(par, FUN = function(p) {
    sub <- filter(summary, variable == p)
    sub <- data.frame(sub$mean, sub$median, sub$sd,
                      sub$`Q_2.5`, sub$`Q_97.5`, sub$`Q_25`,
                      sub$`Q_75`)
    colnames(sub) <- c(paste(p, "mean", sep = "_"),
                       paste(p, "median", sep = "_"),
                       paste(p, "sd", sep = "_"),
                       paste(p, "2.5%", sep = "_"),
                       paste(p, "97.5%", sep = "_"),
                       paste(p, "25%", sep = "_"),
                       paste(p, "75%", sep = "_"))
    sub
  })

  df <- do.call(cbind, list)
  tl <- data.frame(TL = TL)
  cbind(tl, df)
}
