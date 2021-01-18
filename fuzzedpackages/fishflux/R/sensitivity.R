#' A function to check the sensitivity of cnp_model predictions based on the variation of input parameters
#'
#' This function runs the cnp_model fixing all parameters SD's but one to test for sensitivity
#'
#' @param TL total length of a fish in cm
#' @param param list of all parameter means ("_m") and standard deviations ("_sd") Default parameters are set with very low sd's. See \link[fishflux]{cnp_model_mcmc}  for a list of all requested parameters
#' @param iter A positive integer specifying the number of iterations. The default is 1000
#' @param par Charachter vector specifying which input parameter sd's should be used for sensitivity.
#' @param out Charachter vector specifying which output parameter sd's should be returned.
#' @param ... Other arguments that can be used from \link[fishflux]{cnp_model_mcmc}
#' 
#' @details Returns a dataframe with sd's of model predictions. Row names indicate the variable, who's sd was used for the model run.
#' Plots a heatplot with width of the 95%CI of output predictions.
#' 
#' @keywords fish stoichiometry excretion mcmc sensitivity plot
#'
#' @importFrom parallel mcmapply
#' @importFrom tidyr gather
#' @importFrom dplyr %>% summarise left_join group_by mutate
#' @importFrom fishualize scale_fill_fish
#' @importFrom ggplot2 ggplot geom_tile aes geom_text labs theme_bw theme
#' 
#' @examples
#' library(fishflux)
#' sensitivity(TL = 10, param = list(k_sd = 0.2, Dn_sd = 0.2, Dc_sd = 0.1),
#'             par = c("k_sd","Dn_sd","Dc_sd"), out = c("Ic", "In", "Ip", "Gc"))
#' 
#' @export
sensitivity <- function(TL, param, iter = 1000, par,
                        out = c("Ic", "In", "Ip", "Gc",
                                "Gn", "Gp", "Fc", "Fn",
                                "Fp", "Wc", "Wn", "Wp"), ...) {

  #parameter SD's and means
  pm <- c("lt_m", "ac_m", "an_m", "ap_m", "Dc_m",
          "Dn_m", "Dp_m", "linf_m", "k_m", "t0_m",
          "theta_m", "r_m", "h_m", "lwa_m", "lwb_m",
          "mdw_m", "v_m", "F0nz_m", "F0pz_m", "Qc_m",
          "Qn_m", "Qp_m", "a_m", "f0_m")
  psd <- c("lt_sd", "ac_sd", "an_sd", "ap_sd", "Dc_sd",
           "Dn_sd", "Dp_sd", "linf_sd", "k_sd", "t0_sd",
           "theta_sd", "r_sd", "h_sd", "lwa_sd", "lwb_sd",
           "mdw_sd", "v_sd", "F0nz_sd", "F0pz_sd", "Qc_sd",
           "Qn_sd", "Qp_sd", "a_sd", "f0_sd")

  parm <- par[par %in% pm]
  parsd <- par[par %in% psd]

  param_m <- param[parm]
  param_sd <- param[parsd]

  sd_low <- 0.000000001

  param_sdl <- param_sd
  param_sdl[seq_len(length(parsd))] <- sd_low
  param_msdl <- append(param_m, param_sdl)

  #run cnp_model for all sd's with rest very low
  sd <- parsd
  res_sd <- as.data.frame(
    mcmapply(sd, FUN = function(x) {
    param_msdl[x] <- param_sd[x]
    mod <- cnp_model_mcmc(TL, param_msdl, iter, ...)$summary
    ext <- mod[match(out, mod$variable), "Q_97.5"] - mod[match(out, mod$variable), "Q_2.5"]
    ext
  }))

  row.names(res_sd) <- sapply(out, function(x) {
      rn <- paste(x, "_CI", sep = "")
      rn
    })

  res_sd <- as.data.frame(t(res_sd))

  #plot
  res <- res_sd
  res$input_sd <- row.names(res)
  res <- gather(res, "key", "value", - input_sd)
  sum <- summarise(group_by(res, key), sum = sum(value))
  res <- res %>% left_join(sum) %>%
    mutate(scale = value / sum)
  plot <- ggplot(res) +
    geom_tile(aes(x = key, y = input_sd, fill = scale)) +
    scale_fill_fish(option = "Trimma_lantana", end = 0.9) +
    geom_text(aes(x = key, y = input_sd,
                  label = formatC(value, format = "e", digits = 1))) +
    labs(x = "", y = "", fill = "Relative width 95% CI") +
    theme_bw() + theme(legend.position = "bottom")
  print(plot)
  res_sd
}
