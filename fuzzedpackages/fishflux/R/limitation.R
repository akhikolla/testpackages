#' A function to evaluate element limitation of the model
#'
#' This function allows you extract the proportions of the iterations for which c, n and p are the limiting element in the model.
#' 
#' @param mod Model output from cnp_model_mcmc().
#' @param plot Argument to specify if results should be shown in a plot.
#' 
#' @details Returns a data frame with:
#' \describe{
#'   \item{tl}{Total length, in cm}
#'   \item{nutrient}{c, n or p}
#'   \item{prop_lim}{the proportion of iterations for which there is limitation by the element}
#' }
#' 
#' @keywords fish plot limitation
#' 
#' @importFrom dplyr %>% mutate bind_rows
#' @importFrom tidyr gather
#' @importFrom rstan extract
#' @importFrom fishualize scale_color_fish_d
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_bw theme
#' 
#' @examples
#' library(fishflux)
#' mod <- cnp_model_mcmc(TL = 5, param = list(Qc_m = 40, Qn_m = 10, Qp_m = 4,
#'                                            Dc_sd = 0.1, Dn_sd = 0.05, Dp_sd = 0.05))
#' limitation(mod)
#' 
#' @export
limitation <- function(mod, plot = TRUE) {

  if (length(unique(mod$summary$TL)) == 1) {
    ee <- rstan::extract(mod$stanfit,"lim")[[1]]
    c <- length(which(ee == 1)) / length(ee)
    n <- length(which(ee == 2)) / length(ee)
    p <- length(which(ee == 3)) / length(ee)
    lim <- (data.frame(c = c,
                       n = n,
                       p = p)) %>%
      mutate(tl = unique(mod$summary$TL)) %>%
      gather("nutrient", "prop_lim", - tl)
  } else {
    lim <- lapply(mod$stanfit, function(x) {
      ee <- rstan::extract(x, "lim")[[1]]
      c <- length(which(ee == 1)) / length(ee)
      n <- length(which(ee == 2)) / length(ee)
      p <- length(which(ee == 3)) / length(ee)
      data.frame(c = c,
                 n = n,
                 p = p)
    }) %>%
      bind_rows() %>%
      mutate(tl = unique(mod$summary$TL)) %>%
      gather("nutrient", "prop_lim", - tl)
  }

  if (plot) {
    p <- ggplot(lim) +
      geom_point(aes(x = tl, y = prop_lim, color = nutrient)) +
      geom_line(aes(x = tl, y = prop_lim, color = nutrient)) +
      scale_color_fish_d(option = "Hypsypops_rubicundus") +
      labs(x = "TL (cm)", y = "Proportion of iterations", color = "Limiting element") +
      theme_bw() +
      theme(legend.position = "bottom")
    print(p)
  }
  lim
}
