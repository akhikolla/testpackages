#' A function to plot results model
#'
#' This function allows you to plot an overview of the model results in function of the total length of fish
#' 
#' @param mod Model output from cnp_model_mcmc()
#' @param x Variable to be put on x-axis, "biomass" or "tl"
#' @param y Output variable(s) to be plotted. Can be a character or a character vector.
#' @param probs Width of the confidence
#' 
#' @keywords fish plot bioenergetic model stoichiometry
#'
#' @importFrom rstan extract
#' @importFrom dplyr %>% bind_rows mutate group_by
#' @importFrom tidyr gather
#' @importFrom fishualize scale_color_fish_d scale_fill_fish_d
#' @importFrom ggplot2 ggplot aes theme_bw labs scale_fill_brewer
#' @importFrom tidybayes stat_lineribbon
#' 
#' @examples
#' library(fishflux)
#' mod <- cnp_model_mcmc(TL = 5:15, param = list(
#'           Qc_m = 40, Qn_m = 10, Qp_m = 4, Dn_sd = 0.05))
#' plot_cnp(mod = mod, y = c("Fp", "Gp", "Wp", "Ip"),
#'          x = "tl", probs = c(0.5, 0.8))
#' plot_cnp(mod = mod, y = "Fp", x = "tl",
#'          probs = c(0.5, 0.8, 0.95))
#' 
#' @export
plot_cnp <- function(mod, y, x = "tl", probs = c(0.8, 0.95)) {
  vars <- c("lt", "lwa", "lwb", y)

  iter <- lapply(mod$stanfit, FUN = function(x, vars) {
    rstan::extract(x, vars)
  }, vars)
  iter <- iter %>%
    lapply( FUN = get_iter) %>%
    bind_rows() %>%
    mutate(lt = round(lt))

  iter$w <- mean(iter$lwa) * iter$lt^mean(iter$lwb)

  if (length(y) > 1) {
    iter_t <- gather(iter, "output", "value", y)
    if (x == "biomass") {
      plot <-
        ggplot(group_by(iter_t, iter),
               aes(x = w, y = value, color = output,
                   fill = output)) +
        stat_lineribbon(alpha = 0.4, show.legend = FALSE, .width = probs) +
        scale_color_fish_d(option = "Scarus_quoyi") +
        scale_fill_fish_d(option = "Scarus_quoyi") +
        theme_bw() +
        labs(x = "Biomass (g)", y = "Output (g/day)", fill = "output", color = "output")
    } else if (x == "tl") {
      plot <-
        ggplot(group_by(iter_t, iter),
               aes(x = lt, y = value, color = output,
                   fill = output)) +
        stat_lineribbon(alpha = 0.4, .width = probs) +
        scale_color_fish_d(option = "Scarus_quoyi") +
        scale_fill_fish_d(option = "Scarus_quoyi") +
        theme_bw() +
        labs(x = "Total length (cm)", y = "Output (g/day)", fill = "output", color = "output")
    }
  } else if (length(y) == 1) {
    colnames(iter) <- c("lt", "lwa", "lwb", "value", "iter", "w")
    if (x == "biomass") {
      plot <-
        ggplot(group_by(iter, iter),
               aes(x = w, y = value)) +
        stat_lineribbon(alpha = 0.9, show.legend = FALSE, .width = probs) +
        scale_fill_brewer() +
        theme_bw() +
        labs(x = "Biomass (g)", y = paste(y, "(g/day)", sep = " "))
    } else if (x == "tl") {
      plot <-
        ggplot(group_by(iter, iter),
               aes(x = lt, y = value)) +
        stat_lineribbon(alpha = 0.9,  .width = probs ) +
        scale_fill_brewer() +
        theme_bw() +
        labs(x = "Total length (cm)", y = paste(y, "(g/day)", sep = " "))
    }
  }
 plot
}

#' get_iter
#' @param x something
#' @importFrom plyr ldply
get_iter <- function(x) {
  get <- t(ldply(x))
  colnames(get) <- get[1, ]
  get <- data.frame(apply(get[-1, ], 2, as.numeric))
  get$iter <- seq_len(nrow(get))
  get
}
