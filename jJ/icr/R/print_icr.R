# Copyright (C) 2017  Alexander Staudt
#
# This file is part of icr.
#
# icr is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# icr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with icr; if not, see <https://www.gnu.org/licenses/>.


#' @importFrom stats quantile sd
#' @export
print.icr <- function(x, ..., level = 0.95) {
    if (level > 1 || level < 0) {
        stop("level must lie within interval [0,1]")
    }

    # create summary from icr-object
    sig <- 1 - level

    # string representation of metric
    string_metrics <- c("nominal", "ordinal", "interval", "ratio")
    string_metric <- string_metrics[x$metric]

    # coders, units, metric
    h <- data.frame(alpha = round(x$alpha, digits = 3),
                    coders = x$n_coders,
                    units = x$n_units,
                    level = string_metric)
    f_h <- format(h, digits = 3, justify = "right")

    # alpha, quantiles, standard errors
    if (x$bootstrap == TRUE) {
        ci_1 <- quantile(x$bootstraps, c(sig / 2, 1 - sig / 2), na.rm = TRUE)
        sd_1 <- sd(x$bootstraps, na.rm = TRUE)
        nboot <- x$nboot
        b_alpha <- mean(x$bootstraps, na.rm = TRUE)
        if (is.nan(b_alpha)) {
            b_alpha <- NA
        }
    } else {
        ci_1 <- c(NA, NA)
        sd_1 <- NA
        nboot <- NA
        b_alpha <- NA
    }
    if (x$bootnp == TRUE) {
        ci_2 <- quantile(x$bootstrapsNP, c(sig / 2, 1 - sig / 2), na.rm = TRUE)
        sd_2 <- sd(x$bootstrapsNP, na.rm = TRUE)
        nnp <- x$nnp
        b_alphaNP <- mean(x$bootstrapsNP, na.rm = TRUE)
        if (is.nan(b_alphaNP)) {
            b_alphaNP <- NA
        }
    } else {
        ci_2 <- c(NA, NA)
        sd_2 <- NA
        nnp <- NA
        b_alphaNP <- NA
    }

    ll <- format(sig / 2 * 100, scientific = FALSE)
    ul <- format((1 - sig / 2) * 100, scientific = FALSE)

    results <- data.frame(matrix(NA, nrow = 2, ncol = 6), check.names = FALSE)
    results[, 1] <- round(c(b_alpha, b_alphaNP), digits = 3)
    results[, 2] <- round(c(sd_1, sd_2), digits = 3)
    results[, 3] <- round(c(ci_1[1], ci_2[1]), digits = 3)
    results[, 4] <- round(c(ci_1[2], ci_2[2]), digits = 3)
    results[, 5] <- c("Krippendorff", "nonparametric")
    results[, 6] <- c(nboot, nnp)
    colnames(results) <- c("Alpha",
                           "Std. Error",
                           paste0(ll, " %"),
                           paste0(ul, " %"),
                           "Boot. technique",
                           "Bootstraps")
    f_results <- format(results, digits = 3, justify = "right")

    # alpha_min
    alpha_min <- data.frame(alpha_min = c(0.9, 0.8, 0.7, 0.67, 0.6, 0.5),
                            krippendorff = NA,
                            nonparametric = NA)
    if (length(x$bootstraps) > 1) {
        for (i in 1:6) {
            alpha_min[i, 2] <-
                round(sum(x$bootstraps > alpha_min[i, 1]) / x$nboot, digits = 3)
        }
    } else {
        alpha_min[, 2] <- NA
    }
    if (length(x$bootstrapsNP) > 1) {
        for (i in 1:6) {
            alpha_min[i, 3] <-
                round(sum(x$bootstrapsNP > alpha_min[i, 1]) / x$nnp, digits = 3)
        }
    } else {
        alpha_min[, 3] <- NA
    }

    f_alpha_min <- format(alpha_min, digits = 3, justify = "right")

    # print results
    cat("\n", " Krippendorff's alpha ", "\n\n", sep = "")
    print(f_h, row.names = FALSE)

    cat("\n")
    cat(" Bootstrapped alpha", "\n", sep = "")
    print(f_results, row.names = FALSE)

    cat("\n")
    cat(" P(alpha > alpha_min):\n")
    print(f_alpha_min, row.names = FALSE)
}
