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

#' @importFrom stats quantile density
#' @importFrom graphics plot lines polygon title legend
#' @export
plot.icr <- function(x, ..., level = 0.95, return_data = FALSE) {
    if (x$bootstrap == FALSE & x$bootnp == FALSE) {
        stop("No bootstrapped values of alpha available")
    }

    # create data frame for plotting
    density_kh <- NA
    density_np <- NA
    df_kh <- data.frame(x = NULL, y = NULL, type = NULL)
    df_np <- data.frame(x = NULL, y = NULL, type = NULL)

    if (x$bootstrap == TRUE) {
        density_kh <- density(x$bootstraps)
        ci_kh <- quantile(x$bootstraps,
                          c((1 - level) / 2, 1 - (1 - level) / 2))
        df_density_kh <- data.frame(x = density_kh$x,
                                    y = density_kh$y,
                                    ci = FALSE,
                                    ci_limit = FALSE,
                                    type = "Krippendorff",
                                    stringsAsFactors = FALSE)
        df_ci_kh <- data.frame(x = ci_kh,
                               y = 0,
                               ci = TRUE,
                               ci_limit = TRUE,
                               type = "Krippendorff",
                               stringsAsFactors = FALSE)
        df_kh <- rbind(df_density_kh, df_ci_kh)
        df_kh <- df_kh[order(df_kh$x), ]
        df_kh$ci[df_kh$x >= ci_kh[1] & df_kh$x <= ci_kh[2]] <- TRUE
    }

    if (x$bootnp == TRUE) {
        density_np <- density(x$bootstrapsNP)
        ci_np <- quantile(x$bootstrapsNP,
                          c((1 - level) / 2, 1 - (1 - level) / 2))
        df_density_np <- data.frame(x = density_np$x,
                                    y = density_np$y,
                                    ci = FALSE,
                                    ci_limit = FALSE,
                                    type = "nonparametric",
                                    stringsAsFactors = FALSE)
        df_ci_np <- data.frame(x = ci_np,
                               y = 0,
                               ci = TRUE,
                               ci_limit = TRUE,
                               type = "nonparametric",
                               stringsAsFactors = FALSE)
        df_np <- rbind(df_density_np, df_ci_np)
        df_np <- df_np[order(df_np$x), ]
        df_np$ci[df_np$x >= ci_np[1] & df_np$x <= ci_np[2]] <- TRUE
    }

    df <- rbind(df_kh, df_np)
    x_range <- range(df$x)
    y_range <- range(df$y)

    # plot
    plot(x = x_range, y = y_range, type = "n", xlab = "value", ylab = "density")

    lines(x = df$x[df$type == "Krippendorff" & df$ci_limit == FALSE],
          y = df$y[df$type == "Krippendorff" & df$ci_limit == FALSE],
          col = "blue")
    lines(x = df$x[df$type == "nonparametric" & df$ci_limit == FALSE],
          y = df$y[df$type == "nonparametric" & df$ci_limit == FALSE],
          col = "red")

    polygon(x = df$x[df$ci == TRUE & df$type == "Krippendorff"],
            y = df$y[df$ci == TRUE & df$type == "Krippendorff"],
            col = "blue", density = 30, angle = 30, border = FALSE)
    polygon(x = df$x[df$ci == TRUE & df$type == "nonparametric"],
            y = df$y[df$ci == TRUE & df$type == "nonparametric"],
            col = "red", density = 30, angle = 150, border = FALSE)

    if (length(unique(df$type)) == 2) {
        title(expression(paste("Bootstrapped ", alpha)))
        legend(x = 0, y = y_range[2],
               legend = c("Krippendorff", "nonparametric"),
               col = c("blue", "red"), lty = c(1:1),
               box.lty = 0, cex = 0.8)
    } else if (unique(df$type) == "Krippendorff") {
        title(expression(paste("Bootstrapped ", alpha, " (Krippendorff)")))
    } else if (unique(df$type) == "nonparametric") {
        title(expression(paste("Bootstrapped ", alpha, " (nonparametric)")))
    }

    # return plot data frame
    if (return_data == TRUE) {
        return(df)
    }
}
