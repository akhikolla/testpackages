#' Plotting the predicted survival function from an rocTree object
#'
#' Plot the predicted survival function from rocTree objects.
#'
#' @importFrom ggplot2 ggplot geom_step aes xlab ylab geom_smooth
#'
#' @noRd
#' @export
plot.predict.rocTree <- function(x, ...) {
    if (!is.predict.rocTree(x)) stop("Response must be a 'predict.rocTree' object")
    if (names(x$pred)[[2]] == "Survival") {   
        tmp <- data.frame(x = x$pred$Time, y = x$pred$Survival)
        gg <- ggplot(tmp, aes(x = x, y = y)) + geom_step(lwd = I(1.1)) +
            xlab("Time") + ylab("Survival probabilities")
    } 
    if (names(x$pred)[[2]] == "hazard") {
        tmp <- data.frame(x = x$pred$Time, y = x$pred$hazard)
        gg <- ggplot(tmp, aes(x = x, y = y)) + ## geom_step(lwd = I(1.1)) +
            xlab("Time") + ylab("Hazard estimate") +
            geom_smooth(lwd = I(1.1), colour="black", method = "loess", se = FALSE)
    }
    gg
}
