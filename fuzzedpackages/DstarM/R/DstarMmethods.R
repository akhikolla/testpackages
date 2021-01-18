# custom print function for `D*M` output
#' @export
print.DstarM.fitD <- function(x, na.print = "-", ...) {

  dots <- list(...)
  out <- x$Bestvals[c(x$restr.mat)]
  out[duplicated(c(x$restr.mat))] <- NA
  dim(out) <- dim(x$restr.mat)
  colnames(out) <- paste("Condition", seq_len(ncol(out)))
  rmns <- gsub("[0-9]", "", names(x$Bestvals))
  rmns <- rmns[!duplicated(rmns)]
  rownames(out) <- rmns
  print(out, na.print = "-")
}

#' @export
print.DstarM.fitND <- function(x, na.print = "-", ...) {
  print(x$descriptives, na.print = na.print, ...)
}

#' @export
print.DstarM.fitObs <- function(x, na.print = "-", ...) {
  print.default(x, na.print = na.print, ...)
}

#' @export
print.DstarM.Descriptives <- function(x, na.print = "-", ...) {

  # nm = names(x)[1]
  dots <- list(...)
  what <- dots$what
  if (is.null(what)) {
    what <- "cr"
  } else {
    if (!(what %in% c("cr", "c", "r"))) {
      stop(sprintf("Argument what ('%s') must be 'cr', 'c', or 'r'.",
        what))
    }
    dots$what <- NULL
  }
  idx <- switch(what, cr = 1, c = 2, r = 3)
  msg0 <- switch(what, cr = "Condition Response Pairs", c = "Conditions",
    r = "Responses")
  dots$x <- cbind(x$table$props[[idx]], x$table$counts[[idx]])
  dimnames(dots$x) <- list(1:dim(dots$x)[1L], switch(what, cr = rep(x$table$responses,
    2), c = c("Proportion", "Counts"), r = x$table$responses))
  xLen <- nchar(utils::capture.output(print(dots$x[, 1, drop = FALSE],
    digits = NULL)))[1]
  cat(c(" ", msg0, "\n"))
  if (idx == 1) {
    cat(c("  Prop", rep("", max(c(0, xLen - 7))), rep("", xLen * (dim(dots$x)[2L]/2 -
      1)), "Counts\n"))
  }
  do.call(print, dots)
  print(x$graph)
}

#' @export
summary.DstarM.fitD <- function(object, ...) {
  if (names(object)[1L] != "Bestvals") {
    stop("No summary method available for this object.")
  }
  if (!is.null(object$note)) {
    cat(object$note)
  }
  cat("\n")
  cat(sprintf("%s analysis done on %g observations in %s conditions.",
    ifelse(object$DstarM, "D*M", "Chi-Square"), sum(object$n), object$ncondition))
  cat("\n")
  cat("Coefficients: \n \n")
  print(object)
  cat("\n")
  if (!length(object$fixed)) {
    cat("No parameters were fixed. ")
  } else {
    if (dim(object$fixed$fixedMat)[2L] == 1) {
      cat(sprintf("Parameter %s was fixed to %s. ", object$fixed$fixedMat[1L,
        ], object$fixed$fixedMat[2, ]))
    } else if (dim(object$fixed$fixedMat)[2L] == 2L) {
      cat(sprintf("Parameters %s were fixed to %s respectively. ",
        paste(object$fixed$fixedMat[1L, 1L], "and", object$fixed$fixedMat[1L,
          2L]), paste(object$fixed$fixedMat[2L, 1L], "and", object$fixed$fixedMat[2L,
          2L])))
    } else {
      cat(sprintf("Parameters %s were fixed to %s respectively. ",
        paste0(paste(object$fixed$fixedMat[1L, -dim(object$fixed$fixedMat)[2L]],
          collapse = ", "), ", and ", object$fixed$fixedMat[1L, dim(object$fixed$fixedMat)[2L]]),
        paste0(paste(object$fixed$fixedMat[2L, -dim(object$fixed$fixedMat)[2L]],
          collapse = ", "), ", and ", object$fixed$fixedMat[2L, dim(object$fixed$fixedMat)[2L]])))
    }
  }
  cat(sprintf("The value of the objective function was %g after %g iterations.",
    object$GlobalOptimizer$optim$bestval, object$Debug$niter))
  if (object$Debug$niter == object$Debug$itermax) {
    cat("\n")
    cat("WARNING: Maximum iterations was reached. Be aware of potential convergence issues.")
    cat("\n")
  }
}

simplePlot <- function(dots, def.args) {
  # helper function for plot methods

  ind <- unlist(lapply(dots[names(def.args)], is.null))
  dots[names(def.args)[ind]] <- def.args[ind]
  do.call(graphics::matplot, dots)
  if (!is.null(colnames(dots$y))) {

    nc <- NCOL(dots$y)
    if (is.null(dots$lty))
      dots$lty <- rep(1, nc)

    if (is.null(dots$col))
      dots$col <- seq_len(nc)
    graphics::legend("topright", colnames(dots$y), col = dots$col, lty = dots$lty,
      bty = "n")

  }
}

#' @export
plot.DstarM.fitD <- function(x, y = NULL, what = c("model", "data"), ...) {

  dots <- list(...)
  what <- match.arg(what)

  def.args <- list(bty = "n", xlab = "Reaction Time", las = 1, ylab = "density",
    x = x$tt, type = "b", lty = 1, pch = 1)
  if (what == "model") {
    dots$y <- x$modelDist
  } else if (what == "data") {
    dots$y <- x$g.hat
  }
  simplePlot(dots, def.args)
}

#' @export
plot.DstarM.fitND <- function(x, y = NULL, ...) {
  dots <- list(...)
  def.args <- list(bty = "n", xlab = "Reaction Time", las = 1, ylab = "density",
    x = x$tt, type = "b", lty = 1, pch = 1)
  def.args$xlim <- range(x$ttr)
  dots$y <- x$r.hat
  simplePlot(dots, def.args)
}

#' @export
plot.DstarM.fitObs <- function(x, y = NULL, ...) {
  dots <- list(...)
  def.args <- list(bty = "n", xlab = "Reaction Time", las = 1, ylab = "density",
    x = x$tt, type = "b", lty = 1, pch = 1)
  dots$y <- x$obsNorm
  simplePlot(dots, def.args)
}
