# deprecated methods

# custom plot function for `D*M` output
plot.DstarM <- function(x, what = "model", ...) {
  dots <- list(...)
  def.args <- list(bty = "n", xlab = "Reaction Time", las = 1, ylab = "density",
    x = x$tt, type = "b", lty = 1, pch = 1)
  if (!is.null(x$byPp)) {
    # output from byParticipant
    idx <- sapply(x, is.DstarM)
    dots$y <- switch(x$byPp, estDstarM = do.call(cbind, lapply(x[idx],
      `[[`, ifelse(what == "model", "modelDist", "g.hat"))), estND = do.call(cbind,
      lapply(x[idx], `[[`, "r.hat")), estObserved = do.call(cbind,
      lapply(x[idx], `[[`, "obsNorm")))
    dots$x <- x[[idx[1]]]$tt
    if (x$byPp == "estND") {
      def.args$xlim <- range(unlist(lapply(x[idx], `[[`, "ttr")))
    }
  } else if (!is.null(x$ttr)) {
    # output from r.hat()
    def.args$xlim <- range(x$ttr)
    dots$y <- x$r.hat
  } else if (names(x)[1] == "obsNorm") {
    # output from estObserved
    dots$y <- x$obsNorm
  } else if (names(x)[1] == "Bestvals") {
    # output from estDstarM()
    if (what == "model") {
      dots$y <- x$modelDist
    } else if (what == "data") {
      dots$y <- x$g.hat
    }
  } else {
    stop("No plot method available for this object.", call. = FALSE)
  }
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

print.DstarM <- function(x, na.print = "-", ...) {
  # exists for legacy purposes
  other <- FALSE
  if (is.list(x)) {
    # estDstarM, estND, estObserved
    nm <- names(x)[1]
    dots <- list(...)
    if (nm == "Bestvals") {
      # x from estDstarM(pars = ...)
      out <- x$Bestvals[c(x$restr.mat)]
      out[duplicated(c(x$restr.mat))] <- NA
      dim(out) <- dim(x$restr.mat)
      colnames(out) <- paste("Condition", seq_len(ncol(out)))
      rownames(out) <- sapply(strsplit(paste0(names(x$Bestvals)), "(?<=[a-zA-Z])(?=[0-9])",
        perl = TRUE), `[[`, 1)[1:nrow(out)]
      print(out, na.print = "-")
    } else if (nm == "objVals") {
      # x from estDstarM(pars = ...)
      if (isTRUE(dots$sum)) {
        cat("Sum of objective function values:\n")
        print(sum(x$objVals))
      } else {
        cat("Objective function values:\n")
        print(x$objVals)
      }
    } else if (nm == "r.hat") {
      print(x$descriptives)
    } else if (nm == "counts") {
      # x from rtDescriptives
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
      dots$x <- cbind(x$props[[idx]], x$counts[[idx]])
      dimnames(dots$x) <- list(1:dim(dots$x)[1L], switch(what, cr = rep(x$responses,
        2), c = c("Proportion", "Counts"), r = x$responses))
      xLen <- nchar(utils::capture.output(print(dots$x[, 1, drop = FALSE],
        digits = NULL)))[1]
      cat(c(" ", msg0, "\n"))
      if (idx == 1) {
        cat(c("  Prop", rep("", max(c(0, xLen - 7))), rep("", xLen *
          (dim(dots$x)[2L]/2 - 1)), "Counts\n"))
      }
      do.call(print, dots)
    } else {
      other <- TRUE
    }
  } else if (is.matrix(x)) {
    # x from calcIC
    out <- matrix(nrow = dim(x)[1L], ncol = dim(x)[2L] + 1)
    out[, 1] <- format(x[, 1], digits = 7)
    out[, c(2, 4)] <- x[, c(2, 4)]
    out[, 3] <- format(x[, 3], digits = 4)
    out[, 5] <- format(x[, 5], digits = 4)
    out[, 6] <- ifelse(x[, 5] < 0.001, "***", ifelse(x[, 5] < 0.01, "**",
      ifelse(x[, 5] < 0.05, "*", "n.s.")))
    rownames(out) <- rownames(x)
    colnames(out) <- c(colnames(x), "signif")
    out[is.na(x)] <- ""  # as which(is.na(x))
    print(out, na.print = "", quote = FALSE)
    writeLines("\nSignif: |p < .001: *** |p < 0.01: ** |p < .05: * |p > .05: n.s.")
  } else {
    # if all else fails
    other <- TRUE
  }
  if (other) {
    print.default(x)
  }
}
