
conCovOpt <- function(x, outcome = NULL, 
                      type = if (inherits(x, "configTable")) attr(x, "type") else "cs", 
                      maxCombs = 1e7, approx = FALSE, allConCov = FALSE){
  eps = 1e-12
  ct <- configTable(x, type = type, rm.dup.factors = FALSE, rm.const.factors = FALSE)
  cti <- ctInfo(ct)
  responses <- cti$resp_nms
  if (is.null(outcome)){
    outcome <- responses
  } else {
    stopifnot(outcome %in% responses)
  }
  f <- attr(ct, "n")

  out <- vector("list", length(outcome))
  names(out) <- outcome
  for (outc in outcome){
    out[[outc]] <- .conCovOpt1outcome(ct, cti$scores, f, outc, eps = eps, 
                                      maxCombs = maxCombs, approx = approx, allConCov = allConCov)
  }
  attr(out, "configTable") <- ct
  class(out) <- "conCovOpt"
  out
}

.conCovOpt1outcome <- function(ct, sc, f, outcome, eps, maxCombs, approx = FALSE, allConCov = FALSE){
  type <- attr(ct, "type")
  x_df <- as.data.frame(ct)
  y <- sc[, outcome]  ## ctInfo() applied twice!!
  outcomeVar <- if (type == "mv") sub("=.+", "", outcome) else (outcome)
  
  # step 1: grouping wrt lhs-factors --------------------------------------
  ct_without_outcome <- configTable(x_df[-match(outcomeVar, names(ct))], type = type,
  																	rm.dup.factors = FALSE, rm.const.factors = FALSE,
  																	verbose = FALSE)
  sc <- ctInfo(ct_without_outcome)$scores
  noGroups <- nrow(ct_without_outcome) == nrow(ct)
  if (noGroups){
    y_g <- y
    g_freqs <- rep(1L, length(y))
    # u_exoGroups <- seq_along(g_freqs)
    yUnique <- rep(TRUE, nrow(ct_without_outcome))
  } else {
    cases_grouped <- unname(attr(ct_without_outcome, "cases"))
    g_freqs <- attr(ct_without_outcome, "n")
    u_exoGroups <- match(as.integer(unlist(cases_grouped)),
                            seq_along(f))
    exoGroups <- C_relist_Int(u_exoGroups, g_freqs)
    y_g <- relist1(y[u_exoGroups], g_freqs)
    yUnique <- vapply(y_g, dplyr::n_distinct, integer(1)) == 1L
  }

  # step 2: possible lhs-values --------------------------------------
  y1 <- rep(NA, nrow(ct_without_outcome))
  y1[yUnique] <- if (noGroups){
    y 
  } else {
    vapply(y_g[yUnique], "[", 1, FUN.VALUE = y_g[[c(1, 1)]])
  }
  yMatch <- rowAnys(sc == y1)
  yMax <- y1 >= rowMaxs(sc)
  yMin <- y1 <= rowMins(sc)
  #data.frame(yUnique, yMatch, yMax, yMin)
  
  obvious <- yUnique & rowAnys(as.matrix(data.frame(yMatch, yMax, yMin)))
  #cbind(ct_without_outcome, obvious)

  stopifnot(nrow(ct_without_outcome) == length(y_g))
  n1 <- length(y_g)
    
  possible <- vector("list", n1)
  possible[yUnique & yMatch] <- as.list(y_g[yUnique & yMatch])
  possible[yUnique & !yMatch & yMax] <- as.list(rowMaxs(sc)[yUnique & !yMatch & yMax])
  possible[yUnique & !yMatch & yMin] <- as.list(rowMins(sc)[yUnique & !yMatch & yMin])
  possible[lengths(possible)>1] <- lapply(possible[lengths(possible)>1], "[", 1L)
  sel <- vapply(possible, is.null, logical(1))
  if (any(sel) && type %in% c("cs", "mv")){
    possible[sel] <- rep_len(list(0:1), sum(sel))
  } else if (any(sel) && type == "fs"){
    yranges <- vapply(y_g[sel], range, y_g[[1]][c(1, 1)])
    sc_sorted <- t(sc[sel, , drop = FALSE])
    sc_sorted[] <- as.vector(sc_sorted)[order(as.vector(col(sc_sorted)), as.vector(sc_sorted))]
    selected <- which(sel)
    for (i in seq_len(sum(sel))){
      possible[[selected[i]]] <- getPossibleValues(unique.default(sc_sorted[, i]), yranges[, i])
    }
  }

  # Approximate mode: reduce the number of possible values considered
  if (approx){
    ff <- C_relist_Int(if (noGroups) f else f[u_exoGroups], g_freqs)
    wm <- function(x, n) median(rep(x, n))
    findClosest <- function(x, m){
      dist <- abs(x - m)
      x[dist == min(dist)]
    }
    wms <- mapply(wm, y_g, ff, SIMPLIFY = TRUE)
    possible <- mapply(findClosest, possible, wms, SIMPLIFY = FALSE)
  }
  
  #possible
  n_reprodList <- round(product(lengths(possible)))
  if (n_reprodList > maxCombs){
    warning(sprintf("Outcome %s: n_reprodList (=%4.2e) is larger than maxCombs (=%4.2e)",
                    outcome, n_reprodList, maxCombs), call. = FALSE)
    out <- data.frame(con = numeric(0), cov = numeric(0), id = integer(0))
    attr(out, "reprodList") <- possible
    attr(out, "exoGroups") <- if (noGroups) seq_along(g_freqs) else exoGroups
    return(out)
  }

  # step 3: preparing structures for expanding --------------------------------
  ff <- C_relist_Int(if (noGroups) f else f[u_exoGroups], g_freqs)
  s_f <- vapply(ff, sum, integer(1))
  Sx_base <- as.vector(vapply(possible, function(x) as.numeric(x[1]), numeric(1)) %*% s_f)
                         
  Sy <- sum(y*f)
  dx <- Map(function(x, f)(x-x[1])*f, possible, s_f)
  myfn <- function(v, yr, f){
    drop(f %*% (outer(yr, v, pmin) - v[1]))
    }
  dminxy <- Map(myfn, possible, y_g, ff)

  # step 4: expanding --------------------------------
  conCov_allCombs <- C_iterate(dx, dminxy, Sx_base, Sy)
  colnames(conCov_allCombs) <- c("con", "cov")

  # step 5: extract "best" con-cov combinations -----------------------------
  out <- getOptim(conCov_allCombs, eps = eps)
  attr(out, "reprodList") <- possible
  attr(out, "exoGroups") <- if (noGroups) seq_along(ff) else exoGroups
  if (allConCov) attr(out, "allConCov") <- conCov_allCombs
  out
}

# Aux functions
getPossibleValues <- function(possibleValues, yrange){
  csm <- sign(possibleValues - yrange[[1]]) + sign(possibleValues - yrange[[2]])
  from <- if (any(csm == -1)) possibleValues[max(which(csm == -1))] else if (any(csm == -2)) possibleValues[max(which(csm == -2))] else 0
  to <- if (any(csm == 1)) possibleValues[min(which(csm == 1))] else if (any(csm == 2)) possibleValues[min(which(csm == 2))] else 1
  possibleValues[possibleValues >= from & possibleValues <= to]
}  

getOptim <- function(x, eps = 1e-12){
  stopifnot(identical(ncol(x), 2L))
  x <- distinct(data.frame(x, id = seq_len(nrow(x))))
  n <- nrow(x)
  if (n <= 1L) return(x)
  ord <- order(-x[, 1], -x[, 2])
  x <- x[ord, , drop = FALSE]
  cum <- cummax(x[, 2])
  out <- c(TRUE, x[-n, 1] > x[-1, 1] & cum[-1] > cum[-n])
  x <- x[out, , drop = FALSE]
  if (any(eq <- (diff(x$cov) < eps))) x <- x[-(which(eq)+1), , drop = FALSE] # "Fuzzy Dedup" wrt cov
  if (any(eq <- (diff(x$con) > -eps))) x <- x[-which(eq), , drop = FALSE] # "Fuzzy Dedup" wrt con
  rownames(x) <- NULL
  x
}


# print method  
print.conCovOpt <- function(x, ...){
  cat("--- conCovOpt: optimal consistency-coverage pairs ---\n")
	for (outc in names(x)){
		cat("\nOutcome ", outc, ":\n", sep = "")
		print(x[[outc]], ...)
	}
	cat("\n")
	invisible(x)
}


# plot method  
plot.conCovOpt <- function(x, con = 1, cov = 1, ...){
  d <- do.call(rbind, 
    Map(function(cc, outc) data.frame(cc[c("con", "cov")], outcome = rep(outc, nrow(cc))), 
        x, names(x)))
  ggplot(d, aes_string("con", "cov", col = "outcome")) +
    geom_rect(data = data.frame(con, cov),
              aes_string(xmin = "con", ymin = "cov", xmax = 1, ymax = 1),
              col = "lightgray", alpha = 0.05) +
    geom_point() + geom_line()
}

