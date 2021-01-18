
# selectMax
selectMax <- function(x, crit = quote(con*cov), cond = quote(TRUE)){
  stopifnot(inherits(x, "conCovOpt"))
  for (i in seq_along(x)){
    xi <- x[[i]]
    subsetCond <- eval(cond, xi, parent.frame())
    xi <- subset(xi, subsetCond) 
    critCond <- eval(crit, xi, parent.frame())
    xi <- xi[which.max(critCond), , drop = F]
    attributes(xi)[c("reprodList", "exoGroups", "allConCov")] <- 
      attributes(x[[i]])[c("reprodList", "exoGroups", "allConCov")]
    x[[i]] <- xi
  }
  attr(x, "crit") <- sapply(x[vapply(x, nrow, integer(1))>0], 
                            function(d) eval(crit, d))
  attr(x, "parms") <- list(crit = crit, cond = cond)
  class(x) <- "selectMax"
  x
}
print.selectMax <- function(x, ...){
  stopifnot(vapply(x, nrow, integer(1)) <= 1L)
  pr <- do.call(rbind, x)
  pr <- data.frame(outcome = rownames(pr), pr, row.names = NULL,
                   stringsAsFactors = FALSE)
  pr <- pr[order(attr(x, "crit"), decreasing = TRUE), , drop = FALSE]
  rownames(pr) <- NULL
  print(pr, ...)
  invisible(x)
}

# multipleMax
multipleMax <- function(x, outcome){
  stopifnot(inherits(x, "selectMax"), 
            outcome %in% names(x),
            length(outcome) == 1)
  eps <- 1e-12
  acc <- attr(x[[outcome]], "allConCov")
  if (is.null(acc)) stop("Missing attribute 'allConCov'")
  acc <- data.frame(outcome = outcome, 
                    acc,
                    id = seq_len(nrow(acc)))
  crit <- attr(x, "parms")$crit
  cond <- attr(x, "parms")$cond
  acc <- subset(acc, eval(cond, acc))
  out <- do.call(transform, list(acc, crit = crit))
  critMax <- eval(crit, x[[outcome]])
  out <- subset(out, crit >= critMax - eps)
  rownames(out) <- NULL
  out$crit <- NULL
  out
}

# ------------------------------------------------------------------------------

# findOutcomes
findOutcomes <- function(x, con = 1, cov = 1, ...){
  con_threshold <- con
  cov_threshold <- cov
  b <- conCovOpt(x, ...)
  sm <- selectMax(b, cond = quote(con >= con_threshold & cov >= cov_threshold))
  out <- data.frame(Factor = names(b), 
                    outcome = vapply(sm, nrow, integer(1L)) > 0L)
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------

# Auxiliary function getIndices
# Extract indices (applicable to reprodList) from id numbers stored in selctBase
# Note: Does something similiar to cna:::rowID()
getIndices <- function(i, ll) .getInd1(i-1L, ll) + 1L
.getInd1 <- function(i, ll){
  stopifnot(i-1 <= prod(ll))
  n <- length(ll)
  if (n == 1) return(matrix(i))
  pr1 <- prod(ll[-1])
  cbind((i) %/% pr1,
        .getInd1((i) %% pr1, ll[-1]))
}
if (F){
  getIndices(1:12, 12)
  getIndices(1:12, c(2, 6))
  getIndices(1:12, 3:4)
  getIndices(1:12, 4:3)
  getIndices(1:12, c(6, 2))
  getIndices(1:24, 2:4)
  getIndices(1:24, c(3, 2, 4))
}

# ------------------------------------------------------------

# reprodAssign
reprodAssign <- function(x, outcome, id = xi$id){
  stopifnot(inherits(x, "selectMax"), outcome %in% names(x),
            length(outcome) == 1)
  xi <- x[[outcome]]
  stopifnot(length(id) == 1)
  if (nrow(xi) == 0L) 
    stop("There is no solution for outcome ", outcome, " stored in ", deparse(substitute(x)))
  poss <- attr(xi, "reprodList")
  if (!(id %in% seq_len(prod(lengths(poss))))) stop("Invalid 'id' value specified")
  ii <- drop(getIndices(id, lengths(poss)))
  lhsSc_g <- mapply("[", poss, ii)
  exoGroups <- attr(xi, "exoGroups")
  ll <- lengths(exoGroups)
  out <- numeric(sum(ll))
  out[unlist(exoGroups)] <- rep(lhsSc_g, ll)
  out
}

# ------------------------------------------------------------------------------

# DNFbuild
DNFbuild <- function(x, outcome, reduce = c("rreduce", "ereduce", "none"), id = xi$id){
  stopifnot(inherits(x, "selectMax"), outcome %in% names(x))
  # resolve reduce arg
  if (isTRUE(reduce)) reduce <- "rreduce"
  if (isFALSE(reduce) | is.null(reduce)) reduce <- "none"
  reduce <- match.arg(reduce)
  # configTable
  ct <- attr(x, "configTable")
  type <- attr(ct, "type")
  if (type == "fs") stop("DNFbuild() has no implementation of the fs case.")
  xi <- x[[outcome]]
  if (nrow(xi) == 0L) 
    stop("There is no solution for outcome ", outcome, " stored in ", deparse(substitute(x)))
  poss <- attr(xi, "reprodList")
  lhsSc <- reprodAssign(x, outcome, id)
  d <- as.data.frame(ct)
  outcomeVar <- if (type == "mv") sub("=.+", "", outcome) else outcome
  d[[outcomeVar]] <- NULL
  stopifnot(nrow(d) == length(lhsSc))
  dups <- duplicated(d)
  d <- d[!dups, , drop = FALSE]
  lhsSc <- lhsSc[!dups]
  d <- subset(d, lhsSc==1)
  b <- matrix(colnames(d), nrow(d), ncol(d), byrow = TRUE)
  if (type == "cs"){
    b[d == 0] <- tolower(b[d == 0])
  } else if (type == "mv"){
    b <- array(paste0(b, "=", as.matrix(d)), dim(b))
  }
  out <- C_recCharList2char(list(split(b, row(b))), " + ")
  if (reduce == "rreduce") out <- rreduce(out, ct, full = FALSE)
  if (reduce == "ereduce") out <- ereduce(out, ct, full = FALSE)
  out
}
