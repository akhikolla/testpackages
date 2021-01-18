#' Network Degree Distribution Estimated from Labeled Snowball
#' Sample with Multiple Inclusion (LSMI)
#'
#' lsmi_dd computes an empirical network degree distribution and estimates
#' mean degree based on data from an LSMI sample from a network;
#' see Equations 6 and 7 by \insertCite{thompson_etal_2016;textual}{snowboot}
#' and Equation 1 by \insertCite{chen_etal_2018_snowboot;textual}{snowboot}
#' on the details of the calculations.
#'
#' The samples produced with \code{\link{lsmi}} or \code{\link{lsmi_union}} contain
#' just node IDs arranged into lists of seeds and waves (no details on the
#' node degrees or other node features). This information is
#' sufficient to study some properties of a network (e.g., network motifs --
#' not yet implemented in the package). To estimate a degree distribution or
#' mean degree, both the LSMI sample and the original network object are required.
#' If the LSMI object \code{x} is not supplied, the function will attempt
#' sampling an LSMI automatically, using the arguments supplied in "\code{...}"
#' that will be passed to the \code{\link{lsmi}} function.
#'
#' @inheritParams lsmi
#' @param x the LSMI sample obtained from the network \code{net}, for example,
#' with \code{\link{lsmi}} function or as a subset of the output
#' by \code{\link{lsmi_union}}.
#' @param ... arguments passed to the \code{\link{lsmi}} function
#' (ignored if \code{x} is specified, see Details).
#'
#' @return A list object of class "\code{snowboot}" consisting of:
#' \item{fk}{A named numeric vector with estimated probabilities \eqn{\hat{f}(k)}
#' of degrees \eqn{k}, where \eqn{k = 0, 1, \ldots,} \code{max(c(ds, dns))}
#' (i.e., \eqn{k} ranges from 0 to the maximum node degree observed in the LSMI sample).
#' The names of the vector elements are \eqn{k}.}
#' \item{mu}{An estimate of the mean degree.}
#' \item{ds}{An integer vector of degrees of seed nodes.}
#' \item{dns}{An integer vector of degrees of non-seed nodes (i.e., nodes
#' recorded in the waves of neighbors).}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{lsmi}}, \code{\link{lsmi_union}}, \code{\link{boot_dd}}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#'
#' #Obtain an LSMI sample and, at the next step,
#' #use it to estimate the degree distribution:
#' lsmiSample <- lsmi(net, n.seed = 5, n.wave = 3)
#' fkEstimate1 <- lsmi_dd(lsmiSample, net)$fk
#'
#' #Obtain an LSMI sample and estimate the degree
#' #distribution in a single step:
#' fkEstimate2 <- lsmi_dd(net = net, n.seed = 5, n.wave = 3)$fk
#'
#' #Use the output of lsmi_union to get the estimate:
#' lsmiUnionSample <- lsmi_union(net, n.seeds = c(5, 10), n.wave = 3)
#' fkEstimate3 <- lsmi_dd(lsmiUnionSample$lsmi_big, net)$fk
#'
lsmi_dd <- function(x = NULL, net, ...) {
  if(is.null(x)) {
    x <- lsmi(net, ...)
  }
  IDseeds <- sapply(x, function(y) y[[1]][[1]])
  IDnonseeds <- unlist(sapply(x, function(y) unlist(y[-1])))
  IDnonseeds <- na.omit(IDnonseeds)
  ds <- net$degree[IDseeds] #degrees of seeds
  dns <- net$degree[IDnonseeds] #degrees of non-seeds
  ns <- length(ds)
  Es <- mean(ds)
  dmax <- max(ds, dns)
  dds <- table(c(ds, 0:dmax)) - 1 #degree frequency of seeds
  p0 <- dds[1] / ns
  dds <- dds[-1]
  k <- c(1:dmax)
  ddns <- table(c(dns, k)) - 1 #degree frequency of non-seeds
  ddns <- ddns / k
  fk = (dds + (1 - p0) * Es * ddns) / (ns + Es * sum(ddns))
  k <- c(0:dmax)
  fk <- c(p0, fk)
  mu = sum(k * fk)
  res <- list(fk = fk, mu = mu, ds = ds, dns = dns)
  class(res) <- "snowboot"
  return(res)
}
