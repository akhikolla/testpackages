qegp3 <- function(p, kappa=1, sigma, xi, u=0, lower.tail=TRUE, log.p=FALSE){
  if (lower.tail){
    p <- if(log.p) {p / kappa} else {p^(1/kappa)}
    qgpd(p, sigma, xi, u, lower.tail=lower.tail, log.p=log.p)
  } else {
      ## transform, very carefully, to the appropriate thing
      ## for qgpd
      if (log.p) {
          p <- .log1mexp(.log1mexp(p)/kappa)
      } else {
          p <- -expm1(log1p(-p)/kappa)
      }
      qgpd(p, sigma, xi, u, lower.tail=lower.tail, log.p=log.p)
  }
}
