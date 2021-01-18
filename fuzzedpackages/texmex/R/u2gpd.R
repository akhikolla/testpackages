u2gpd <-
function(u, p = 1, th=0, sigma, xi) {
  qgpd((1 - u) / p, sigma=sigma, u=th, xi=xi, lower.tail=FALSE)
}
