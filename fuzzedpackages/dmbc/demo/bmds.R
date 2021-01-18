library(dmbc)

# BMDS Example 1: Airline Distances Between Cities
airline <- read.csv(file = system.file("extdata", "airline.csv",
  package = "dmbc"))
airline.nm <- airline[, 1]
airline <- airline[, 2:31]
colnames(airline) <- airline.nm
airline <- as.dist(airline)

min_p <- 1
max_p <- 6
burnin <- 200
nsim <- 1000
totiter <- burnin + nsim

airline.mds <- cmdscale(airline, max_p)
airline.bmds <- bmds(airline, min_p, max_p, burnin, nsim)

opar <- par(mfrow = c(1, 2))
plot(min_p:max_p, airline.bmds$mdsIC$mdsic, type = "b",
  main = "MDS Information Criterion", xlab = "p", ylab = "MDSIC")
MDSICmin <- which.min(airline.bmds$mdsIC$mdsic)
points((min_p:max_p)[MDSICmin], airline.bmds$mdsIC$mdsic[MDSICmin],
  col = "red", pch = 10, cex = 1.75, lwd = 1.5)

airline.bmds.x.mode <- bmds_get_x_mode(airline, airline.bmds, MDSICmin, min_p,
  max_p, start = (burnin + 1), end = totiter)
airline.bmds.d <- dist(airline.bmds.x.mode)
airline.mds.d <- dist(airline.mds[, 1:((min_p:max_p)[MDSICmin])])
plot(airline, airline.bmds.d, type = "n", xlab = "observed",
  ylab = "estimated", main = "Airline Distances \n Between Cities",
  xlim = c(0, max(airline, airline.bmds.d)),
  ylim = c(0, max(airline, airline.bmds.d)))
abline(0, 1, lty = 2, col = "gray")
points(airline, airline.mds.d, pch = 19, col = "cyan", cex = .5)
points(airline, airline.bmds.d, pch = 19, col = "magenta", cex = .5)
legend(x = "bottomright", legend = c("Classical MDS", "Bayesian MDS"),
  pch = c(19, 19), col = c("cyan", "magenta"))
par(opar)

# BMDS Example 2: Careers of Lloyds Bank Employees, 1905-1950
lloyds <- read.csv(file = system.file("extdata", "lloyds.csv",
  package = "dmbc"))
lloyds.nm <- lloyds[, 1]
lloyds <- lloyds[, 2:81]
colnames(lloyds) <- lloyds.nm
lloyds <- as.dist(lloyds)

min_p <- 1
max_p <- 12
burnin <- 200
nsim <- 1000
totiter <- burnin + nsim

lloyds.mds <- cmdscale(lloyds, max_p)
lloyds.bmds <- bmds(lloyds, min_p, max_p, burnin, nsim)

opar <- par(mfrow = c(1, 2))
plot((min_p:max_p), lloyds.bmds$mdsIC$mdsic, type = "b",
  main = "MDS Information Criterion", xlab = "p", ylab = "MDSIC")
MDSICmin <- which.min(lloyds.bmds$mdsIC$mdsic)
points((min_p:max_p)[MDSICmin], lloyds.bmds$mdsIC$mdsic[MDSICmin],
  col = "red", pch = 10, cex = 1.75, lwd = 1.5)

lloyds.bmds.x.mode <- bmds_get_x_mode(lloyds, lloyds.bmds, MDSICmin,
  min_p, max_p, start = (burnin + 1), end = totiter)
lloyds.bmds.d <- dist(lloyds.bmds.x.mode)
lloyds.mds.d <- dist(lloyds.mds[, 1:((min_p:max_p)[MDSICmin])])
plot(lloyds, lloyds.bmds.d, type = "n", xlab = "observed",
  ylab = "estimated", main = "Careers of Lloyds \n Bank Employees, 1905-1950",
  xlim = c(0, max(lloyds, lloyds.bmds.d)),
  ylim = c(0, max(lloyds, lloyds.bmds.d)))
abline(0, 1, lty = 2, col = "gray")
points(lloyds, lloyds.mds.d, pch = 19, col = "cyan", cex = .5)
points(lloyds, lloyds.bmds.d, pch = 19, col = "magenta", cex = .5)
legend(x = "topleft", legend = c("Classical MDS", "Bayesian MDS"),
  pch = c(19, 19), col = c("cyan", "magenta"))
par(opar)

# BMDS Example 3: Road distances (in km) between 21 cities in Europe
data(eurodist, package = "datasets")

min_p <- 1
max_p <- 10
burnin <- 200
nsim <- 1000
totiter <- burnin + nsim

eurodist.mds <- cmdscale(eurodist, max_p)
eurodist.bmds <- bmds(eurodist, min_p, max_p, burnin, nsim)

opar <- par(mfrow = c(1, 2))
plot((min_p:max_p), eurodist.bmds$mdsIC$mdsic, type = "b",
  main = "MDS Information Criterion", xlab = "p", ylab = "MDSIC")
MDSICmin <- which.min(eurodist.bmds$mdsIC$mdsic)
points((min_p:max_p)[MDSICmin], eurodist.bmds$mdsIC$mdsic[MDSICmin],
  col = "red", pch = 10, cex = 1.75, lwd = 1.5)

eurodist.bmds.x.mode <- bmds_get_x_mode(eurodist, eurodist.bmds,
  MDSICmin, min_p, max_p, start = (burnin + 1), end = totiter)
eurodist.bmds.d <- dist(eurodist.bmds.x.mode)
eurodist.mds.d <- dist(eurodist.mds[, 1:((min_p:max_p)[MDSICmin])])
plot(eurodist, eurodist.bmds.d, type = "n", xlab = "observed",
  ylab = "estimated", main = "Road distances (in km) \n between 21 cities in Europe",
  xlim = c(0, max(eurodist, eurodist.bmds.d)),
  ylim = c(0, max(eurodist, eurodist.bmds.d)))
abline(0, 1, lty = 2, col = "gray")
points(eurodist, eurodist.mds.d, pch = 19, col = "cyan", cex = .5)
points(eurodist, eurodist.bmds.d, pch = 19, col = "magenta", cex = .5)
legend(x = "topleft", legend = c("Classical MDS", "Bayesian MDS"),
  pch = c(19, 19), col = c("cyan", "magenta"))
par(opar)
