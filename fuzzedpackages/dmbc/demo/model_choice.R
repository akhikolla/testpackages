library(dmbc)

# load data
data(simdiss, package = "dmbc")

# model selection
pmax <- 2
Gmax <- 2
prm.prop <- list(z = 1.5, alpha = .75)
burnin <- 20000
nsim <- 10000
seed <- 1809

set.seed(seed)

control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
  alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
  thin = 10, store.burnin = TRUE, procrustes = TRUE, relabel = FALSE)
sim.ic <- dmbc_IC(data = simdiss, pmax = pmax, Gmax = Gmax, control = control,
  est = "mean")

# update the dmbc_ic object with additional values of p and G
pmax <- pmax + 1
Gmax <- Gmax + 2
new.ic <- update(sim.ic, pmax = pmax, Gmax = Gmax)

# summarize results
summary(new.ic)

# plot the results
library(bayesplot)
color_scheme_set("mix-yellow-blue")
p <- plot(new.ic, size = c(4, 1.5))
p + panel_bg(fill = "gray90", color = NA)
