library(dmbc)

# load data
data(simdiss, package = "dmbc")

# plotting the data
library(bayesplot)
cols <- color_scheme_set("brightblue")
plot(simdiss, colors = unlist(cols)[c(1, 6)], font = NA, cex.font = 0.75)

G <- 3
p <- 2
prm.prop <- list(z = 1.5, alpha = .75)
burnin <- 20000
nsim <- 10000
seed <- 2301

set.seed(seed)

control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
  alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
  nchains = 1, thin = 10, store.burnin = TRUE, threads = 1,
  parallel = "snow")
sim.dmbc <- dmbc(simdiss, p, G, control)

# summarize and plot the results
summary(sim.dmbc, include.burnin = FALSE)

color_scheme_set("mix-red-blue")
plot(sim.dmbc, what = "trace", regex_pars = "eta")

###

# Parallel version
# set.seed(101)
control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
  alpha.prop = prm.prop[["alpha"]], random.start = TRUE, threads = 2,
  parallel = "snow", nchains = 2, seed = 101, verbose = TRUE)
sim.dmbc <- dmbc(simdiss, p, G, control)

###

# an example plot using bayesplot
library(bayesplot)
sim_list <- dmbc_fit_list_to_list(sim.dmbc, include.burnin = TRUE)
mcmc_trace(sim_list, regex_pars = "lambda")

# other graphs using mcmcplots
library(mcmcplots)
sim_list_mcmc <- dmbc_fit_list_to_mcmc.list(sim.dmbc, include.burnin = TRUE)
traplot(sim_list_mcmc, regex = "lambda", greek = TRUE, style = "plain")
traplot(sim_list_mcmc, parms = c("z_p[1, 1, 1]", "z_p[5, 1, 1]",
  "z_p[10, 1, 1]", "z_p[3, 2, 1]", "z_p[9, 2, 1]", "z_p[16, 2, 1]"))
traplot(sim_list_mcmc, regex = "log", greek = TRUE, style = "plain")
autplot1(sim_list_mcmc[, "alpha[1]", drop = FALSE])

###

# parameter estimates
dmbc_get_postmean(sim.dmbc, chain = 1)
dmbc_get_ml(sim.dmbc, chain = 1)
dmbc_get_map(sim.dmbc, chain = 1)

###

# coda analyses
library(coda)
sim_dmbc_sub <- subset(sim.dmbc, regex_pars = c("alpha", "lambda"))
plot(sim_dmbc_sub)
cumuplot(sim_dmbc_sub)
gelman.plot(sim_dmbc_sub)
geweke.plot(sim_dmbc_sub)
raftery.diag(sim_dmbc_sub)
HPDinterval(sim_dmbc_sub)
heidel.diag(sim_dmbc_sub)
densplot(sim_dmbc_sub)

###

# extracting estimated latent configuration
library(bayesplot)
library(ggplot2)
z <- dmbc_get_configuration(sim.dmbc, chain = 1, est = "mean",
  labels = 1:16)
# summary(z)
color_scheme_set("mix-pink-blue")
graph <- plot(z, size = 2, size_lbl = 3)
graph + panel_bg(fill = "gray90", color = NA)
