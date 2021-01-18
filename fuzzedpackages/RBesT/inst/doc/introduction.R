## ---- include=FALSE------------------------------------------------------
library(RBesT)
library(knitr)
knitr::opts_chunk$set(
    fig.width = 1.62*4,
    fig.height = 4
    )
## setup up fast sampling when run on CRAN
is_CRAN <- !identical(Sys.getenv("NOT_CRAN"), "true")
## NOTE: for running this vignette locally, please uncomment the
## following line:
## is_CRAN <- FALSE
.user_mc_options <- list()
if (is_CRAN) {
    .user_mc_options <- options(RBesT.MC.warmup=50, RBesT.MC.iter=100, RBesT.MC.chains=2, RBesT.MC.thin=1)
}

## ----results="asis",echo=FALSE-------------------------------------------
kable(AS)

## ------------------------------------------------------------------------
library(RBesT)
set.seed(34563)
map_mcmc <- gMAP(cbind(r, n-r) ~ 1 | study,
                 data=AS,
                 tau.dist="HalfNormal",
                 tau.prior=1,
                 beta.prior=2,
                 family=binomial)
print(map_mcmc)

## a graphical representation of model checks is available
pl <- plot(map_mcmc)

## a number of plots are immediately defined
names(pl)

## forest plot with model estimates
print(pl$forest_model)

## ------------------------------------------------------------------------
set.seed(36546)
map_mcmc_sens <- update(map_mcmc, tau.prior=1/2)
print(map_mcmc_sens)

## ------------------------------------------------------------------------
map <- automixfit(map_mcmc)
print(map)
plot(map)$mix

## ------------------------------------------------------------------------
round(ess(map, method="elir"))   ## default method
round(ess(map, method="moment"))
round(ess(map, method="morita"))

## ------------------------------------------------------------------------
## add a 20% non-informative mixture component
map_robust <- robustify(map, weight=0.2, mean=1/2)
print(map_robust)

## ------------------------------------------------------------------------
round(ess(map_robust))

## ------------------------------------------------------------------------
p_truth       <- seq(0.1,0.95,by=0.01)
uniform_prior <- mixbeta(c(1,1,1))
treat_prior   <- mixbeta(c(1,0.5,1)) # prior for treatment used in trial
lancet_prior  <- mixbeta(c(1,11,32)) # prior for control   used in trial
decision      <- decision2S(0.95, 0, lower.tail=FALSE)

design_uniform   <- oc2S(uniform_prior, uniform_prior, 24, 6, decision)
design_nonrobust <- oc2S(treat_prior,   map          , 24, 6, decision)
design_robust    <- oc2S(treat_prior,   map_robust   , 24, 6, decision)

typeI_uniform   <- design_uniform(  p_truth, p_truth)
typeI_nonrobust <- design_nonrobust(p_truth, p_truth)
typeI_robust    <- design_robust(   p_truth, p_truth)

ocI <- rbind(data.frame(p_truth=p_truth, typeI=typeI_robust,    prior="robust"),
             data.frame(p_truth=p_truth, typeI=typeI_nonrobust, prior="non-robust"),
             data.frame(p_truth=p_truth, typeI=typeI_uniform,   prior="uniform")
             )

library(ggplot2)
theme_set(theme_bw()) # nice plotting theme
qplot(p_truth, typeI, data=ocI, colour=prior, geom="line", main="Type I Error")

## ------------------------------------------------------------------------
delta <- seq(0,0.7,by=0.01)
m <- summary(map)["mean"]
p_truth1 <- m +   delta
p_truth2 <- m + 0*delta

power_uniform   <- design_uniform(  p_truth1, p_truth2)
power_nonrobust <- design_nonrobust(p_truth1, p_truth2)
power_robust    <- design_robust(   p_truth1, p_truth2)

ocP <- rbind(data.frame(p_truth1=p_truth1, p_truth2=p_truth2, delta=delta, power=power_robust,    prior="robust"),
             data.frame(p_truth1=p_truth1, p_truth2=p_truth2, delta=delta, power=power_nonrobust, prior="non-robust"),
             data.frame(p_truth1=p_truth1, p_truth2=p_truth2, delta=delta, power=power_uniform,   prior="uniform")
             )

qplot(delta, power, data=ocP, colour=prior, geom="line", main="Power")


## ------------------------------------------------------------------------
## Critical values at which the decision flips are given conditional
## on the outcome of the second read-out; as we like to have this as a
## function of the treatment group outcome, we flip label 1 and 2
decision_flipped <- decision2S(0.95, 0, lower.tail=TRUE)
crit_uniform     <- decision2S_boundary(uniform_prior, uniform_prior, 6, 24, decision_flipped)
crit_nonrobust   <- decision2S_boundary(map          , treat_prior  , 6, 24, decision_flipped)
crit_robust      <- decision2S_boundary(map_robust   , treat_prior  , 6, 24, decision_flipped)
treat_y2 <- 0:24
## Note that -1 is returned to indicated that the decision is never 1
ocC <- rbind(data.frame(y2=treat_y2, y1_crit=crit_robust(treat_y2),    prior="robust"),
             data.frame(y2=treat_y2, y1_crit=crit_nonrobust(treat_y2), prior="non-robust"),
             data.frame(y2=treat_y2, y1_crit=crit_uniform(treat_y2),   prior="uniform")
             )

qplot(y2, y1_crit, data=ocC, colour=prior, geom="step", main="Critical values y1(y2)")

## ------------------------------------------------------------------------
## just positive
decision(postmix(treat_prior, n=24, r=15), postmix(map, n=6, r=3))
## negative
decision(postmix(treat_prior, n=24, r=14), postmix(map, n=6, r=4))

## ------------------------------------------------------------------------
r_placebo <- 1
r_treat   <- 14

## first obtain posterior distributions...
post_placebo <- postmix(map_robust,  r=r_placebo, n=6)
post_treat   <- postmix(treat_prior, r=r_treat  , n=24)

## ...then calculate probability that the difference is smaller than
## zero
prob_smaller <- pmixdiff(post_treat, post_placebo,  0, lower.tail=FALSE)

prob_smaller

prob_smaller > 0.95

## alternativley we can use the decision object
decision(post_treat, post_placebo)

## ------------------------------------------------------------------------
sessionInfo()

## ----include=FALSE-------------------------------------------------------
options(.user_mc_options)

