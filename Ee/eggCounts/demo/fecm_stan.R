##################################################################
## Illustrate faecal egg count reduction modelling
##################################################################
# load libraries
library('eggCounts')
require('coda',quietly = TRUE)
require('rstan',quietly = TRUE)

# simulate faecal egg counts
set.seed(1)
counts <- simData2s(n = 15, preMean = 500, delta = 0.07, kappa = 1, 
                   f = 15, paired = TRUE)

# look at simulated counts: dataframe with columns 
head(counts, 5)

# plot FECs
plotCounts(counts[,c("obsPre","obsPost")])

# run a paired model with individual efficacy
model <- fecr_stan(counts$obsPre, counts$obsPost, rawCounts = FALSE, preCF = 15, 
                   paired = TRUE, zeroInflation = FALSE, indEfficacy = TRUE)

# compute the probability that the reduction is below 95%
fecr_probs(model$stan.samples, threshold = 0.95, plot = TRUE)

# extract posterior samples 
model_mcmc <- stan2mcmc(model$stan.samples)

# check the traceplots and densities of posterior samples
plot(model_mcmc[,c("kappa","delta_mu","delta_shape")])
