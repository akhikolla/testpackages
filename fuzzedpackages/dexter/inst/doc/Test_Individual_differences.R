## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(echo = TRUE, dev='CairoPNG')

par_hook = function(before, options, envir)
{
  if(before)
  {
    do.call(par, options$par)
  }
}
knit_hooks$set(par = par_hook)

library(dexter)
library(dplyr)

## -----------------------------------------------------------------------------
items = tibble(item_id=paste0('i',1:20), item_score=1, delta=runif(20, -2, 2))
sim_Rasch = r_score(items)

theta = rep(0.5, 2000)
simulated = sim_Rasch(theta)

## ---- fig.align='center', fig.width=7,par=list(mfrow=c(1,2))------------------
hist(rowSums(simulated), main='', xlab='sumScore')
plot(ecdf(rowSums(simulated)), bty='l', main='ecdf', xlab='sumScore' )

## ---- fig.align='center', fig.width=7,par=list(mfrow=c(1,2))------------------
mm = fit_inter(simulated)

plot(mm, show.observed = TRUE, 
     items = c('i1','i2'))

## ---- fig.align='center', fig.height=4, fig.width=4---------------------------
dd = individual_differences(simulated)
plot(dd)

## -----------------------------------------------------------------------------
dd

## ---- fig.align='center', fig.height=4, fig.width=4,results='hide',message=FALSE----
db = start_new_project(verbAggrRules, ":memory:")
add_booklet(db, verbAggrData, "data")

dd = individual_differences(db, booklet_id=="data")
plot(dd)

## -----------------------------------------------------------------------------
dd

## ---- include=FALSE-----------------------------------------------------------
close_project(db)

