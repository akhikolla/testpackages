library(blockcluster)

## load data
data(gaussiandata)
#####################################################################
## cemInitStep + BEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BEM")
## usage of coclusterContinuous function in its most simplest form
message("Continuous Data, cemInit, BEM")
set.seed(42); out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## cemInitStep + BCEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BCEM")
## usage of coclusterContinuous function in its most simplest form
message("Continuous Data, cemInit, BCEM")
set.seed(42); out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## cemInitStep + BSEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BSEM")
## usage of coclusterContinuous function in its most simplest form
message("Continuous Data, cemInit, BSEM")
set.seed(42); out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)
#####################################################################

#####################################################################
## emInitStep + BEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BEM")
## usage of coclusterContinuous function in its most simplest form
message("Continuous Data, emInit, BEM")
set.seed(42); out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## emInitStep + BCEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BCEM")
## usage of coclusterContinuous function in its most simplest form
message("Continuous Data, emInit, BCEM")
set.seed(42); out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## emInitStep + BSEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BSEM")
## usage of coclusterContinuous function in its most simplest form
message("Continuous Data, emInit, BSEM")
set.seed(42); out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)
#####################################################################

#####################################################################
## Test with wrong nbcocluster
## cemInitStep + BEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BEM")
## usage of coclusterContinuous function in its most simplest form
message("Continuous Data, cemInit, BEM, wrong nbcocluster")
set.seed(42); out<-coclusterContinuous(gaussiandata,nbcocluster=c(3,3), strategy = strat)
