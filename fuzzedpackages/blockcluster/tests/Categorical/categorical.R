library(blockcluster)

## load data
data(categoricaldata)
#####################################################################
## cemInitStep + BEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, cemInit, BEM")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## cemInitStep + BCEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BCEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, cemInit, BCEM")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## cemInitStep + BSEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BSEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, cemInit, BSEM")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)
#####################################################################

#####################################################################
## emInitStep + BEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, emInit, BEM")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## emInitStep + BCEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BCEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, emInit, BCEM")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## emInitStep + BSEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BSEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, emInit, BSEM")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)
#####################################################################

#####################################################################
## randomInit + BEM
strat = coclusterStrategy(initmethod = "randomInit", algo = "BEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, randomInit, BEM")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## randomInit + BCEM
strat = coclusterStrategy(initmethod = "randomInit", algo = "BCEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, randomInit, BCEM")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## randomInit + BSEM
strat = coclusterStrategy(initmethod = "randomInit", algo = "BSEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, randomInit, BSEM")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)
#####################################################################

#####################################################################
## Test with wrong nbcocluster
## cemInitStep + BEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BEM")
## usage of coclusterCategorical function in its most simplest form
message("Categorical Data, cemInit, BEM, wrong nbcocluster")
set.seed(42); out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,3), strategy = strat)

