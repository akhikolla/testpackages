library(blockcluster)

## load data
data(binarydata)
#####################################################################
## cemInitStep + BEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BEM")
## usage of coclusterBinary function in its most simplest form
message("Binary Data, cemInit, BEM")
set.seed(42); out<-coclusterBinary(binarydata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## cemInitStep + BCEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BCEM")
## usage of coclusterBinary function in its most simplest form
message("Binary Data, cemInit, BCEM")
set.seed(42); out<-coclusterBinary(binarydata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## cemInitStep + BSEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BSEM")
## usage of coclusterBinary function in its most simplest form
message("Binary Data, cemInit, BSEM")
set.seed(42); out<-coclusterBinary(binarydata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)
#####################################################################

#####################################################################
## emInitStep + BEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BEM")
## usage of coclusterBinary function in its most simplest form
message("Binary Data, emInit, BEM")
set.seed(42); out<-coclusterBinary(binarydata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## emInitStep + BCEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BCEM")
## usage of coclusterBinary function in its most simplest form
message("Binary Data, emInit, BCEM")
set.seed(42); out<-coclusterBinary(binarydata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## emInitStep + BSEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BSEM")
## usage of coclusterBinary function in its most simplest form
message("Binary Data, emInit, BSEM")
set.seed(42); out<-coclusterBinary(binarydata,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)
#####################################################################

#####################################################################
## Test with wrong nbcocluster
## cemInitStep + BEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BEM")
## usage of coclusterBinary function in its most simplest form
message("Binary Data, cemInit, BEM, wrong nbcocluster")
set.seed(42); out<-coclusterBinary(binarydata,nbcocluster=c(3,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

