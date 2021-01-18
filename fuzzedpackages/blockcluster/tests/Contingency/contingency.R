library(blockcluster)

## load data
#####################################################################
data(contingencydataunknown)
## cemInitStep + BEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BEM")
## usage of coclusterContingency function in its most simplest form
message("Contingency Data, cemInit, BEM")
set.seed(42); out<-coclusterContingency(contingencydataunknown,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## cemInitStep + BCEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BCEM")
## usage of coclusterContingency function in its most simplest form
message("Contingency Data, cemInit, BCEM")
set.seed(42); out<-coclusterContingency(contingencydataunknown,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

## cemInitStep + BSEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BSEM")
## usage of coclusterContingency function in its most simplest form
message("Contingency Data, cemInit, BSEM")
set.seed(42); out<-coclusterContingency(contingencydataunknown,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)
#####################################################################

#####################################################################
data(contingencydataunknown)
## emInitStep + BEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BEM")
## usage of coclusterContingency function in its most simplest form
message("Contingency Data, emInit, BEM")
set.seed(42); out<-coclusterContingency(contingencydataunknown,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

data(contingencydataunknown)
## emInitStep + BCEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BCEM")
## usage of coclusterContingency function in its most simplest form
message("Contingency Data, emInit, BCEM")
set.seed(42); out<-coclusterContingency(contingencydataunknown,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)

data(contingencydataunknown)
## emInitStep + BSEM
strat = coclusterStrategy(initmethod = "emInitStep", algo = "BSEM")
## usage of coclusterContingency function in its most simplest form
message("Contingency Data, emInit, BSEM")
set.seed(42); out<-coclusterContingency(contingencydataunknown,nbcocluster=c(2,3), strategy = strat)
sort(out@rowproportions)
sort(out@columnproportions)
#####################################################################

#####################################################################
## Test with wrong nbcocluster
## cemInitStep + BEM
strat = coclusterStrategy(initmethod = "cemInitStep", algo = "BEM")
## usage of coclusterContingency function in its most simplest form
message("Contingency Data, cemInit, BEM, wrong nbcocluster")
set.seed(42); out<-coclusterContingency(contingencydataunknown,nbcocluster=c(3,3), strategy = strat)

