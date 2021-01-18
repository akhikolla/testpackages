require(etm)

## Print with a bit less precision to avoid lots of notes in the comparison
old <- options(digits = 4)

### Simple test

time <- id <- 1:10
from <- rep(0, 10)
to <- rep(1, 10)

data1 <- data.frame(id, from, to, time)
tra1 <- matrix(FALSE, 2, 2)
tra1[1, 2] <- TRUE

etm1 <- etm(data1, c("0", "1"), tra1, NULL, 0)

all.equal(as.vector(trprob(etm1, "0 0")), cumprod((10:1 - 1) / (10:1)))

etm1$cov["0 0", "0 0", ]

all.equal(etm1$cov["0 0", "0 0",], trcov(etm1, "0 0"))

### A simple test from AHR's author, where the first time is censored
if (!require(survival)) {
    stop("This test requires the survival package")
}

data <- data.frame(id=1:10, time=1:10, from=0, to=1, status=TRUE)

tra <- matrix(FALSE, nrow=2, ncol=2)
tra[1, 2] <- TRUE

data$to[1] <- "cens"
data$status[1] <- FALSE

fit.km <- survfit(Surv(time, status) ~ 1, data=data)
fit.etm <- etm(data, c("0","1"), tra, "cens", s=0, t="last", covariance=FALSE)

all.equal(fit.km$surv[data$status], fit.etm$est[1,1,], check.attributes = FALSE)

data$to[2] <- "cens"
data$status[2] <- FALSE

fit.km <- survfit(Surv(time, status) ~ 1, data=data)
fit.etm <- etm(data, c("0","1"), tra, "cens", s=0, t="last", covariance=FALSE)

all.equal(fit.km$surv[data$status], fit.etm$est[1,1,], check.attributes = FALSE)

### a bit more complicated

time <- id <- 1:10
from <- rep(0, 10)
to <- rep(c(1, 2), 5)
data2 <- data.frame(id, from, to, time)

tra2 <- matrix(FALSE, 3, 3)
tra2[1, 2:3] <- TRUE

etm2 <- etm(data2, c("0", "1", "2"), tra2, NULL, 0)

aa <- table(time, to)

cif1 <- cumsum(aa[, 1] / 10)
cif2 <- cumsum(aa[, 2] / 10)
surv <- cumprod((10:1 - 1) / (10:1))

all.equal(trprob(etm2, "0 1"), cif1)
all.equal(trprob(etm2, "0 2"), cif2)
all.equal(as.vector(trprob(etm2, "0 0")), surv)

## a test on id
data2$id <- letters[1:10]

etm3 <- etm(data2, c("0", "1", "2"), tra2, NULL, 0)

all.equal(trprob(etm2, "0 1"), trprob(etm3, "0 1"))
all.equal(trprob(etm2, "0 2"), trprob(etm3, "0 2"))
all.equal(trprob(etm2, "0 0"), trprob(etm3, "0 0"))


### Test on sir.cont

data(sir.cont)

## Modification for patients entering and leaving a state
## at the same date
## Change on ventilation status is considered
## to happen before end of hospital stay
sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time), ]
for (i in 2:nrow(sir.cont)) {
    if (sir.cont$id[i]==sir.cont$id[i-1]) {
        if (sir.cont$time[i]==sir.cont$time[i-1]) {
            sir.cont$time[i-1] <- sir.cont$time[i-1] - 0.5
        }
    }
}

### Computation of the transition probabilities
## Possible transitions.
tra <- matrix(ncol=3,nrow=3,FALSE)
tra[1, 2:3] <- TRUE
tra[2, c(1, 3)] <- TRUE

## etm
prob.sir <- etm(sir.cont, c("0", "1", "2"), tra, "cens", 1)

prob.sir

summ.sir <- summary(prob.sir)
all.equal(summ.sir[['0 1']]$P, as.vector(trprob(prob.sir, "0 1")))
subset(summ.sir[[3]],time<183) # issue with precision on different architectures

## gonna play a bit with the state names
dd <- sir.cont
dd$from <- ifelse(dd$from == 0, "initial state", "ventilation")
dd$to <- as.character(dd$to)
for (i in seq_len(nrow(dd))) {
    dd$to[i] <- switch(dd$to[i],
                    "0" = "initial state",
                    "1" = "ventilation",
                    "2" = "end of story",
                    "cens" = "cens"
                    )
}

test <- etm(dd, c("initial state", "ventilation", "end of story"), tra, "cens", 1)

all.equal(test$est["initial state", "initial state", ],
          prob.sir$est["0", "0", ])
all.equal(trprob(test, "initial state initial state"), trprob(prob.sir, "0 0"))
all.equal(trprob(test, "initial state ventilation"), trprob(prob.sir, "0 1"))
all.equal(trprob(test, "initial state end of story"), trprob(prob.sir, "0 2"))

all.equal(trcov(test, "initial state end of story"), trcov(prob.sir, "0 2"))

aa <- summary(test)
all.equal(summ.sir[[6]], aa[[6]])
all.equal(summ.sir[[4]], aa[[4]])

### Test on abortion data

data(abortion)

from <- rep(0, nrow(abortion))
to <- abortion$cause
entry <- abortion$entry
exit <- abortion$exit
id <- 1:nrow(abortion)
data <- data.frame(id, from, to, entry, exit, group = abortion$group)

## Computation of the CIFs
tra <- matrix(FALSE, 4, 4)
tra[1, 2:4] <- TRUE

cif.control <- etm(data[data$group == 0, ], c("0", "1", "2", "3"),
                        tra, NULL, 0)
cif.exposed <- etm(data[data$group == 1, ], c("0", "1", "2", "3"),
                        tra, NULL, 0)

all.equal(trprob(cif.control, "0 1"), cif.control$est["0", "1", ])
all.equal(trcov(cif.control, c("0 1", "0 2")), cif.control$cov["0 1", "0 2", ])

trprob(cif.control, "0 1")
trprob(cif.control, "0 2")
trprob(cif.control, "0 0")

trcov(cif.control, "0 1")
trcov(cif.control, "0 2")
trcov(cif.control, "0 0")

aa <- summary(cif.control)
aa$"0 1"
all.equal(aa$"0 1"$P, as.vector(trprob(cif.control, "0 1")))

### test on los data

data(los.data) # in package changeLOS

## putting los.data in the long format (see changeLOS)
my.observ <- prepare.los.data(x=los.data)

tra <- matrix(FALSE, 4, 4)
tra[1, 2:4] <- TRUE
tra[2, 3:4] <- TRUE

tr.prob <- etm(my.observ, c("0","1","2","3"), tra, NULL, 0)

tr.prob
summary(tr.prob)

cLOS <- etm::clos(tr.prob, aw = TRUE)

cLOS


### Tests on pseudo values
t_pseudo <- closPseudo(my.observ, c("0","1","2","3"), tra, NULL,
                       formula = ~ 1, aw = TRUE)

cLOS$e.phi == t_pseudo$theta[, "e.phi"]
cLOS$e.phi.weights.1 == t_pseudo$theta[, "e.phi.weights.1"]
cLOS$e.phi.weights.other == t_pseudo$theta[, "e.phi.weights.other"]

mean(t_pseudo$pseudoData$ps.e.phi)

### tests on etmprep

### creation of fake data in the wild format, following an illness-death model
## transition times
tdisease <- c(3, 4, 3, 6, 8, 9)
tdeath <- c(6, 9, 8, 6, 8, 9)

## transition status
stat.disease <- c(1, 1, 1, 0, 0, 0)
stat.death <- c(1, 1, 1, 1, 1, 0)

## a covariate that we want to keep in the new data
set.seed(1313)
cova <- rbinom(6, 1, 0.5)

dat <- data.frame(tdisease, tdeath,
                  stat.disease, stat.death,
                  cova)

## Possible transitions
tra <- matrix(FALSE, 3, 3)
tra[1, 2:3] <- TRUE
tra[2, 3] <- TRUE

## data preparation
newdat <- etmprep(c(NA, "tdisease", "tdeath"),
                  c(NA, "stat.disease", "stat.death"),
                  data = dat, tra = tra,
                  cens.name = "cens", keep = "cova")

newdat

ref <- data.frame(id = c(1, 1, 2, 2, 3, 3, 4, 5, 6),
                  entry = c(0, 3, 0, 4, 0, 3, 0, 0, 0),
                  exit = c(3, 6, 4, 9, 3, 8, 6, 8, 9),
                  from = c(0, 1, 0, 1, 0, 1, 0, 0, 0),
                  to = c(rep(c(1, 2), 3), 2, 2, "cens"),
                  cova = c(1, 1, 0, 0, 1, 1, 0, 1, 1))
ref$from <- factor(as.character(ref$from), levels = c("0", "1", "2", "cens"))
ref$to <- factor(as.character(ref$to), levels = c("0", "1", "2", "cens"))

all.equal(ref, newdat)


######################################
### Test the stratified calls
######################################

if (require("kmi", quietly = TRUE)) {
    library(etm)

    data(icu.pneu)
    my.icu.pneu <- icu.pneu

    my.icu.pneu <- my.icu.pneu[order(my.icu.pneu$id, my.icu.pneu$start), ]
    masque <- diff(my.icu.pneu$id)

    my.icu.pneu$from <- 0
    my.icu.pneu$from[c(1, masque) == 0] <- 1

    my.icu.pneu$to2 <- my.icu.pneu$event
    my.icu.pneu$to2[my.icu.pneu$status == 0] <- "cens"
    my.icu.pneu$to2[c(masque, 1) == 0] <- 1


    my.icu.pneu$to <- ifelse(my.icu.pneu$to2 %in% c(2, 3), 2,
                             my.icu.pneu$to2)

    my.icu.pneu <- my.icu.pneu[, c("id", "start", "stop", "from", "to",
                                   "to2", "age", "sex")]
    names(my.icu.pneu)[c(2, 3)] <- c("entry", "exit")

    bouh_strat <- etm(my.icu.pneu, c("0", "1", "2"), tra_ill(), "cens", 0, strata = "sex")

    bouh_female <- etm(my.icu.pneu[my.icu.pneu$sex == "F", ],
                       c("0", "1", "2"), tra_ill(), "cens", 0)

    all(bouh_strat[[1]]$est == bouh_female$est)

    ## Test the summary
    the_summary <- summary(bouh_strat)
    the_summary

    ## Test trprob
    all(trprob(bouh_strat, "0 1")[[1]] == trprob(bouh_female, "0 1"))
    all(trprob(bouh_strat, "0 1", c(0, 5, 10, 15))[[1]] == trprob(bouh_female, "0 1", c(0, 5, 10, 15)))

    ## Test trcov
    all(trcov(bouh_strat, "0 1")[[1]] == trcov(bouh_female, "0 1"))
    all(trcov(bouh_strat, c("0 1", "0 2"))[[1]] == trcov(bouh_female, c("0 1", "0 2")))
    all(trcov(bouh_strat, "0 1", c(0, 5, 10, 15))[[1]] == trcov(bouh_female, "0 1", c(0, 5, 10, 15)))
} else {
    print("These tests require the kmi package")
}

options(old)
