Sys.setenv("R_TESTS" = "")
library(eggCounts)
library(testthat)
# simulate egg counts
set.seed(1)
count <- simData1s()
counts <- simData2s()
if (Sys.info()["sysname"] != "SunOS"){
  
fecrtCI(counts$masterPre, counts$masterPost, paired=TRUE)
  
# test raw counts and calculated epg result same mcmc samples
set.seed(1)
t1 <- stan2mcmc(fec_stan(count[,"obs"], rawCounts=FALSE, CF=50, zeroInflation = TRUE, nburnin=500,nsamples=1000,thinning=1)$stan.samples)
set.seed(1)
t2 <- stan2mcmc(fec_stan(count[,"obs"]/50, rawCounts=TRUE, CF=rep(50,10), zeroInflation = TRUE, nburnin=500,nsamples=1000,thinning=1)$stan.samples)
expect_that(t1, equals(t2))

# test 2-sample default models compiles okay
t6<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,paired=TRUE,zeroInflation = TRUE,indEfficacy = FALSE,nsamples=1000,
              nburnin=500)$stan.samples
expect_that(t6,is_a("stanfit"))

t7<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,paired=FALSE,zeroInflation = TRUE,indEfficacy = FALSE,nsamples=1000, nburnin=500)$stan.samples
expect_that(t7,is_a("stanfit"))

t8<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,paired=TRUE,zeroInflation = FALSE,indEfficacy = FALSE,nsamples=1000,
              nburnin=500)$stan.samples
expect_that(t8,is_a("stanfit"))

t9<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"][-1],rawCounts=TRUE,paired=FALSE,zeroInflation = FALSE,indEfficacy = FALSE,nsamples=1000, nburnin=500)$stan.samples
expect_that(t9,is_a("stanfit"))

t10<-fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,paired=TRUE,zeroInflation = FALSE, indEfficacy = TRUE, nsamples=1000, nburnin=500)$stan.samples
expect_that(t10,is_a("stanfit"))

expect_that(fecr_probs(t6, plot = FALSE), is_a("numeric")) 
expect_that(fecr_probs(t10, plot = FALSE), is_a("numeric")) 

# warning checks
expect_warning(fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE,paired=TRUE,zeroInflation = FALSE, indEfficacy = FALSE, nsamples=600, nburnin=300))

expect_warning(fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=TRUE, paired=TRUE, zeroInflation = TRUE, indEfficacy = FALSE, nsamples=1000, nburnin=500, adaptDelta = 0.3))

# check incorrect inputs
expect_that(fecr_stan(counts[,"masterPre"],counts[,"masterPost"],rawCounts=FALSE,preCF=50),throws_error())
# expect_that(fecr_stan(counts[,"obsPre"],counts[,"obsPost"],preCF=50.5), throws_error())
expect_that(fecr_stan(counts[,"obsPre"],counts[,"obsPost"],preCF=50, paired=FALSE), throws_error())}

# check functions as of version 2.0

expect_that(getPrior_delta(0.4,0.7,p=0.6, plot=FALSE), is_a("matrix")) 
expect_that(getPrior_mu(200,0.3,500,0.7, plot=FALSE), is_a("matrix")) 

# t11 <- fecr_stanExtra(counts[,"masterPre"],counts[,"masterPost"], rawCounts=TRUE, modelName = "Po")$stan.samples
# expect_that(t11,is_a("stanfit"))

