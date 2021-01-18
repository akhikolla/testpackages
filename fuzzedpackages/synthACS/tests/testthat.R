library(data.table)
library(acs)
library(synthACS)
library(testthat)

# test_check("synthACS")

# LOAD TEST DATA FOR LOCAL TESTS
#------------------------------------------------
# load("./tests/testthat/dat-test_micro.xz")
# load("./tests/testthat/dat-par_sim_anneal.Rdata")

# MAKE TEST DATA
#------------------------------------------------
# library(synthACS)
# ca_geo <- geo.make(state= 'CA', county= '*')
# ca_dat <- pull_synth_data(2012, 5, ca_geo)
# save.image("./synthACS/tests/testthat/acsdat.Rdata")
# towork.Rda == object transit_work from paper

