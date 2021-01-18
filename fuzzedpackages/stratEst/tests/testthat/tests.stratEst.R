library(stratEst)

test_that("core systematic tests",  {
  # unsorted data
  set.seed(1)
  sorted_data <- data.DF2011[data.DF2011$treatment=="D5R48",]
  unsorted_data <- sorted_data[sample(nrow(sorted_data)),]
  unsorted <- stratEst.model(unsorted_data,strategies=strategies.DF2011,inner.runs = 5, outer.runs = 2,sample.id = "treatment", verbose = F)
  set.seed(1)
  sorted <- stratEst.model(sorted_data,strategies=strategies.DF2011,inner.runs = 5, outer.runs = 2,sample.id = "treatment", verbose = F)
  unsorted_shares <- round(as.numeric(unsorted$shares),3)
  sorted_shares <- round(as.numeric(sorted$shares),3)
  expect_equal(sorted_shares[1],unsorted_shares[1])
  expect_equal(sorted_shares[2],unsorted_shares[2])
  expect_equal(sorted_shares[4],unsorted_shares[4])
  expect_equal(sorted_shares[5],unsorted_shares[5])

  # fix shares
  data <- data.DF2011[data.DF2011$treatment=="D5R40",]
  strats <- lapply( 1:3, function(x) stratEst.strategy(choices = c("c","d"), inputs = c("cc","cd","dc","dd")) )
  not_fixed <- stratEst.model(data,strats,inner.runs = 10, outer.runs = 2, verbose = F)
  not_fixed_samples <- stratEst.model(data,strats,sample.id="treatment",inner.runs = 10, outer.runs = 2,verbose = F)
  one_fixed_share <- stratEst.model(data,strats,shares=c(0.5,NA,NA),inner.runs = 10, outer.runs = 5,verbose = F)
  expect_equal(round(as.numeric(one_fixed_share$shares[1]),1),0.5)
  two_fixed_shares <- stratEst.model(data,strats,shares=c(0.5,0.2,NA),inner.runs = 10, outer.runs = 5,verbose = F)
  expect_equal(round(as.numeric(two_fixed_shares$shares[c(1,2)]),1),c(0.5,0.2))
  all_fixed_shares <- stratEst.model(data,strats,shares=c(0.5,0.2,0.3),inner.runs = 10, outer.runs = 5,verbose = F)
  expect_equal(round(as.numeric(all_fixed_shares$shares),1),c(0.5,0.2,0.3))
  # some_fixed_samples <- stratEst.model(data.DF2011[data.DF2011$treatment %in% c("D5R32","D5R40","D5R48"),],strats,shares=list("treatment.D5R32" = c(0.5,NA,NA),"treatment.D5R40" = c(NA,0.2,NA), "treatment.D5R48" = c(NA,NA,0.1)),sample.id="treatment",inner.runs = 10, outer.runs = 2,verbose = F)
  # expect_equal(round(c(as.numeric(some_fixed_samples$shares$treatment.D5R32[1]),as.numeric(some_fixed_samples$shares$treatment.D5R40[2]),as.numeric(some_fixed_samples$shares$treatment.D5R48[3])),1),c(0.5,0.2,0.1))
  # all_fixed_samples <- stratEst.model(data.DF2011[data.DF2011$treatment %in% c("D5R32","D5R40","D5R48"),],strats,shares=list("treatment.D5R32" = c(0.5,0.2,0.3),"treatment.D5R40" = c(0.6,0.2,0.2), "treatment.D5R48" = c(0.3,0.6,0.1)),sample.id="treatment",inner.runs = 10, outer.runs = 2,verbose = F)
  # expect_equal(round(as.numeric(unlist(all_fixed_samples$shares)),1),c(0.5,0.2,0.3,0.6,0.2,0.2,0.3,0.6,0.1))

  # fix probs

  # fix trembles
#
#   # restrict trembles
#   data <- data.DF2011[data.DF2011$treatment=="D5R40",]
#   restrict_trembles_no <- stratEst.model(data,strategies.PD[c("ALLD","TFT","WSLS")],r.trembles="no",verbose = F)
#   expect_equal(5,length(restrict_trembles_no$trembles.par))
#   restrict_trembles_global <- stratEst.model(data,strategies.PD[c("ALLD","TFT","WSLS")],r.trembles="global",verbose = F)
#   expect_equal(1,length(restrict_trembles_global$trembles.par))
#   restrict_trembles_strats <- stratEst.model(data,strategies.PD[c("ALLD","TFT","WSLS")],r.trembles="strategies",verbose = F)
#   expect_equal(3,length(restrict_trembles_strats$trembles.par))
#   restrict_trembles_states <- stratEst.model(data,strategies.PD[c("ALLD","TFT","WSLS")],r.trembles="states",verbose = F)
#   expect_equal(2,length(restrict_trembles_states$trembles.par))
#
#   # restrict probs
#   data <- data.DF2011[data.DF2011$treatment=="D5R40",]
#   strats <- lapply( 1:3, function(x) stratEst.strategy(choices = c("c","d"), inputs = c("cc","cd","dc","dd")) )
#   restrict_probs_no <- stratEst.model(data,strats,r.probs="no",verbose = F)
#   expect_equal(24,length(restrict_probs_no$probs.par))
#   restrict_probs_global <- stratEst.model(data,strats,r.probs="global",verbose = F)
#   expect_equal(2,length(restrict_probs_global$probs.par))
#   restrict_probs_strategies <- stratEst.model(data,strats,r.probs="strategies",verbose = F)
#   expect_equal(6,length(restrict_probs_strategies$probs.par))
#   restrict_probs_states <- stratEst.model(data,strats,r.probs="states",verbose = F)
#   expect_equal(8,length(restrict_probs_states$probs.par))
#
#   # select
#   set.seed(1)
#   data <- data.DF2011[data.DF2011$treatment=="D75R48",]
#   strats <- lapply( 1:3, function(x) stratEst.strategy(choices = c("c","d"), inputs = c("cc","cd","dc","dd")) )
#   select_strategies_icl <- stratEst.model(data,strats,select = "strategies", crit ="icl",verbose = F)
#   expect_equal(2,length(select_strategies_icl$shares))
#   set.seed(1)
#   select_trembles_global <- stratEst.model(data,strats,response="pure",r.trembles = "global",select = "trembles", crit ="icl",verbose = F)
#   expect_equal(3,length(select_trembles_global$trembles.par))
#   set.seed(1)
#   select_all_global <- stratEst.model(data,strats,r.trembles = "global",r.probs = "global",select = c("strategies","trembles","probs"), crit ="icl",verbose = F)
#   expect_equal(14,length(select_all_global$probs.par))
#   expect_equal(3,length(select_all_global$shares.par))
})


# test_that("additional systematic tests",  {
#   skip_on_cran()
#   # response
#   data <- data.DF2011[data.DF2011$treatment=="D5R32",]
#   set.seed(1)
#   strats <- lapply( 1:3, function(x) stratEst.strategy(choices = c("c","d"), inputs = c("cc","cd","dc","dd")) )
#   response_pure <- stratEst.model(data,strats,response= "pure",verbose = F)
#   set.seed(1)
#   response_mixed <- stratEst.model(data,strats,response= "mixed",verbose = F)
#   expect_equal(1,round(as.numeric(all( response_pure$probs.par == 0 | response_pure$probs.par == 1 ))))
#
#
#   #select
#   set.seed(1)
#   data <- data.DF2011[data.DF2011$treatment=="D75R48",]
#   strats <- lapply( 1:3, function(x) stratEst.strategy(choices = c("c","d"), inputs = c("cc","cd","dc","dd")) )
#   select_strategies_aic <- stratEst.model(data,strats,select = "strategies", crit ="aic",verbose = F)
#   expect_equal(2,round(as.numeric(length(select_strategies_aic$shares))))
#   set.seed(1)
#   select_strategies_bic <- stratEst.model(data,strats,select = "strategies", crit ="bic",verbose = F)
#   expect_equal(2,round(as.numeric(length(select_strategies_bic$shares))))
#   set.seed(1)
#   select_trembles_global <- stratEst.model(data,strats,response="pure",r.trembles = "global",select = "trembles", crit ="icl",verbose = F)
#   expect_equal(3,round(as.numeric(length(select_trembles_global$trembles.par))))
#   set.seed(1)
#   select_probs_global <- stratEst.model(data,strats,r.probs = "global",select = "probs", crit ="icl",verbose = F)
#   expect_equal(14,round(as.numeric(length(select_probs_global$probs.par))))
#   set.seed(1)
#   select_both_global <- stratEst.model(data,strats,r.trembles = "global",r.probs = "global",select = c("probs","trembles"), crit ="icl",verbose = F)
#   expect_equal(14,round(as.numeric(length(select_both_global$probs.par))))
# })

test_that("multivariate output test",  {
  skip_on_cran()
  N = 32                                       # number of subjects
  num_probs = c(2,3,4,5)                       # number of distinct probs of the machines
  tremble = 0.2
  shares = runif(5)                            # generate 5 shares
  shares = shares/sum(shares)                  # normalize shares
  for( m in 1:length(num_probs)){
    strats <- list()
    for( s in 1:5 ){
      check = 0
      while( check == 0 ){
        strategy_mat = matrix(runif(5*num_probs[m]),5,num_probs[m])
        strategy_mat = matrix(as.numeric(strategy_mat == apply(strategy_mat,1,max)),5,num_probs[m])
        if( sum(as.numeric(apply(strategy_mat,1,sum) == 1)) == 5 ){ check = 1 }
      }
      strats[[s]] <- stratEst.strategy(choices = as.character(c(1:num_probs[m])), inputs = as.character(c(1:5)), prob.choices = c(t(strategy_mat)) )
    }

    data <- stratEst.simulate(strategies = strats,shares = shares)

    P = stratEst.model(data,strats,response="pure",outer.runs = 2, verbose = F)
    M = stratEst.model(data,strats,response="mixed",outer.runs = 2, verbose = F)

    P_l = nrow( P$probs )
    M_l = nrow( M$probs )

    expect_equal( 25 , round(as.numeric(P_l)) )
    expect_equal( 25  , round(as.numeric(M_l)) )

  }
})

