library(indelmiss)
context("Model fit output testing")
suppressWarnings(RNGversion("3.5.0"))
#right value, right class
test_that("Example 1", {
  skip_on_cran()
  loadNamespace("phangorn")
  set.seed(1)
  usertree <- ape::rtree(n = 7, br = rbeta(n = 7, shape1 = 1, shape2 = 7))
  data <- phangorn::simSeq(usertree, l = 5000, type = "USER", levels = c(0, 1), bf = c(1/(1 + 5), 5/(1 + 5)), Q = 1)
  datab <- matrix(as.numeric(as.character(data)), nrow = 7)
  userphyl <- t(datab)
  #Run the models.
  indel_user <- indelrates(datasource = "user", usertree = usertree, userphyl = userphyl, toi = 1, zerocorrection = TRUE, rootprob = "stationary", modelnames = c("M3", "M4"), optmethod = "nlminb")
  expect_equal(indel_user$results$M3$par, c(0.618884279006853, 2.95133445015348), tolerance = 0.1, scale = 1  )
  expect_equal(indel_user$results$M3$se, c(0.0147842484621341, 0.073182063204754), tolerance = 0.1, scale = 1  )
  expect_equal(indel_user$bg, list(1:13) )
  expect_equal(indel_user$results$M3$objective, 13491.1432157728, tolerance = 0.1, scale = 1  )
  expect_equal(indel_user$results$M3$df, 2 )
  expect_equal(indel_user$results$M3$AIC, -26986.286, tolerance = 0.1, scale = 1  )
  expect_equal(indel_user$results$M3$BIC, -26999.3203863828, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel_user$results$M4$par, c(0.614960843592059, 2.94050736262886, 0.00322102309478782), tolerance = 0.1, scale = 1  )
  expect_equal(indel_user$results$M4$se, c(0.0154105867104352, 0.0740736756527944, 0.00391236290876835
  ), tolerance = 0.1, scale = 1  )
  expect_equal(indel_user$bg, list(1:13) )
  expect_equal(indel_user$results$M4$objective, 13490.7946836198, tolerance = 0.1, scale = 1  )
  expect_equal(indel_user$results$M4$df, 3 )
  expect_equal(indel_user$results$M4$AIC, -26987.59, tolerance = 0.1, scale = 1  )
  expect_equal(indel_user$results$M4$BIC, -27007.1415795743, tolerance = 0.1, scale = 1  )
})

test_that("Example 2", {
  set.seed(1)
  indel1 <- indelrates(verbose = FALSE, datasource = "simulation")
  expect_equal(indel1$results$M1$par, 0.96672156017793, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M1$se, 0.021292353919985, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$bg, list(1:9) )
  expect_equal(indel1$results$M1$objective, 13991.8370991254, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M1$df, 1 )
  expect_equal(indel1$results$M1$AIC, -27985.674, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M1$BIC, -27992.1911931914, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel1$results$M2$par, c(0.966721576550205, 0), tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M2$se, c(0.0224627233157562, 0.0123046116115758) , tolerance = 0.1, scale = 1 )
  expect_equal(indel1$bg, list(1:9) )
  expect_equal(indel1$results$M2$objective, 13991.8370991254, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M2$df, 2 )
  expect_equal(indel1$results$M2$AIC, -27987.674, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M2$BIC, -28000.7083863828, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel1$results$M3$par, c(0.967369808843821, 1.01031309753709), tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M3$se, c(0.021459632871592, 0.0348671725136782) , tolerance = 0.1, scale = 1 )
  expect_equal(indel1$bg, list(1:9) )
  expect_equal(indel1$results$M3$objective, 13990.5236336781, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M3$df, 2 )
  expect_equal(indel1$results$M3$AIC, -27985.048, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M3$BIC, -27998.0823863828, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel1$results$M4$par, c(0.967365060558715, 1.01032125984796, 0), tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M4$se, c(0.0226014774923262, 0.0349802677500273, 0.0124005375170138), tolerance = 0.1, scale = 1  )
  expect_equal(indel1$bg, list(1:9) )
  expect_equal(indel1$results$M4$objective, 13990.5236336028, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M4$df, 3 )
  expect_equal(indel1$results$M4$AIC, -27987.048, tolerance = 0.1, scale = 1  )
  expect_equal(indel1$results$M4$BIC, -28006.5995795742, tolerance = 0.1, scale = 1  )
})

test_that("Example 3", {
  set.seed(1)
  indel2 <- indelrates(datasource = "simulation", seed = 1, taxa = 5, brlensh = c(1, 8), mu = 1, nu = 5, phyl = 5000, pmiss = 0, toi = 1,  zerocorrection = TRUE, rootprob = "stationary",  modelnames = c("M1", "M2", "M3", "M4"), optmethod = "nlminb")
  expect_equal(indel2$results$M1$par, 0.454999255951111, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M1$se, 0.0103789340113083, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$bg, list(1:9) )
  expect_equal(indel2$results$M1$objective, 10558.7944879097, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M1$df, 1 )
  expect_equal(indel2$results$M1$AIC, -21119.588, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M1$BIC, -21126.1051931914, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel2$results$M2$par, c(0.441387304237004, 0.0151123429435782), tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M2$se, c(0.0110135126957606, 0.00512141501934544), tolerance = 0.1, scale = 1  )
  expect_equal(indel2$bg, list(1:9) )
  expect_equal(indel2$results$M2$objective, 10554.1484594603, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M2$df, 2 )
  expect_equal(indel2$results$M2$AIC, -21112.296, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M2$BIC, -21125.3303863828, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel2$results$M3$par, c(0.593450786172035, 2.88917740005654) , tolerance = 0.1, scale = 1 )
  expect_equal(indel2$results$M3$se, c(0.016164220296224, 0.0969372866570588), tolerance = 0.1, scale = 1  )
  expect_equal(indel2$bg, list(1:9) )
  expect_equal(indel2$results$M3$objective, 9914.17116230853, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M3$df, 2 )
  expect_equal(indel2$results$M3$AIC, -19832.342, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M3$BIC, -19845.3763863828, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel2$results$M4$par, c(0.593450890187949, 2.88918238481361, 0), tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M4$se, c(0.0173477787482988, 0.0976009558697179, 0.00555730933531838
  ), tolerance = 0.1, scale = 1  )
  expect_equal(indel2$bg, list(1:9) )
  expect_equal(indel2$results$M4$objective, 9914.17116231024, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M4$df, 3 )
  expect_equal(indel2$results$M4$AIC, -19834.342, tolerance = 0.1, scale = 1  )
  expect_equal(indel2$results$M4$BIC, -19853.8935795743, tolerance = 0.1, scale = 1  )
})

test_that("Example 4", {
  skip_on_cran()
  set.seed(1)
  indel3 <- indelrates(datasource = "simulation", seed = 1, taxa = 5, brlensh = c(1, 8), mu = 1, nu = 5, phyl = 5000, pmiss = c(0, 0.300, 0.500, 0, 0), toi = "all",zerocorrection = TRUE, rootprob = "maxlik", modelnames = c("M3", "M4"),optmethod = "nlminb")
    
  expect_equal(indel3$results$M3$par, c(0.385727045370546, 5.53344820031837, 0.99), tolerance = 0.1, scale = 1  )
  expect_equal(indel3$results$M3$se, c(0.0663898915945263, 0.0867488350152869, 0.0174962850264404), tolerance = 0.1, scale = 1  )
  expect_equal(indel3$bg, list(1:9) )
  expect_equal(indel3$results$M3$objective, 14058.3306274577, tolerance = 0.1, scale = 1  )
  expect_equal(indel3$results$M3$df, 3 )
  expect_equal(indel3$results$M3$AIC, -28122.662, tolerance = 0.1, scale = 1  )
  expect_equal(indel3$results$M3$BIC, -28142.2135795743, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel3$results$M4$par, c(0.58235445899531, 2.98700863991035, 0, 0.300044233462258, 0.502901240329879, 
                                        0.00607446977611889, 0, 0.174223678015219), tolerance = 0.1, scale = 1  )
  expect_equal(indel3$results$M4$se, c(0.0850557906702689, 0.417786391060708, 0.0139481332118715, 
                                       0.0129600900226938, 0.0104402686345968, 0.0128008393638187, 0.0127712676335168, 
                                       0.0286346341960002), tolerance = 0.1, scale = 1  )
  expect_equal(indel3$bg, list(1:9) )
  expect_equal(indel3$results$M4$objective, 12474.437350036, tolerance = 0.1, scale = 1  )
  expect_equal(indel3$results$M4$df, 8 )
  expect_equal(indel3$results$M4$AIC, -24964.874, tolerance = 0.1, scale = 1  )
  expect_equal(indel3$results$M4$BIC, -25017.0115455313, tolerance = 0.1, scale = 1  )
})

test_that("Example 5", {
  skip_on_cran()
  set.seed(1)
  indel4 <- indelrates(datasource = "simulation", seed = 1, taxa = 10, brlensh = c(1, 8), mu = 1, nu = 5, phyl = 5000, pmiss = 0, toi = 1, bgtype = "ancestornodes", bg = c(15), zerocorrection = TRUE, rootprob = "maxlik", modelnames = c("M3", "M4"), optmethod = "nlminb")
  
  expect_equal(indel4$results$M3$par, c(0.630059740399771, 2.90570383259757, 0.611904192624238, 3.00737387388619, 
                                        0.170585645072776), tolerance = 0.1, scale = 1  )
  expect_equal(indel4$results$M3$se, c(0.0216284595065479, 0.106063142230079, 0.0158841922980097, 
                                       0.086643019353326, 0.0100851481629415), tolerance = 0.1, scale = 1  )
  expect_equal(indel4$bg, list(c(15, 16, 7, 4, 17, 5, 6), c(1, 2, 3, 8, 9, 10, 11, 12, 
                                                            13, 14, 18, 19, 15)) )
  expect_equal(indel4$results$M3$objective, 18250.2891816449, tolerance = 0.1, scale = 1  )
  expect_equal(indel4$results$M3$df, 5 )
  expect_equal(indel4$results$M3$AIC, -36510.578, tolerance = 0.1, scale = 1  )
  expect_equal(indel4$results$M3$BIC, -36543.1639659571, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel4$results$M4$par, c(0.630449653689202, 2.90123485477509, 0.608373649646646, 2.93396857472047, 
                                        0.0068001492467262, 0.162664613107343), tolerance = 0.1, scale = 1  )
  expect_equal(indel4$results$M4$se, c(0.0216309098112965, 0.106012999147008, 0.015800653079229, 0.0927215170316333, 
                                       0.00353695076961093, 0.0105055622525034), tolerance = 0.1, scale = 1  )
  expect_equal(indel4$bg, list(c(15, 16, 7, 4, 17, 5, 6), c(1, 2, 3, 8, 9, 10, 11, 12, 
                                                            13, 14, 18, 19, 15)) )
  expect_equal(indel4$results$M4$objective, 18248.2994453113, tolerance = 0.1, scale = 1  )
  expect_equal(indel4$results$M4$df, 6 )
  expect_equal(indel4$results$M4$AIC, -36508.598, tolerance = 0.1, scale = 1  )
  expect_equal(indel4$results$M4$BIC, -36547.7011591485, tolerance = 0.1, scale = 1  )
})

test_that("Example 6", {
  skip_on_cran()
  set.seed(1)
  indel5 <- indelrates(verbose = FALSE, datasource = "simulation", seed = 1, taxa = 5, brlensh = c(1, 8), mu = 1, nu = 3, phyl = 5000, pmiss = 0, toi = 1, bgtype = "listofnodes", bg = list(c(7, 1, 2), c(6, 8, 3, 7, 9, 5, 4, 9)), zerocorrection = TRUE, rootprob = "maxlik", modelnames = c("M1", "M2", "M3", "M4"), optmethod = "nlminb")
  
  expect_equal(indel5$results$M3$par, c(0.685285939054463, 1.92027513917066, 0.647763452426702, 1.93943852743902, 
                                        0.258745026475232), tolerance = 0.1, scale = 1  )
  expect_equal(indel5$results$M3$se, c(0.0378587136068766, 0.164512061702596, 0.0245080040440439, 
                                       0.123733788878302, 0.0104720200517138), tolerance = 0.1, scale = 1  )
  expect_equal(indel5$bg, list(c(7, 1, 2), c(6, 8, 3, 7, 9, 5, 4, 9)) )
  expect_equal(indel5$results$M3$objective, 11052.4732831761, tolerance = 0.1, scale = 1  )
  expect_equal(indel5$results$M3$df, 5 )
  expect_equal(indel5$results$M3$AIC, -22114.946, tolerance = 0.1, scale = 1  )
  expect_equal(indel5$results$M3$BIC, -22147.5319659571, tolerance = 0.1, scale = 1  )
  
  expect_equal(indel5$results$M4$par, c(0.660373035600555, 1.94728626056224, 0.646525771837654, 1.95562651404274, 
                                        0.00508680454636529, 0.259933203537567), tolerance = 0.1, scale = 1 )
  expect_equal(indel5$results$M4$se, c(0.0604097978739725, 0.173727345813147, 0.0246631298194466, 
                                       0.128060207393939, 0.00964958823756083, 0.0107900386057528), tolerance = 0.1, scale = 1  )
  expect_equal(indel5$bg, list(c(7, 1, 2), c(6, 8, 3, 7, 9, 5, 4, 9)) )
  expect_equal(indel5$results$M4$objective, 11052.3361028594, tolerance = 0.1, scale = 1  )
  expect_equal(indel5$results$M4$df, 6 )
  expect_equal(indel5$results$M4$AIC, -22116.672, tolerance = 0.1, scale = 1  )
  expect_equal(indel5$results$M4$BIC, -22155.7751591485, tolerance = 0.1, scale = 1  )
})

# com <- indelrates(verbose = TRUE, usertree = mycobacteriumdata$tree, userphyl = mycobacteriumdata$phyl, datasource = "user", toi = c(3:9), 
#                   bgtype = "listofnodes", bg = list(c(2, 3, 15), c(4, 14, 15), 
#                                                     c(17:19, 5:8), c(9, 16, 17), c(10, 11, 12, 13, 16, 14, 
#                                                                                    1)), zerocorrection = TRUE, rootprob = "maxlik", 
#                   optmethod = "nlminb", numhessian = TRUE, modelnames = c("M4"), 
#                   control = list(eval.max = 50000, iter.max = 50000))