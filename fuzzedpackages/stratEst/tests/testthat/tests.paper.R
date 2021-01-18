library(stratEst)

test_that("Introductory example paper" , {
  skip_on_cran()
  set.seed(1)
  rps = c("r", "p", "s")
  mixed = stratEst.strategy(choices = rps)
  nash = stratEst.strategy(choices = rps, prob.choices = rep(1/3, 3))

  expect_equal(1,as.numeric(ncol(mixed)==3 & nrow(mixed==1)))
  expect_equal(1,as.numeric(all(colnames(mixed) == c("prob.r","prob.p","prob.s")) & all(is.na(mixed))))
  expect_equal(1,as.numeric(ncol(nash)==3 & nrow(nash==1)))
  expect_equal(1,as.numeric(all(colnames(nash) == c("prob.r","prob.p","prob.s")) & all(nash==1/3)))

  last.choice = c(NA, rps)
  imitate =  stratEst.strategy(choices = rps, inputs = last.choice,
                               num.states = 4,
                               prob.choices = c(rep(1/3, 3), 1, 0, 0,
                                                0, 1, 0, 0, 0, 1),
                               tr.inputs = rep(c(2, 3, 4), 4))

  expect_equal(1,as.numeric(ncol(imitate)==7 & nrow(imitate==4)))
  expect_equal(1,as.numeric(all(colnames(imitate) == c("prob.r","prob.p","prob.s","tremble","tr(r)","tr(p)","tr(s)")) & all(unlist(imitate[,1:3])==c(1/3,1,0,0,1/3,0,1,0,1/3,0,0,1)) & all(is.na(imitate[,4])) & all(unlist(imitate[,5:7])==rep(c(2,3,4),each=4)) ))

  data.WXZ2014 <- stratEst.data(data = WXZ2014, choice = "choice",
                                input = c("choice"), input.lag = 1,
                                id = "id", game = "game",
                                period = "period")


  model.nash <- stratEst.model(data = data.WXZ2014,
                               strategies = list("nash" = nash))
  model.mixed <- stratEst.model(data = data.WXZ2014,
                                strategies = list("mixed" = mixed))
  model.imitate <- stratEst.model(data = data.WXZ2014,
                                  strategies = list("imitate" = imitate))
  model.mixture <- stratEst.model(data = data.WXZ2014,
                                  strategies = list("nash" = nash,
                                                    "imitate" = imitate))

  models <- list(model.nash, model.mixed, model.imitate, model.mixture)
  compare <- round(do.call(rbind, unlist(lapply(models, stratEst.check), recursive = F)))
  rownames(compare) <- c("model.nash", "model.mixed", "model.imitate",
                         "model.mixture")

  expect_equal(1,as.numeric(all( c(compare) == c(-23730,-23704,-23206,-22358,0,2,1,2,47460,47412,46414,44721,47460,47417,46416,
                                             44725,47460,47417,46416,44728) ) ) )

  t.probs <- stratEst.test(model = model.mixed, par = "probs", values = 1/3)
  expect_equal(1,as.numeric( all( unlist(t.probs) == c(0.3223,0.3566,0.3212,-0.0111,0.0232,-0.0122,0.0014,0.0013,0.0012,-8.0838,17.6404,-10.3417,70,70,70,0,0,0))))

  expect_equal(1,as.numeric(all(round(model.mixture$shares,2) == c(0.58,0.42))))
  expect_equal(1,as.numeric(round(model.mixture$trembles.par,3) == 0.391))

})

test_that("Simulated data" , {
  skip_on_cran()
  set.seed(1)
  lr <- c("left","right")
  mixed <- stratEst.strategy( choices = lr, inputs = lr, num.states = 1 )
  pure <- stratEst.strategy( choices = lr, inputs = lr, prob.choices = c(1,0,0,1),
                             tr.inputs = c(1,2,1,2) )
  strategies <- list( "mixed" = mixed, "pure" = pure )
  p <- runif(1)
  t <- runif(1)/4
  beta <- rnorm(1)
  s <- exp(beta)/sum( 1 + exp(beta) )
  sim.shares <- c(s,1-s)
  mixed$prob.left <- p
  mixed$prob.right <- 1-p
  pure$tremble <- t
  sim.strategies <- list( "mixed" = mixed, "pure" = pure )
  sim.data <- stratEst.simulate( strategies = sim.strategies, shares = sim.shares,
                                 num.ids = 100, num.games = 10, num.periods = rep(5,10) )
  model <- stratEst.model( data = sim.data, strategies = strategies, verbose = F )
  sim.data$intercept <- rep(1,nrow(sim.data))
  model.lcr <- stratEst.model( data = sim.data, strategies = strategies,
                               covariates = "intercept", verbose = F )
  pars <- c(s,1-s,p,1-p,t)
  test.pars <- stratEst.test( model, values = pars )
  expect_equal(1,as.numeric(all(round(unlist(test.pars),4)==c(0.5400,0.4600,0.2711,0.7289,0.0939,-0.0058,0.0058, 0.0056,-0.0056,0.0009,0.0498,0.0498,0.0088,0.0088,0.0061,-0.1160,0.1160,0.6335,-0.6335,0.1443,97.0000,97.0000,97.0000,97.0000,97.0000,0.9079,0.9079,0.5279,0.5279,0.8856))))

  strategy <- sim.data$strategy
  choice <- sim.data$choice
  input <- sim.data$input
  s.sample <- mean(strategy == "mixed")
  p.sample <- mean( choice[strategy == "mixed"] == "left" )
  t.sample <- mean( choice[strategy == "pure"] != input[strategy == "pure"] )
  expect_equal(1,as.numeric(all(round(c(s.sample,1-s.sample,p.sample,1-p.sample,t.sample),4)==c(0.5400,0.4600,0.2711,0.7289,0.0939))))

})

test_that("Replication example DalBo and Frechette, 2011" , {
  skip_on_cran()
  expect_equal(1,as.numeric(1,all( colnames(strategies.DF2011[["TFT"]]) == c("prob.d","prob.c","tremble","tr(cc)","tr(cd)","tr(dc)","tr(dd)") )))
  expect_equal(1,as.numeric(all(strategies.DF2011[["TFT"]]$prob.d == c(0,1))))
  expect_equal(1,as.numeric(all(strategies.DF2011[["TFT"]]$prob.c == c(1,0))))
  expect_equal(1,as.numeric(all(is.na(strategies.DF2011[["TFT"]]$tremble))))
  expect_equal(1,as.numeric(all(strategies.DF2011[["TFT"]]$'tr(cc)' == c(1,1))))
  expect_equal(1,as.numeric(all(strategies.DF2011[["TFT"]]$'tr(cd)' == c(2,2))))
  expect_equal(1,as.numeric(all(strategies.DF2011[["TFT"]]$'tr(dc)' == c(1,1))))
  expect_equal(1,as.numeric(all(strategies.DF2011[["TFT"]]$'tr(dd)' == c(2,2))))
  data.DF2011 <- stratEst.data( data = DF2011, choice = "choice",
                                input = c("choice","other.choice"), input.lag = 1 )
  model.DF2011 <- stratEst.model( data = data.DF2011,
                                  strategies = strategies.DF2011,
                                  sample.id="treatment" , verbose = F )
  expect_equal(1,as.numeric(all(unlist(c(round(do.call(rbind,model.DF2011$shares),2))) == c(0.92,0.78,0.53,0.65,0.11,0.00,0.00,0.08,0.07,0.00,0.30,0.08,0.00,0.04,0.00,0.00,0.27,0.12,0.08,
  0.10,0.38,0.35,0.33,0.56,0.00,0.00,0.02,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.24))))
})
