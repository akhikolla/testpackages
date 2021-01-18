library(stratEst)

test_that("behavior strategies",  {
  set.seed(1)
  N = 200
  Obs = 20

  intercept = 0.2
  intercept2 = 0
  dummy = 1
  dummy2 = 0
  response = 0.25
  tremble = 0.10

  state <- c(1,2)
  r1 <- c(1,0)
  t1 <- c(1,1)
  t2 <- c(2,2)

  TFT1 = stratEst.strategy( inputs = c("1","2") , choices = c("0","1") , num.states = 2 , prob.choices = c(0,1,1,0) , tr.inputs = c(1,2,1,2))
  TFT2 = stratEst.strategy( inputs = c("1","2") , choices = c("0","1") , num.states = 2  , tr.inputs = c(1,2,1,2))

  strategies = list("TFT1"=TFT1,"TFT2"=TFT2)
  covariates = matrix(c(rep(1,N),rep(0,N/2),rep(1,N/2)),N,2)

  coefficient_mat = matrix(c(intercept,dummy),2,1)
  strategy_selection = matrix(NA,N,2)
  strategy_selection[,1] <- 1/( 1 + exp(covariates %*% coefficient_mat) )
  strategy_selection[,2] <- exp(covariates %*% coefficient_mat)/( 1 + exp(covariates %*% coefficient_mat) )

  expected_shares = apply( strategy_selection,2,mean)

  rand_vec = runif(N)
  strategy_vec = ifelse( rand_vec > strategy_selection[,2] , 1 , 2 )
  strat_id = rep( strategy_vec , each = (Obs*2))

  real_shares = rep(NA,2)
  real_shares[1] = sum(as.numeric(strategy_vec==1))/length(strategy_vec)
  real_shares[2] = sum(as.numeric(strategy_vec==2))/length(strategy_vec)

  TFT1 =  matrix(c(1,2,(1-tremble),(0+tremble),1,1,2,2),2,4)
  TFT2 =  matrix(c(1,2,response,response,1,1,2,2),2,4)

  strategy_mat = rbind(TFT1,TFT2)

  id = rep(NA,Obs*2*N)
  game = rep(1,Obs*2*N)
  period = rep(c(1:(Obs*2)),N)
  input = rep(c(rep(1,Obs),rep(2,Obs)),N)
  output = rep(NA,Obs*2*N)
  rand_vec = runif(Obs*2*N)
  pr_vec = rep(NA,Obs*2*N)

  for( i in 1:N ){
    id[ ((i-1)*Obs*2 + 1) : ((i-1)*Obs*2 + Obs*2) ] = rep(i,Obs*2)
    pr_vec[ ((i-1)*Obs*2 + 1) : ((i-1)*Obs*2 + Obs*2) ] = c(rep(strategy_mat[ 1 + (strategy_vec[i]-1)*2 , 2 ],Obs),rep(strategy_mat[ 2 + (strategy_vec[i]-1)*2 , 2 ],Obs))
  }

  choice <- ifelse( rand_vec >= pr_vec , 0 , 1 )

  covar = as.numeric(id>100)

  data <- as.data.frame(cbind(id,game,period,input,choice,covar))
  data <- stratEst.data(data,input = "input")

  model = stratEst.model( data, strategies, covariates = "covar" , select=c("strategies","probs","trembles"), crit = "icl", inner.runs = 100, inner.max = 10, outer.runs = 1,outer.max = 200,outer.tol = 0,lcr.runs = 100, verbose = F)

  expect_equal( 0.350, round(as.numeric(model$shares[1]),3) )
  expect_equal( 0.231, round(as.numeric(model$probs.par[2]),3) )
  expect_equal( 0.110, round(as.numeric(model$trembles.par[1]),3) )
  expect_equal( 1.386, round(as.numeric(model$coefficients.par[1]),3) )
  expect_equal( 2 , length( as.numeric(model$shares ) ) )
  expect_equal( 4 , length( as.numeric(model$probs.par ) ) )
  expect_equal( 1 , length( as.numeric(model$trembles.par ) ) )

})
