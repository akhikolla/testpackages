context('test polytomous')

library(dplyr)
library(dexter)
library(tidyr)


test_that('discriminations', {
  
  # set a seed for determinism in cran check
  # expect this check to succeed with > 99% of seeds
  set.seed(123)
  
  
  # weird a's, rendering scores impossible
  items = tibble(item_id = sprintf('i%02i',1:15), 
                 item_score = c(5,1,1,2,2,sample(1:5,10,TRUE)),
                 beta = c(runif(5,-.5,.5),runif(5,-1,0),runif(5,0,1)),
                 module_id = c(rep('M',5),rep('E',5),rep('H',5)))
  
  items$beta = items$beta - mean(items$beta)
  
  scoring_rules = bind_rows(tibble(item_id = sprintf('i%02i',1:15), item_score=0,response=0),
                            mutate(select(items, item_id, item_score), response=item_score))
  
  routing_rules = mst_rules(bk1 = M[0:4] --+ E, bk2 = M[5:Inf] --+ H)
  
  theta = rnorm(10000)
  
  sim_dat = r_score(items)(theta)
  rownames(sim_dat) = 1:length(theta)
  
  rsp_data = bind_rows(
    pivot_longer(as_tibble(sim_dat[rowSums(sim_dat[,1:5])<=4,1:10], rownames='person_id'),
                 -person_id, names_to="item_id", values_to='response') %>%
      mutate(booklet_id='bk1'),
    pivot_longer(as_tibble(sim_dat[rowSums(sim_dat[,1:5])>4,c(1:5,11:15)], rownames='person_id'),
                 -person_id, names_to="item_id", values_to='response') %>%
      mutate(booklet_id='bk2'))
  
  rsp_data$test_id = 'weird_a'
  

  db = create_mst_project(":memory:")
  add_scoring_rules_mst(db, scoring_rules)
  
  create_mst_test(db,
                  test_design = items,
                  routing_rules = routing_rules,
                  test_id = 'weird_a')
  
  add_response_data_mst(db, rsp_data)
  
  f=fit_enorm_mst(db)
  
  expect_gt(cor(items$beta,coef(f)$beta), 0.98)
  

})



test_that('poly NR problem', {

  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")
  
  items = coef(fit_enorm(db))
  
  close_project(db)
  
  db = create_mst_project(":memory:")
  add_scoring_rules_mst(db, verbAggrRules)
  
  
  dt = verbAggrData[, unique(as.character(verbAggrRules$item_id))]
  
  dsg = tibble(item_id=colnames(dt)[order(colSums(dt),decreasing=TRUE)],
               module_id=c(rep('A', 9), rep('S',6),rep('B',9)))

  
  routing_rules = mst_rules(bk1 = S[0:4] --+ A, bk2 = S[5:Inf] --+ B)
  
  create_mst_test(db,
                  test_design = dsg,
                  routing_rules = routing_rules,
                  test_id = 'poly',routing='last')
  
  dt$person_id = 1:nrow(dt)
  mod1_score = rowSums(dt[,dsg$item_id[10:15]])
  
  bk1 = dt[mod1_score<=4,c(dsg$item_id[1:15],'person_id')]
  bk2 = dt[mod1_score>4,c(dsg$item_id[10:24],'person_id')]
  
  add_booklet_mst(db,bk1,'poly','bk1')
  add_booklet_mst(db,bk2,'poly','bk2')
  
  f=fit_enorm_mst(db)
  
  tst = items %>%
    inner_join(coef(f),by=c('item_id','item_score')) %>%
    inner_join(dsg,by='item_id') %>%
    mutate(col=as.integer(as.factor(module_id)))
  
  expect_gt(cor(tst$beta.x,tst$beta.y),0.97)
  
  #plot(tst$beta.x,tst$beta.y,col=tst$col)
  #abline(0,1)
  
  
  
})







test_that('poly normal2', {
  
  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")
  
  items = coef(fit_enorm(db))
  
  close_project(db)
  
  db = create_mst_project(":memory:")
  add_scoring_rules_mst(db, verbAggrRules)
  
  
  dt = verbAggrData[, unique(as.character(verbAggrRules$item_id))]
  
  dsg = tibble(item_id=colnames(dt)[order(colSums(dt),decreasing=TRUE)],
               module_id=c(rep('A', 8), rep('S',8),rep('B',8)))
  
  routing_rules = mst_rules(bk1 = S[0:5] --+ A, bk2 = S[6:Inf] --+ B)
  
  create_mst_test(db,
                  test_design = dsg,
                  routing_rules = routing_rules,
                  test_id = 'poly')
  
  dt$person_id = 1:nrow(dt)
  mod1_score = rowSums(dt[,dsg$item_id[9:16]])
  
  bk1 = dt[mod1_score<=5,c(dsg$item_id[1:16],'person_id')]
  bk2 = dt[mod1_score>5,c(dsg$item_id[9:24],'person_id')]
  
  add_booklet_mst(db,bk1,'poly','bk1')
  add_booklet_mst(db,bk2,'poly','bk2')
  
  f=fit_enorm_mst(db)
  
  tst = items %>%
    inner_join(coef(f),by=c('item_id','item_score')) %>%
    inner_join(dsg,by='item_id') %>%
    mutate(col=as.integer(as.factor(module_id)))
    
  #plot(tst$beta.y,tst$beta.x,col=tst$col)
  #abline(0,1)
  
  expect_gt(cor(tst$beta.x,tst$beta.y),0.97)
  
  g=fit_enorm_mst(db,method='Bayes')
  
  tst = items %>%
    inner_join(coef(g),by=c('item_id','item_score')) %>%
    inner_join(dsg,by='item_id') %>%
    mutate(col=as.integer(as.factor(module_id)))
  
  expect_gt(cor(tst$beta,tst$mean_beta),0.97)
  
  #plot(tst$beta,tst$mean_beta,col=tst$col)
  #abline(0,1)
  
  
  
})

