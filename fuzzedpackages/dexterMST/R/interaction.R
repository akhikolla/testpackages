

#' Fit the interaction model on a single multi-stage booklet
#' 
#' @param db a db handle
#' @param test_id id of the test as defined in \code{\link{create_mst_test}}
#' @param booklet_id id of the booklet as defined in \code{\link{create_mst_test}}
#' 
fit_inter_mst = function(db, test_id, booklet_id)
{
  prd = list(test_id = test_id, booklet_id = booklet_id)
  design = dbGetQuery(db, 
                      'SELECT * FROM Booklet_design WHERE test_id=:test_id AND booklet_id=:booklet_id 
                      ORDER BY module_nbr;', prd)
  
  routing = dbGetQuery(db,'SELECT routing FROM Tests WHERE test_id=:test_id;',list(test_id=test_id))$routing
  
  md = dbReadTable(db,'module_design') %>%
    inner_join(design, by=c('test_id','module_id')) 
  
  md$item_id = ffactor(md$item_id)
  
  
  # to~do: what about missing scorecat in a booklet which exists in others?
  
  rsp = dbGetQuery(db,
                    'SELECT person_id, item_id, item_score
                      FROM Responses 
                    INNER JOIN Scoring_rules USING(item_id, response)
                    WHERE test_id=:test_id AND booklet_id=:booklet_id;',
                    prd)
  if(NROW(rsp)==0)
    stop(paste('no data for test:', test_id,'and booklet:',booklet_id))
  
  # to~do: take levels from db/module_design
  rsp$person_id = ffactor(rsp$person_id, as_int=TRUE)
  rsp$item_id = ffactor(rsp$item_id)
  
  if(is.unsorted(rsp$person_id))
    rsp = arrange(rsp, .data$person_id)
  

  rsp$booklet_score = im_booklet_score(rsp$person_id, rsp$item_score)
  
  # to~do:
  # een predicaat hier toestaan (zolang het niet leidt tot meer design complexiteit)
  
  max_score = dbGetQuery(db,
                          'SELECT MAX(item_score) 
                            FROM Scoring_rules
                              INNER JOIN module_design USING(item_id)
                                INNER JOIN Booklet_design USING(test_id, module_id)
                            WHERE test_id=:test_id AND booklet_id=:booklet_id;', prd)[1,1]
  
  suf = sufstats_im(rsp,max_score)
  
  ssIS = suf$ssIS %>%
    inner_join(md,by='item_id') %>%
    arrange(.data$module_nbr, .data$item_id)
  
  ssI = ssIS %>%
    mutate(indx = row_number()) %>%
    group_by(.data$item_id,.data$module_nbr) %>%
    summarise(first = min(.data$indx), last = max(.data$indx), nCat = n(), sufC=sum(.data$sufC)) %>%
    ungroup() %>%
    arrange(.data$first)
  
  if(routing=='all')
    bkl_max = last(design$module_exit_score_max)
  else
    bkl_max = sum(design$module_exit_score_max)
  
  plt = suf$plt
  
  scoretab = suf$scoretab[1:(bkl_max+1)]

  res = Estim_MST(a = pull(ssIS, .data$item_score), 
                  first = split(ssI$first, ssI$module_nbr), 
                  last = split(ssI$last, ssI$module_nbr), 
                  min_scores = pull(design, .data$module_exit_score_min), 
                  max_scores = pull(design, .data$module_exit_score_max), 
                  sufI = pull(ssIS, .data$sufI), 
                  sufC = pull(ssI, .data$sufC),
                  scoretab = scoretab, 
                  routing = routing)
  
  regs = res$regs
  res$regs = NULL
  
  # regressions for plots
  regs$observed = plt
  
  regs$itrRM = rowsum(regs$ctrRM * ssIS$item_score, ssIS$item_id, reorder=FALSE)
  regs$itrIM = rowsum(regs$ctrIM * ssIS$item_score, ssIS$item_id, reorder=FALSE)
  

  out = list(inputs = list(ssI = ssI, ssIS = ssIS, scoretab = tibble(booklet_score=0:bkl_max,n=scoretab), 
                           design = design, routing = routing,
                           test_id=test_id, booklet_id=booklet_id), 
             pars = res, regs = regs)
  class(out) = append('im_mst', class(out))
  out
}

#' plots for the interaction model
#' 
#' @param x output of \code{\link{fit_inter_mst}}
#' @param item_id id of the item to plot
#' @param show.observed plot the observed mean item scores for each test score
#' @param curtains percentage of most extreme values to cover with curtains, 0 to omit curtains
#' @param zoom if TRUE, limits the plot area to the test score range allowed by the routing rules
#' @param ... further arguments to plot
#' 
plot.im_mst = function(x, item_id=NULL, show.observed = TRUE, curtains = 10, zoom = FALSE, ...){
  
  user.args = list(...)
  zI = x$regs$itrIM
  zR = x$regs$itrRM
  
  if(x$inputs$routing=='all')
  {
    route_max=last(x$inputs$design$module_exit_score_max)
    route_min=last(x$inputs$design$module_exit_score_min)
  } else
  {
    route_max=sum(x$inputs$design$module_exit_score_max)
    route_min=sum(x$inputs$design$module_exit_score_min)
  }
  booklet_max = sum(x$inputs$ssIS$item_score[x$inputs$ssI$last])
  

  if(is.null(item_id))
    item_id=x$inputs$ssI$item_id
  
  for(itm in item_id)
  {
    if(! itm %in% x$inputs$ssI$item_id)
      stop(paste0('item "',itm,'" not found'))
    
    max_item = max(filter(x$inputs$ssIS, .data$item_id == itm) %>% pull(.data$item_score))
    
    mdl_nbr = filter(x$inputs$ssI, .data$item_id == itm) %>% pull(.data$module_nbr)
    mdl_id = filter(x$inputs$design, .data$module_nbr == mdl_nbr) %>% pull(.data$module_id)
    
    
    plt.rng = c(0, booklet_max)
    if(zoom)
      plt.rng = c(route_min, route_max)
    
    do.call(plot, modifyList(list(x = plt.rng,y = c(0, max_item), type="n",
                                  main = paste0(itm, ', module ', mdl_nbr, " (", mdl_id, ")") , 
                                  xlab = 'test score', ylab = 'item score', bty='l'),
                             user.args))
    
    
    if(curtains > 0)
    {
      qnt = x$inputs$scoretab %>%
        filter((cumsum(.data$n) <= sum(.data$n) * (100 - 0.5 * curtains)/100) & (cumsum(.data$n) >= sum(.data$n) * 0.5 * curtains/100))  %>%
        summarise(mn = min(.data$booklet_score), mx = max(.data$booklet_score))
      
      usr = graphics::par('usr')
      graphics::rect(usr[1], usr[3], pull(qnt, .data$mn), usr[2], col="#EEEEEE", border=NA,xpd=FALSE,lwd=0)
      graphics::rect(pull(qnt, .data$mx), usr[3], usr[2], usr[4], col="#EEEEEE", border=NA,xpd=FALSE,lwd=0)
      abline(h=usr[3])
      abline(v=usr[1])
    }
    
    if(show.observed) 
    {
      obs = filter(x$regs$observed, .data$item_id == itm)
      points(pull(obs, .data$booklet_score), pull(obs, .data$mean_item_score),col="coral",pch=20)
    }
    if(zoom)
    {
      lines(route_min:route_max, zI[row.names(zI)==itm,(route_min:route_max)+1 ], col="gray80", lwd=3)
      lines(route_min:route_max, zR[row.names(zR)==itm,(route_min:route_max)+1])
    } else
    {
      i = zI[row.names(zI)==itm,]
      r = zR[row.names(zI)==itm,]
      rng = route_min:route_max
      if(route_min>0)
      {
        i=c(0,i)
        r=c(0,r)
        rng=c(0:route_min, rng)
      }
      if(route_max<booklet_max)
      {
        i=c(i,rep(0,booklet_max-route_max+1))
        r=c(r,rep(0,booklet_max-route_max+1))
        rng = c(rng,route_max:booklet_max)
        # this could be done outside of loop
      }
      
      lines(rng,i,col="gray80", lwd=3)
      lines(rng,r)
       
    }

    
  }
  invisible(NULL)
}


print.im_mst = function(x, ...)
{
  routing = paste0(x$inputs$design$module_id, 
                   paste0('[',x$inputs$design$module_exit_score_min,':',x$inputs$designmodule_exit_score_max,']'),
                  collapse=' --+ ')
  
  cat(paste0('Interaction model parameters for booklet "',x$inputs$booklet_id,'" in test "',x$inputs$test_id,
             '"\n\ndesign: ',routing,
             '\n\nn items: ',nrow(x$inputs$ssI),'\nn persons: ',sum(x$inputs$scoretab$n),'\n'))
  invisible(x)
}


coef.im_mst = function(object, ...)
{
  object$inputs$ssIS %>%
    select(.data$item_id, .data$item_score) %>%
    add_column(bRM = object$pars$bRM, cRM = object$pars$cRM, bIM = object$pars$bIM, cIM = object$pars$cIM,
               se.c = object$pars$se.c, fit.IM = object$pars$fit.stats)
}


