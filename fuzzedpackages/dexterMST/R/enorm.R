##########################################
#' Fit the extended nominal response model on MST data
#'
#' Fits an Extended NOminal Response Model (ENORM) using conditional maximum likelihood (CML)
#' or a Gibbs sampler for Bayesian estimation; both adapted for MST data
#' 
#' 
#' @param db an dextermst db handle
#' @param predicate logical predicate to select data to include in the analysis, see details
#' @param fixed_parameters data.frame with columns `item_id`, `item_score` and `beta`
#' @param method If CML, the estimation method will be Conditional Maximum Likelihood. 
#' If Bayes, a Gibbs sampler will be used to produce a sample from the posterior.
#' @param nDraws Number of Gibbs samples when estimation method is Bayes.
#' 
#' @return
#' object of type 'mst_enorm'. Can be cast to a data.frame of item parameters 
#' using function `coef` or used in dexter's \code{\link[dexter]{ability}} functions
#' 
#' @details 
#'  You can use the predicate to include or omit responses from the analysis, e.g.
#'  `p = fit_enorm_mst(db, item_id != 'some_item' & student_birthdate > '2005-01-01')`
#'  
#' DexterMST will automatically correct the routing rules for the purpose of the current analysis. 
#' There are some caveats though. Predicates that lead to many different designs, e.g. a predicate like
#' \code{response != 'NA'} (which is perfectly valid but can potentially create 
#' almost as many tests as there are students) might take very long to compute. 
#' 
#' Predicates that remove complete modules from a test, e.g. \code{module_nbr !=2} or \code{module_id != 'RU4'} 
#' will cause an error and should be avoided. 
#' 
#' 
#' @references 
#' Zwitser, R. J. and Maris, G (2015). Conditional statistical inference with multistage testing designs. 
#' Psychometrika. Vol. 80, no. 1, 65-84.
#' 
#' Koops, J. and Bechger, T. and Maris, G. (in press); Bayesian inference for multistage and other 
#' incomplete designs. In Research for Practical Issues and Solutions in Computerized Multistage Testing.
#' Routledge, London. 
#' 
fit_enorm_mst = function(db, predicate = NULL, fixed_parameters = NULL, method=c("CML", "Bayes"), nDraws=1000)
{
  method = match.arg(method)
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env() 
  fit_enorm_mst_(db, qtpredicate, env, fixed_parameters, method, nDraws)
}

fit_enorm_mst_ = function(db, qtpredicate, env, fixed_parameters=NULL, method=c("CML", "Bayes"), nDraws=1000)
{
  method = match.arg(method)
  respData = get_mst_data(db, qtpredicate, env)
  
  if(NROW(respData$x) == 0)
    stop('selection is empty, no data to analyze')
  
  plt = respData$suf_stats$plt
  ssIS = respData$suf_stats$ssIS
  
  ssI = ssIS %>%
    mutate(rn = row_number()) %>%
    group_by(.data$item_id) %>%
    summarise(first = min(.data$rn), last = max(.data$rn), item_max_score = max(.data$item_score)) %>%
    ungroup() %>%
    mutate_if(is.double, as.integer) %>%
    arrange(.data$first) 
  

  routing = respData$routing %>%
    select(-.data$test_id,-.data$booklet_id) %>%
    arrange(.data$bid, .data$module_nbr)
  
  # calibration_design gives wrong results
  # and weirdly appears to be completely unnecessary
  #design = calibration_design(respData, ssI, ssIS, routing) %>%
  #  arrange(bid, module_nbr, first)

  design = respData$booklet_design %>%
    inner_join(ssI, by='item_id') %>%
    inner_join(routing, by=c('bid','module_nbr')) %>%
    select(.data$bid, .data$module_nbr, .data$first, .data$last, .data$item_id) %>%
    arrange(.data$bid, .data$module_nbr, .data$first)
  
  if(!is_connected(design$bid, design$first, design$last))
  {
    # check is not perfect yet, I think actual responses within booklets are needed according to Fischer
    stop('Your design is not connected',call.=FALSE)
  }
  
  modules = design %>%
    count(.data$bid, .data$module_nbr, name = 'nit') %>%
    inner_join(routing, by=c('bid','module_nbr')) %>%
    arrange(.data$bid, .data$module_nbr)
  
  booklets = routing %>%
    group_by(.data$bid, .data$routing) %>%
    arrange(.data$module_nbr) %>%
    summarise(last_max = last(.data$module_exit_score_max), sum_max = sum(.data$module_exit_score_max), 
              last_min = last(.data$module_exit_score_min), sum_min = sum(.data$module_exit_score_min), 
              nmod=n()) %>%
    ungroup() %>%
    mutate(max_score = if_else(routing=='all', .data$last_max, .data$sum_max),
           min_score = if_else(routing=='all', .data$last_min, .data$sum_min)) %>%
    select(.data$bid, .data$routing, .data$nmod, .data$max_score, .data$min_score) %>%
    arrange(.data$bid)
  
  scoretab = plt %>%
    select(.data$bid, .data$booklet_score,.data$N) %>%
    distinct(.data$bid, .data$booklet_score, .keep_all=TRUE)  %>%
    right_join(tibble(bid = rep(booklets$bid, booklets$max_score+1),
                      booklet_score = unlist(lapply(booklets$max_score, function(s) 0:s))),
               by = c('bid','booklet_score')) %>%
    mutate(N = coalesce(.data$N,0L)) %>%
    arrange(.data$bid, .data$booklet_score)

  
  # to~do: accept range of different inputs, like in dexter
  if (!is.null(fixed_parameters))
  {
    if(length(setdiff(c('item_id','item_score','beta'), colnames(fixed_parameters))) > 0)
      stop('fixed_parameters needs to have the following columns: (item_id, item_score, beta)')
    
    fixed_not_found = setdiff(fixed_parameters$item_id, as.character(ssI$item_id) )
    
    if(length(fixed_not_found) > 0)
    {
      message('some of your fixed parameters have no match in the data, they will be ignored')
      print(fixed_not_found)
    }
    fixed_parameters = fixed_parameters %>%
      arrange(.data$item_id, .data$item_score) %>%
      mutate(rn = row_number())
    
    ffl = fixed_parameters %>%
      group_by(.data$item_id) %>%
      summarise(first = min(.data$rn), last=max(.data$rn))
    
    dx = dexter.toDexter(fixed_parameters$beta, fixed_parameters$item_score, ffl$first, ffl$last, 
                    re_normalize=FALSE)
    dx_b = dx$est$b[-dx$inputs$ssI$first]
    
    fixed_parameters$b = dx_b
    
    fixed_parameters = ssIS %>%
      select(.data$item_id, .data$item_score) %>%
      mutate(item_id = as.character(.data$item_id)) %>%
      left_join(fixed_parameters, by=c('item_id','item_score')) %>%
      arrange(.data$item_id, .data$item_score) %>%
      pull(.data$b)
    
    if(!anyNA(fixed_parameters))
      stop('Nothing to calibrate, all item parameters have been fixed')
    if(all(is.na(fixed_parameters)))
      stop('none of your fixed parameters belong to items in the data')
  }
  

  calibrate = if.else(method=='CML', Calibrate_MST, Calibrate_Bayes_MST) 
  
  res = calibrate( first = ssI$first, last = ssI$last,
                   sufI = ssIS$sufI, a = ssIS$item_score, 
                   scoretab=scoretab,
                   booklets=booklets,
                   modules = modules,
                   design = design,
                   fixed_b = fixed_parameters,
                   nIter=nDraws)
  
  bids = distinct(respData$routing, .data$bid, .data$test_id) %>%
    group_by(.data$test_id) %>%
    mutate(booklet_id=gsub(paste0('^',.data$test_id[1],'-?'),'',as.character(.data$bid),perl=TRUE))

  
  out = list(mst_inputs = list(ssI=ssI, ssIS=ssIS, design=inner_join(design,bids,by='bid'), 
                               booklets=booklets, modules=modules),
             abl_tables = list(),
             inputs = list(method = method,
                           ssI = mutate(ssI,
                                        first = as.integer(first + rank(first) - 1L), 
                                        last = as.integer(last + rank(last))),
                           ssIS = bind_rows(ssIS, 
                                            tibble(item_id=pull(ssI, 'item_id'), item_score = 0L, 
                                                   sufI = as.integer(NA))) %>% 
                             arrange(.data$item_id, .data$item_score),
                           plt = rename(plt, booklet_id='bid'),
                           has_fixed_parms = !is.null(fixed_parameters)))
             
  out$inputs$design = design %>%
    select(booklet_id=.data$bid,.data$item_id) %>%
    inner_join(out$inputs$ssI[,c('item_id','first','last')],by='item_id')

  
  
  if(method=='Bayes')
  {
    colnames(res$beta) = paste(ssIS$item_id, ssIS$item_score)
    out$mst_est = list(beta=res$beta)
    
  } else
  {
    out$mst_est = bind_cols(select(ssIS, .data$item_id, .data$item_score), 
                            res[names(res) != 'acov.beta'])
  }
  renorm = is.null(fixed_parameters) || all(is.na(fixed_parameters))
  # for compatibility with dexter
  out$est = dexter.toDexter(out$mst_est$beta, out$mst_inputs$ssIS$item_score, 
                            out$mst_inputs$ssI$first, out$mst_inputs$ssI$last,
                            re_normalize=renorm)$est
  
  if(method=='CML')
  {
    out$est$acov.beta = res$acov.beta
    out$est$b = drop(out$est$b)
  }
  class(out) = append(c('mst_enorm', 'prms') ,class(out))
  
  out$abl_tables$mle = ability_tables(out, standard_errors=FALSE) %>%
    filter(is.finite(.data$theta)) 
  
  out$abl_tables$mle$booklet_id = ffactor(out$abl_tables$mle$booklet_id, levels=levels(routing$bid))
  
  out
}


print.mst_enorm = function(x, ...)
{
  msg = paste0('Parameters for the Extended Nominal Response model based on Multi Stage calibration\n\n',
               nrow(x$inputs$ssI), ' items\n\n',
               'Use coef() to retrieve the item parameters\n')
  cat(msg)
  invisible(msg)
}






