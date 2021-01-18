
# restructures a,b,first,last for functions that work on a single booklet 
# does NOT change order of items
single_booklet = function(b,a,first,last)
{
  indx = unlist(mapply(':',first, last,SIMPLIFY=FALSE))
  
  ncat = last-first+1L
  last = cumsum(ncat)
  first = last-ncat+1L
  
  
  list(b=b[indx],a=a[indx], first=first, last=last)

}

# expected scores on domains given an mst design

E_profile_MS_enorm = function(b, a, first,last, mod_min, mod_max, mnit, AB,routing) 
{
  max_score = sum(a[last])
  
  first = as.integer(first-1L)
  last = as.integer(last-1L)
  
  g = elsym_C(ROUTING[[routing]], b, a, first, last, mod_min, mod_max, mnit, max_score)
  
  sapply(sort(unique(AB)), function(A)
  {
    ab = if_else(AB==A,0L,1L)
    
    prof_enorm(b, a, first,last,ROUTING[[routing]], mnit, mod_min, mod_max, max_score, ab)/g
  })
  
}

# to~do (volgende versie), allow a design as a list of routing, design, etc. 

##########################################
#' Profile analysis
#'
#' Expected and observed domain scores per booklet and test score
#'
#' @param parms An object returned by \code{\link{fit_enorm_mst}}
#' @param item_property the name of the item property used to define the domains.
#' @param domains data.frame with column item_id and a column whose name matches `item_property` 
#' 
profile_tables_mst = function(parms, domains, item_property)
{
  if(!item_property %in% colnames(domains))
    stop(paste('column', item_property, 'not found in domains'))
  
 domains[[item_property]] = factor(domains[[item_property]])
 domains$item_id = factor(domains$item_id,levels=levels(parms$mst_inputs$design$item_id))
  
  design = parms$mst_inputs$design %>% 
    left_join(domains[,c('item_id',item_property)], by='item_id') %>%
    arrange(.data$bid,.data$module_nbr,.data$first)
  
  if(anyNA(design[[item_property]]))
    stop('all items need to be categorized')
  
  
  b = parms$mst_est$b
  a = parms$mst_inputs$ssIS$item_score
  
  out = design %>%
    group_by(.data$bid, .data$test_id, .data$booklet_id) %>%
    do({
      sbk = single_booklet(b,a, .$first, .$last)
      
      mod = filter(parms$mst_inputs$modules, .data$bid==.$bid[1])
      
      p = E_profile_MS_enorm(sbk$b, sbk$a,sbk$first, sbk$last, 
               mod$module_exit_score_min, mod$module_exit_score_max, 
               count(., .data$module_nbr)$n, as.integer(.[[item_property]]),
               mod$routing[1])
        
      if(mod$routing[1]=='all')
        mod = mod[nrow(mod),]
      
      mn =  sum(mod$module_exit_score_min)
      mx = sum(mod$module_exit_score_max)
      
      tibble(booklet_score=rep(mn:mx, n_distinct(.[[item_property]])),
             item_domain = rep(sort(unique(.[[item_property]])),each=mx-mn+1),
             expected_domain_score=as.vector(p[1+(mn:mx),])) %>%
        arrange(.data$booklet_score, .data$item_domain)

  }) %>%
    ungroup() %>%
    select(-.data$bid) %>%
    mutate_if(is.factor,as.character)
  
  names(out)[names(out)=='item_domain'] = item_property
  out
}
