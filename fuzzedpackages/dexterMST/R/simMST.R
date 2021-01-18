




#' Simulate multistage testing data
#' 
#' Simulates data from an extended nominal response model according to an mst design
#' 
#' @param pars item parameters, can be either: 
#' a data.frame with columns item_id, item_score, beta or a dexter or dexterMST parameters object
#' @param theta vector of person abilities
#' @param test_design data.frame with columns item_id, module_id, item_position
#' @param routing_rules output pf \code{\link{mst_rules}}
#' @param routing 'all' or 'last' routing
#' 
#' 
#' 
sim_mst = function(pars, theta, test_design, routing_rules, routing=c('last','all'))
{
  routing_type = match.arg(routing)
  dat = dexter::r_score(pars)(theta)
  
  nmod=max(routing_rules$module_nbr)
  if(nmod <=2)
    routing_type='last'
  
  test_design$item_id = as.character(test_design$item_id)
  test_design$module_id = as.character(test_design$module_id)
  
  routing_rules$booklet_id = as.character(routing_rules$booklet_id )
  routing_rules$module_id = as.character(routing_rules$module_id)
  
  stopifnot(setequal(test_design$module_id,routing_rules$module_id))
  
  mdlist = split(test_design$item_id, test_design$module_id)
  msum = matrix(0L,length(theta), length(mdlist))
  mdl = names(mdlist)
  colnames(msum) = mdl
  
  for(module_id in mdl)
    msum[,module_id] = rowSums(dat[,mdlist[[module_id]]])

  if(routing_type=='last')
  {
    routing_rules$exit_min = coalesce(routing_rules$exit_min ,0L)
    routing_rules$exit_max = coalesce(routing_rules$exit_max ,as.integer(1e9))
  } else
  {
    routing_rules = routing_rules %>%
      group_by(.data$booklet_id) %>%
      arrange(.data$module_nbr) %>%
      mutate(exit_min = coalesce(.data$exit_min,lag(.data$exit_min,default=0L)),
             exit_max=c(.data$exit_max[-n()],coalesce(.data$exit_max[n()],as.integer(1e9)))) %>%
      ungroup()
  }
  
  lapply(split(routing_rules, routing_rules$booklet_id), function(rl){
    indx = rep(TRUE,length(theta))
    
    if(routing_type == 'last')
    {
      for(i in 1:nrow(rl))
        indx[indx] = between(msum[indx,rl$module_id[i]],rl$exit_min[i],rl$exit_max[i])
    } else
    {
      sm=integer(length(theta))
      for(i in 1:nrow(rl))
      {
        indx[indx] = between(msum[indx,rl$module_id[i]] + sm[indx],rl$exit_min[i],rl$exit_max[i])
        sm[indx] = sm[indx] + msum[indx,rl$module_id[i]]
      }
    }
    
    items = unlist(mdlist[rl$module_id])
    persons = which(indx)
    tibble(person_id = rep(persons,length(items)), 
           item_id = rep(items,each=length(persons)),
           item_score = as.integer(dat[persons,items]))
  }) %>%
    bind_rows(.id='booklet_id') %>%
    mutate_if(is.factor, as.character)
}