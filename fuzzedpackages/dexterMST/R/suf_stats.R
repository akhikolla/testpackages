

suf_stats_nrm = function(respData, max_score)
{
  suf_stats = suf_stats_nrm_c(respData$bid,respData$booklet_score,
                            respData$item_id,respData$item_score, 
                            nlevels(respData$item_id), max_score)
  
  class(suf_stats$ssIS$item_id) = 'factor'
  levels(suf_stats$ssIS$item_id) = levels(respData$item_id)
  
  nit = n_distinct(suf_stats$ssIS$item_id)
  
  if(NROW(filter(suf_stats$ssIS, .data$item_score==0))<nit)
    stop("items need to have minimum score of 0")
  
  suf_stats$ssIS = filter(suf_stats$ssIS, .data$item_score>0)
  
  if(n_distinct(suf_stats$ssIS$item_id)<nit)
    stop("items may not have a maximum score of 0")
  
  
  suf_stats
}



sufstats_im = function(rsp, max_score)
{
  suf_stats =  suf_stats_im_c(rsp$booklet_score, rsp$item_id, rsp$item_score, nlevels(rsp$item_id),max_score)
  
  class(suf_stats$ssIS$item_id) = 'factor'
  levels(suf_stats$ssIS$item_id) = levels(rsp$item_id)
  
  class(suf_stats$plt$item_id) = 'factor'
  levels(suf_stats$plt$item_id) = levels(rsp$item_id)
  
  suf_stats
}


# calibration_design = function(respData, ssI, ssIS, routing)
# {
#   min_max = function(bid, module_nbr)
#   {
#     w = which(routing$bid==bid & routing$module_nbr==module_nbr)
#     
#     if(module_nbr==1 || routing$routing[w] == 'last')
#       return(c(routing$module_exit_score_min[w], routing$module_exit_score_max[w]))
#     
#     
#     c(routing$module_exit_score_min[w] - routing$module_exit_score_min[w-1],
#       routing$module_exit_score_max[w] - routing$module_exit_score_min[w-1])
#   }
#   
#   
#   respData$booklet_design %>%
#     inner_join(ssI, by='item_id') %>%
#     mutate(item_min_score=0L) %>%
#     group_by(bid) %>%
#     do({
#       bk = .
# 
#       out = bk %>% 
#           group_by(module_nbr) %>%
#           do({
#             md = .
#             mm = min_max(md$bid[1],md$module_nbr[1])
#             exit_min = mm[1]
#             exit_max = mm[2]
#             
#             repeat
#             {
#               mmin = sum(md$item_min_score)
#               
#               decr = which(md$item_max_score  - md$item_min_score > exit_max - mmin)
#               for(i in decr)
#               {
#                 for(j in md$last[i]:md$first[i])
#                 {
#                   if(ssIS$item_score[j] <= exit_max - mmin)
#                   {
#                     md$last[i] = j
#                     md$item_max_score[i] = ssIS$item_score[j]
#                     break
#                   }
#                   if(j==md$first[i]) # dit item in deze module kan enkel 0 scores krijgen
#                     md$last[i] = NA_integer_
#                 }
#               }
#               md = filter(md, !is.na(last))
#               
#               msum = sum(md$item_max_score)
#               incr = which(md$item_min_score + msum - md$item_max_score < exit_min)
#               for(i in incr)
#               {
#                 for(j in md$first[i]:md$last[i])
#                 {
#                   if(ssIS$item_score[j] + msum - md$item_max_score[i] >= exit_min )
#                   {
#                     md$item_min_score[i] = ssIS$item_score[j]
#                     md$first[i] = j
#                     break;
#                   }
#                   # else hoeft hier denk ik niet, er zouden geen leerlingen kunnen zijn in dit boekje
#                 }
#               }
#               if(length(incr)==0 && length(decr)==0)
#                 break
#             }
#             md
#           })
#       if(!(nrow(bk) == nrow(out) && all(bk$first == out$first & bk$last == out$last)))
#       {
#         ## to#do: temp remove temp warning
#         warn('temp warning: design adjusted for impossible item scores')
#       }
#       out
#     }) %>%
#     select(bid, module_nbr, item_id, first, last) %>%
#     arrange(bid, module_nbr, first)
# }
# 
# 
# 
