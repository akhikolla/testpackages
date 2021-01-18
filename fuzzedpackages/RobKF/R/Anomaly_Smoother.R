Get_ancestor = function(x,x_list,kk){
  
  current_id = x$id
  
  while (x$id > current_id - kk){
    
    old_x = x
    
    x = x_list[[x$id+1-x$horizon]][[x$ancestor_id]]
    
  }
  
  if (x$id == current_id - kk){
    
    return(x)
    
  } 
  
  if (x$id == current_id - kk - 1){
    
    return(list(which_type = old_x$which_type, component = old_x$component, horizon = 1))
    
  } else {
    
    return(list(which_type = "None", component = 0,horizon = 1))
    
  }
  
}

Get_k_ancestors = function(jj,x,kk){
  
  out = lapply(x[[jj]], Get_ancestor, kk=kk, x_list=x)
  
  return(out)
  
}

Anomaly_Smoother = function(x,time=NULL,horizon=NULL){
  
  if (is.null(time)){
    time = length( x[["particles"]])-1
  }
  
  if (is.null(horizon)){
    horizon = x$horizon - 1
  }
  
  smoothed_particle_list = lapply(as.list((horizon+1):(time+1)),x=x[["particles"]],Get_k_ancestors,kk = horizon)
  
  for (ii in 1:horizon){

    add   = Get_k_ancestors(time+1, x[["particles"]], horizon-ii)

    smoothed_particle_list = c(smoothed_particle_list,list(add))

  }
  
  x[["particles"]] = smoothed_particle_list
  
  return(x)
  
}