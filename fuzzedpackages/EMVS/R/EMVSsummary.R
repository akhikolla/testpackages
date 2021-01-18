EMVSsummary<-function(result){
        
  if(result$independent == F) {
        posts=result$log_g_function
        posts[is.finite(posts)==F]=NaN
        models=apply(result$prob_inclusion,2,function(x){as.numeric(x>0.5)})
        list(log_g_function=posts,models=models)
  } else {
    models=apply(result$prob_inclusion,2,function(x){as.numeric(x>0.5)})
    list(models=models)
  }

}
