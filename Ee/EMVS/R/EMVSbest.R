
EMVSbest<-function(result){

  if(result$independent == F){
        posts=result$log_g_function
        posts[is.finite(posts)==F]=NaN
	which<-which.max(posts)
        logpost<-posts[which]
        print("Best Model Found")
        list(log_g_function=logpost,indices=(1:dim(result$betas)[2])[result$prob_inclusion[which,]>0.5])
  } else {
  print("Log posterior not available for independent prior case")
}
}


