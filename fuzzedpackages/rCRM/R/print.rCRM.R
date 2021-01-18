

###################
#####  Print  #####
###################

print.rCRM=function(x, digits=3, ...) {

  tem=switch(x$family,"2P"="2-parameter")
  cat("\nRegularized CRM: ",tem)

  cat("\n\nThe target probability is",x$tp,"\n")

  cat("\nThe estimated probability of DLT:\n\n")
  temi=rbind(as.character(x$dose0), sprintf(paste0("%.",digits,"f"),x$prob), "")
  temi[3,x$dose.closest]="*"
  colnames(temi)=paste0("D",1:ncol(temi))
  rownames(temi)=c("Dose","Prob","Closest")
  temi=noquote(temi)
  print(temi)
}


