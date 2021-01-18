regmed_dat_check <-
function(x,y,mediator){
  
  ###  mediator names ###
  mediator.names<-colnames(mediator)
  if(is.null(mediator.names)) mediator.names<-paste0("med.",1:ncol(mediator))
  
  ### check agreement of object dimensions ###
  if(any(c(length(x),length(y))!=nrow(mediator))) stop("dimensions of x,y and mediator do not agree")

  ### deal with missing data ###

  if(any(c(is.na(x),is.na(mediator),is.na(y)))) {

     if(options()$na.action=="na.pass") {

	stop("na.pass not allowed with missing data")

     }else{

	if(options()$na.action=="na.fail"){

		stop("missing data in x,y or mediator")

	}else{

		keep.obs <- !(is.na(y) | is.na(x) | rowSums(is.na(mediator)))
		x<-x[keep.obs]
		mediator<-mediator[keep.obs,]
		y<-y[keep.obs]

        }
      }
  }

  return(list(x=x,y=y,mediator=mediator,mediator.names=mediator.names))

}
