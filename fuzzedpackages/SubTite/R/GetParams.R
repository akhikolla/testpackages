#'Obtains true simulation parameters for each supported distribution function to correspond to a probability of the truth by time T1.
#' @param Family What distribution Family to simulate from. Options include: Exponential,Gamma, Lognormal, Uniform, Weibull.
#' @param ParamNum Parameter index for user set value.
#' @param Param #Groups X #Doses Matrix containing one parameter for each subgroup and dose.
#' @param GroupProb #Groups X #Doses Matrix containing the true toxicity probability by time T1.
#' @param T1 Toxicity observation window.
#' @return A list containing the hyperparameter matrices to input into the SimTrial function. Also plots the hazard of toxicity for each subgroup and dose.
#' @importFrom graphics legend lines par par plot
#' @importFrom stats dlnorm dweibull qnorm
#' @examples
#' GroupProb =matrix(c(.05,.3,.6,.7,.8,.01,.02,.13,.27,.5),nrow=2,byrow=TRUE)
#' ##True Simulation distribution
#' Family="Weibull"
#' T1=6
#' Param = GroupProb*0 + 4 ##Late onset weibull
#' SimTruth = GetParams("Weibull",1,Param,GroupProb,T1)
#'@export
GetParams = function(Family,ParamNum,Param,GroupProb,T1){

  SimTruth = as.list(c(0,0))

  ##Error spot
  if(!(ParamNum %in% c(1,2))){
    message("ParamNum must be 1 or 2")
  }

  if(!(Family == "Uniform")){
    if((!(nrow(Param) == nrow(GroupProb))) || (!(ncol(Param) == ncol(GroupProb)))){
      message("Param and GroupProb are not the same dimension")
    }

  }

  stop=0


  if(Family=="Uniform"){
    SimTruth[[1]]=GroupProb
    SimTruth[[2]]=GroupProb
    ###Plot hazards....
    x=seq(.01,T1,.01)
    if(3>5){
      ##Disable plotting for this
    plot(x,rep(1/(T1-0),length(x)),type="l",main="Uniform Density",ylab="Density",xlab="Time" )
  }
}

  if(Family=="Weibull"){
    ##Check if param1 or param2 is specified
    if(ParamNum==1){
      ##Alpha is specified
      a=Param
      b=T1*(-log(1-GroupProb))^{-1/a}

      if(sum(b<0)>0){
        message("Resulting scale parameter is negative.")
        stop=1
      }


    }else{
      ##Beta is specified
      b=Param
      a=log(-log(1-GroupProb))/(log(T1/b))

      if(sum(b<0)>0){
        message("Resulting shape parameter is negative.")
        stop=1
      }


    }

    if(stop==0){
      ###Plot the densities...
      x=seq(.1,T1,.01)
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      COL1 = 1:15
      COL1 = COL1[-7]


      if(3>5){
      ##Lets not do plotting for now...
      ##Do this for each subgroup...
      for(g in 1:nrow(GroupProb)){
        plot(x,dweibull(x,a[g,1],b[g,1]),type="l",main=paste0("Weibull Density Subgroup ", g),
             ylim=c(0,max(dweibull(x,a[g,],b[g,]))),ylab="Density")
        for(j in 2:ncol(GroupProb)){
          lines(x,dweibull(x,a[g,j],b[g,j]),type="l",col=COL1[j])
        }

        leglabel = "Dose 1"
        for(j in 2:ncol(GroupProb)){
          leglabel = c(leglabel,paste0("Dose ",j))
        }



        legend("topleft",inset=c(1,max(max(dweibull(x,a[g,],b[g,])))),col = COL1,legend =leglabel,lty=1)



      }
}


      SimTruth[[1]]=a
      SimTruth[[2]]=b

      return(SimTruth)

    }

  }





  if(Family=="Lognormal"){
    ##Check if param1 or param2 is specified
    if(ParamNum==1){
      ##mean is specified
      a=Param
      b=(log(T1)-a)/qnorm(GroupProb)

      if(sum(b<0)>0){
        message("Resulting standard deviation parameter is negative.")
        stop=1
      }


    }else{
      ##sd  is specified
      b=Param
      a=log(T1)-b*qnorm(GroupProb)




    }

    if(stop==0){
      ###Plot the densities...
      x=seq(.1,T1,.01)


      if(3>5){

      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      COL1 = 1:15
      COL1 = COL1[-7]
      ##Do this for each subgroup...
      for(g in 1:nrow(GroupProb)){
        plot(x,dlnorm(x,a[g,1],b[g,1]),type="l",
             main=paste0("Lognormal Density Subgroup ", g),
             ylim=c(0,max(dlnorm(x,a[g,],b[g,]))),ylab="Density")
        for(j in 2:ncol(GroupProb)){
          lines(x,dlnorm(x,a[g,j],b[g,j]),type="l",col=COL1[j])
        }

        leglabel = "Dose 1"
        for(j in 2:ncol(GroupProb)){
          leglabel = c(leglabel,paste0("Dose ",j))
        }



        legend("topleft",inset=c(1,max(max(dlnorm(x,a[g,],b[g,])))),col = COL1,legend =leglabel,lty=1)



      }

}

      SimTruth[[1]]=a
      SimTruth[[2]]=b

      return(SimTruth)

    }

  }



  if(Family=="Gamma"){

    if(ParamNum==1){
      a=Param
      b=Param ##Storage for B output
      ##Find the needed parameter using a while loop
      ##Use While loop to get the other probability

      MAX = 0


      for(k in 1:nrow(Param)){
        for(j in 1:ncol(Param)){
          prob=0
          b1=.01
          while(prob<GroupProb[k,j]){
            b1=b1+.01
            prob=pgamma(T1,a[k,j],rate=b1)
          }



          b[k,j]=b1

        }
      }


    }else{

      a=Param
      b=Param ##Storage for B output
      ##Find the needed parameter using a while loop
      ##Use While loop to get the other probability
      for(k in 1:nrow(Param)){
        for(j in 1:ncol(Param)){
          prob=0
          b1=.01
          while(prob<GroupProb[k,j]){
            b1=b1+.01
            prob=pgamma(T1,b1,rate=b[k,j])
          }


          a[k,j]=b1

        }


      }


    }


  }

  if(!(Family %in% c("Gamma","Lognormal","Weibull","Uniform"))){
    message("Not a supported distribution family.")
  }

}
