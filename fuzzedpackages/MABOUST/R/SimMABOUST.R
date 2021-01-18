#' Simulate the MABOUST Trial design.
#'
#' Simulates trial replicates of the MABOUST trial and reports Operating Characteristics (OCs).
#' @param nTreat Number of treatments in consideration, i.e. K.
#' @param nCat Number of ordinal outcome categories, i.e. J.
#' @param UT Vector of numerical utility scores to give outcomes 1,...,J.
#' @param DeltaVEC Vector of \eqn{\bf{\Delta}} values to test.
#' @param gamma Length 3 vector of cutoff parameters.
#' @param PSPIKE Prior probability of a pairwise null effect.
#' @param B Number of MCMC iterations to perform.
#' @param nSims Number of trial replications to complete.
#' @param NLOOK Vector containing how many patients should be evaluated before each interim decision.
#' @param PROBS K-list of J-vectors containing ordinal outcome probabilities.
#' @param XPROB List of matrices containing discrete values various covariates can take, along with their probabilities.
#' @param Beta Covariate Effect Vector on Outcome.
#' @return The set of active treatments to continue, an optimal treatment, or a set of equally optimal treatments. Also reports posterior mean utilities and ordinal outcome probabilities as well as pairwise comparisons of utility similarity, when appropriate.
#' @references
#'  Chapple, A.G., Bennani, Y., Clement, M. (2020). "MABOUST: A Multi-Armed Bayesian Ordinal Outcome Utility-Based Sequential Trial". Submitted.
#' @examples
#' ##Clinical Parameters
#' nCat = 6
#' nTreat = 3
#' UT = c(0,10,20,80,90,100)  ###Utilities
#' DeltaVEC  = c(5,10)   ###Vector of deltas to try
#' NLOOK = c(20,50)  ###Interim Looks
#' ###Which treatments are active?
#' ACTIVE = c(1,0,1) ###Treatments 1, 3 are active
#' FUTILITY = 1 ###Futility look is allowed.
#' ###Design parameters
#' gamma= c(.5, .05, .05)
#' PSPIKE = .9
#' set.seed(1)
#' ##Generate Random Data
#' n=300
#' Y=sample(1:nCat,n,replace=TRUE)
#' T1 = sample(1:nTreat,n,replace=TRUE)
#' XPROB = as.list(rep(NA,3))
#' XPROB[[1]]=rbind(0:10,round(dpois(0:10,2),2)) ###CCI
#' XPROB[[2]]=rbind(c(-1,0,1),c(.5,.4,.1)) ###O2 Status
#' XPROB[[3]]=rbind(c(-2,-1,0,1),c(.27,.38,.18,.17))
#' Beta =
#' ###Number of iterations
#' B=100
#' ##Get Simulation Parameters
#'  #' ##Get Simulation Parameters
#'  PROBS = as.list(rep(NA,3))
#'  PROBS[[1]]=c(.33,.11,.42,.02,.11,.01)
#'  PROBS[[2]]=c(.24,.11,.48,.05,.11,.01)
#'  PROBS[[3]]=c(.14, .20, .48, .03, .12, .03)
#'  Beta=c(-.13, -.07, -.10)
#'  nSims=1 ##Number of sims to run
#'  SimMABOUST(nSims,NLOOK, nTreat,nCat, UT, DeltaVEC,gamma,PSPIKE, B, PROBS, Beta, XPROB)
#' @export
SimMABOUST=function(nSims,
                    NLOOK, ##Sample sizes required for interim looks
                    nTreat, ##Number of treatments in consideration, i.e. K.
                    nCat, ##Number of ordinal outcome categories, i.e. J.
                    UT, ##Vector of numerical utility scores to give outcomes 1,...,J.
                    DeltaVEC, ##Vector of \Delta values to test.
                    gamma, ##Vector of cutoff parameters.
                    PSPIKE, ###Prior probability of a pairwise null effect.
                    B, ##Number of MCMC iterations to perform.
                    PROBS, ###List of outcome probabilities.
                    Beta, ##Vector of covariate X relationship on outcome.
                    XPROB ##List of probability distributions of each X outcome.

){

  PROBSTRUE = PROBS
  UTTRUE=rep(NA,nTreat)
  for(j in 1:nTreat){
    UTTRUE[j]=PROBS[[j]]%*%UT
  }





  ##Storage matrices
  SIZE=rep(NA,nSims)
  INF=rep(1,nTreat)
  INFPROBS=as.list(rep(NA,nSims))
  WORST = matrix(nrow=nSims,ncol=nTreat)
  TREATMENTSIZES = WORST
  UTSTORE=WORST
  CD=rep(NA,nSims)

  DATAHOLD=as.list(rep(NA,nSims))

  for(rep1 in 1:nSims){

    if(rep1%%100==0){
      cat(paste0("Replication # ", rep1,". Time: ", Sys.time(),"
"))
    }

    SKIP=0
    SKIPFUT=0


    n=NLOOK[1] ##Sample size of first look

    ###################
    ###Generate Data###
    ###################
    ##Treatment
    T1=sample(1:nTreat,n,replace=T)
    X = matrix(nrow=n,ncol=length(XPROB))

    for(j in 1:ncol(X)){
      X[,j]=sample(XPROB[[j]][1,],n,prob=XPROB[[j]][2,], replace=TRUE)
    }



    VALS = X%*%Beta

    Y=rep(NA,n)
    ##Generate outcome for each i....
    for(i in 1:n){
      prob1 = PROBS[[T1[[i]]]]
      cumprob = cumsum(prob1)[-length(prob1)]
      ETA = log(cumprob/(1-cumprob))
      ETA = ETA + VALS[i]
      cumprob = c(0,exp(ETA)/(1+exp(ETA)),1)
      prob1=diff(cumprob)
      ###Generate Y
      Y[i]=sample(1:nCat,1,prob=prob1)
    }




    ###########
    ###MCMC####
    ###########
    ####Hyperparameters
    B=2000 ##Number of iterations in MCMC
    RESULTS=MCMC_MABOUST(Y-1, T1-1,   X,   B, nTreat,nCat,PSPIKE)





    ###Get the Utilities###
    UTMAT = matrix(ncol=nTreat,nrow=nrow(RESULTS[[1]]))
    for(b in 1:nrow(UTMAT)){
      HOLD=GetProbs(nCat,RESULTS[[1]][b,])
      for(j in 1:nTreat){
        UTMAT[b,j]=HOLD[[j]]%*%UT
      }
    }


    ###Remove first half
    UTMAT=UTMAT[(round(nrow(UTMAT)/2):nrow(UTMAT)),]



    INF=rep(1,length(INF))
    KEEP=rep(1,length(INF))


    ################################
    ###Should we drop a treatment?##
    ################################
    for(j in 1:nTreat){
      IND=1:nTreat
      IND=IND[-j]
      KEEP1=KEEP[-j]
      IND = IND[KEEP1==1]
      ##Sequentially check the INF...
      if(length(IND)>0){
        for(m2 in 1:length(DeltaVEC)){
          INF[j]=mean(UTMAT[,j]+DeltaVEC[m2]<UTMAT[,IND[1]])
          if(length(IND)>1){
            for(m1 in 2:length(IND)){
              INF[j]=max(INF[j],mean(UTMAT[,j]+DeltaVEC[m2]<UTMAT[,IND[m1]]))
            }
          }
          if(INF[j]>CUTOFF(DeltaVEC[m2],length(Y),sum(KEEP),6,gamma)){
            KEEP[j]=0
            break
          }
        }
      }
    }


    if(sum(KEEP)==1){
      SKIP=1
    }







    if(SKIP==0){
      for(k in 2:length(NLOOK)){
        #################
        T2=0
        Y1=0
        CCI=0
        AGE=0
        X1=0
        ###Next 100...
        T2=sample(which(KEEP==1),NLOOK[k]-NLOOK[k-1],replace=T)
        T1=c(T1,T2)
        Y1=T2



        X1 = matrix(nrow=length(T2),ncol=length(XPROB))
        for(j in 1:ncol(X1)){
          X1[,j]=sample(XPROB[[j]][1,],length(T2),prob=XPROB[[j]][2,],replace=TRUE)
        }


        VALS = X1%*%Beta

        ##Generate outcome for each i....
        for(i in 1:length(T2)){
          prob1 = PROBS[[T2[[i]]]]
          cumprob = cumsum(prob1)[-length(prob1)]
          ETA = log(cumprob/(1-cumprob))
          ETA = ETA + VALS[i]
          cumprob = c(0,exp(ETA)/(1+exp(ETA)),1)
          prob1=diff(cumprob)
          ###Generate Y
          Y1[i]=sample(1:nCat,1,prob=prob1)
        }





        Y=c(Y,Y1)
        X=rbind(X,X1)


        RESULTS=MCMC_MABOUST(Y-1,   ####Outcome vector (minus 1 for c++, i.e. 0,...,J-1)
                             T1-1,  ##Treatment indicator vector
                             X, ##Covariate MAtrix
                             B,nTreat,nCat,PSPIKE) ##Number of iterations for the MCMC




        UTMAT = matrix(ncol=nTreat,nrow=nrow(RESULTS[[1]]))
        PROBSOUT=as.list(rep(NA,nTreat))
        for(j in 1:length(PROBS)){
          PROBSOUT[[j]]=matrix(ncol=nCat,nrow=nrow(RESULTS[[1]]))
        }





        for(b in 1:nrow(UTMAT)){
          HOLD=GetProbs(nCat,RESULTS[[1]][b,])
          for(j in 1:nTreat){
            UTMAT[b,j]=HOLD[[j]]%*%UT
            PROBSOUT[[j]][b,]=HOLD[[j]]
          }
        }



        UTMAT=UTMAT[(round(nrow(UTMAT)/2):nrow(UTMAT)),]


        ################################
        ###Should we drop a treatment?##
        ################################



        INF=1-KEEP

        ################################
        ###Should we drop a treatment?##
        ################################




        for(j in 1:nTreat){
          if(KEEP[j]==1){
            IND=1:nTreat
            IND=IND[-j]
            KEEP1=KEEP[-j]
            IND = IND[KEEP1==1]
            if(length(IND)>0){
              ##Sequentially check the INF...
              for(m2 in 1:length(DeltaVEC)){
                INF[j]=mean(UTMAT[,j]+DeltaVEC[m2]<UTMAT[,IND[1]])
                if(length(IND)>1){
                  for(m1 in 2:length(IND)){
                    INF[j]=max(INF[j],mean(UTMAT[,j]+DeltaVEC[m2]<UTMAT[,IND[m1]]))
                  }
                }
                if(INF[j]>CUTOFF(DeltaVEC[m2],length(Y),sum(KEEP),6,gamma)){
                  KEEP[j]=0
                  break
                }
              }
            }
          }
        }







        if(sum(KEEP)==1){


          break
        }


        ###Don't allow futility until half patients are enrolled

        if(length(Y)>=(0*max(NLOOK))){
          if(sum(KEEP)>1){
            ##Check if we have a futility stop
            ##Matrix of comparisons...
            COMPKEEP = matrix(rep(1,sum(KEEP)^2),ncol=sum(KEEP),nrow=sum(KEEP))
            for(j in 1:sum(KEEP)){
              for(k1 in 1:sum(KEEP)){
                if(j>k1){
                  COMPKEEP[j,k1]=mean(abs(UTMAT[,which(KEEP==1)[j]]- UTMAT[,which(KEEP==1)[k1]])<DeltaVEC[1])
                }
              }
            }
            if(mean(COMPKEEP>CUTOFF(0,length(Y),sum(KEEP),6,gamma))==1){
              SKIP=1
              SKIPFUT=1
              break
            }


          }

        }




      }
    }

    if(SKIP==1){
      k=0
    }


    DATAHOLD[[rep1]]=cbind(Y,T1)

    SIZE[rep1]=length(Y)


    for(j in 1:nTreat){
      TREATMENTSIZES[rep1,j]=sum(T1==j)
    }


    ###Which arms were dropped?
    WORST[rep1,]=1-KEEP




    ##Correct decision....
    ##What is the set of equally optimal treatments?
    MAXUT = max(UTTRUE)

    KEEP2 = rep(0,length(KEEP))
    KEEP2[ which(abs(UTTRUE-MAXUT)<DeltaVEC[1])]=1
    CD[rep1]=mean(KEEP==KEEP2)




  }



  ###Store results



  TREATPROBS  = matrix(unlist(PROBS),byrow=TRUE,nrow=nTreat)
  rownames(TREATPROBS)=paste("Treatment ", 1:nTreat)
  colnames(TREATPROBS)=paste("Outcome ", 1:nCat)

  SIMTRUTH = as.list(rep(NA,3))
  names(SIMTRUTH)=c("True Treatment Probabilities",
                    "True Covariate - Outcome Relationship",
                    "True Covariate Distributions")
  SIMTRUTH[[1]]=TREATPROBS
  SIMTRUTH[[2]]=Beta
  SIMTRUTH[[3]]=XPROB

  RESULTS = as.list(rep(NA,4))
  names(RESULTS)=c("OCs Summary", "Treatment Specific Summary","Trial Parameters","Simulation Parameters")
  VEC=c(mean(CD),mean(CD==1),mean(SIZE))
  names(VEC)=c("Correct Decision %", "Generalized Power", "Average Sample Size")
  RESULTS[[1]]=VEC

  HOLD = matrix(nrow=3,ncol=nTreat)
  HOLD[1,]= round(UTTRUE,2)
  HOLD[2,]=round(colMeans(WORST),2)
  HOLD[3,]=round(colMeans(TREATMENTSIZES,na.rm=TRUE),2)

  rownames(HOLD)=c("True Utilities: ",
                   "Probability of Stopping: ",
                   "Average Trial Size: ")
  colnames(HOLD)=paste("Treatment ", 1:nTreat)

  RESULTS[[2]]=HOLD




  ###Design parameters...
  DESIGN =as.list(rep(NA,5))
  DESIGN[[1]]=c(nTreat,nCat)
  names(DESIGN[[1]])=c("# Treatments","# Outcomes")
  DESIGN[[2]]=NLOOK
  DESIGN[[3]]=UT
  DESIGN[[4]]=DeltaVEC
  DESIGN[[5]]=gamma
  DESIGN[[6]]=PSPIKE
  names(DESIGN)=c("# Treatments and Outcomes",
                  "Sequential Decision Looks",
                  "Elicited Outcome Utilities",
                  "Vector of Delta Values",
                  "Gamma - Cutoff Vector",
                  "Prior Probability of Pairwise Clustering")

  RESULTS[[3]]=DESIGN


  RESULTS[[4]]=SIMTRUTH


  return(RESULTS)

}
