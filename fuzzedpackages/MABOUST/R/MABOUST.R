#' Conduct the MABOUST Trial design.
#'
#' Performs posterior sampling for the MABOUST design and determines whether the trial should continue and what treatment(s) are optimal.
#' @param Y Ordinal Outcome Vector, labeled 1,...,J
#' @param T1 Treatment Indicator, labeled 1,...,K.
#' @param X Matrix of patient covariates.
#' @param ACTIVE Binary indicator of active treatments. This vector must be length K, and have a 1 for each entry corresponding to an active treatment and 0 otherwise.
#' @param FUTILITY Binary indicator of whether a futility decision will be allowed.
#' @param nTreat Number of treatments in consideration, i.e. K.
#' @param nCat Number of ordinal outcome categories, i.e. J.
#' @param UT Vector of numerical utility scores to give outcomes 1,...,J.
#' @param DeltaVEC Vector of \eqn{\Delta} values to test.
#' @param gamma Length 3 vector of cutoff parameters.
#' @param PSPIKE Prior probability of a pairwise null effect.
#' @param B Number of MCMC iterations to perform.
#' @importFrom graphics boxplot
#' @return The set of active treatments to continue, an optimal treatment, or a set of equally optimal treatments. Also reports posterior mean utilities and ordinal outcome probabilities as well as pairwise comparisons of utility similarity, when appropriate.
#' @references
#' [1] Chapple and Clement (2020), MABOUST: A Multi-Armed Bayesian Ordinal Outcome Utility-Based Sequential Trial. Submitted.
#' @examples
#' ##Clinical Parameters
#' nCat = 6
#' nTreat = 3
#' UT = c(0,10,20,80,90,100)
#' DeltaVEC  = c(5,10)
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
#' X=matrix(rnorm(n*2),ncol=2)
#' ###Number of iterations
#' B=100
#' MABOUST(Y, T1, X, ACTIVE, FUTILITY, nTreat, nCat, UT, DeltaVEC, gamma, PSPIKE,B )
#' @export
MABOUST=function(Y, ##Ordinal Outcome Vector, labeled 1,...,J
                 T1, ##Treatment Indicator, labeled 1,...,K.
                 X, ##Matrix of patient covariates.
                 ACTIVE, ##Vector of Active Treatments.
                 FUTILITY, ##Is a Futility Look allowed here?
                 nTreat, ##Number of treatments in consideration, i.e. K.
                 nCat, ##Number of ordinal outcome categories, i.e. J.
                 UT, ##Vector of numerical utility scores to give outcomes 1,...,J.
                 DeltaVEC, ##Vector of \Delta values to test.
                 gamma, ##Vector of cutoff parameters.
                 PSPIKE, ###Prior probability of a pairwise null effect.
                 B ##Number of MCMC iterations to perform
){


  STOPFUT=0

  ###Check if any Y, or T1 labeled a 0...
  if(min(Y)<=0 | min(T1)<=0){
    warning("Either Y or T1 is mispecified.")
  }else{



    RESULTS=MCMC_MABOUST(Y-1,   ####Outcome vector (minus 1 for c++, i.e. 0,...,J-1)
                         T1-1,  ##Treatment indicator vector
                         X, ##Covariate MAtrix
                         B,
                         nTreat,
                         nCat,
                         PSPIKE) ##Number of iterations for the MCMC







    UTMAT = matrix(ncol=nTreat,nrow=nrow(RESULTS[[1]]))
    PROBSOUT=as.list(rep(NA,nTreat))
    for(j in 1:length(PROBSOUT)){
      PROBSOUT[[j]]=matrix(ncol=nCat,nrow=nrow(RESULTS[[1]]))
    }





    for(b in 1:nrow(UTMAT)){
      HOLD=GetProbs(nCat,RESULTS[[1]][b,])
      for(j in 1:nTreat){
        UTMAT[b,j]=HOLD[[j]]%*%UT
        PROBSOUT[[j]][b,]=HOLD[[j]]
      }
    }


    MEANLIST = lapply(PROBSOUT,function(z){colMeans(z[(round(nrow(z)/2)):nrow(z),])})

    POSTMEAN = round(matrix(unlist(MEANLIST),ncol=nCat,nrow=nTreat,byrow=TRUE),3)

    colnames(POSTMEAN)=paste0("Outcome",1:nCat)
    rownames(POSTMEAN)=paste0("Treatment",1:nTreat)

    OUTLIST = as.list(rep(NA,3))




    UTMAT=UTMAT[(round(nrow(UTMAT)/2):nrow(UTMAT)),]

    ###MAke a boxplot of posterior utilities...
    UTHOLD = NA
    TRTHOLD = NA
    for(j in 1:nTreat){
      UTHOLD= c(UTHOLD,UTMAT[,j])
      TRTHOLD = c(TRTHOLD,rep(j,length(UTMAT[,j])))
    }

    UTHOLD=UTHOLD[-1]
    TRTHOLD=TRTHOLD[-1]

    boxplot(UTHOLD~TRTHOLD,main="Posterior Distributions of Mean Utility",ylab="Mean Utility",xlab="Treatment")


    UTHOLD = colMeans(UTMAT)
    names(UTHOLD)=paste0("Treatment ",1:nTreat)

    OUTLIST[[2]]=UTHOLD
    OUTLIST[[3]]=POSTMEAN

    names(OUTLIST)[1:3]=c("Treatment #s to continue: ","Posterior Mean Utilities","Posterior Mean Outcome Probabilities")




    INF=1-ACTIVE
    KEEP = ACTIVE
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








    ###Don't allow futility until half patients are enrolled


    if(sum(KEEP)>1){
      if(FUTILITY==1){
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
        STOPFUT=1
      }


      if(STOPFUT==1){

        OUTLIST[[1]]=which(KEEP==1)
        names(OUTLIST)[1]="Set of Optimal Treatment #s: "


      }else{
        OUTLIST[[1]]=which(KEEP==1)

      }


      colnames(COMPKEEP)=paste0("Treatment", which(KEEP==1))
      rownames(COMPKEEP)=paste0("Treatment", which(KEEP==1))

      OUTLIST=c(OUTLIST,NA)
      OUTLIST[[4]]=COMPKEEP
      names(OUTLIST)[4]=paste0("Pairwise Posterior Probabilities of Utilities being within ", DeltaVEC[1])



      }else{

        OUTLIST[[1]]=which(KEEP==1)

      }


    }else{

      OUTLIST[[1]]=which(KEEP==1)
      names(OUTLIST)[1]="Optimal Treatment #: "

    }


    return(OUTLIST)

  }

}
