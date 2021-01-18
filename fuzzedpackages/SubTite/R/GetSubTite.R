#' Gives the subgroup specific optimal dose vector.
#'  Returns a list containing the optimal doses to enroll each subgroup at and the subgroups that should have their accrual suspended temporarily.
#' @param Y Vector containing observed event or censoring times.
#' @param I Vector containing event indicators (1 if patient experiences an event for a patient).
#' @param Doses Vector containing numerical doses assigned to patients in the trial.
#' @param Groups Vector containing group assignment of patients, 1 is baseline group.
#' @param Include Binary vector indicating whether each patient record should be included in the decision making process.
#' @param ID Vector of patient IDs. Can be numeric or character valued.
#' @param T1 Reference time for toxicity.
#' @param Target Target cumulative toxicity probability vector at time T1.
#' @param Dose Vector containing the standardized doses considered.
#' @param Upper Cutoff values used to determine if accrual in a subgroup should be suspended.
#' @param cohort Number of patients needed to be assigned at a dose level prior to escalation.
#' @param Conservative Binary Indicator of Whether conservative escalation, i.e. not allowing escalation until cohort patients have been fully evaluated at the highest tried dose level.
#' @param meanmu Prior mean for baseline intercept.
#' @param meanslope Prior mean for baseline slope.
#' @param MeanInts Vector of prior means for the group specific intercept parameters.
#' @param MeanSlopes Vector of prior means for the group specific slope parameters.
#' @param varint Prior variance for the intercept parameters.
#' @param varbeta Prior variance for the slope parameters.
#' @param phetero Prior probability of heterogeneous subgroups.
#' @param Borrow Parameter to specify subgroup borrowing/clustering. 0=No borrowing, 1=Borrowing but no clustering, 2=Borrowing and clustering.
#' @param B Number of Iterations to run for MCMC
#' @return Returns a list with two objects, a vector of optimal doses for each subgroup and matrix of posterior toxicity probabilities at each dose level within each subgroup.
#' @references
#' [1] Chapple and Thall (2017), Subgroup Specific Dose Finding in Phase I Clinical Trials Based on Time to Toxicity Within a Fixed Follow Up Period.
#' @examples
#' T1=28 ##Reference time for toxicity
#' Target=rep(.3,2) ##Target toxicity probability
#' Upper=rep(.95,2) ##Upper cutoffs for excessive toxicity
#' ##How many patients in each subgroup have been assigned at each dose level?
#' cohort=3 ##Cohort size required for escalation
#' Conservative = 1 ##Conservative escalation
#' ##Only can escalate with a fully evaluated cohort at the highest dose level.
#' ##Matrix of umber of patients tried or fully evaluated at each dose level.
#' ##Hyperparameters
#' meanmu=-0.4467184 ##Common Intercept hypermean
#' meanslope= 0.8861634 ##Common slope hypermean
#' MeanInts = -0.5205379 ##Group Intercept hypermeans
#' MeanSlopes = 0.1888923 ##Group slope hyperneabs
#' varint=5 #Prior Variance of the intercept betas
#' varbeta=1 ##Prior Variance of slope betas
#' phetero=.9 ##Prior Probability of hetergeneity
#' Borrow=0 ##Borrowing specification, 0=none, 1=some, 2=clustering.
#' B=5000 ##Number of iterations
#' Borrow=2
#' Y=c(28,26,29,28,29,5,1)
#' RawDose=c(350,420,530,660,825)
#' Dose=(RawDose-mean(RawDose))/sd(RawDose)
#' I <- c(0,0,0,0,0,0,0)
#' Doses <- rep(2,7)
#' Groups <- c(0,1,1,0,0,1,1)
#' Include <- rep(1,7)
#' ID=1:length(Y)
#' Z=GetSubTite(Y, I,Doses, Groups, Include,ID,cohort, Conservative,
#' T1,Target,  Upper, Dose,  meanmu, meanslope,
#'  MeanInts,  MeanSlopes ,varint,varbeta,phetero, Borrow,B)
#' Z
#'@export
GetSubTite=function(Y, I,Doses, Groups,  Include, ID, cohort,Conservative,T1, Target,
                    Upper, Dose,  meanmu, meanslope,
                    MeanInts,  MeanSlopes ,varint,varbeta,phetero,Borrow,B){



  Doses2=Doses

  if(sum(Doses %in% Dose)>0){
    warning("Doses assigned to patients should be numbered.")
  }


  ###



  DATAHOLD = data.frame(cbind(ID,round(Y,3),I,Doses,Groups,Include))
  colnames(DATAHOLD)=c("Patient ID","Toxicity/Censoring Time","Toxicity Indicator","Assigned Dose #","Subgroup","Included?")



  ###Re-write Y so that it is within the appropriate follow up window.
  Y[Y>T1]=T1  ##




  ERRHOLD=c(length(Target), length(Upper), length(MeanInts)+1, length(MeanSlopes)+1)

  HOLD=0
  ##Check for errors in dimension specification
  for(k in 1:length(ERRHOLD)){
    for(m in 1:length(ERRHOLD)){
      if(ERRHOLD[k] != ERRHOLD[m]){
        HOLD=1
      }
    }
  }

  if(HOLD==1){
    message("Target toxicity vector, toxicity threshold, or subgroup hyperparameter vector has incorrect dimensions")
  }else{


    ###Contains Design parameters
    DESIGN = as.list(rep(NA,14))

    names(DESIGN) = c("Standardized dose levels:",
                      "Target toxicity probabilities:",
                      "Posterior probability thresholds for overly toxic subgroups:",
                      "Escalation scheme",
                      "Prior mean for the baseline intercept = ",
                      "Prior mean for the baseline slope = ",
                      "Prior means for other subgroup intercepts = ",
                      "Prior means for other subgroup slopes = ",
                      "Prior intercept variance = ",
                      "Prior slope variance = ",
                      "Borrow/Clustering Settings",
                      "Borrow Indicator",
                      "DLT observation time: ",
                      "Number of MCMC iterations = ")

    DESIGN[[5]]=meanmu
    DESIGN[[6]]=meanslope
    DESIGN[[7]]=MeanInts
    DESIGN[[8]]=MeanSlopes
    DESIGN[[9]]=varint
    DESIGN[[10]]=varbeta
    DESIGN[[12]]=Borrow
    DESIGN[[13]]=T1
    DESIGN[[14]]=B

    DESIGN[[1]]=Dose

    names(Target)=paste0("Subgroup ",1:length(Target))
    DESIGN[[2]]=Target

    names(Upper)=paste0("Subgroup ",1:length(Target))
    DESIGN[[3]]=Upper






    if(Conservative==1){
      DESIGN[[4]]=paste0("Conservative escalation requiring ", cohort," patients to be fully evaluated.")

    }else{
      DESIGN[[4]]=paste0("Aggressive escalation requiring ", cohort," patients to be treated, but not fully evaluated.")
    }

    if(Borrow==2){
      DESIGN[[11]]=paste0("Clustering and borrowing with a prior probability of heterogeneity of: ",phetero)
    }


    if(Borrow==1){
      DESIGN[[11]]=paste0("Borrowing between the subgroups but no clustering")
    }


    if(Borrow==0){
      DESIGN[[11]]=paste0("No borrowing or clustering between subgroups")
    }


    TIME = paste0("Model run on: ",Sys.time())

    DESIGN=c(TIME,DESIGN)
    names(DESIGN)[[1]]="Date/Time of escalation decision"

    DESIGN = c("Sub-TITE Package Version: 4.0.0",DESIGN)




    ###Check if the Groups are labeled incorrectly
    ##Groups should be labeled 0, ..., G-1
    ##If they are not labeled right, re-write them
    if(min(Groups)!=0 | max(Groups) != (length(Target)-1)){

      warning("Subgroup vector is not labeled correctly! Please re-format from 0,...,G-1.")


    }else{








      ###Let's write our dose-tried matrix...
      ###The vector Doses will now have the numeric values
      DoseTried = matrix(nrow=length(Upper),ncol=length(Dose))
      NumTreated = DoseTried
      ##Let's fill this in manually....
      for(j in 1:length(Dose)){
        for(k in 1:length(Upper)){
          NumTreated[k,j]=sum(Groups==(k-1) & Doses==j)
        }
      }


      ###Holder for number of patients tried
      NUMTREAT=NumTreated

      NUMTOX = NUMTREAT

      for(j in 1:length(Dose)){
        for(k in 1:length(Upper)){
          NUMTOX[k,j]=sum(Groups==(k-1) & I==1 & Doses==j)
        }
      }




      ###re-write lower doses....
      for(k in 1:nrow(NumTreated)){
        for(j in 1:(ncol(NumTreated)-1)){
          if(NumTreated[k,j]==0){
            ###Have any doses above this been used?
            if(sum(NumTreated[k,(j+1):ncol(NumTreated)])>0){
              ##Some doses above this have been tried already
              NumTreated[k,j]=cohort
            }

          }
        }
      }



      ###Output the number of patients treated and dlt at each dose
      HOLD = NUMTREAT

      for(k in 1:nrow(NUMTREAT)){
        for(j in 1:ncol(NUMTREAT)){
          HOLD[k,j]=paste0(NUMTOX[k,j],"/",NUMTREAT[k,j])
        }
      }


      colnames(HOLD)=1:ncol(NUMTREAT)
      rownames(HOLD)=paste0("Subgroup ",1:nrow(HOLD))

      HOLDTHIS = noquote(HOLD)


      ###Do we have conservative escalation?
      if(Conservative==1){
        ###We have conservative escalation, so let's re-write the NumTreated (DoesTried) Matrix
        ###So that we cannot escalate unless we have fully evaluated the largest dose level
        for(k in 1:nrow(NumTreated)){
          ###What is the highest TRIED dose level?
          which1=max(which(NumTreated[k,]>0))

          ##Are we at the highest dose level?
          if(which1<ncol(NumTreated)){
            ###We can still escalate
            ###Let's check to make sure that this dose level has the right number of fully evaluated patients
            NUM  = sum( (Y[Doses==which1 & Groups==(k-1)]>=T1 & I[Doses==which1 & Groups==(k-1)]==0) | I[Doses==which1 & Groups==(k-1)]==1)
            ###Now re-write this with NUM
            NumTreated[k,which1]=NUM


          }


        }

      }





      DoseTried=NumTreated




      ##Repackage MeanInts and MeanSlopes
      MeanInts=c(0,MeanInts)
      MeanSlopes=c(0,MeanSlopes)
      Stopped=rep(0,nrow(DoseTried))
      ##This matrix contains posterior mean toxicity probabilities at each dose for each subgroup.
      ##The last column in the matrix has whether or not each group should be stopped.
      RESULTS=MCMC( Y,I,  Dose[Doses],  Groups,  T1,  Target,  Upper, Dose, meanmu,  meanslope,
                    MeanInts,  MeanSlopes, varint,  varbeta, phetero, Stopped,  length(Y),  Borrow,B)

      CLUST=RESULTS[[2]]


      RESULTS1=RESULTS ##HOLDER FOR STOPPED GROUPS
      RESULTS=RESULTS[[1]]
      Stopped= RESULTS[,ncol(DoseTried)+1]
      ##Get optimal Dose
      OptDose= Stopped
      ##Check if ALL groups are stopped, if so run a separate trial in each.
      if(Borrow>0){
        if(sum(Stopped)==length(Stopped)){
          message("Borrowing has caused all groups to stop due to excessive toxicity.
Separate models will be fit to each group to ensure the design is not stopping too early due to borrowing.")
          RESULTS=MCMC( Y,I,  Dose[Doses],  Groups,  T1,  Target,  Upper, Dose, meanmu,  meanslope,
                        MeanInts,  MeanSlopes, varint,  varbeta, phetero,    Stopped=rep(0,nrow(DoseTried)),  length(Y),  0,B)
          RESULTS1=RESULTS ##HOLDER FOR STOPPED GROUPS
          RESULTS=RESULTS[[1]]
        }
      }



      PROBS = data.frame(RESULTS[,1:ncol(DoseTried)])  ##Used for looping over probs
      ##Can only choose a dose level if we've escalated correctly
      PROBS1=PROBS
      Y1=DoseTried<cohort ##Flag which dose levels haven't been treated enough.
      ##Cant escalate past that dose

      for(k in 1:length(Stopped)){
        j=1
        if(sum(1-Y1[k,])<ncol(DoseTried)){
          ##Checks if all doses meet cohort criteria
          while(Y1[k,j]==FALSE){
            j=j+1
          }
          ##We can go up one more
          j=j+1
          ##Are we at the highest dose level? If not, we need to check if we can escalate.
          if(j<ncol(DoseTried)){
            ##Reset PROB1 with -1000, so now we can't pick it
            PROBS1[k,j:ncol(DoseTried)]=-10000
          }
        }
      }




      ##Now get optimal doses
      for(k in 1:length(Stopped)){
        if(Stopped[k]==0){
          a1 = abs(PROBS1[k,]-Target[k])
          OptDose[k]=which(a1==min(a1)) ##Minimum distance from Target probability
        }else{
          OptDose[k]=NA ##Na for stopped groups
        }
      }





      Z=as.list(c(0,0,0,0,0,0))

      cat("
Decisions for next patient or optimal dose.

")

      for(k in 1:length(Stopped)){
        if(!is.na(OptDose[k])){
          cat(paste0("Next reccomended dose level for subgroup ",k,": Dose ",OptDose[k],"
"))
        }else{
          cat(paste0("Subgroup ",k," is too toxic, do not enroll these patients at this time.
"))
        }
      }


      Z[[1]]=OptDose
      for(k in 1:nrow(PROBS)){
        rownames(PROBS)[k]=paste0("Subgroup ", k)
      }
      ##Columns by dose levels
      colnames(PROBS)=1:ncol(DoseTried)

      Z[[2]]=PROBS
      ##Posterior probability of stopping
      Z[[3]]=RESULTS[,ncol(RESULTS)]


      Z[[5]]=HOLDTHIS

      names(Z)=c("Optimal Dose","Posterior Mean Toxicity Probability",
                 "Posterior Probability of Overly Toxic Subgroup",
                 "Clustering Parameters", "Number of Treated and DLTs",
                 "Data")

      if(Borrow==2){

        ###Cluster membership...
        CLUST=CLUST+1
        ##Cluster assignment
        CLUST1 = as.data.frame(matrix(nrow=(nrow(PROBS)),ncol=nrow(PROBS)))


        for(k in 1:nrow(PROBS)){
          for(j in 1:nrow(PROBS)){
            CLUST1[k,j]=mean(CLUST[,k]==j)
          }
        }

        for(k in 1:nrow(PROBS)){
          rownames(CLUST1)[k]=paste0("Subgroup ", k)
          colnames(CLUST1)[k]=paste0("Latent Subgroup",k)

        }

        NCLUST=rep(0, nrow(CLUST))

        for(j in 1:nrow(CLUST)){
          NCLUST[j] = length(unique(CLUST[j,]))
        }


        ## print(RESULTS1[[2]])

        Z2=as.list(c(0,0,0))
        Z2[[1]]=CLUST1
        Z2[[2]]=table(NCLUST)/nrow(CLUST)


        ##Finally, let's get the clustering breakdown
        ##For G groups, there are (G choose 2) + (G choose 3) + .... (G choose G-1) + 2 cluster configurations

        G=nrow(DoseTried)

        if(G==2){
          HOLD = c(0,0)
          HOLD[1]=mean(CLUST[,1]==CLUST[,2])
          HOLD[2]=mean(CLUST[,1]!=CLUST[,2])

          names(HOLD)=c("1-2","1,2")
        }else{

          if(G<=4){
            NClust=2
            for(j in 2:(G-1)){
              NClust = NClust + choose(G,j)
            }



            ##Make Matrix
            HOLD = rep(NA,NClust)
            ##First do the pairs
            if(G==3){
              ##5 clusters
              ##1-2,3
              ##2-3, 1
              ##1-3, 2
              ##1,2,3
              ##1-2-3

              HOLD[1]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,3]!=CLUST[,1]))
              HOLD[2]=mean((CLUST[,2]==CLUST[,3])*(CLUST[,2]!=CLUST[,1]))
              HOLD[3]=mean((CLUST[,1]==CLUST[,3])*(CLUST[,2]!=CLUST[,1]))
              HOLD[4]=mean(NCLUST==G)
              HOLD[5]=mean(NCLUST==1)
              names(HOLD)=c("1-2,3", "2-3,1","1-3,2","1,2,3","1-2-3")

            }



            if(G==4){
              HOLD=rep(NA,15)
              names(HOLD)=c("1-2,3-4", "1-3,2-4","1-4,2-3",
                            "1-2,3,4", "1-3,2,4", "1-4,2,3",
                            "2-3,1,4","2-4,1,3","3-4,1,3",
                            "1-2-3,4", "1-2-4,3","1-3-4,2",
                            "2-3-4,1","1,2,3,4","1-2-3-4")
              ##15 clusters
              ##1-2,3-4
              ##1-3, 2-4
              ##1-4, 2-3 *
              HOLD[1]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,3]==CLUST[,4])*(CLUST[,3] != CLUST[,1]))
              HOLD[2]=mean((CLUST[,1]==CLUST[,3])*(CLUST[,2]==CLUST[,4])*(CLUST[,4] != CLUST[,1]))
              HOLD[3]=mean((CLUST[,1]==CLUST[,4])*(CLUST[,2]==CLUST[,3])*(CLUST[,3] != CLUST[,1]))
              ##1-2, 3,4
              #1-3, 2,4
              #1-4, 2,3
              #2-3,1,4
              #2-4,1,3
              #3-4, 1, 2
              HOLD[4]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,3]!=CLUST[,1])*(CLUST[,3] != CLUST[,4])*(CLUST[,4] != CLUST[,1]))
              HOLD[5]=mean((CLUST[,1]==CLUST[,3])*(CLUST[,2]!=CLUST[,4])*(CLUST[,4] != CLUST[,1])*(CLUST[,4] != CLUST[,2]))
              HOLD[6]=mean((CLUST[,1]==CLUST[,4])*(CLUST[,2]!=CLUST[,3])*(CLUST[,3] != CLUST[,1])*(CLUST[,2] != CLUST[,1]))
              HOLD[7]=mean((CLUST[,2]==CLUST[,3])*(CLUST[,1]!=CLUST[,2])*(CLUST[,1] != CLUST[,4])*(CLUST[,4] != CLUST[,2]))
              HOLD[8]=mean((CLUST[,2]==CLUST[,4])*(CLUST[,1]!=CLUST[,4])*(CLUST[,1] != CLUST[,2])*(CLUST[,3] != CLUST[,2]))
              HOLD[9]=mean((CLUST[,2]==CLUST[,3])*(CLUST[,1]!=CLUST[,2])*(CLUST[,4] != CLUST[,2])*(CLUST[,1] != CLUST[,4]))
              ##1-2-3,4
              ##1-2-4,3
              ##1-3-4,2*
              ##2-3-4,1
              ##1,2,3,4
              ##1-2-3-4

              ##3 pairs

              ##4 3 ways
              HOLD[10]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,1]==CLUST[,3])*(CLUST[,1] != CLUST[,4]))
              HOLD[11]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,1]==CLUST[,4])*(CLUST[,1] != CLUST[,3]))
              HOLD[12]=mean((CLUST[,1]==CLUST[,3])*(CLUST[,1]==CLUST[,4])*(CLUST[,1] != CLUST[,2]))
              HOLD[13]=mean((CLUST[,2]==CLUST[,3])*(CLUST[,3]==CLUST[,4])*(CLUST[,1] != CLUST[,4]))
              HOLD[14]=mean(NCLUST==G)
              HOLD[15]=mean(NCLUST==1)

            }




          }else{
            HOLD = "Cannot Return Unique Clusters for more than 4 subgroups"
          }
        }






        Z2[[3]]=HOLD
        names(Z2)= c("Latent Subgroup Posterior","Posterior # of Clusters","Posterior Probability Cluster Configuration")
        Z[[4]]=Z2

      }else{
        Z[[4]]="No Clustering"
      }


      ###Write the dataframe into the last item of the list
      Z[[length(Z)]]=DATAHOLD
      Z=c(1,Z)
      Z[[1]]=DESIGN
      names(Z)[[1]]="Design Parameters"


      ###Add last... # completely evaluated at each dose...
      Z=c(Z,1)
      names(Z)[[length(Z)]]="Number Fully Evaluated At Each Dose"
      MAT = matrix(nrow=length(Target),ncol=length(Dose))
      rownames(MAT)=paste0("Subgroup ",1:length(Target))
      colnames(MAT)=1:length(Dose)


      ##Fill in....
      for(j in 1:nrow(MAT)){
        for(k in 1:ncol(MAT)){
          MAT[j,k]=sum(I[Groups==(j-1) & Doses==k]) + sum((1-I[Groups==(j-1) & Doses==k]) & Y[Groups==(j-1) & Doses==k]>=T1)
        }
      }


      Z[[length(Z)]]=MAT


      return(Z)

    }
  }

}
