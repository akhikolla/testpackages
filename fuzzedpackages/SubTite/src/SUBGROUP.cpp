#include "RcppArmadillo.h"
#include <time.h>
#include <Rmath.h>
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#include <numeric>
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <iostream>     // std::cout
#include <cmath>
#include <cfloat>
#include <stdio.h>
#include <float.h>
#define PI 3.14159265


// [[Rcpp::depends(RcppArmadillo)]]



using namespace Rcpp;


int Sample2(arma::vec groupprob){
  arma::vec cumprob=groupprob;
  int m=0;

  for(m=1;m<groupprob.n_rows;m++){
    cumprob[m]=cumprob[m]+cumprob[m-1];
  }
  //Now we have the vector of cumulative probabilities, let's draw a random unif
  double U=as_scalar(arma::randu(1));

  int Which=0;


  if(U<cumprob[0]){
    Which=0;
  }else{

    for(m=0;m<(groupprob.n_rows-2);m++){
      if( (U>cumprob[m]) && (U<cumprob[m+1]) ){
        Which=m+1;
      }

    }

    if(U>cumprob[groupprob.n_rows-2]){
      Which=groupprob.n_rows -1;
    }


  }


  return(Which);

}


int GetFullyFollowed(arma::vec Y, //Survival Times
                     arma::vec I, //Censoring Indicators
                     arma::vec Doses, //Doses given to patients
                     arma::vec Dose, //Standardized doses in consideration in trial
                     double T1){
  int sum1=0;

  int m=0;

  for(m=0;m<Y.n_rows;m++){
    if(Doses[m]==Dose[0]){
      if(I[m]==1){
        sum1++;
      }

      if(Y[m]==T1){
        sum1++;
      }




    }


  }





  return(sum1);


}



double MaxVec(arma::vec Y){
  int J1=Y.n_rows;
  int j=0;
  double max=Y[0];
  for(j=1;j<J1;j++){
    if(Y[j]>max){
      max=Y[j];
    }
  }

  return(max);

}




double MinVec(arma::vec Y){
  int J1=Y.n_rows;
  int j=0;
  double max=Y[0];
  for(j=1;j<J1;j++){
    if(Y[j]<max){
      max=Y[j];
    }
  }

  return(max);

}


double min1(double a, double b){
  double z=0;
  if(a>=b){
    z=b;
  }else{
    z=a;
  }

  return(z);
}



double max1(double a, double b){
  double z=0;
  if(a>=b){
    z=a;
  }else{
    z=b;
  }

  return(z);
}





///Likelihood for the cluster enabled version, we can also use this for the NON cluster enabled version
//[[Rcpp::export]]
double LikeStoppedCluster(arma::vec Y,  //Vector of Toxicity Times
                          arma::vec I,  // Vector of Censoring Indi
                          arma::vec Dose, //Vector of numerical Doses given to Patients
                          arma::vec Group, //Vector of Group Membership, 0 is baseline
                          double mu, // Shared intercept
                          double slope, //Shared slope
                          arma::vec a, //Vector of intercepts for  groups
                          arma::vec b, //vector of slopes for  groups
                          double T1, // Reference Time
                          int nPats, //Number of patients currently
                          arma::vec GroupMem1, //Current clustering membership
                          arma::vec Stopped //Stopping vector
){


  double LogL = 0; // Will store likelihood to return

  //Baseline vector
  arma::vec eta(nPats);

  //Needed for For loops
  int m=0;
  int k=1;
  int j=0;

  for(m=0;m<eta.n_rows;m++){
    for(k=0;k<a.n_rows;k++){
      if(GroupMem1(Group[m])==k){

        eta(m) = a(k)+exp(b(k)+slope)*Dose(m)+mu;

      }
    }
  }


  //Now we have the linear predictors in the vector \eta for all the patients
  //Now let's compute the likelihood.
  arma::vec Y1=Y/T1; //Fraction of followup time needed to get the likelihood
  arma::vec eta1=exp(eta);

  for(m=0;m<eta.n_rows;m++){
    //Check if this group is stopped. If it is, we don't use it in the likelihood AT ALL
    if(Stopped(GroupMem1(Group(m)))==0){
      if(I(m)==1){
        //Not censored so we use the PDF
        LogL = LogL + eta[m] - log(1+eta1[m]);

      }else{
        LogL = LogL + log(1+(1-Y1[m])*eta1[m])-log(1+eta1[m]);


      }
    }
  }






  return(LogL);


}


//[[Rcpp::export]]
double LikeStoppedSeparate(arma::vec Y,  //Vector of Toxicity Times
                           arma::vec I,  // Vector of Censoring Indi
                           arma::vec Dose, //Vector of numerical Doses given to Patients
                           arma::vec Group, //Vector of Group Membership, 0 is baseline
                           arma::vec a, //Vector of intercepts for  groups
                           arma::vec b, //vector of slopes for  groups
                           double T1, // Reference Time
                           int nPats, //Number of patients currently
                           arma::vec Stopped //Stopping vector
){


  double LogL = 0; // Will store likelihood to return
  //Baseline vector
  arma::vec eta(nPats);

  //Needed for For loops
  int m=0;
  int k=1;
  int j=0;


  for(m=0;m<eta.n_rows;m++){
    for(k=0;k<a.n_rows;k++){
      if(Group[m]==k){
        //Now we have one total parameter that we use for borrowing.
        eta(m) = a(k)+exp(b(k))*Dose(m);
      }
    }
  }



  //Now we have the linear predictors in the vector \eta for all the patients
  //Now let's compute the likelihood.
  arma::vec Y1=Y/T1; //Fraction of followup time needed to get the likelihood
  arma::vec eta1=exp(eta);

  for(m=0;m<eta.n_rows;m++){
    if(I[m]==1){
      //Not censored so we use the PDF
      LogL = LogL + eta[m] - log(1+exp(eta[m]));

    }else{
      LogL = LogL + log(1+(1-Y1[m])*exp(eta[m]))-log(1+exp(eta[m]));

    }

  }







  return(LogL);


}

















//Clusters on a random group. TEST THIS
//[[Rcpp::export]]
int GetRandGroup(arma::vec INC,arma::vec Stopped){ //Vector containing stopped indicators)


  //How many are UNCLUSTERED AND not stopped
  double NumGroup = sum(INC.t()*(1-Stopped));
  //Get probabilities for uniform assignment
  double prob = 1/NumGroup; //Can assign with equal probability
  arma::vec z9(2);

  //Random uniform
  double U=as_scalar(arma::randu(1));




  //While loop until we get to prob
  int m=0;
  double prob1=prob;

  if(U>prob){
    while(U>prob1){
      m++; //Increment m
      prob1=(m+1)*prob;
    }
  }



  //m holds the number of the stopped and correctly binned group.

  //Let's find where this is...
  int k=-1;
  int  sum1=-1;
  while(sum1<m){
    k++;
    sum1 = sum1+INC(k)*(1-Stopped(k));
  }





  return(k);

}





//' Performs MCMC and returns needed values for dose-finding in a list.
//' @param Y Vector containing observed event or censoring times.
//' @param I Vector containing event indicators (1 if patient experiences an event for a patient).
//' @param Doses Vector containing Doses of patients in trial.
//' @param Groups Vector containing group assignment of patients, 0 is baseline group.
//' @param T1 Reference time for toxicity.
//' @param Target Target cumulative toxicity probability vector at time T1.
//' @param Dose Vector containing the standardized doses considered.
//' @param Upper Cutoff values used to determine if accrual in a subgroup should be suspended.
//' @param meanmu Prior mean for baseline intercept.
//' @param meanslope Prior mean for baseline slope.
//' @param MeanInts Vector of prior means for the group specific intercept parameters.
//' @param MeanSlopes Vector of prior means for the group specific slope parameters.
//' @param varint Prior variance for the intercept parameters.
//' @param varbeta Prior variance for the slope parameters.
//' @param phetero Prior probability of heterogeneous subgroups.
//' @param SubRout Parameter to specify subgroup borrowing/clustering. 0=No borrowing, 1=Borrowing but no clustering, 2=Borrowing and clustering.
//' @param B Number of Iterations to run for MCMC
//' @param Stopped Current vector of STOPPED groups
//' @param NumPat Number of patients
//'@importFrom Rcpp evalCpp
//'@useDynLib SubTite
//'@return A list of quantities needed for determining the next dose to enroll each subgroup.
//'@export
//[[Rcpp::export]]
List MCMC( arma::vec Y, //Vector of Times
           arma::vec I, //Vector of Toxicity Indicators
           arma::vec Doses, //Standardized Dose Values given to each patient
           arma::vec Groups, //SubGroup Membership
           double T1, //Reference Time
           arma::vec Target, //Target Toxicity Probability for each group
           arma::vec Upper, //Vector of Upper Cutoffs for each groups
           arma::vec Dose, //Standardized values of Doses considered
           double meanmu, //This is the prior mean of the baseline int
           double meanslope, //Prior mean of baseline slope
           arma::vec MeanInts, //Prior Means of Intercept Parameters for all G groups
           arma::vec MeanSlopes, //Prior Means of Slope Parameters for all G groups
           double varint, //Prior Variance of Intercept Parameters
           double varbeta, //Prior Variance of Slope Parameters
           double phetero, //Prior probability of heterogeneity
           arma::vec Stopped, //Current vector of STOPPED groups
           int NumPat, //Number of patients
           int SubRout, // Contains a 0 for simple model with no borrowing, 1 for borrowing model, 2 for borrow-clustering model.
           double B //Number of iterations for MCMC
){


  int g=0;

  int StoreInx=0;

  //Change this to be an input
  arma::vec PROBIN(MeanInts.n_rows);
  for(g=0;g<PROBIN.n_rows;g++){
    //Prior probability of NO clustering
    PROBIN(g)=phetero;
  }

  int MEANS=0;

  int IntIN=0;
  int IntOUT=0;

  //Group Slope Prior Var
  double varbeta1=varbeta;
  //Group Int Prior Var
  double varint1=varint;



  int Group=0;
  //Important For loop integer
  int m =0; //For the inner MCMC looprep
  int i=0; //For the patient index
  int rep=0; //For the simulation repetition.


  //First Fill out Target Vector
  int nDose=Dose.n_rows;



  double B1=B/2;

  //Make List objects we'll use





  //Innitialize parameters for MCMC
  double mu=meanmu;
  double slope=meanslope;
  double sig=T1;

  double NewMean=0;
  double NewSlope=0;
  int Which1=0;
  int J=Dose.n_rows; //Number of Doses


  //Important quantities we will need in the MCMC
  int  G=MeanInts.n_rows; //Number of groups
  double alpha=0;
  double U=0;
  double signew=0;
  double slopenew=0;
  double Munew=0;

  //For loop integers
  int j=0;
  int k=0;


  //Vector Of group intercepts
  arma::vec a(G);
  a.zeros();
  //Vector of group slopes
  arma::vec b=a;
  //Vectors for Group MCMC



  //Vectors for proposals
  arma::vec aprop = a;
  arma::vec bprop=b;





  arma::vec NumA(G);
  arma::vec NumB(G);
  arma::vec IntA(G);
  arma::vec IntB(G);
  IntB.zeros();
  IntA.zeros();
  NumA.zeros();
  NumB.zeros();

  double NumMu = 2;
  double NumSlope=2;
  double IntMu=1;
  double IntSlope=1;

  arma::vec avar(G);
  arma::vec bvar(G);
  avar.zeros();
  avar=avar+1;
  bvar.zeros();
  bvar=bvar+1;


  double muvar=1;  ///Proposal variance for shared intercept
  double slopevar=1; //Proposal variance for shared slope
  double trialtime=0;
  NumericVector z9(2);
  NumericVector zprop(5);

  double NumSig=2;
  double IntSig=1;
  double sigvar=1;

  arma::vec etaG(G);
  arma::mat DoseProb(B1,G*nDose);
  arma::mat MeanDose(G,nDose);
  arma::vec sigstore(B1); //What is sigstore??
  arma::vec mustore(B1); //Storage for shared intercept
  arma::vec slopestore(B1); //Storage for shared slope
  arma::mat astore(B1,G); //Storage for group specific intercepts
  arma::mat bstore=astore; //Storage for group specific slopes
  double eta1=0;
  double DEC=0;

  arma::vec GroupMem(G);
  GroupMem.zeros();
  arma::vec GroupMemProp=GroupMem;

  arma::vec INVEC(G);
  INVEC.zeros();

  double eta2=0;

  arma::mat MeanVec(G,nDose);
  arma::vec NTox(G);
  arma::vec OptDose(G);
  MeanVec.zeros();
  NTox.zeros();
  OptDose.zeros();
  arma::mat CLUST(B1,G);
  CLUST.zeros();
  arma::mat INSTORE=CLUST;







  arma::vec INVECNEW=INVEC;
  arma::vec Y2;
  arma::vec I2;
  arma::vec Groups2;
  arma::vec Doses2;

  Y2.zeros();
  I2.zeros();
  Groups2.zeros();
  Doses2.zeros();


  //Start all parameters at their prior distributions
  mu=meanmu;
  slope=meanslope;
  INVEC.zeros();
  INVECNEW.zeros();
  INVEC = INVEC+1;
  a=MeanInts;
  b=MeanSlopes;

  //Start counter
  for(m=0;m<GroupMem.n_rows;m++){
    GroupMem[m]=m;
  }
  //Proposal for combinations
  GroupMemProp=GroupMem;


  //No Borrowing, No Clustering

  if(SubRout==0){
    //Reset our vectors
    MeanInts = MeanInts + mu;
    MeanSlopes = MeanInts + slope;
    a=MeanInts;
    b=MeanSlopes;

    for(m=0;m<B;m++){

      if(m<(B/2 + 2)){
        if(m%100==0){



          for(k=0;k<G;k++){
            if((IntA[k]/NumA[k])>.5){
              avar[k]=avar[k]*2;
            }

            if((IntA[k]/NumA[k])<.2){
              avar[k]=avar[k]/2;
            }








            if((IntB[k]/NumB[k])>.5){
              bvar[k]=bvar[k]*2;
            }

            if((IntB[k]/NumB[k])<.2){
              bvar[k]=bvar[k]/2;
            }



          }




          //Reset our counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB=IntB.zeros()+1;






        }
      }



      for(k=0;k<aprop.n_rows;k++){
        //We must have the right group AND it must NOT be stopped.
        if(Stopped(k)==0){
          //Propose new value
          aprop=a;
          aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);
          //Get our Alphas
          //Prior Ratio
          alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 ;
          //Loglikelihood ratio
          alpha=alpha+   LikeStoppedSeparate(Y, I,  Doses, Groups,  aprop, b, sig,NumPat,Stopped) -  LikeStoppedSeparate(Y, I,  Doses, Groups,  a, b, sig,NumPat, Stopped);
          //Metropolis Hastings
          U=log(as_scalar(arma::randu(1)));
          z9(0)=LikeStoppedSeparate(Y, I,  Doses, Groups,  a, b, sig,NumPat, Stopped);
          z9(1)=LikeStoppedSeparate(Y, I,  Doses, Groups,  aprop, b, sig,NumPat,Stopped) ;


          //          Rf_PrintValue(wrap(z9));

          if(U<alpha){
            a(k)=aprop(k);
            IntA[k]=IntA[k]+1;
          }

          NumA[k]=NumA[k]+1;


        }
      }
      //Now do a MH step for the slope

      for(k=0;k<aprop.size();k++){

        if(Stopped(k)==0){
          bprop=b;
          bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];
          //Prior Ratio
          alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1;
          //Likelihood ratio
          alpha=alpha+   LikeStoppedSeparate(Y, I,  Doses, Groups,   a, bprop, sig,NumPat,Stopped) -  LikeStoppedSeparate(Y, I,  Doses, Groups,  a, b, sig,NumPat, Stopped);
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            b[k]=bprop[k];
            IntB[k]=IntB[k]+1;
          }
          NumB[k]=NumB[k]+1;

        }

      }







      if(m>(B-B1-1)){
        //Make OptDose and StoppedGroups
        StoreInx=m-B1;

        for(j=0;j<G;j++){;
          astore(StoreInx,j)=a(j);
          bstore(StoreInx,j)=b(j);
        }



        //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
        //First nDose



        for(k=0;k<G;k++){

          for(j=0;j<nDose;j++){
            eta1=0;

            eta1 = a(k)+exp(b(k))*Dose[j];

            z9(0)=eta1;


            //For all groups add the intercept and slope*Groups

            eta1=exp(eta1);




            if(eta1>100000){
              DoseProb(StoreInx, k*nDose+j) = 1;

            }else{
              DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

            }



          }


        }
        //Now we have a matrix containing our sampled dose probabilities









      }





      //End of MCMC

    }

  }


  //Borrow strength, No cluster

  if(SubRout==1){


    for(m=0;m<B;m++){

      if(m<(B/2 + 2)){
        if(m%100==0){


          for(k=0;k<G;k++){
            if((IntA[k]/NumA[k])>.5){
              avar[k]=avar[k]*2;
            }

            if((IntA[k]/NumA[k])<.2){
              avar[k]=avar[k]/2;
            }








            if((IntB[k]/NumB[k])>.5){
              bvar[k]=bvar[k]*2;
            }

            if((IntB[k]/NumB[k])<.2){
              bvar[k]=bvar[k]/2;
            }



          }


          if((IntMu/NumMu)>.5){
            muvar=muvar*2;
          }

          if((IntMu/NumMu)<.2){
            muvar=muvar/2;
          }


          if((IntSlope/NumSlope)>.5){
            slopevar=slopevar*2;
          }

          if((IntSlope/NumSlope)<.2){
            slopevar=slopevar/2;
          }




          //Reset our counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB=IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          IntSig=1;
          NumSig=2;




        }
      }



      for(k=0;k<aprop.n_rows;k++){
        //We must have the right group AND it must NOT be stopped.
        if((INVEC(k)>0) && (Stopped(k)==0)){
          //Propose new value
          aprop=a;
          aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);
          //Get our Alphas
          //Prior Ratio
          alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 ;
          //Loglikelihood ratio
          alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);
          //Metropolis Hastings
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            a(k)=aprop(k);
            IntA[k]=IntA[k]+1;
          }

          NumA[k]=NumA[k]+1;


        }
      }
      //Now do a MH step for the slope

      for(k=0;k<aprop.size();k++){

        if((INVEC(k)>0) && (Stopped(k)==0)){
          bprop=b;
          bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];
          //Prior Ratio
          alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1;
          //Likelihood ratio
          alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            b[k]=bprop[k];
            IntB[k]=IntB[k]+1;
          }
          NumB[k]=NumB[k]+1;

        }

      }









      Munew = mu + as_scalar(arma::randn(1))*muvar;



      alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint;

      alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);


      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        mu=Munew;
        IntMu=IntMu+1;
      }

      NumMu=NumMu+1;




      //Sample the new slope and the new \beta_g coefficients jointly




      slopenew = slope + as_scalar(arma::randn(1))*slopevar;




      alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


      alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);


      U=log(as_scalar(arma::randu(1)));




      if(U<alpha){
        slope=slopenew;
        IntSlope=IntSlope+1;
      }
      NumSlope=NumSlope+1;


      if(m>(B-B1-1)){
        //Make OptDose and StoppedGroups
        StoreInx=m-B1;


        for(j=0;j<G;j++){;
          astore(StoreInx,j)=a(j);
          bstore(StoreInx,j)=b(j);
        }

        mustore[StoreInx]=mu;
        slopestore[StoreInx]=slope;
        sigstore[StoreInx]=sig;


        //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
        //First nDose

        aprop.zeros();
        bprop.zeros();
        for(j=0;j<G;j++){
          for(k=0;k<G;k++){
            if(GroupMem[j]==k){
              aprop[j]=a[k];
              bprop[j]=b[k];
            }
          }
        }




        for(k=0;k<G;k++){

          for(j=0;j<nDose;j++){

            eta1 = aprop(k)+exp(bprop(k)+slope)*Dose[j]+mu;




            //For all groups add the intercept and slope*Groups

            eta1=exp(eta1);

            if(eta1>100000){
              DoseProb(StoreInx, k*nDose+j) = 1;

            }else{
              DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

            }



          }


        }
        //Now we have a matrix containing our sampled dose probabilities









      }





      //End of MCMC

    }

  }


  //Borrow strength, AND cluster


  if(SubRout==2){



    for(m=0;m<B;m++){

      if(m<(B/2 + 2)){
        if(m%100==0){


          for(k=0;k<G;k++){
            if((IntA[k]/NumA[k])>.5){
              avar[k]=avar[k]*2;
            }

            if((IntA[k]/NumA[k])<.2){
              avar[k]=avar[k]/2;
            }








            if((IntB[k]/NumB[k])>.5){
              bvar[k]=bvar[k]*2;
            }

            if((IntB[k]/NumB[k])<.2){
              bvar[k]=bvar[k]/2;
            }



          }


          if((IntMu/NumMu)>.5){
            muvar=muvar*2;
          }

          if((IntMu/NumMu)<.2){
            muvar=muvar/2;
          }


          if((IntSlope/NumSlope)>.5){
            slopevar=slopevar*2;
          }

          if((IntSlope/NumSlope)<.2){
            slopevar=slopevar/2;
          }




          //Reset our counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB=IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          IntSig=1;
          NumSig=2;




        }
      }



      for(k=0;k<aprop.n_rows;k++){
        //We must have the right group AND it must NOT be stopped.
        if((INVEC(k)>0) && (Stopped(k)==0)){
          //Propose new value
          aprop=a;
          aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);
          //Get our Alphas
          //Prior Ratio
          alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 ;
          //Loglikelihood ratio
          alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);
          //Metropolis Hastings
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            a(k)=aprop(k);
            IntA[k]=IntA[k]+1;
          }

          NumA[k]=NumA[k]+1;


        }
      }
      //Now do a MH step for the slope

      for(k=0;k<aprop.size();k++){

        if((INVEC(k)>0) && (Stopped(k)==0)){
          bprop=b;
          bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];
          //Prior Ratio
          alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1;
          //Likelihood ratio
          alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            b[k]=bprop[k];
            IntB[k]=IntB[k]+1;
          }
          NumB[k]=NumB[k]+1;

        }

      }




      //Clustering moves









      Munew = mu + as_scalar(arma::randn(1))*muvar;



      alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint;

      alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);


      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        mu=Munew;
        IntMu=IntMu+1;
      }

      NumMu=NumMu+1;







      U=as_scalar(arma::randu(1));
      //Global Clustering move
      if(U<.333333){
        //Cluster or uncluster ALL groups
        if(sum(INVEC)==G){
          //All groups are unclustered.
          //Auto Delete Move

          GroupMemProp=GroupMem;
          //Now Randomly Draw Our New Group Assignment
          INVECNEW.zeros();
          //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
          GroupMemProp.zeros();
          g=GetRandGroup(INVEC,Stopped); //Holds the group spot
          GroupMemProp=GroupMemProp+g;
          //Alphaprop and betaprop are equal to the coefficients for group g
          aprop.zeros();
          bprop.zeros();
          aprop =aprop+a(g);
          bprop=bprop+b(g);

          alpha=0;
          //Here are the old priors ALL are UNCLUSTERED
          for(k=0;k<bprop.n_rows;k++){
            if(INVEC[k]==1){
              alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);
            }
          }


          //Density when 0 is automatically 1 it's spiked
          alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



          U=log(as_scalar(arma::randu(1)));



          if(U<alpha){
            b=bprop;
            a=aprop;
            GroupMem=GroupMemProp;
            INVEC.zeros();
            //Only one entry is in it's home bin
            INVEC(g)=1; //This is the only one included
          }






        }else{


          if(sum(INVEC)==1){
            //Right now EVERYTHING is clustered in one bin
            //Let's propose unclustering EVERYTHING


            for(k=0;k<aprop.n_rows;k++){
              if(INVEC[k]==0){
                //Let's only do this move on those that are clustered, i.e. in other bins
                aprop[k]=a[k]+as_scalar(arma::randn(1));
                bprop[k]=b[k]+as_scalar(arma::randn(1));
                GroupMemProp[k]=k;
              }
            }



            alpha=0;
            //Here are the new priors ALL are NEWLY UNCLUSTERED
            for(k=0;k<bprop.n_rows;k++){
              if(INVEC[k]==0){
                alpha=alpha -  .5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);
              }
            }


            //Density when 0 is automatically 1 it's spiked
            alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);






            U=log(as_scalar(arma::randu(1)));



            if(U<alpha){

              for(k=0;k<aprop.n_rows;k++){
                b[k]=bprop[k];
                a[k]=aprop[k];
                INVEC[k]=1;
                GroupMem[k]=k;
              }


            }





          }else{
            //Randomly Decide to add or Delete ALL

            U=as_scalar(arma::randu(1));


            if(U<.5){
              //Add ALL


              //We only remove those that are NOT in it's home bin
              for(k=0;k<aprop.n_rows;k++){
                if(INVEC[k]==0){
                  //Let's only do this move on those that are clustered, i.e. in other bins
                  aprop[k]=a[k]+as_scalar(arma::randn(1));
                  bprop[k]=b[k]+as_scalar(arma::randn(1));
                  GroupMemProp[k]=k;
                }
              }




              alpha=0;
              //Here are the new priors ALL are NEWLY UNCLUSTERED
              for(k=0;k<bprop.n_rows;k++){
                if(INVEC[k]==0){
                  alpha=alpha -  .5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);
                }
              }


              //Density when 0 is automatically 1 it's spiked
              alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);






              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){

                for(k=0;k<aprop.n_rows;k++){
                  b[k]=bprop[k];
                  a[k]=aprop[k];
                  INVEC[k]=1;
                  GroupMem[k]=k;
                }


              }



            }else{
              //Cluster everything in a bin
              //All groups are unclustered.
              //Auto Delete Move
              //Now Randomly Draw Our New Group Assignment
              INVECNEW.zeros();
              //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
              GroupMemProp.zeros();
              g=GetRandGroup(INVEC,Stopped); //Holds the group spot
              GroupMemProp=GroupMemProp+g;
              //Alphaprop and betaprop are equal to the coefficients for group g
              aprop.zeros();
              bprop.zeros();
              aprop =aprop+a(g);
              bprop=bprop+b(g);

              alpha=0;
              //Here are the old priors ALL are UNCLUSTERED
              for(k=0;k<bprop.n_rows;k++){
                if(INVEC[k]==1){
                  alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);
                }
              }


              //Density when 0 is automatically 1 it's spiked
              alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){
                b=bprop;
                a=aprop;
                GroupMem=GroupMemProp;
                INVEC.zeros();
                //Only one entry is in it's home bin
                INVEC(g)=1; //This is the only one included
              }




            }





          }




        }









      }else{


        //Must do delete move
        if(sum(INVEC)==G){
          //Delete randomly

          GroupMemProp=GroupMem;
          //Now Randomly Draw Our New Group Assignment
          //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
          GroupMemProp=GroupMem;
          g=GetRandGroup(INVEC,Stopped); //Holds the group spot to DELETE
          //Where should we put it??
          INVECNEW=INVEC;
          INVECNEW[g]=0;
          j=GetRandGroup(INVECNEW,Stopped); //Holds the group to cluster with!!
          GroupMemProp[g]=j;

          //Set up proposals
          aprop=a;
          bprop=b;
          //Cluster
          aprop(g) =a(j);
          bprop(g)=b(j);

          alpha=0;
          //Here are the old priors ALL are UNCLUSTERED

          //I just need to delete entry g
          alpha=alpha +  .5*pow((b[g]-MeanSlopes[g]),2)/varbeta1 +.5*pow((a[g]-MeanInts[g]),2)/varint1 + log(1-PROBIN[g])-log(PROBIN[g]);


          //Density when 0 is automatically 1 it's spiked
          alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



          U=log(as_scalar(arma::randu(1)));



          if(U<alpha){
            b=bprop;
            a=aprop;
            GroupMem=GroupMemProp;
            INVEC=INVECNEW;

          }






        }else{


          if(sum(INVEC)==1){
            //Must do add move

            //Add move



            GroupMemProp=GroupMem;
            //Now Randomly Draw Our New Group Assignment
            //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
            GroupMemProp=GroupMem;
            g=GetRandGroup(1-INVEC,Stopped); //Holds the group spot to Return to it's cluster
            //Reset GorupMemProp
            GroupMemProp[g]=g;
            INVECNEW=INVEC;
            INVECNEW[g]=1;
            //Set up proposals
            aprop=a;
            bprop=b;
            //Cluster
            aprop(g) =a(g)+as_scalar(arma::randn(1));
            bprop(g)=b(g)+as_scalar(arma::randn(1));

            alpha=0;
            //Here are the old priors ALL are UNCLUSTERED

            //I just need to delete entry g
            alpha=alpha -  .5*pow((bprop[g]-MeanSlopes[g]),2)/varbeta1 -.5*pow((aprop[g]-MeanInts[g]),2)/varint1 - log(1-PROBIN[g])+log(PROBIN[g]);


            //Density when 0 is automatically 1 it's spiked
            alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



            U=log(as_scalar(arma::randu(1)));



            if(U<alpha){
              b=bprop;
              a=aprop;
              GroupMem=GroupMemProp;
              INVEC=INVECNEW;

            }





          }else{

            //Add or delete move here

            if(U<.66666){
              //Delete randomly

              GroupMemProp=GroupMem;
              //Now Randomly Draw Our New Group Assignment
              //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
              GroupMemProp=GroupMem;
              g=GetRandGroup(INVEC,Stopped); //Holds the group spot to DELETE
              //Where should we put it??
              INVECNEW=INVEC;
              INVECNEW[g]=0;
              j=GetRandGroup(INVECNEW,Stopped); //Holds the group to cluster with!!
              GroupMemProp[g]=j;

              //Set up proposals
              aprop=a;
              bprop=b;
              //Cluster
              aprop(g) =a(j);
              bprop(g)=b(j);

              alpha=0;
              //Here are the old priors ALL are UNCLUSTERED

              //I just need to delete entry g
              alpha=alpha +  .5*pow((b[g]-MeanSlopes[g]),2)/varbeta1 +.5*pow((a[g]-MeanInts[g]),2)/varint1 + log(1-PROBIN[g])-log(PROBIN[g]);


              //Density when 0 is automatically 1 it's spiked
              alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){
                b=bprop;
                a=aprop;
                GroupMem=GroupMemProp;
                INVEC=INVECNEW;

              }





            }else{
              //Add move



              GroupMemProp=GroupMem;
              //Now Randomly Draw Our New Group Assignment
              //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
              GroupMemProp=GroupMem;
              g=GetRandGroup(1-INVEC,Stopped); //Holds the group spot to Return to it's cluster
              //Reset GroupMemProp
              GroupMemProp[g]=g;


              //Where should we put it??
              INVECNEW=INVEC;
              INVECNEW[g]=1;
              //Set up proposals
              aprop=a;
              bprop=b;
              //Cluster
              aprop(g) =a(g)+as_scalar(arma::randn(1));
              bprop(g)=b(g)+as_scalar(arma::randn(1));

              alpha=0;
              //Here are the old priors ALL are UNCLUSTERED

              //I just need to delete entry g
              alpha=alpha -  .5*pow((bprop[g]-MeanSlopes[g]),2)/varbeta1 -.5*pow((aprop[g]-MeanInts[g]),2)/varint1 - log(1-PROBIN[g])+log(PROBIN[g]);


              //Density when 0 is automatically 1 it's spiked
              alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){
                b=bprop;
                a=aprop;
                GroupMem=GroupMemProp;
                INVEC=INVECNEW;

              }



            }











          }

        }
      }

      // Rf_PrintValue(wrap(GroupMem));


      //Sample the new slope and the new \beta_g coefficients jointly




      slopenew = slope + as_scalar(arma::randn(1))*slopevar;




      alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


      alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);


      U=log(as_scalar(arma::randu(1)));




      if(U<alpha){
        slope=slopenew;
        IntSlope=IntSlope+1;
      }
      NumSlope=NumSlope+1;


      if(m>(B-B1-1)){
        //Make OptDose and StoppedGroups
        StoreInx=m-B1;


        for(j=0;j<G;j++){
          CLUST(StoreInx,j) = GroupMem[j];
          INSTORE(StoreInx,j)=INVEC[j];
        }
        //GroupCLUST.row(StoreInx) = GroupMem.t();
        // Rf_PrintValue(wrap(CLUST.row(StoreInx)));



        //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
        //First nDose


        aprop.zeros();
        bprop.zeros();
        for(j=0;j<G;j++){
          for(k=0;k<G;k++){
            if(GroupMem[j]==k){
              aprop[j]=a[k];
              bprop[j]=b[k];
            }
          }
        }




        for(k=0;k<G;k++){

          for(j=0;j<nDose;j++){

            eta1 = aprop(k)+exp(bprop(k)+slope)*Dose[j]+mu;




            //For all groups add the intercept and slope*Groups

            eta1=exp(eta1);

            if(eta1>100000){
              DoseProb(StoreInx, k*nDose+j) = 1;

            }else{
              DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

            }



          }


        }
        //Now we have a matrix containing our sampled dose probabilities









      }





      //End of MCMC

    }

  }




  arma::vec STOPPROB(G);
  STOPPROB.zeros();

  //Do we have any new stopped groups?
  for(k=0;k<G;k++){
    //Don't check this for already stopped groups
    if(Stopped(k)==0){
      eta2=0;

      for(j=0;j<DoseProb.n_rows;j++){
        if(DoseProb(j,k*nDose)>Target[k]){
          eta2=eta2+1;
        }
      }

      STOPPROB(k)=eta2/DoseProb.n_rows;

      //Is this bigger than our posterior probability threshold???
      if(eta2>(Upper[k]*DoseProb.n_rows)){
        Stopped(k)=1;
      }

    }

  }


  //Let's return a storage matrix with the probs and also
  //Stopped groups
  //LAST VALUE contains stopped indicator
  arma::mat RETURNMAT(G,nDose+2); //Number of groups

  //Get Mean Toxicity probabilities
  for(k=0;k<G;k++){
    for(j=0;j<nDose;j++){
      RETURNMAT(k,j)=sum(DoseProb.col(k*nDose+j))/DoseProb.n_rows;
    }

    //Fill in the stopped value and the stopped
    RETURNMAT(k,nDose)=Stopped(k);

    //Fill in the stopping probability
    RETURNMAT(k,nDose+1)=STOPPROB(k);

  }



  //Rf_PrintValue(wrap(CLUST));



  List z1 = List::create(RETURNMAT,CLUST,INSTORE);




  return(z1);

}











//' Performs MCMC and returns needed values for dose-finding in a list.
//' @param Y Vector containing observed event or censoring times.
//' @param I Vector containing event indicators (1 if patient experiences an event for a patient).
//' @param Doses Vector containing Doses of patients in trial.
//' @param Groups Vector containing group assignment of patients, 0 is baseline group.
//' @param T1 Reference time for toxicity.
//' @param Target Target cumulative toxicity probability vector at time T1.
//' @param Dose Vector containing the standardized doses considered.
//' @param Upper Cutoff values used to determine if accrual in a subgroup should be suspended.
//' @param meanmu Prior mean for baseline intercept.
//' @param meanslope Prior mean for baseline slope.
//' @param MeanInts Vector of prior means for the group specific intercept parameters.
//' @param MeanSlopes Vector of prior means for the group specific slope parameters.
//' @param varint Prior variance for the intercept parameters.
//' @param varbeta Prior variance for the slope parameters.
//' @param phetero Prior probability of heterogeneous subgroups.
//' @param SubRout Parameter to specify subgroup borrowing/clustering. 0=No borrowing, 1=Borrowing but no clustering, 2=Borrowing and clustering.
//' @param B Number of Iterations to run for MCMC
//' @param Stopped Current vector of STOPPED groups
//' @param NumPat Number of patients
//'@importFrom Rcpp evalCpp
//'@useDynLib SubTite
//'@return A matrix of quantities needed for determining the next dose to enroll each subgroup while using the SimTrial function.
//'@export
//[[Rcpp::export]]
arma::mat MCMCSIM( arma::vec Y, //Vector of Times
                   arma::vec I, //Vector of Toxicity Indicators
                   arma::vec Doses, //Standardized Dose Values given to each patient
                   arma::vec Groups, //SubGroup Membership
                   double T1, //Reference Time
                   arma::vec Target, //Target Toxicity Probability for each group
                   arma::vec Upper, //Vector of Upper Cutoffs for each groups
                   arma::vec Dose, //Standardized values of Doses considered
                   double meanmu, //This is the prior mean of the baseline int
                   double meanslope, //Prior mean of baseline slope
                   arma::vec MeanInts, //Prior Means of Intercept Parameters for all G groups
                   arma::vec MeanSlopes, //Prior Means of Slope Parameters for all G groups
                   double varint, //Prior Variance of Intercept Parameters
                   double varbeta, //Prior Variance of Slope Parameters
                   double phetero, //prior probability of heterogeneity
                   arma::vec Stopped, //Current vector of STOPPED groups
                   int NumPat, //Number of patients
                   int SubRout, // Contains a 0 for simple model with no borrowing, 1 for borrowing model, 2 for borrow-clustering model.
                   double B
){


  int g=0;

  int StoreInx=0;

  //Change this to be an input
  arma::vec PROBIN(MeanInts.n_rows);
  for(g=0;g<PROBIN.n_rows;g++){
    //Prior probability of NO clustering
    PROBIN(g)=phetero;
  }

  int MEANS=0;

  int IntIN=0;
  int IntOUT=0;

  //Group Slope Prior Var
  double varbeta1=varbeta;
  //Group Int Prior Var
  double varint1=varint;



  int Group=0;
  //Important For loop integer
  int m =0; //For the inner MCMC looprep
  int i=0; //For the patient index
  int rep=0; //For the simulation repetition.


  //First Fill out Target Vector
  int nDose=Dose.n_rows;



  double B1=B/2;

  //Make List objects we'll use





  //Innitialize parameters for MCMC
  double mu=meanmu;
  double slope=meanslope;
  double sig=T1;

  double NewMean=0;
  double NewSlope=0;
  int Which1=0;
  int J=Dose.n_rows; //Number of Doses


  //Important quantities we will need in the MCMC
  int  G=MeanInts.n_rows; //Number of groups
  double alpha=0;
  double U=0;
  double signew=0;
  double slopenew=0;
  double Munew=0;

  //For loop integers
  int j=0;
  int k=0;


  //Vector Of group intercepts
  arma::vec a(G);
  a.zeros();
  //Vector of group slopes
  arma::vec b=a;
  //Vectors for Group MCMC



  //Vectors for proposals
  arma::vec aprop = a;
  arma::vec bprop=b;





  arma::vec NumA(G);
  arma::vec NumB(G);
  arma::vec IntA(G);
  arma::vec IntB(G);

  NumA.zeros();
  NumB.zeros();
  IntA.zeros();
  IntB.zeros();

  double NumMu = 2;
  double NumSlope=2;
  double IntMu=1;
  double IntSlope=1;

  arma::vec avar(G);
  arma::vec bvar(G);
  avar.zeros();
  avar=avar+1;
  bvar.zeros();
  bvar=bvar+1;


  double muvar=1;  ///Proposal variance for shared intercept
  double slopevar=1; //Proposal variance for shared slope
  double trialtime=0;
  NumericVector z9(2);
  NumericVector zprop(5);

  double NumSig=2;
  double IntSig=1;
  double sigvar=1;

  arma::vec etaG(G);
  arma::mat DoseProb(B1,G*nDose);
  arma::mat MeanDose(G,nDose);
  arma::vec sigstore(B1); //What is sigstore??
  arma::vec mustore(B1); //Storage for shared intercept
  arma::vec slopestore(B1); //Storage for shared slope
  arma::mat astore(B1,G); //Storage for group specific intercepts
  arma::mat bstore=astore; //Storage for group specific slopes
  double eta1=0;
  double DEC=0;

  arma::vec GroupMem(G);
  GroupMem.zeros();
  arma::vec GroupMemProp=GroupMem;

  arma::vec INVEC(G);
  INVEC.zeros();

  double eta2=0;

  arma::mat MeanVec(G,nDose);
  arma::vec NTox(G);
  arma::vec OptDose(G);
  MeanVec.zeros();
  NTox.zeros();
  OptDose.zeros();
  arma::mat CLUST(B1,G);
  CLUST.zeros();
  arma::mat INSTORE=CLUST;







  arma::vec INVECNEW=INVEC;
  arma::vec Y2;
  arma::vec I2;
  arma::vec Groups2;
  arma::vec Doses2;

  Y2.zeros();
  I2.zeros();
  Groups2.zeros();
  Doses2.zeros();

  //Start all parameters at their prior distributions
  mu=meanmu;
  slope=meanslope;
  INVEC.zeros();
  INVECNEW.zeros();
  INVEC = INVEC+1;
  a=MeanInts;
  b=MeanSlopes;

  //Start counter
  for(m=0;m<GroupMem.n_rows;m++){
    GroupMem[m]=m;
  }
  //Proposal for combinations
  GroupMemProp=GroupMem;

  //No Borrowing, No Clustering

  if(SubRout==0){
    //Reset our vectors
    MeanInts = MeanInts + mu;
    MeanSlopes = MeanInts + slope;
    a=MeanInts;
    b=MeanSlopes;

    for(m=0;m<B;m++){

      if(m<(B/2 + 2)){
        if(m%100==0){



          for(k=0;k<G;k++){
            if((IntA[k]/NumA[k])>.5){
              avar[k]=avar[k]*2;
            }

            if((IntA[k]/NumA[k])<.2){
              avar[k]=avar[k]/2;
            }








            if((IntB[k]/NumB[k])>.5){
              bvar[k]=bvar[k]*2;
            }

            if((IntB[k]/NumB[k])<.2){
              bvar[k]=bvar[k]/2;
            }



          }




          //Reset our counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB=IntB.zeros()+1;






        }
      }



      for(k=0;k<aprop.n_rows;k++){
        //We must have the right group AND it must NOT be stopped.
        if(Stopped(k)==0){
          //Propose new value
          aprop=a;
          aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);
          //Get our Alphas
          //Prior Ratio
          alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 ;
          //Loglikelihood ratio
          alpha=alpha+   LikeStoppedSeparate(Y, I,  Doses, Groups,  aprop, b, sig,NumPat,Stopped) -  LikeStoppedSeparate(Y, I,  Doses, Groups,  a, b, sig,NumPat, Stopped);
          //Metropolis Hastings
          U=log(as_scalar(arma::randu(1)));
          z9(0)=LikeStoppedSeparate(Y, I,  Doses, Groups,  a, b, sig,NumPat, Stopped);
          z9(1)=LikeStoppedSeparate(Y, I,  Doses, Groups,  aprop, b, sig,NumPat,Stopped) ;


          //          Rf_PrintValue(wrap(z9));

          if(U<alpha){
            a(k)=aprop(k);
            IntA[k]=IntA[k]+1;
          }

          NumA[k]=NumA[k]+1;


        }
      }
      //Now do a MH step for the slope

      for(k=0;k<aprop.size();k++){

        if(Stopped(k)==0){
          bprop=b;
          bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];
          //Prior Ratio
          alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1;
          //Likelihood ratio
          alpha=alpha+   LikeStoppedSeparate(Y, I,  Doses, Groups,   a, bprop, sig,NumPat,Stopped) -  LikeStoppedSeparate(Y, I,  Doses, Groups,  a, b, sig,NumPat, Stopped);
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            b[k]=bprop[k];
            IntB[k]=IntB[k]+1;
          }
          NumB[k]=NumB[k]+1;

        }

      }







      if(m>(B-B1-1)){
        //Make OptDose and StoppedGroups
        StoreInx=m-B1;

        for(j=0;j<G;j++){;
          astore(StoreInx,j)=a(j);
          bstore(StoreInx,j)=b(j);
        }



        //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
        //First nDose



        for(k=0;k<G;k++){

          for(j=0;j<nDose;j++){
            eta1=0;

            eta1 = a(k)+exp(b(k))*Dose[j];

            z9(0)=eta1;


            //For all groups add the intercept and slope*Groups

            eta1=exp(eta1);




            if(eta1>100000){
              DoseProb(StoreInx, k*nDose+j) = 1;

            }else{
              DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

            }



          }


        }
        //Now we have a matrix containing our sampled dose probabilities









      }





      //End of MCMC

    }

  }


  //Borrow strength, No cluster

  if(SubRout==1){


    for(m=0;m<B;m++){

      if(m<(B/2 + 2)){
        if(m%100==0){


          for(k=0;k<G;k++){
            if((IntA[k]/NumA[k])>.5){
              avar[k]=avar[k]*2;
            }

            if((IntA[k]/NumA[k])<.2){
              avar[k]=avar[k]/2;
            }








            if((IntB[k]/NumB[k])>.5){
              bvar[k]=bvar[k]*2;
            }

            if((IntB[k]/NumB[k])<.2){
              bvar[k]=bvar[k]/2;
            }



          }


          if((IntMu/NumMu)>.5){
            muvar=muvar*2;
          }

          if((IntMu/NumMu)<.2){
            muvar=muvar/2;
          }


          if((IntSlope/NumSlope)>.5){
            slopevar=slopevar*2;
          }

          if((IntSlope/NumSlope)<.2){
            slopevar=slopevar/2;
          }




          //Reset our counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB=IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          IntSig=1;
          NumSig=2;




        }
      }



      for(k=0;k<aprop.n_rows;k++){
        //We must have the right group AND it must NOT be stopped.
        if((INVEC(k)>0) && (Stopped(k)==0)){
          //Propose new value
          aprop=a;
          aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);
          //Get our Alphas
          //Prior Ratio
          alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 ;
          //Loglikelihood ratio
          alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);
          //Metropolis Hastings
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            a(k)=aprop(k);
            IntA[k]=IntA[k]+1;
          }

          NumA[k]=NumA[k]+1;


        }
      }
      //Now do a MH step for the slope

      for(k=0;k<aprop.size();k++){

        if((INVEC(k)>0) && (Stopped(k)==0)){
          bprop=b;
          bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];
          //Prior Ratio
          alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1;
          //Likelihood ratio
          alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            b[k]=bprop[k];
            IntB[k]=IntB[k]+1;
          }
          NumB[k]=NumB[k]+1;

        }

      }









      Munew = mu + as_scalar(arma::randn(1))*muvar;



      alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint;

      alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);


      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        mu=Munew;
        IntMu=IntMu+1;
      }

      NumMu=NumMu+1;




      //Sample the new slope and the new \beta_g coefficients jointly




      slopenew = slope + as_scalar(arma::randn(1))*slopevar;




      alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


      alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);


      U=log(as_scalar(arma::randu(1)));




      if(U<alpha){
        slope=slopenew;
        IntSlope=IntSlope+1;
      }
      NumSlope=NumSlope+1;


      if(m>(B-B1-1)){
        //Make OptDose and StoppedGroups
        StoreInx=m-B1;


        for(j=0;j<G;j++){;
          astore(StoreInx,j)=a(j);
          bstore(StoreInx,j)=b(j);
        }

        mustore[StoreInx]=mu;
        slopestore[StoreInx]=slope;
        sigstore[StoreInx]=sig;


        //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
        //First nDose

        aprop.zeros();
        bprop.zeros();
        for(j=0;j<G;j++){
          for(k=0;k<G;k++){
            if(GroupMem[j]==k){
              aprop[j]=a[k];
              bprop[j]=b[k];
            }
          }
        }




        for(k=0;k<G;k++){

          for(j=0;j<nDose;j++){

            eta1 = aprop(k)+exp(bprop(k)+slope)*Dose[j]+mu;




            //For all groups add the intercept and slope*Groups

            eta1=exp(eta1);

            if(eta1>100000){
              DoseProb(StoreInx, k*nDose+j) = 1;

            }else{
              DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

            }



          }


        }
        //Now we have a matrix containing our sampled dose probabilities









      }





      //End of MCMC

    }

  }

  int m7=0; //For clustering

  //Borrow strength, AND cluster
  if(SubRout==2){



    for(m=0;m<B;m++){

      if(m<(B/2 + 2)){
        if(m%100==0){


          for(k=0;k<G;k++){
            if((IntA[k]/NumA[k])>.5){
              avar[k]=avar[k]*2;
            }

            if((IntA[k]/NumA[k])<.2){
              avar[k]=avar[k]/2;
            }








            if((IntB[k]/NumB[k])>.5){
              bvar[k]=bvar[k]*2;
            }

            if((IntB[k]/NumB[k])<.2){
              bvar[k]=bvar[k]/2;
            }



          }


          if((IntMu/NumMu)>.5){
            muvar=muvar*2;
          }

          if((IntMu/NumMu)<.2){
            muvar=muvar/2;
          }


          if((IntSlope/NumSlope)>.5){
            slopevar=slopevar*2;
          }

          if((IntSlope/NumSlope)<.2){
            slopevar=slopevar/2;
          }




          //Reset our counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB=IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          IntSig=1;
          NumSig=2;




        }
      }



      for(k=0;k<aprop.n_rows;k++){
        //We must have the right group AND it must NOT be stopped.
        if((INVEC(k)>0) && (Stopped(k)==0)){
          //Propose new value
          aprop=a;
          aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);
          //Get our Alphas
          //Prior Ratio
          alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 ;
          //Loglikelihood ratio
          alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);
          //Metropolis Hastings
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            a(k)=aprop(k);
            IntA[k]=IntA[k]+1;
          }

          NumA[k]=NumA[k]+1;


        }
      }
      //Now do a MH step for the slope

      for(k=0;k<aprop.size();k++){

        if((INVEC(k)>0) && (Stopped(k)==0)){
          bprop=b;
          bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];
          //Prior Ratio
          alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1;
          //Likelihood ratio
          alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);
          U=log(as_scalar(arma::randu(1)));
          if(U<alpha){
            b[k]=bprop[k];
            IntB[k]=IntB[k]+1;
          }
          NumB[k]=NumB[k]+1;

        }

      }




      //Clustering moves









      Munew = mu + as_scalar(arma::randn(1))*muvar;



      alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint;

      alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);


      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        mu=Munew;
        IntMu=IntMu+1;
      }

      NumMu=NumMu+1;







      U=as_scalar(arma::randu(1));
      //Global Clustering move
      if(U<.333333){
        //Cluster or uncluster ALL groups
        if(sum(INVEC)==G){
          //All groups are unclustered.
          //Auto Delete Move

          GroupMemProp=GroupMem;
          //Now Randomly Draw Our New Group Assignment
          INVECNEW.zeros();
          //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
          GroupMemProp.zeros();
          g=GetRandGroup(INVEC,Stopped); //Holds the group spot
          GroupMemProp=GroupMemProp+g;
          //Alphaprop and betaprop are equal to the coefficients for group g
          aprop.zeros();
          bprop.zeros();
          aprop =aprop+a(g);
          bprop=bprop+b(g);

          alpha=0;
          //Here are the old priors ALL are UNCLUSTERED
          for(k=0;k<bprop.n_rows;k++){
            if(INVEC[k]==1){
              alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);
            }
          }


          //Density when 0 is automatically 1 it's spiked
          alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



          U=log(as_scalar(arma::randu(1)));



          if(U<alpha){
            b=bprop;
            a=aprop;
            GroupMem=GroupMemProp;
            INVEC.zeros();
            //Only one entry is in it's home bin
            INVEC(g)=1; //This is the only one included
          }






        }else{


          if(sum(INVEC)==1){
            //Right now EVERYTHING is clustered in one bin
            //Let's propose unclustering EVERYTHING


            for(k=0;k<aprop.n_rows;k++){
              if(INVEC[k]==0){
                //Let's only do this move on those that are clustered, i.e. in other bins
                aprop[k]=a[k]+as_scalar(arma::randn(1));
                bprop[k]=b[k]+as_scalar(arma::randn(1));
                GroupMemProp[k]=k;
              }
            }



            alpha=0;
            //Here are the new priors ALL are NEWLY UNCLUSTERED
            for(k=0;k<bprop.n_rows;k++){
              if(INVEC[k]==0){
                alpha=alpha -  .5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);
              }
            }


            //Density when 0 is automatically 1 it's spiked
            alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);






            U=log(as_scalar(arma::randu(1)));



            if(U<alpha){

              for(k=0;k<aprop.n_rows;k++){
                b[k]=bprop[k];
                a[k]=aprop[k];
                INVEC[k]=1;
                GroupMem[k]=k;
              }


            }





          }else{
            //Randomly Decide to add or Delete ALL

            U=as_scalar(arma::randu(1));


            if(U<.5){
              //Add ALL


              //We only remove those that are NOT in it's home bin
              for(k=0;k<aprop.n_rows;k++){
                if(INVEC[k]==0){
                  //Let's only do this move on those that are clustered, i.e. in other bins
                  aprop[k]=a[k]+as_scalar(arma::randn(1));
                  bprop[k]=b[k]+as_scalar(arma::randn(1));
                  GroupMemProp[k]=k;
                }
              }




              alpha=0;
              //Here are the new priors ALL are NEWLY UNCLUSTERED
              for(k=0;k<bprop.n_rows;k++){
                if(INVEC[k]==0){
                  alpha=alpha -  .5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);
                }
              }


              //Density when 0 is automatically 1 it's spiked
              alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);






              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){

                for(k=0;k<aprop.n_rows;k++){
                  b[k]=bprop[k];
                  a[k]=aprop[k];
                  INVEC[k]=1;
                  GroupMem[k]=k;
                }


              }



            }else{
              //Cluster everything in a bin
              //All groups are unclustered.
              //Auto Delete Move
              //Now Randomly Draw Our New Group Assignment
              INVECNEW.zeros();
              //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
              GroupMemProp.zeros();
              g=GetRandGroup(INVEC,Stopped); //Holds the group spot
              GroupMemProp=GroupMemProp+g;
              //Alphaprop and betaprop are equal to the coefficients for group g
              aprop.zeros();
              bprop.zeros();
              aprop =aprop+a(g);
              bprop=bprop+b(g);

              alpha=0;
              //Here are the old priors ALL are UNCLUSTERED
              for(k=0;k<bprop.n_rows;k++){
                if(INVEC[k]==1){
                  alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);
                }
              }


              //Density when 0 is automatically 1 it's spiked
              alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){
                b=bprop;
                a=aprop;
                GroupMem=GroupMemProp;
                INVEC.zeros();
                //Only one entry is in it's home bin
                INVEC(g)=1; //This is the only one included
              }




            }





          }




        }









      }else{


        //Must do delete move
        if(sum(INVEC)==G){
          //Delete randomly

          GroupMemProp=GroupMem;
          //Now Randomly Draw Our New Group Assignment
          //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
          GroupMemProp=GroupMem;
          g=GetRandGroup(INVEC,Stopped); //Holds the group spot to DELETE
          //Where should we put it??
          INVECNEW=INVEC;
          INVECNEW[g]=0;
          j=GetRandGroup(INVECNEW,Stopped); //Holds the group to cluster with!!
          GroupMemProp[g]=j;

          //Set up proposals
          aprop=a;
          bprop=b;
          //Cluster
          aprop(g) =a(j);
          bprop(g)=b(j);

          alpha=0;
          //Here are the old priors ALL are UNCLUSTERED

          //I just need to delete entry g
          alpha=alpha +  .5*pow((b[g]-MeanSlopes[g]),2)/varbeta1 +.5*pow((a[g]-MeanInts[g]),2)/varint1 + log(1-PROBIN[g])-log(PROBIN[g]);


          //Density when 0 is automatically 1 it's spiked
          alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



          U=log(as_scalar(arma::randu(1)));



          if(U<alpha){
            b=bprop;
            a=aprop;
            GroupMem=GroupMemProp;
            INVEC=INVECNEW;

          }






        }else{


          if(sum(INVEC)==1){
            //Must do add move

            //Add move



            GroupMemProp=GroupMem;
            //Now Randomly Draw Our New Group Assignment
            //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
            GroupMemProp=GroupMem;
            g=GetRandGroup(1-INVEC,Stopped); //Holds the group spot to Return to it's cluster
            //Reset GorupMemProp
            GroupMemProp[g]=g;
            INVECNEW=INVEC;
            INVECNEW[g]=1;
            //Set up proposals
            aprop=a;
            bprop=b;
            //Cluster
            aprop(g) =a(g)+as_scalar(arma::randn(1));
            bprop(g)=b(g)+as_scalar(arma::randn(1));

            alpha=0;
            //Here are the old priors ALL are UNCLUSTERED

            //I just need to delete entry g
            alpha=alpha -  .5*pow((bprop[g]-MeanSlopes[g]),2)/varbeta1 -.5*pow((aprop[g]-MeanInts[g]),2)/varint1 - log(1-PROBIN[g])+log(PROBIN[g]);


            //Density when 0 is automatically 1 it's spiked
            alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



            U=log(as_scalar(arma::randu(1)));



            if(U<alpha){
              b=bprop;
              a=aprop;
              GroupMem=GroupMemProp;
              INVEC=INVECNEW;

            }





          }else{

            //Add or delete move here

            if(U<.66666){
              //Delete randomly

              GroupMemProp=GroupMem;
              //Now Randomly Draw Our New Group Assignment
              //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
              GroupMemProp=GroupMem;
              g=GetRandGroup(INVEC,Stopped); //Holds the group spot to DELETE
              //Where should we put it??
              INVECNEW=INVEC;
              INVECNEW[g]=0;
              j=GetRandGroup(INVECNEW,Stopped); //Holds the group to cluster with!!
              GroupMemProp[g]=j;

              //Set up proposals
              aprop=a;
              bprop=b;
              //Cluster
              aprop(g) =a(j);
              bprop(g)=b(j);


              //Make sure that the other groups also move wiith this one
              //i.e. for all m: \zeta_m = g, set \zeta_m = j.
              //Also... Set a_m = a_j and b_m = a_j for all these.
              for(m7=0;m7<aprop.n_rows;m7++){
if(GroupMem(m7)==g){
  //Used to be clustered with group g, everything should go with it!
  aprop(m7)=a(j);
  bprop(m7)=b(j);
  GroupMemProp(m7)=j;

}
              }


              alpha=0;
              //Here are the old priors ALL are UNCLUSTERED

              //I just need to delete entry g
              alpha=alpha +  .5*pow((b[g]-MeanSlopes[g]),2)/varbeta1 +.5*pow((a[g]-MeanInts[g]),2)/varint1 + log(1-PROBIN[g])-log(PROBIN[g]);


              //Density when 0 is automatically 1 it's spiked
              alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){
                b=bprop;
                a=aprop;
                GroupMem=GroupMemProp;
                INVEC=INVECNEW;

              }





            }else{
              //Add move



              GroupMemProp=GroupMem;
              //Now Randomly Draw Our New Group Assignment
              //EVERYTHING IS CLUSTERED ON ONE RANDOM GROUP
              GroupMemProp=GroupMem;
              g=GetRandGroup(1-INVEC,Stopped); //Holds the group spot to Return to it's cluster
              //Reset GroupMemProp
              GroupMemProp[g]=g;


              //Where should we put it??
              INVECNEW=INVEC;
              INVECNEW[g]=1;
              //Set up proposals
              aprop=a;
              bprop=b;
              //Cluster
              aprop(g) =a(g)+as_scalar(arma::randn(1));
              bprop(g)=b(g)+as_scalar(arma::randn(1));

              //Nothing else should move with group g parameters. Group g is now being unclustered,
              //So nothing else can be clustered on top of group g at this iteration.


              alpha=0;
              //Here are the old priors ALL are UNCLUSTERED

              //I just need to delete entry g
              alpha=alpha -  .5*pow((bprop[g]-MeanSlopes[g]),2)/varbeta1 -.5*pow((aprop[g]-MeanInts[g]),2)/varint1 - log(1-PROBIN[g])+log(PROBIN[g]);


              //Density when 0 is automatically 1 it's spiked
              alpha =  alpha + LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,NumPat,GroupMemProp,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem,Stopped);



              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){
                b=bprop;
                a=aprop;
                GroupMem=GroupMemProp;
                INVEC=INVECNEW;

              }



            }











          }

        }
      }

      // Rf_PrintValue(wrap(GroupMem));


      //Sample the new slope and the new \beta_g coefficients jointly




      slopenew = slope + as_scalar(arma::randn(1))*slopevar;




      alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


      alpha=alpha+   LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,NumPat,GroupMem,Stopped) -  LikeStoppedCluster(Y, I,  Doses, Groups,  mu,slope, a, b, sig,NumPat,GroupMem, Stopped);


      U=log(as_scalar(arma::randu(1)));




      if(U<alpha){
        slope=slopenew;
        IntSlope=IntSlope+1;
      }
      NumSlope=NumSlope+1;


      if(m>(B-B1-1)){
        //Make OptDose and StoppedGroups
        StoreInx=m-B1;


        for(j=0;j<G;j++){
          CLUST(StoreInx,j) = GroupMem[j];
          INSTORE(StoreInx,j)=INVEC[j];
        }
        //GroupCLUST.row(StoreInx) = GroupMem.t();
        // Rf_PrintValue(wrap(CLUST.row(StoreInx)));



        //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
        //First nDose


        aprop.zeros();
        bprop.zeros();
        for(j=0;j<G;j++){
          for(k=0;k<G;k++){
            if(GroupMem[j]==k){
              aprop[j]=a[k];
              bprop[j]=b[k];
            }
          }
        }




        for(k=0;k<G;k++){

          for(j=0;j<nDose;j++){

            eta1 = aprop(k)+exp(bprop(k)+slope)*Dose[j]+mu;




            //For all groups add the intercept and slope*Groups

            eta1=exp(eta1);

            if(eta1>100000){
              DoseProb(StoreInx, k*nDose+j) = 1;

            }else{
              DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

            }



          }


        }
        //Now we have a matrix containing our sampled dose probabilities









      }





      //End of MCMC

    }

  }




  arma::vec STOPPROB(G);
  STOPPROB.zeros();

  //Do we have any new stopped groups?
  for(k=0;k<G;k++){
    //Don't check this for already stopped groups
    if(Stopped(k)==0){
      eta2=0;

      for(j=0;j<DoseProb.n_rows;j++){
        if(DoseProb(j,k*nDose)>Target[k]){
          eta2=eta2+1;
        }
      }

      STOPPROB(k)=eta2/DoseProb.n_rows;

      //Is this bigger than our posterior probability threshold???
      if(eta2>(Upper[k]*DoseProb.n_rows)){
        Stopped(k)=1;
      }

    }

  }


  //Let's return a storage matrix with the probs and also
  //Stopped groups
  //LAST VALUE contains stopped indicator
  arma::mat RETURNMAT(G,nDose+1); //Number of groups

  //Get Mean Toxicity probabilities
  for(k=0;k<G;k++){
    for(j=0;j<nDose;j++){
      RETURNMAT(k,j)=sum(DoseProb.col(k*nDose+j))/DoseProb.n_rows;
    }

    //Fill in the stopped value and the stopped
    RETURNMAT(k,nDose)=Stopped(k);



  }



  //Rf_PrintValue(wrap(CLUST));


  // Rprintf("OUTMCMC");


  return(RETURNMAT);

}








//' Simulates a Sub-TITE trial design
//'
//' Simulates replicates from a Sub-TITE trial with user specified true toxicity time distributions for different doses and subgroups and returns average summary statistics of the trial.
//' @param nSims Number of Trials to Simulate.
//' @param Nmax Maximum Number of Patients to enroll in the trial.
//' @param T1 Reference time for toxicity.
//' @param Target Target cumulative toxicity probability (or subgroup specific vector) at time T1.
//' @param Dose Standardized vector of doses to try.
//' @param DoseStart Dose (or vector of Doses) to enroll the first patient in each subgroup at.
//' @param Accrue Expected montly patient accrual rate.
//' @param groupprob Probability vector of subgroup assignment.
//' @param Upper Cutoff values used to determine if accrual in a subgroup should be suspended.
//' @param meanmu Prior mean of the baseline intercept parameter.
//' @param meanslope Prior mean of the baseline slope parameter.
//' @param MeanInts G-1 length vector of subgroup specific prior intercept means.
//' @param MeanSlopes G-1 length vector of subgroup specific prior slope means.
//' @param Family What distribution Family to simulate from. Options include: Exponential,Gamma, Lognormal, Uniform, Weibull.
//' @param Param1 nGroups X nDose matrix of first parameter values.
//' @param Param2 NGroups X nDose matrix of second parameter values.
//' @param varint Prior Variance of Intercept Parameters.
//' @param varbeta Prior Variance of Slope Parameters.
//' @param phetero Prior prob of heterogeneity.
//' @param NSep Number of patients to assign based on no borrowing.
//' @param NBorrow Number of patients to assign based on no clustering
//' @param cohort Number of patients to enroll before escalating.
//' @param FULLY Do we have to fully evaluate a cohort before escalating?
//' @return A list of simulation outputs to be processed in R.
//[[Rcpp::export]]
List SimTrial1(int nSims, //Number of Simulations to Run
               int Nmax, //Max patients to enroll in trial
               double  T1, //Reference Time
               arma::vec Target, //Reference Probability Target (or vector)
               arma::vec Dose, //Standardized Vector of Doses
               arma::vec DoseStart, //Vector (or just one) Starting Dose Level
               arma::vec Upper, //Thresholds for subgroup suspension
               double Accrue, //Expected Monthly accrual rate
               arma::vec groupprob, //Subgroup Enrollment Probabilities
               int Family, //Indicator of Distribution family to simulate from. Alphabetical order
               arma::mat Param1, //Group Specific Parameter matrix for each dose and subgroup to simulate from
               arma::mat Param2,
               double meanmu, //Prior mean of intercept
               double meanslope, //Prior mean of slope
               arma::vec MeanInts, //Prior means of group intercepts
               arma::vec MeanSlopes,
               double varint,
               double varbeta,
               double phetero, //prior probability of heterogeneity
               double NSep, //Number of patients to enroll before borrowing.
               double NBorrow, //Number of patients to enroll before clustering.
               int cohort, //number of patients to fully evaluate before escalating.
               int FULLY ){ // Do we fully evaluate a cohort before escalating?

  int g=0;

  //Change this to be an input
  arma::vec PROBIN(Param1.n_rows-1);
  for(g=0;g<PROBIN.n_rows;g++){
    PROBIN(g)=phetero;
  }

  int MEANS=0;

  int IntIN=0;
  int IntOUT=0;

  //Group Slope Prior Var
  double varbeta1=varbeta;
  //Group Int Prior Var
  double varint1=varint;



  int Group=0;
  //Important For loop integer
  int m =0; //For the inner MCMC looprep
  int i=0; //For the patient index
  int rep=0; //For the simulation repetition.


  //First Fill out Target Vector
  int nDose=Dose.n_rows;


  //Important Integer Quantities
  int B=2000; //Number of iterations for MCMC


  int B1=B/2;

  //Make List objects we'll use

  //Trial Vectors
  arma::vec Y(Nmax);
  arma::vec I(Nmax);
  arma::vec Groups(Nmax);
  arma::vec Doses(Nmax);
  arma::vec Times(Nmax);
  arma::vec ACC(Nmax);
  Y.zeros();
  I.zeros();
  Groups.zeros();
  Doses.zeros();
  Times.zeros();
  ACC.zeros();


  //Innitialize parameters for MCMC
  double mu=meanmu;
  double slope=meanslope;
  double sig=T1;

  double NewMean=0;
  double NewSlope=0;
  int Which1=0;



  //Important quantities we will need in the MCMC
  int  J=MeanInts.n_rows;
  double alpha=0;
  double U=0;
  double signew=0;
  double slopenew=0;
  double Munew=0;

  //For loop integers
  int j=0;
  int k=0;


  //Check to see if Dosestart = 0 otherwise rearrange for c++
  DoseStart=DoseStart-1;


  int StoreInx=0;
  //Vector Of group intercepts
  arma::vec a(J);
  //Vector of group slopes
  arma::vec b=a;
  //Vectors for Group MCMC



  //Vectors for proposals
  arma::vec aprop = a;
  arma::vec bprop=b;





  arma::vec NumA(J);
  arma::vec NumB(J);
  arma::vec IntA(J);
  arma::vec IntB(J);

  double NumMu = 2;
  double NumSlope=2;
  double IntMu=1;
  double IntSlope=1;

  arma::vec avar(J);
  arma::vec bvar(J);

  NumA.zeros();
  NumB.zeros();
  IntA.zeros();
  IntB.zeros();
  avar.zeros();
  bvar.zeros();


  double muvar=1;
  double slopevar=1;
  double trialtime=0;
  NumericVector z9(2);
  NumericVector zprop(5);

  double NumSig=2;
  double IntSig=1;
  double sigvar=1;

  arma::vec etaG(J);
  arma::mat DoseProb(B1,J*nDose);
  arma::mat MeanDose(J,nDose);
  arma::vec sigstore(B1);
  arma::vec mustore(B1);
  arma::vec slopestore(B1);
  arma::mat astore(B1,J);
  arma::mat bstore=astore;
  double eta1=0;
  double DEC=0;

  arma::vec GroupMem(J);
  GroupMem.zeros();
  etaG.zeros();
  DoseProb.zeros();
  sigstore.zeros();
  mustore.zeros();
  astore.zeros();
  bstore.zeros();

  arma::vec GroupMemProp=GroupMem;

  arma::vec INVEC(J);
  INVEC.zeros();

  int eta2=0;

  arma::mat MeanVec(J,nDose);
  arma::vec NTox(J);
  arma::vec OptDose(J);
  arma::vec StoppedGroups(J);
  arma::vec SuspendGroups(J);
  arma::vec GLast(J);
  arma::vec nTreated(J);
  double stopped=0;
  arma::vec TrialTimes(nSims);
  arma::mat OptimalDoses(nSims,J);
  arma::mat NTOX = OptimalDoses;
  arma::vec GroupVec(nDose);
  arma::vec TriedGroups(J);

  MeanVec.zeros();
  NTox.zeros();
  OptDose.zeros();
  StoppedGroups.zeros();
  SuspendGroups.zeros();
  GLast.zeros();
  nTreated.zeros();
  TrialTimes.zeros();
  OptimalDoses.zeros();
  NTox.zeros();
  GroupVec.zeros();
  TriedGroups.zeros();

  arma::mat DoseStore(nSims,Nmax);
  arma::mat GroupStore(nSims,Nmax);
  DoseStore.zeros();
  GroupStore.zeros();
  arma::mat IStore=GroupStore;
  arma::mat DoseTried(J,nDose);
  DoseTried.zeros();
  arma::vec nTreated1(J);
  nTreated1.zeros();
  nTreated1.zeros();
  INVEC.zeros();

  arma::vec INVECNEW=INVEC;
  arma::vec Y2;
  arma::vec I2;
  arma::vec Groups2;
  arma::vec Doses2;
  arma::vec GetGroup=StoppedGroups;
  Y2.zeros();
  I2.zeros();
  Groups2.zeros();
  Doses2.zeros();


  //Need these for dose assignment
  arma::vec HOLD1(nDose);
  HOLD1.zeros();
  arma::mat HOLD(J,nDose+1);
  HOLD.zeros();
  arma::vec DOSETRIED1(nDose);
  DOSETRIED1.zeros();
  int count = 0; //For Counting full followup


  //Wrap nreps around
  for(rep=0;rep<nSims;rep++){


    if(rep==0){
      Rprintf("Starting Simulations. Every 50 done will be notified. \n");
    }else{

      if(rep%50==0){

        Rf_PrintValue(wrap(rep));
        Rprintf("Simulations Finished. \n");

      }
    }

    //Wrap i around

    MeanVec.zeros();
    //Reset all trial vectors
    trialtime=0;
    stopped=0;
    OptDose.zeros();
    StoppedGroups.zeros();
    SuspendGroups.zeros();
    TriedGroups.zeros();
    Y.zeros();
    I.zeros();
    Times.zeros();
    ACC.zeros();
    Groups.zeros();
    Doses.zeros();
    GLast.zeros()+10000;
    nTreated.zeros();
    DoseTried.zeros();

    NTox.zeros();


    //Enroll first patient

    i=0;

    //Trial is starting
    trialtime=trialtime+R::rexp(1/Accrue);
    Group=Sample2(groupprob);


    nTreated[Group]=nTreated[Group]+1;

    if(DoseStart[Group]==0){
      nTreated1[Group]=nTreated1[Group]+1;
    }

    Groups[i]=Group;

    Doses[i]=Dose[DoseStart[Group]];

    ACC[i]=trialtime;

    if(Family==0){
      Times(i)=R::rexp(Param1(Group,DoseStart[Group]));
    }



    if(Family==1){
      Times[i]=R::rgamma(Param1(Group,DoseStart[Group]),1/Param2(Group,DoseStart[Group]));
    }

    if(Family==2){
      Times[i]=R::rlnorm(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
    }


    if(Family==3){
      if(as_scalar(arma::randu(1))<Param1(Group,DoseStart[Group])){
        Times[i]=as_scalar(arma::randu(1))*T1;
      }else{
        Times[i]=T1+1;
      }
    }

    if(Family==4){
      Times[i]=R::rweibull(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
    }


    DoseTried(Group,DoseStart(Group))=1;
    TriedGroups[Group]=1;




    for(i=0;i<Nmax;i++){
      //Start the Trial, Let's do this

      DEC=0; //This will be our counter for stopped group

      while(DEC<1){

        StoppedGroups.zeros();

        trialtime=trialtime+R::rexp(1/Accrue);
        Group=Sample2(groupprob);




        //Enrolling in a new subgroup
        if(TriedGroups[Group]==0){

          DEC=1;
          nTreated[Group]=nTreated[Group]+1;

          if(DoseStart[Group]==0){
            nTreated1[Group]=nTreated1[Group]+1;
          }

          Groups[i]=Group;

          Doses[i]=Dose[DoseStart[Group]];

          ACC[i]=trialtime;

          if(Family==0){
            Times(i)=R::rexp(Param1(Group,DoseStart[Group]));
          }



          if(Family==1){
            Times[i]=R::rgamma(Param1(Group,DoseStart[Group]),1/Param2(Group,DoseStart[Group]));
          }

          if(Family==2){
            Times[i]=R::rlnorm(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
          }


          if(Family==3){
            if(as_scalar(arma::randu(1))<Param1(Group,DoseStart[Group])){
              Times[i]=as_scalar(arma::randu(1))*T1;
            }else{
              Times[i]=T1+1;
            }
          }

          if(Family==4){
            Times[i]=R::rweibull(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
          }


          DoseTried(Group,DoseStart(Group))=1;
          TriedGroups[Group]=1;

        }else{
          //Not in a new subgroup, let's do our MCMC




          //Setup our Data
          I.zeros();
          Y.zeros();

          //Setup Y
          for(k=0;k<i;k++){
            Y[k]=min1(min1(Times[k],trialtime-ACC[k]),T1);
            if(Times[k]==Y[k]){
              I[k]=1;
            }
          }



          if(i<=NSep){
            //No Borrowing
            HOLD=  MCMCSIM(  Y,  I,  Doses,  Groups,  T1,  Target,  Upper,  Dose,  meanmu,
                             meanslope,  MeanInts, MeanSlopes,  varint,  varbeta, phetero,  StoppedGroups,  i+1,  0,  B );
            StoppedGroups = HOLD.col(Dose.n_rows);
          }else{



            if(i<=NBorrow){
              //No Borrowing
              HOLD=  MCMCSIM(  Y,  I,  Doses,  Groups,  T1,  Target,  Upper,  Dose,  meanmu,
                               meanslope,  MeanInts, MeanSlopes,  varint,  varbeta, phetero,  StoppedGroups,  i+1,  1,  B );

              StoppedGroups = HOLD.col(Dose.n_rows);


            }else{
              HOLD=  MCMCSIM(  Y,  I,  Doses,  Groups,  T1,  Target,  Upper,  Dose,  meanmu,
                               meanslope,  MeanInts, MeanSlopes,  varint,  varbeta, phetero,  StoppedGroups,  i+1,  2,  B );

              StoppedGroups = HOLD.col(Dose.n_rows);


            }


            //Check if we need to re-run this separately due to complete stopping

            if(sum(StoppedGroups)==J){
              StoppedGroups.zeros();
              HOLD=  MCMCSIM(  Y,  I,  Doses,  Groups,  T1,  Target,  Upper,  Dose,  meanmu,
                               meanslope,  MeanInts, MeanSlopes,  varint,  varbeta, phetero,  StoppedGroups,  i+1,  0,  B );
              StoppedGroups = HOLD.col(Dose.n_rows);
            }

          }








          if(sum(StoppedGroups)==J){
            stopped=1;


            break;
          }


          //Let's check and make sure this patient is NOT in a stopped group
          if(StoppedGroups[Group]==0){
            DEC=1;
          }

        }


        //Ok now let's enroll this patient, Whats the Optimal Dose




      }
      //End of while, we have a patient enrolling
      //From a non-stopped subgroup

      //What dose should we assign this patient to?



      //Determine Optimal Dose
      k=Group;

      //Are we doing FULL follow up or not?
      if(FULLY==1){
        DOSETRIED1=DoseTried.row(k).t();
        for(j=DoseStart(k);j<Dose.n_rows;j++){
          //Need a counter...
          count=0;
          for(m=0;m<i;m++){
            //Check if observation m is from group k
            if(Groups(m)==k){
              //Check if Y(m)==T1 OR I(m)==1 Both CANNOT be true
              if(I(m)==1){
                count++;
              }else{
                if(Y(m)==T1){
                  count++;
                }
              }



            }

          }

          //Reset this value
          DOSETRIED1(j)=count;


        }

      }else{
        //Not FULL followup
        DOSETRIED1=DoseTried.row(k).t();
      }



      //Fill in hold 1
      for(j=0;j<Dose.n_rows;j++){
        HOLD1(j)=HOLD(k,j);
      }

      for(j=DoseStart(k);j<(Dose.n_rows-1);j++){
        if(DOSETRIED1(j)<cohort){
          //We can't escalate past this dose level
          HOLD1(j+1)=-1000;
        }
      }


      //Rewrite HOLD1 based on dosetried vector


      eta1=MinVec(abs(HOLD1-Target(k)));
      j=0;




      GroupVec = abs(HOLD1-Target(k));

      while(GroupVec[j]!=eta1 ){
        j++;
        //Loop Proceeds until the minimum difference is reached
      }

      OptDose[k]=j;




      //Now we have the subgroup specific optimal doses in the vector OptDose


      //Now let's see if the dose is bigger than the largest dose.
      if(OptDose[Group]>=nDose){
        OptDose[Group]=nDose-1;
      }


      //Assign Patient to subgroup, optdose, add accrual time.

      Groups[i]=Group;
      if(OptDose[Group]==0){
        //If we assign one at the lowest dose add this
        nTreated1[Group]=nTreated1[Group]+1;
      }
      nTreated[Group]=nTreated[Group]+1;


      Doses[i]=Dose[OptDose[Group]];


      ACC[i]=trialtime;
      //Based on the Group and Optimal dose, randomly generate the next patient Data


      if(Family==0){
        Times(i)=R::rexp(Param1(Group,OptDose[Group]));
      }



      if(Family==1){
        Times[i]=R::rgamma(Param1(Group,OptDose[Group]),Param2(Group,OptDose[Group]));
      }

      if(Family==2){
        Times[i]=R::rlnorm(Param1(Group,OptDose[Group]),Param2(Group,OptDose[Group]));
      }


      if(Family==3){
        if(as_scalar(arma::randu(1))<Param1(Group,OptDose[Group])){
          Times[i]=as_scalar(arma::randu(1))*T1;
        }else{
          Times[i]=T1+1;
        }
      }

      if(Family==4){
        Times[i]=R::rweibull(Param1(Group,OptDose[Group]),Param2(Group,OptDose[Group]));
      }


      DoseTried(Group,OptDose(Group))=DoseTried(Group,OptDose(Group))+1;








      //End of the i loop, now we have ALL the patients
    }




    //Now let's make the final decision.
    if(stopped==1){
      //The trial Ended Early, i-1 is the index of the last patient enrolled.


      for(k=0;k<i;k++){
        for(j=0;j<J;j++){
          if(Groups[k]==j){
            NTox[j]=NTox[j]+I[k];
          }
        }
      }


      //All Optimal Doses are -1
      for(k=0;k<J;k++){
        OptDose[k]=-1;
      }


    }else{
      //The trial did NOT end early. Follow patients out until end of trial and get the optimal doses
      trialtime=trialtime+T1;



      StoppedGroups.zeros();

      Y.zeros();
      I.zeros();

      //Setup Y
      for(k=0;k<Nmax;k++){

        Y[k]=min1(Times[k],T1);


        if(Times[k]<T1){
          I[k]=1;
        }
      }



      if(i<=NSep){
        //No Borrowing
        HOLD=  MCMCSIM(  Y,  I,  Doses,  Groups,  T1,  Target,  Upper,  Dose,  meanmu,
                         meanslope,  MeanInts, MeanSlopes,  varint,  varbeta, phetero,  StoppedGroups,  Nmax,  0,  B );
        StoppedGroups = HOLD.col(Dose.n_rows);
      }else{



        if(i<=NBorrow){
          //No Borrowing
          HOLD=  MCMCSIM(  Y,  I,  Doses,  Groups,  T1,  Target,  Upper,  Dose,  meanmu,
                           meanslope,  MeanInts, MeanSlopes,  varint,  varbeta, phetero,  StoppedGroups,  Nmax,  1,  B );

          StoppedGroups = HOLD.col(Dose.n_rows);


        }else{
          HOLD=  MCMCSIM(  Y,  I,  Doses,  Groups,  T1,  Target,  Upper,  Dose,  meanmu,
                           meanslope,  MeanInts, MeanSlopes,  varint,  varbeta, phetero, StoppedGroups,  Nmax,  2,  B );

          StoppedGroups = HOLD.col(Dose.n_rows);


        }


        //Check if we need to re-run this separately due to complete stopping

        if(sum(StoppedGroups)==J){
          StoppedGroups.zeros();
          HOLD=  MCMCSIM(  Y,  I,  Doses,  Groups,  T1,  Target,  Upper,  Dose,  meanmu,
                           meanslope,  MeanInts, MeanSlopes,  varint,  varbeta, phetero,  StoppedGroups, Nmax,  0,  B );
          StoppedGroups = HOLD.col(Dose.n_rows);
        }

      }






      if(sum(StoppedGroups)==J){
        stopped=1;
        //We should declare NO DOSE OPTIMAL

      }












      if(stopped==1){
        //No optimal doses at all, all groups are not good.
        OptDose.zeros();
        OptDose=OptDose-1;
      }else{

        for(k=0;k<J; k++){
          if(StoppedGroups(k)==0){

            //Figure out the dose here
            DOSETRIED1=DoseTried.row(k).t();

            //We CANNOT escalate at all here

            //Fill in hold 1
            for(j=0;j<Dose.n_rows;j++){
              HOLD1(j)=HOLD(k,j);
            }

            for(j=DoseStart(k);j<Dose.n_rows;j++){
              if(DOSETRIED1(j)<cohort){
                //We can't escalate to this dose level
                HOLD1(j)=-1000;
              }
            }


            //Rewrite HOLD1 based on dosetried vector


            eta1=MinVec(abs(HOLD1-Target(k)));
            j=0;




            GroupVec = abs(HOLD1-Target(k));

            while(GroupVec[j]!=eta1 ){
              j++;
              //Loop Proceeds until the minimum difference is reached
            }

            OptDose[k]=j;




            //Now we have the subgroup specific optimal doses in the vector OptDose



            //Now let's see if the dose is bigger than the largest dose.
            if(OptDose[k]>=nDose){
              OptDose[k]=nDose-1;
            }




          }else{
            //This subgroup is too toxic
            OptDose(k)=-1;

          }



        }


      }


      //Now we have our optimal dose-vector for this trial.
      //Let's get other needed OC quantities.

      //Num Toxicities
      for(k=0;k<Nmax;k++){
        for(j=0;j<J;j++){
          if(Groups[k]==j){
            NTox[j]=NTox[j]+I[k];
          }
        }
      }


    }

    OptDose=OptDose+1;


    TrialTimes(rep)=trialtime;
    for(j=0;j<J;j++){
      OptimalDoses(rep,j)=OptDose(j);
      NTOX(rep,j)=NTox(j);
    }




    //Fill in Doses and Groups Treated
    for(j=0;j<Nmax;j++){
      DoseStore(rep,j)=Doses(j);
      GroupStore(rep,j)=Groups(j);
      IStore(rep,j)=I(j);
    }



    //End Trial simulation rep
  }



  List z1 = List::create(OptimalDoses,NTOX,TrialTimes,IStore,DoseStore,GroupStore);





  return(z1);




}



