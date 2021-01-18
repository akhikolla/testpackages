#include <time.h>
#include <Rmath.h>
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <iostream>     // std::cout
#include <cmath>
#include <cfloat>

// [[Rcpp::depends(RcppArmadillo)]]



using namespace Rcpp;




double abs1(double a){
  if(a<0){
    return(-a);
  }else{
    return(a);
  }


}


//Returns the maximum from a vector
double MaxVec(arma::vec Y){
  int J1=Y.n_rows;
  int j;
  double max=Y[0];
  for(j=1;j<J1;j++){
    if(Y[j]>max){
      max=Y[j];
    }
  }

  return(max);

}


double logINT(int z){
  double z1=z;
  return(log(z1));
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





double Like(arma::vec YE1, //Binary indicators of efficacy
            arma::vec YT1, //Binary indicators of Toxicity
            arma::vec Doses1, //Vector of Doses given to patients
            arma::vec beta, //Parameter vector containing linear parameters, length 5
            double psi, //  association between toxicity and efficacy
            int Num
){

  double LogL =0;


  //Make YE and YT and Doses
  arma::vec YE(Num);
  arma::vec YT(Num);
  arma::vec Dose(Num);

  //Loop over the doses for the likelihood
  //Doses
  int j=0;
  //Patients
  int i=0;
  //4 EFF-TOX outcomes
  int m=0;

  for(j=0;j<Num;j++){
    YE[j]=YE1[j];
    YT[j]=YT1[j];
    Dose[j]=Doses1[j];

  }


  //Holds tox and eff prob for a given dose
  double piE=0;
  double piT=0;
  //Stores probabilities of the 4 outcomes
  // 1 Neither
  // 2 Eff - No Tox
  // 3 Tox - No Eff
  // 4 Both
  arma::vec probs(4);

  //Association Constand
  double Con = (exp(psi)-1)/(exp(psi)+1);
  double Con1=0;

  for(i=0;i<YE.n_rows;i++){
    //Compute PiE and PiT
    piE = beta[0]+beta[1]*Dose[i]+beta[2]*pow(Dose[i],2);
    piT = beta[3]+beta[4]*Dose[i];
    piE = exp(piE);
    piT=exp(piT);
    piE = piE/(1+piE);
    piT = piT/(1+piT);
    //Now we have the probabilities of toxictity and efficacy for this dose.
    Con1 = piE*(1-piE)*piT*(1-piT)*Con;

    probs[0]=log(Con1 + (1-piE)*(1-piT));
    probs[1]=log((1-piT)*piE - Con1);
    probs[2]= log((1-piE)*piT - Con1);
    probs[3]= log(piE*piT + Con1);

    //Loop over each patient
    if( (YE[i]==0) && (YT[i]==0)){
      LogL = LogL + probs[0];
    }
    if( (YE[i]==1) && (YT[i]==0)){
      LogL = LogL + probs[1];
    }
    if( (YE[i]==0) && (YT[i]==1) ){
      LogL = LogL + probs[2];
    }
    if( (YE[i]==1) && (YT[i]==1)){
      LogL = LogL + probs[3];
    }















  }





  return(LogL);



}


//Function for Likelihood for survival model.
double Like2( arma::vec Y, //Survival Times
              arma::vec I, //Censoring Indicators
              arma::vec YE, //Indicators of short term patient Efficacy
              arma::vec YT, //Indicators of short term patient Toxicity
              arma::vec Doses,  //Vector of Standardized Dose values given to patients
              arma::vec beta, //Parameter vector for (x_j, Y_e, Y_T) in baseline hazard
              arma::vec s, //Vector containing the split point locations
              arma::vec lam, //Vector containing the log hazard heights
              int J //Numer of split points
){




  arma::vec eta = beta[0]*Doses-exp(beta[1])*YE+exp(beta[2])*YT + beta[3]*pow(Doses,2);



  int m=0;
  int l=0;
  double LogL=0;

  //  Rf_PrintValue(wrap(eta));

  for(m=0;m<Y.n_rows;m++){
    if(I[m]==1){
      LogL = LogL + eta[m];
    }

  }



  //Cycle through each interval and obtain survival estimates and add hazard heights
  for(l=0; l<(J+1);l++){
    for(m=0; m<Y.n_rows;m++){
      LogL = LogL - max1(0,min1(s(l+1),Y(m))-s(l))*exp(eta[m]+lam[l]);


      if( (Y(m)>s[l]) && (Y(m)<=s[l+1]) && (I[m]==1)){
        LogL = LogL + lam[l];
      }




    }

  }



  return(LogL);



}




double GetDesire(double PE,
                 double PT,
                 arma::vec Contour){


  double p=0;
  double pu=0;
  double pl=0;
  double tol=.005;


  arma::vec z9(2);

  //Now lets Solve for p to get the desirability
  while(abs1(pow((Contour[2]-1)/(Contour[0]-1),p) + pow(Contour[3]/Contour[1],p)-1) >tol){
    pu = p+.005;
    pl=p-.005;

    if( abs1(pow((Contour[2]-1)/(Contour[0]-1),pu) + pow((Contour[3])/(Contour[1]),pu)-1)<
      abs1(pow((Contour[2]-1)/(Contour[0]-1),pl) + pow((Contour[3])/(Contour[1]),pl)-1)){
      p=pu;
    }else{
      p=pl;
    }




  }



  double desire = 1-pow((pow(((PE-1)/(Contour[0]-1)),p)+pow((PT/Contour[1]),p)),1/p);



  return(desire);


}







double SampBern(double prob){
  double U=as_scalar(arma::randu(1));


  if(U<prob){
    return(1);
  }else{
    return(0);
  }


}


double GetMin(arma::vec DESIRE){
  double min=DESIRE[0];
  int j=0;
  for(j=1;j<DESIRE.n_rows;j++){
    if(DESIRE[j]>-200){
      if(DESIRE[j]<min){
        min=DESIRE[j];
      }
    }
  }
  return(min);



}




int GetDose(arma::vec DESIRE){
  //First off, we don't want ANY desire that's negative.
  int NumD = sum(DESIRE>0);
  arma::vec DESIRE1(NumD);
  arma::vec WHICH1(NumD);

  int j=0;
  int k=0;

  for(j=0;j<NumD;j++){
    while(DESIRE[k]<0){
      k++;
    }
    DESIRE1[j]=DESIRE[k];
    WHICH1[j]=k;
    k=k+1;


  }

  DESIRE1=DESIRE1/sum(DESIRE1);
  arma::vec cumprob=DESIRE1;


  for(j=1;j<NumD;j++){
    cumprob[j]=cumprob[j]+cumprob[j-1];
  }



  double U=as_scalar(arma::randu(1));



  int which1=0;

  if(U>cumprob[0]){
    //Not the lowest acceptable dose considered
    while(U>cumprob[which1]){
      which1++;
    }


  }




  return(WHICH1[which1]);


}


double GetSd(arma::vec X){
  double Mean = sum(X)/X.n_rows;

  return(pow(sum(pow(X-Mean,2))/(X.n_rows-1),.5));


}




int GetDose1(arma::vec DESIRE){
  //First off, we don't want ANY desire that's negative.
  int j=0;
  int k=0;
  arma::vec z9(2);


  int count=0;

  for(j=0;j<DESIRE.n_rows;j++){
    if(DESIRE[j]>-200){
      count++;
    }
  }
  int NumD=count;





  if(count>1){


    arma::vec DESIRE1(NumD);
    arma::vec WHICH1(NumD);



    for(j=0;j<NumD;j++){
      while(DESIRE[k]< (-200)){
        k++;
      }
      DESIRE1[j]=DESIRE[k];
      WHICH1[j]=k;
      k=k+1;


    }




    DESIRE1=(DESIRE1-sum(DESIRE1)/DESIRE1.n_rows)/GetSd(DESIRE1);


    DESIRE1=exp(DESIRE1)/sum(exp(DESIRE1));


    arma::vec cumprob=DESIRE1;



    for(j=1;j<NumD;j++){
      cumprob[j]=cumprob[j]+cumprob[j-1];
    }



    double U=as_scalar(arma::randu(1));




    int which1=0;

    if(U>cumprob[0]){
      //Not the lowest acceptable dose considered
      while(U>cumprob[which1]){
        which1++;

        if(which1==(cumprob.n_rows-1)){
          break;
        }

      }


    }



    return(WHICH1[which1]);

  }else{



    k=0;
    while(DESIRE[k]< (-200)){
      k++;
    }






    return(k);
  }

}





//' Simulates repitions of an Adaptive Eff-Tox Trial.
//'
//' This function simulates repititions of an adaptive Eff-Tox Trial and returns a list containing the optimal dose chosen
//' @param DoseStart Dose to start enrolling cohorts of patients at.
//' @param Dose Vector containing the standardized doses considered.
//' @param Hypermeans Vector containing prior hypermeans of length 6 for Eff-Tox parameters.
//' @param Hypervars Vector containing prior hypervariances of length 6 for Eff-Tox parameters.
//' @param Contour Vector containing 4 entries used to make the desireability function. Contour[1] contains a desired toxicity probability given efficacy, Countour[2] contains a desired efficacy probability given toxicity, and (Contour[3],Contour[4]) is an equally desireable pair of efficacy and toxicity probabilities that are non-zero or one.
//' @param PiLim Vector of length two with PiLim[1] containing the acceptable lower limit on efficacy probability and PiLim[2] containing the acceptable upper limit on toxicity probability.
//' @param ProbLim Vector of length two with ProbLim[1] containing the probability cutoff for acceptable efficacy probability and ProbLim[2] containing the probability cutoff for acceptable toxicity probability.
//' @param B Number of iterations to perform in the MCMC.
//' @param cohort Size of each patient cohort.
//' @param NET Maximum sample size for phase I/II.
//' @param NF Number of patients to assign optimal doses prior to adaptive randomization.
//' @param nSims Number of simulated trials to run.
//' @param PETrue True vector of efficacy probabilities for each dose.
//' @param PTTrue True vector of toxicity probabilities for each dose.
//' @return List containing the vector of optimal doses chosen, a matrix of posterior desireability scores for each trial, and a matrix consisting of patient dose assignments and Toxicity and Efficacy indicators, with each Nmax rows corresponding to a separate trial. Trials that are stopped due to excessive toxicity probabilty or small efficacy probabilities are not included in the final results.
//' @examples
//' ##Doses, YE,YT
//' ##Starting Dose
//' DoseStart=1
//'##Vector of Numerical Doses
//'Dose = c(1,2,3,3.5,5)
//'Dose=(Dose-mean(Dose))/sd(Dose)
//'## Contour Vector
//'Contour = c(.35, .75,.7,.4)
//'##Hypermeans
//'Hypermeans = c(.022,3.45,0,-4.23,3.1,0)
//'Hypervars = c(2.6761, 2.6852, .2, 3.1304, 3.1165, 1)
//'Hypervars=Hypervars^2
//'##Acceptability Criteria
//'PiLim = c(.3,.4)
//'ProbLim=c(.1,.1)
//'##Cohort Size, N^F and N_ET
//'cohort=3
//'NF=15
//'NET=30
//'PTTrue = c(.1,.15,.25,.35,.5)
//'PETrue=c(.2,.4,.6,.65,.7)
//'##Number of iterations for MCMC
//'B=2000
//'### Number of Simulations
//'nSims=1
//'RunAdaptiveEffToxTrial(DoseStart,Dose, Hypermeans,  Hypervars,
//'Contour, PiLim, ProbLim,  cohort, NET,  NF, B, nSims, PETrue, PTTrue )
//'@export
//[[Rcpp::export]]
List RunAdaptiveEffToxTrial( int DoseStart, //Dose level to START the EFF-Tox Trial
                             arma::vec Dose, //Vector of Doses considered in the trial
                             arma::vec Hypermeans, //6 vector of prior means
                             arma::vec Hypervars, //6 vector of prior standard deviations
                             arma::vec Contour, //4 vector of Contour parameter
                             arma::vec PiLim, //2 vector of acceptable limits
                             arma::vec ProbLim, //2 vector of cutoff for acceptabilities
                             int cohort, //Cohort size
                             int NET, //Maximum Sample size to run the EFFtox Trial, must be divisible by cohort
                             int NF, //Minimum sample size to begin adaptive randomization, must be divisible by count
                             int B, // Number of reps to perform in MCMC
                             int nSims, //Number of Simulated trials to run.
                             arma::vec PETrue, //True vector of efficacy probabilities for each dose
                             arma::vec PTTrue ){
int Nmax=NET;
int Nmin=NF;
arma::vec Accept=PiLim;
arma::vec Lower=ProbLim;
  arma::vec MeanUt(Dose.n_rows);
  arma::vec PiEff(Dose.n_rows);
  arma::vec PiTox(Dose.n_rows);
  arma::vec ACCEPTTOX(Dose.n_rows);
  arma::vec ACCEPTEFF(Dose.n_rows);
  arma::vec DESIRE(Dose.n_rows);
  arma::vec ACCEPT(Dose.n_rows);
  ACCEPTTOX.zeros();
  ACCEPTEFF.zeros();

  double p=0;
  double pu=0;
  double pl=0;
  double tol=.005;



  //Now lets Solve for p to get the desirability
  while(abs1(pow((Contour[2]-1)/(Contour[0]-1),p) + pow(Contour[3]/Contour[1],p)-1) >tol){
    pu = p+.005;
    pl=p-.005;

    if( abs1(pow((Contour[2]-1)/(Contour[0]-1),pu) + pow((Contour[3])/(Contour[1]),pu)-1)<
      abs1(pow((Contour[2]-1)/(Contour[0]-1),pl) + pow((Contour[3])/(Contour[1]),pl)-1)){
      p=pu;
    }else{
      p=pl;
    }




  }





  int NCohort = Nmax/cohort;
  int MinCohort=Nmin/cohort;

  //What Doses are Acceptable

  arma::vec ACCTOX(Dose.n_rows);
  arma::vec ACCEFF=ACCTOX;


  //Storage for each trial outputs
  arma::vec Doses(Nmax);
  arma::vec YE(Nmax);
  arma::vec YT(Nmax);


  //Storage matrices for trial output
  arma::vec CHOSENDOSE(nSims);
  arma::mat UTStore(nSims,Dose.n_rows);
  arma::mat TrialStore(nSims*Nmax,3);


  int Count=0;




  int B1 = B/2;

  //Used for for loop
  int b=0;
  //Used in Sim look
  int nREP=0;

  //Storage Matrices
  arma::mat piEStore(B1,Dose.n_rows);
  arma::mat piTStore(B1,Dose.n_rows);
  arma::mat AcceptTox(B1,Dose.n_rows);
  arma::mat AcceptEff(B1,Dose.n_rows);
  AcceptTox.zeros();
  AcceptEff.zeros();
  arma::vec DESIRE2(Dose.n_rows);
  DESIRE2.zeros();


  arma::mat ParamStore(B1,6);
  //Setup innitial
  arma::vec beta(5);
  beta.zeros();
  arma::vec betaprop=beta;
  double psi=0;
  double psiprop=0;
  //Quantities for MCMC
  double U=0;
  double alpha=0;
  int m=0;
  int j=0;
  int StoreInx=0;
  double PE=0;
  double PT=0;

  arma::vec bvar(6);
  bvar.zeros();
  bvar=bvar+1;
  arma::vec Intb = bvar;
  arma::vec Numb = bvar + 1;

  arma::vec Doses2=Doses;

  NumericVector z9(2);
  int CoNum=0;

  int STOPPED=0;
  int i=0;
  int OptDose=0;
  double min1=0;
  int DoseStart1=DoseStart;



  int nRep=0;

  arma::vec DESIRE1=DESIRE;

  arma::vec z5(2);
  arma::vec DoseTried(Dose.n_rows);




  while(nRep<nSims){

    STOPPED=0;
    if(nRep>0){
      if(nRep%100==0){



        Rprintf("Iteration:");
        z9[0]=nRep;
        z9[1]=nRep;
        Rf_PrintValue(wrap(z9));


      }
    }
    Doses2.zeros();
    //Wrap Trial right here

    Doses.zeros();
    DoseTried.zeros();

    for(m=0;m<DoseStart;m++){
    DoseTried[m]=1;
  }

    DoseStart1=DoseStart-1;
  OptDose=DoseStart1;


    //Wrap around each cohort
    for(i=0;i<NCohort;i++){

      //Fill in next cohort with the data size
      for(b=0;b<cohort;b++){
        Doses[i*cohort+b] = Dose[OptDose];
        Doses2[i*cohort+b] = OptDose;

        YE[i*cohort+b]=SampBern(PETrue[OptDose]);
        YT[i*cohort+b]=SampBern(PTTrue[OptDose]);
        DoseTried(OptDose)=1;

      }

      CoNum=(i+1)*cohort;

      Intb=Intb.zeros()+1;
      Numb=Numb.zeros()+2;
      bvar=bvar.zeros()+1;

      for(j=0;j<5;j++){
        beta[j]=Hypermeans[j];
      }
      psi=Hypermeans[5];

      //Perform MCMC
      for(b=0;b<B;b++){


        if(b<(B/2 + 2)){
          if(b%250==0){

            for(j=0;j<6;j++){
              if((Intb[j]/Numb[j])>.5){
                bvar[j]=bvar[j]*2;
              }

              if((Intb[j]/Numb[j])<.2){
                bvar[j]=bvar[j]/2;
              }


              Intb[j]=1;
              Numb[j]=2;

            }





          }
        }




        //All but beta_T can be less than 0, sample those
        for(j=0;j<4;j++){

          betaprop=beta;
          betaprop[j]=as_scalar(arma::randn(1))*bvar[j] + betaprop[j];

          alpha =   Like(YE,YT,Doses,betaprop,psi,CoNum) - Like(YE,YT,Doses, beta, psi,CoNum);
          alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

          U=log(as_scalar(arma::randu(1)));

          if(U<alpha){
            beta[j]=betaprop[j];
            Intb[j]=Intb[j]+1;
          }

          Numb[j] = Numb[j]+1;





        }


        //toxicity slope
        j=4;
        betaprop=beta;
        betaprop[j]=exp(as_scalar(arma::randn(1))*bvar[j] + log(betaprop[j]));

        alpha =   Like(YE,YT,Doses,betaprop,psi,CoNum) - Like(YE,YT,Doses, beta, psi,CoNum);
        alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

        alpha = alpha + log(betaprop[j]) - log(beta[j]);


        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          beta[j]=betaprop[j];
          Intb[j]=Intb[j]+1;
        }

        Numb[j] = Numb[j]+1;




        //Get new Psi

        //toxicity slope
        j=5;

        psiprop=as_scalar(arma::randn(1))*bvar[j] + psi;

        alpha =   Like(YE,YT,Doses,beta,psiprop,CoNum) - Like(YE,YT,Doses, beta, psi,CoNum);
        alpha = alpha - .5*pow(psiprop-Hypermeans[j],2)/Hypervars[j] + .5*pow(psi-Hypermeans[j],2)/Hypervars[j];

        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          psi=psiprop;
          Intb[j]=Intb[j]+1;
        }

        Numb[j]=Numb[j]+1;









        if(b>(B-B1-1)){

          StoreInx = b-B1;

          for(j=0;j<5;j++){
            ParamStore(StoreInx,j)=beta[j];
          }
          ParamStore(StoreInx,5)=psi;

          //Fill in Toxcity and Efficacy probabilities


          for(j=0;j<Dose.n_rows;j++){
            PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
            PT = beta[3]+beta[4]*Dose[j];

            PE = exp(PE);
            PT=exp(PT);
            PE=PE/(1+PE);
            PT = PT/(1+PT);
            piEStore(StoreInx,j)=PE;
            piTStore(StoreInx,j)=PT;


            if(PE>Accept[0]){
              AcceptEff(StoreInx,j)=1;
            }else{
              AcceptEff(StoreInx,j)=0;

            }

            if(PT<Accept[1]){
              AcceptTox(StoreInx,j)=1;
            }else{
              AcceptTox(StoreInx,j)=0;

            }






          }



        }



      }

      //Compute posterior mean



      for(j=0;j<Dose.n_rows;j++){
        ACCEPTEFF[j]=sum(AcceptEff.col(j))>(Lower(0)*B1)  ;
        ACCEPTTOX[j]=sum(AcceptTox.col(j))>(Lower(1)*B1) ;
      }




      //If ==1, for both EFF and TOX then dose j is acceptable

      //Determine Optimal Dose

      //Get Posterior mean of parameters
      for(j=0;j<5;j++){
        beta[j]=sum( ParamStore.col(j))/ParamStore.n_rows;
      }


      for(j=0;j<Dose.n_rows;j++){
        PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
        PT = beta[3]+beta[4]*Dose[j];

        PE = exp(PE);
        PT=exp(PT);
        PE=PE/(1+PE);
        PT = PT/(1+PT);
        PiEff(j)=PE;
        PiTox(j)=PT;

      }






      //Now get the mean desireability for each dose if they are acceptable
      for(j=0;j<Dose.n_rows;j++){
        if((ACCEPTEFF[j]+ACCEPTTOX[j])==2){
          //Dose is ACCEPTABLE
          DESIRE[j]=1-pow(pow(((PiEff(j)-1)/(Contour[0]-1)),p)+pow((PiTox(j)/Contour[1]),p),1/p);
          ACCEPT[j]=1;

        }else{
          ACCEPT[j]=0;
          DESIRE[j]=-1000;
        }



      }



      //Check if we should stop the trial

      if(sum(ACCEPT)>0){


        //Now pick the dose





        if((i<MinCohort) || (sum(ACCEPT)==1)){
          OptDose=0;

          for(j=1;j<Dose.n_rows;j++){
            if(DESIRE[j]>DESIRE[j-1]){
              OptDose=j;
            }
          }

          //Is this dose bigger than the dose higher than what's been tried?
          if(DoseTried(OptDose)==0){
            j=0;

            while(DoseTried(j)==1){
              j++;
            }

            OptDose=j;

          }

        }else{
          //Ok let's adaptively randomize patients

          for(j=0;j<Dose.n_rows;j++){
            if(DoseTried(j)==0){
              ACCEPT(j)=0;
            }
          }


          for(j=1;j<Dose.n_rows;j++){
            if((DoseTried(j)==0) && (DoseTried(j-1)==0)){
              ACCEPT(j)=0;
            }
          }



          if(sum(ACCEPT)==0){
            //Do regular


            OptDose=0;

            for(j=1;j<Dose.n_rows;j++){
              if(DESIRE[j]>DESIRE[j-1]){
                OptDose=j;
              }
            }

            //Is this dose bigger than the dose higher than what's been tried?
            if(DoseTried(OptDose)==0){
              j=0;

              while(DoseTried(j)==1){
                j++;
              }

              OptDose=j;

            }






          }else{





            OptDose=  GetDose1(DESIRE);

            //Is this dose bigger than the dose higher than what's been tried?
            if(DoseTried(OptDose)==0){
              j=0;

              while(DoseTried(j)==1){
                j++;

              }

              OptDose=j;

            }








          }



        }









      }else{
        //Trial Stopped
        OptDose=-1;
        STOPPED=1;
        break;

      }





    }

    if(STOPPED==0){


      Intb=Intb.zeros()+1;
      Numb=Numb.zeros()+2;
      bvar=bvar.zeros()+1;


      for(j=0;j<5;j++){
        beta[j]=Hypermeans[j];
      }
      psi=Hypermeans[5];

      //Perform MCMC
      for(b=0;b<B;b++){


        if(b<(B/2 + 2)){
          if(b%250==0){

            for(j=0;j<6;j++){
              if((Intb[j]/Numb[j])>.5){
                bvar[j]=bvar[j]*2;
              }

              if((Intb[j]/Numb[j])<.2){
                bvar[j]=bvar[j]/2;
              }


              Intb[j]=1;
              Numb[j]=2;

            }





          }
        }




        //All but beta_T can be less than 0, sample those
        for(j=0;j<4;j++){

          betaprop=beta;
          betaprop[j]=as_scalar(arma::randn(1))*bvar[j] + betaprop[j];

          alpha =   Like(YE,YT,Doses,betaprop,psi,Nmax) - Like(YE,YT,Doses, beta, psi,Nmax);
          alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

          U=log(as_scalar(arma::randu(1)));

          if(U<alpha){
            beta[j]=betaprop[j];
            Intb[j]=Intb[j]+1;
          }

          Numb[j] = Numb[j]+1;





        }


        //toxicity slope
        j=4;
        betaprop=beta;
        betaprop[j]=exp(as_scalar(arma::randn(1))*bvar[j] + log(betaprop[j]));

        alpha =   Like(YE,YT,Doses,betaprop,psi,Nmax) - Like(YE,YT,Doses, beta, psi,Nmax);
        alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

        alpha = alpha + log(betaprop[j]) - log(beta[j]);


        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          beta[j]=betaprop[j];
          Intb[j]=Intb[j]+1;
        }

        Numb[j] = Numb[j]+1;




        //Get new Psi

        //toxicity slope
        j=5;

        psiprop=as_scalar(arma::randn(1))*bvar[j] + psi;

        alpha =   Like(YE,YT,Doses,beta,psiprop,Nmax) - Like(YE,YT,Doses, beta, psi,Nmax);
        alpha = alpha - .5*pow(psiprop-Hypermeans[j],2)/Hypervars[j] + .5*pow(psi-Hypermeans[j],2)/Hypervars[j];

        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          psi=psiprop;
          Intb[j]=Intb[j]+1;
        }

        Numb[j]=Numb[j]+1;









        if(b>(B-B1-1)){

          StoreInx = b-B1;

          for(j=0;j<5;j++){
            ParamStore(StoreInx,j)=beta[j];
          }
          ParamStore(StoreInx,5)=psi;

          //Fill in Toxcity and Efficacy probabilities


          for(j=0;j<Dose.n_rows;j++){
            PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
            PT = beta[3]+beta[4]*Dose[j];

            PE = exp(PE);
            PT=exp(PT);
            PE=PE/(1+PE);
            PT = PT/(1+PT);
            piEStore(StoreInx,j)=PE;
            piTStore(StoreInx,j)=PT;


            if(PE>Accept[0]){
              AcceptEff(StoreInx,j)=1;
            }else{
              AcceptEff(StoreInx,j)=0;

            }

            if(PT<Accept[1]){
              AcceptTox(StoreInx,j)=1;
            }else{
              AcceptTox(StoreInx,j)=0;

            }







          }



        }



      }




      for(j=0;j<Dose.n_rows;j++){
        ACCEPTEFF[j]=sum(AcceptEff.col(j))>(Lower(0)*B1)  ;
        ACCEPTTOX[j]=sum(AcceptTox.col(j))>(Lower(1)*B1) ;
      }





      //If ==1, for both EFF and TOX then dose j is acceptable

      //Determine Optimal Dose

      //Get Posterior mean of parameters
      for(j=0;j<5;j++){
        beta[j]=sum( ParamStore.col(j))/ParamStore.n_rows;
      }


      for(j=0;j<Dose.n_rows;j++){
        PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
        PT = beta[3]+beta[4]*Dose[j];

        PE = exp(PE);
        PT=exp(PT);
        PE=PE/(1+PE);
        PT = PT/(1+PT);
        PiEff(j)=PE;
        PiTox(j)=PT;

      }








      //Now get the mean desireability for each dose if they are acceptable
      for(j=0;j<Dose.n_rows;j++){
        if((ACCEPTEFF[j]+ACCEPTTOX[j])==2){
          //Dose is ACCEPTABLE
          DESIRE[j]=1-pow(pow(((PiEff(j)-1)/(Contour[0]-1)),p)+pow((PiTox(j)/Contour[1]),p),1/p);
          ACCEPT[j]=1;

        }else{
          ACCEPT[j]=0;
          DESIRE[j]=-1000;
        }



      }





      if(sum(ACCEPT)>0){
        OptDose=0;

        for(j=1;j<Dose.n_rows;j++){
          if(DESIRE[j]>DESIRE[j-1]){
            OptDose=j;
          }
        }




        //Is this dose bigger than the dose higher than what's been tried?
        if(DoseTried(OptDose)==0){
          j=0;

          while(DoseTried(j)==1){
            j++;
          }

          OptDose=j-1;

        }




        CHOSENDOSE[nRep]=OptDose+1;
        for(j=0;j<Dose.n_rows;j++){
          UTStore(nRep,j)= DESIRE[j];
        }


        for(j=0;j<Nmax;j++){
          TrialStore(nRep*Nmax+j,0)=Doses2(j);
          TrialStore(nRep*Nmax+j,1)=YE(j);
          TrialStore(nRep*Nmax+j,2)=YT(j);



        }


        nRep++;



      }else{
   //     Rprintf("Skipped Simulation      ");
      }


    }else{
  //    Rprintf("Skipped Simulation       ");
    }






  }


  List z1 = List::create(CHOSENDOSE,UTStore,TrialStore,ParamStore);

  return(z1);
}











//' Determines the optimal dose to assign the next patient cohort.
//'
//' This function returns the optimal acceptable dose number to assign the next patient cohort or stops the trial if no dose is deemed acceptable.
//' @param YE   Vector containing observed efficacy indicators.
//' @param YT   Vector containing observed toxicity indicators.
//' @param Doses Vector containing numbered Doses of patients in trial.
//' @param Dose Vector containing the standardized doses considered.
//' @param DosesTried Binary vector corresponding to which doses have been tried.
//' @param Hypermeans Vector containing prior hypermeans of length 6 for Eff-Tox parameters.
//' @param Hypervars Vector containing prior hypervariances of length 6 for Eff-Tox parameters.
//' @param Contour Vector containing 4 entries used to make the desireability function. Contour[1] contains a desired toxicity probability given efficacy, Countour[2] contains a desired efficacy probability given toxicity, and (Contour[3],Contour[4]) is an equally desireable pair of efficacy and toxicity probabilities that are non-zero or one.
//' @param PiLim Vector of length two with PiLim[1] containing the acceptable lower limit on efficacy probability and PiLim[2] containing the acceptable upper limit on toxicity probability.
//' @param ProbLim Vector of length two with ProbLim[1] containing the probability cutoff for acceptable efficacy probability and ProbLim[2] containing the probability cutoff for acceptable toxicity probability.
//' @param B Number of iterations to perform in the MCMC.
//' @return The optimal dose level to administer the next patient cohort.
//'@examples
//'##Doses, YE,YT
//'Doses= c(1,1,1,2,2,2,1,1,1,3,3,3,1,1,1,2,2,2)
//'YE = c(0,0,1,1,1,0,0,0,0,1,1,1,0,0,1,1,1,0)
//'YT=c(0,0,0,1,1,0,1,0,0,1,1,1,0,0,0,1,0,0)
//'##Vector of Numerical Doses
//'Dose = c(1,2,3,3.5,5)
//'Dose=(Dose-mean(Dose))/sd(Dose)
//'##Five doses, but only 3 tried so we have
//'DosesTried=c(1,1,1,0,0)
//'## Contour Vector
//'Contour = c(.35, .75,.7,.4)
//'##Hypermeans
//'Hypermeans = c(.022,3.45,0,-4.23,3.1,0)
//'Hypervars = c(2.6761, 2.6852, .2, 3.1304, 3.1165, 1)
//'Hypervars=Hypervars^2
//'##Acceptability Criteria
//'PiLim = c(.3,.4)
//'ProbLim=c(.1,.1)
//'##Number of iterations
//'B=2000
//'AssignEffTox(YE,YT, Doses, Dose, DosesTried, Hypermeans,  Hypervars, Contour, PiLim,  ProbLim, B )
//'@export
//[[Rcpp::export]]
int AssignEffTox( arma::vec YE, //Observed Efficacy Indicator Vector
                  arma::vec YT, //Observed toxicity indicator vector.
       arma::vec Doses, //Vector of numeric doses assigned to patients
    arma::vec Dose, //Vector of Doses considered in the trial
    arma::vec DosesTried, //Vector of whether or not each dose has been tried
                             arma::vec Hypermeans, //6 vector of prior means
                             arma::vec Hypervars, //6 vector of prior standard deviations
                             arma::vec Contour, //4 vector of Contour parameter
                             arma::vec PiLim, //2 vector of acceptable limits
                             arma::vec ProbLim, //2 vector of cutoff for acceptabilities
                             int B // Number of reps to perform in MCMC
                               ){

  arma::vec Accept=PiLim;
  arma::vec Lower=ProbLim;

int Stopped=0;


  arma::vec MeanUt(Dose.n_rows);
  arma::vec PiEff(Dose.n_rows);
  arma::vec PiTox(Dose.n_rows);
  arma::vec ACCEPTTOX(Dose.n_rows);
  arma::vec ACCEPTEFF(Dose.n_rows);
  arma::vec DESIRE(Dose.n_rows);
  arma::vec ACCEPT(Dose.n_rows);
  ACCEPTTOX.zeros();
  ACCEPTEFF.zeros();

  double p=0;
  double pu=0;
  double pl=0;
  double tol=.005;



  //Now lets Solve for p to get the desirability
  while(abs1(pow((Contour[2]-1)/(Contour[0]-1),p) + pow(Contour[3]/Contour[1],p)-1) >tol){
    pu = p+.005;
    pl=p-.005;

    if( abs1(pow((Contour[2]-1)/(Contour[0]-1),pu) + pow((Contour[3])/(Contour[1]),pu)-1)<
      abs1(pow((Contour[2]-1)/(Contour[0]-1),pl) + pow((Contour[3])/(Contour[1]),pl)-1)){
      p=pu;
    }else{
      p=pl;
    }




  }



  //What Doses are Acceptable

  arma::vec ACCTOX(Dose.n_rows);
  arma::vec ACCEFF=ACCTOX;



  int Count=0;




  int B1 = B/2;

  //Used for for loop
  int b=0;
  //Used in Sim look
  int nREP=0;

  //Storage Matrices
  arma::mat piEStore(B1,Dose.n_rows);
  arma::mat piTStore(B1,Dose.n_rows);
  arma::mat AcceptTox(B1,Dose.n_rows);
  arma::mat AcceptEff(B1,Dose.n_rows);
  AcceptTox.zeros();
  AcceptEff.zeros();
  arma::vec DESIRE2(Dose.n_rows);
  DESIRE2.zeros();


  arma::mat ParamStore(B1,6);
  //Setup innitial
  arma::vec beta(5);
  beta.zeros();
  arma::vec betaprop=beta;
  double psi=0;
  double psiprop=0;
  //Quantities for MCMC
  double U=0;
  double alpha=0;
  int m=0;
  int j=0;
  int StoreInx=0;
  double PE=0;
  double PT=0;

  arma::vec bvar(6);
  bvar.zeros();
  bvar=bvar+1;
  arma::vec Intb = bvar;
  arma::vec Numb = bvar + 1;

  arma::vec Doses2=Doses;

  NumericVector z9(2);
  int CoNum=YE.n_rows;

  int STOPPED=0;
  int i=0;
  int OptDose=0;
  double min1=0;



  int nRep=0;

  arma::vec DESIRE1=DESIRE;

  arma::vec z5(2);
  arma::vec DoseTried(Dose.n_rows);





    STOPPED=0;



    Doses2.zeros();
    //Wrap Trial right here
    OptDose=0;
    Doses.zeros();
    DoseTried.zeros();
    DoseTried[0]=1;
    //Wrap around each cohort

      //Fill in next cohort with the data size


      Intb=Intb.zeros()+1;
      Numb=Numb.zeros()+2;
      bvar=bvar.zeros()+1;

      for(j=0;j<5;j++){
        beta[j]=Hypermeans[j];
      }
      psi=Hypermeans[5];

      //Perform MCMC
      for(b=0;b<B;b++){


        if(b<(B/2 + 2)){
          if(b%250==0){

            for(j=0;j<6;j++){
              if((Intb[j]/Numb[j])>.5){
                bvar[j]=bvar[j]*2;
              }

              if((Intb[j]/Numb[j])<.2){
                bvar[j]=bvar[j]/2;
              }


              Intb[j]=1;
              Numb[j]=2;

            }





          }
        }




        //All but beta_T can be less than 0, sample those
        for(j=0;j<4;j++){

          betaprop=beta;
          betaprop[j]=as_scalar(arma::randn(1))*bvar[j] + betaprop[j];

          alpha =   Like(YE,YT,Doses,betaprop,psi,CoNum) - Like(YE,YT,Doses, beta, psi,CoNum);
          alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

          U=log(as_scalar(arma::randu(1)));

          if(U<alpha){
            beta[j]=betaprop[j];
            Intb[j]=Intb[j]+1;
          }

          Numb[j] = Numb[j]+1;





        }


        //toxicity slope
        j=4;
        betaprop=beta;
        betaprop[j]=exp(as_scalar(arma::randn(1))*bvar[j] + log(betaprop[j]));

        alpha =   Like(YE,YT,Doses,betaprop,psi,CoNum) - Like(YE,YT,Doses, beta, psi,CoNum);
        alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

        alpha = alpha + log(betaprop[j]) - log(beta[j]);


        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          beta[j]=betaprop[j];
          Intb[j]=Intb[j]+1;
        }

        Numb[j] = Numb[j]+1;




        //Get new Psi

        //toxicity slope
        j=5;

        psiprop=as_scalar(arma::randn(1))*bvar[j] + psi;

        alpha =   Like(YE,YT,Doses,beta,psiprop,CoNum) - Like(YE,YT,Doses, beta, psi,CoNum);
        alpha = alpha - .5*pow(psiprop-Hypermeans[j],2)/Hypervars[j] + .5*pow(psi-Hypermeans[j],2)/Hypervars[j];

        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          psi=psiprop;
          Intb[j]=Intb[j]+1;
        }

        Numb[j]=Numb[j]+1;









        if(b>(B-B1-1)){

          StoreInx = b-B1;

          for(j=0;j<5;j++){
            ParamStore(StoreInx,j)=beta[j];
          }
          ParamStore(StoreInx,5)=psi;

          //Fill in Toxcity and Efficacy probabilities


          for(j=0;j<Dose.n_rows;j++){
            PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
            PT = beta[3]+beta[4]*Dose[j];

            PE = exp(PE);
            PT=exp(PT);
            PE=PE/(1+PE);
            PT = PT/(1+PT);
            piEStore(StoreInx,j)=PE;
            piTStore(StoreInx,j)=PT;


            if(PE>Accept[0]){
              AcceptEff(StoreInx,j)=1;
            }else{
              AcceptEff(StoreInx,j)=0;

            }

            if(PT<Accept[1]){
              AcceptTox(StoreInx,j)=1;
            }else{
              AcceptTox(StoreInx,j)=0;

            }






          }



        }



      }

      //Compute posterior mean



      for(j=0;j<Dose.n_rows;j++){
        ACCEPTEFF[j]=sum(AcceptEff.col(j))>(Lower(0)*B1)  ;
        ACCEPTTOX[j]=sum(AcceptTox.col(j))>(Lower(1)*B1) ;
      }




      //If ==1, for both EFF and TOX then dose j is acceptable

      //Determine Optimal Dose

      //Get Posterior mean of parameters
      for(j=0;j<5;j++){
        beta[j]=sum( ParamStore.col(j))/ParamStore.n_rows;
      }


      for(j=0;j<Dose.n_rows;j++){
        PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
        PT = beta[3]+beta[4]*Dose[j];

        PE = exp(PE);
        PT=exp(PT);
        PE=PE/(1+PE);
        PT = PT/(1+PT);
        PiEff(j)=PE;
        PiTox(j)=PT;

      }






      //Now get the mean desireability for each dose if they are acceptable
      for(j=0;j<Dose.n_rows;j++){
        if((ACCEPTEFF[j]+ACCEPTTOX[j])==2){
          //Dose is ACCEPTABLE
          DESIRE[j]=1-pow(pow(((PiEff(j)-1)/(Contour[0]-1)),p)+pow((PiTox(j)/Contour[1]),p),1/p);
          ACCEPT[j]=1;

        }else{
          ACCEPT[j]=0;
          DESIRE[j]=-1000;
        }



      }



      //Check if we should stop the trial

      if(sum(ACCEPT)>0){


        //Now pick the dose





          OptDose=0;

          for(j=1;j<Dose.n_rows;j++){
            if(DESIRE[j]>DESIRE[j-1]){
              OptDose=j;
            }
          }

          //Is this dose bigger than the dose higher than what's been tried?
          if(DoseTried(OptDose)==0){
            j=0;

            while(DoseTried(j)==1){
              j++;
            }

            OptDose=j;

          }

















      }else{
        //Trial Stopped
        OptDose=-1;
        STOPPED=1;

      }








  if(Stopped==1){
    return(-1);
  }else{
    return(OptDose);
  }


}




//' Randomizes Eff-Tox dose proportional to posterior desireability scores.
//'
//' This function returns a random acceptable dose number to assign the next patient cohort or stops the trial if no dose is deemed acceptable.
//' @param YE   Vector containing observed efficacy indicators.
//' @param YT   Vector containing observed toxicity indicators.
//' @param Doses Vector containing numbered Doses of patients in trial.
//' @param Dose Vector containing the standardized doses considered.
//' @param DosesTried Binary vector corresponding to which doses have been tried.
//' @param Hypermeans Vector containing prior hypermeans of length 6 for Eff-Tox parameters.
//' @param Hypervars Vector containing prior hypervariances of length 6 for Eff-Tox parameters.
//' @param Contour Vector containing 4 entries used to make the desireability function. Contour[1] contains a desired toxicity probability given efficacy, Countour[2] contains a desired efficacy probability given toxicity, and (Contour[3],Contour[4]) is an equally desireable pair of efficacy and toxicity probabilities that are non-zero or one.
//' @param PiLim Vector of length two with PiLim[1] containing the acceptable lower limit on efficacy probability and PiLim[2] containing the acceptable upper limit on toxicity probability.
//' @param ProbLim Vector of length two with ProbLim[1] containing the probability cutoff for acceptable efficacy probability and ProbLim[2] containing the probability cutoff for acceptable toxicity probability.
//' @param B Number of iterations to perform in the MCMC.
//' @return A random dose level to administer the next patient cohort.
//'@examples
//'##Doses, YE,YT
//'Doses= c(1,1,1,2,2,2,1,1,1,3,3,3,1,1,1,2,2,2)
//'YE = c(0,0,1,1,1,0,0,0,0,1,1,1,0,0,1,1,1,0)
//'YT=c(0,0,0,1,1,0,1,0,0,1,1,1,0,0,0,1,0,0)
//'##Vector of Numerical Doses
//'Dose = c(1,2,3,3.5,5)
//'Dose=(Dose-mean(Dose))/sd(Dose)
//'##Five doses, but only 3 tried so we have
//'DosesTried=c(1,1,1,0,0)
//'## Contour Vector
//'Contour = c(.35, .75,.7,.4)
//'##Hypermeans
//'Hypermeans = c(.022,3.45,0,-4.23,3.1,0)
//'Hypervars = c(2.6761, 2.6852, .2, 3.1304, 3.1165, 1)
//'Hypervars=Hypervars^2
//'##Acceptability Criteria
//'PiLim = c(.3,.4)
//'ProbLim=c(.1,.1)
//'##Number of iterations
//'B=2000
//'RandomEffTox(YE,YT, Doses, Dose, DosesTried, Hypermeans,  Hypervars, Contour, PiLim,  ProbLim, B )
//'@export
//[[Rcpp::export]]
int RandomEffTox( arma::vec YE, //Observed Efficacy Indicator Vector
                  arma::vec YT, //Observed toxicity indicator vector
                  arma::vec Doses, //Vector of numeric doses assigned to patients
                  arma::vec Dose, //Vector of Doses considered in the trial
                  arma::vec DosesTried, //Vector of whether or not each dose has been tried
                  arma::vec Hypermeans, //6 vector of prior means
                  arma::vec Hypervars, //6 vector of prior standard deviations
                  arma::vec Contour, //4 vector of Contour parameter
                  arma::vec PiLim, //2 vector of acceptable limits
                  arma::vec ProbLim, //2 vector of cutoff for acceptabilities
                  int B // Number of reps to perform in MCMC
){

  int Stopped=0;

  arma::vec Accept=PiLim;
  arma::vec Lower=ProbLim;


  arma::vec MeanUt(Dose.n_rows);
  arma::vec PiEff(Dose.n_rows);
  arma::vec PiTox(Dose.n_rows);
  arma::vec ACCEPTTOX(Dose.n_rows);
  arma::vec ACCEPTEFF(Dose.n_rows);
  arma::vec DESIRE(Dose.n_rows);
  arma::vec ACCEPT(Dose.n_rows);
  ACCEPTTOX.zeros();
  ACCEPTEFF.zeros();

  double p=0;
  double pu=0;
  double pl=0;
  double tol=.005;



  //Now lets Solve for p to get the desirability
  while(abs1(pow((Contour[2]-1)/(Contour[0]-1),p) + pow(Contour[3]/Contour[1],p)-1) >tol){
    pu = p+.005;
    pl=p-.005;

    if( abs1(pow((Contour[2]-1)/(Contour[0]-1),pu) + pow((Contour[3])/(Contour[1]),pu)-1)<
      abs1(pow((Contour[2]-1)/(Contour[0]-1),pl) + pow((Contour[3])/(Contour[1]),pl)-1)){
      p=pu;
    }else{
      p=pl;
    }




  }



  //What Doses are Acceptable

  arma::vec ACCTOX(Dose.n_rows);
  arma::vec ACCEFF=ACCTOX;



  int Count=0;




  int B1 = B/2;

  //Used for for loop
  int b=0;
  //Used in Sim look
  int nREP=0;

  //Storage Matrices
  arma::mat piEStore(B1,Dose.n_rows);
  arma::mat piTStore(B1,Dose.n_rows);
  arma::mat AcceptTox(B1,Dose.n_rows);
  arma::mat AcceptEff(B1,Dose.n_rows);
  AcceptTox.zeros();
  AcceptEff.zeros();
  arma::vec DESIRE2(Dose.n_rows);
  DESIRE2.zeros();


  arma::mat ParamStore(B1,6);
  //Setup innitial
  arma::vec beta(5);
  beta.zeros();
  arma::vec betaprop=beta;
  double psi=0;
  double psiprop=0;
  //Quantities for MCMC
  double U=0;
  double alpha=0;
  int m=0;
  int j=0;
  int StoreInx=0;
  double PE=0;
  double PT=0;

  arma::vec bvar(6);
  bvar.zeros();
  bvar=bvar+1;
  arma::vec Intb = bvar;
  arma::vec Numb = bvar + 1;

  arma::vec Doses2=Doses;

  NumericVector z9(2);
  int CoNum=YE.n_rows;

  int STOPPED=0;
  int i=0;
  int OptDose=0;
  double min1=0;



  int nRep=0;

  arma::vec DESIRE1=DESIRE;

  arma::vec z5(2);
  arma::vec DoseTried(Dose.n_rows);





  STOPPED=0;



  Doses2.zeros();
  //Wrap Trial right here
  OptDose=0;
  Doses.zeros();
  DoseTried.zeros();
  DoseTried[0]=1;
  //Wrap around each cohort

  //Fill in next cohort with the data size


  Intb=Intb.zeros()+1;
  Numb=Numb.zeros()+2;
  bvar=bvar.zeros()+1;

  for(j=0;j<5;j++){
    beta[j]=Hypermeans[j];
  }
  psi=Hypermeans[5];

  //Perform MCMC
  for(b=0;b<B;b++){


    if(b<(B/2 + 2)){
      if(b%250==0){

        for(j=0;j<6;j++){
          if((Intb[j]/Numb[j])>.5){
            bvar[j]=bvar[j]*2;
          }

          if((Intb[j]/Numb[j])<.2){
            bvar[j]=bvar[j]/2;
          }


          Intb[j]=1;
          Numb[j]=2;

        }





      }
    }




    //All but beta_T can be less than 0, sample those
    for(j=0;j<4;j++){

      betaprop=beta;
      betaprop[j]=as_scalar(arma::randn(1))*bvar[j] + betaprop[j];

      alpha =   Like(YE,YT,Doses,betaprop,psi,CoNum) - Like(YE,YT,Doses, beta, psi,CoNum);
      alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        beta[j]=betaprop[j];
        Intb[j]=Intb[j]+1;
      }

      Numb[j] = Numb[j]+1;





    }


    //toxicity slope
    j=4;
    betaprop=beta;
    betaprop[j]=exp(as_scalar(arma::randn(1))*bvar[j] + log(betaprop[j]));

    alpha =   Like(YE,YT,Doses,betaprop,psi,CoNum) - Like(YE,YT,Doses, beta, psi,CoNum);
    alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

    alpha = alpha + log(betaprop[j]) - log(beta[j]);


    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      beta[j]=betaprop[j];
      Intb[j]=Intb[j]+1;
    }

    Numb[j] = Numb[j]+1;




    //Get new Psi

    //toxicity slope
    j=5;

    psiprop=as_scalar(arma::randn(1))*bvar[j] + psi;

    alpha =   Like(YE,YT,Doses,beta,psiprop,CoNum) - Like(YE,YT,Doses, beta, psi,CoNum);
    alpha = alpha - .5*pow(psiprop-Hypermeans[j],2)/Hypervars[j] + .5*pow(psi-Hypermeans[j],2)/Hypervars[j];

    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      psi=psiprop;
      Intb[j]=Intb[j]+1;
    }

    Numb[j]=Numb[j]+1;









    if(b>(B-B1-1)){

      StoreInx = b-B1;

      for(j=0;j<5;j++){
        ParamStore(StoreInx,j)=beta[j];
      }
      ParamStore(StoreInx,5)=psi;

      //Fill in Toxcity and Efficacy probabilities


      for(j=0;j<Dose.n_rows;j++){
        PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
        PT = beta[3]+beta[4]*Dose[j];

        PE = exp(PE);
        PT=exp(PT);
        PE=PE/(1+PE);
        PT = PT/(1+PT);
        piEStore(StoreInx,j)=PE;
        piTStore(StoreInx,j)=PT;


        if(PE>Accept[0]){
          AcceptEff(StoreInx,j)=1;
        }else{
          AcceptEff(StoreInx,j)=0;

        }

        if(PT<Accept[1]){
          AcceptTox(StoreInx,j)=1;
        }else{
          AcceptTox(StoreInx,j)=0;

        }






      }



    }



  }

  //Compute posterior mean



  for(j=0;j<Dose.n_rows;j++){
    ACCEPTEFF[j]=sum(AcceptEff.col(j))>(Lower(0)*B1)  ;
    ACCEPTTOX[j]=sum(AcceptTox.col(j))>(Lower(1)*B1) ;
  }




  //If ==1, for both EFF and TOX then dose j is acceptable

  //Determine Optimal Dose

  //Get Posterior mean of parameters
  for(j=0;j<5;j++){
    beta[j]=sum( ParamStore.col(j))/ParamStore.n_rows;
  }


  for(j=0;j<Dose.n_rows;j++){
    PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
    PT = beta[3]+beta[4]*Dose[j];

    PE = exp(PE);
    PT=exp(PT);
    PE=PE/(1+PE);
    PT = PT/(1+PT);
    PiEff(j)=PE;
    PiTox(j)=PT;

  }






  //Now get the mean desireability for each dose if they are acceptable
  for(j=0;j<Dose.n_rows;j++){
    if((ACCEPTEFF[j]+ACCEPTTOX[j])==2){
      //Dose is ACCEPTABLE
      DESIRE[j]=1-pow(pow(((PiEff(j)-1)/(Contour[0]-1)),p)+pow((PiTox(j)/Contour[1]),p),1/p);
      ACCEPT[j]=1;

    }else{
      ACCEPT[j]=0;
      DESIRE[j]=-1000;
    }



  }



  //Check if we should stop the trial

  if(sum(ACCEPT)>0){


    //Now pick the dose



    OptDose=  GetDose1(DESIRE);

    //Is this dose bigger than the dose higher than what's been tried?
    if(DoseTried(OptDose)==0){
      j=0;

      while(DoseTried(j)==1){
        j++;

      }

      OptDose=j;

    }













  }else{
    //Trial Stopped
    OptDose=-1;
    STOPPED=1;

  }








  if(Stopped==1){
    return(-1);
  }else{
    return(OptDose);
  }


}




double MinVec(arma::vec Y){
  int J1=Y.n_rows;
  int j;
  double max=Y[0];
  for(j=1;j<J1;j++){
    if(Y[j]<max){
      max=Y[j];
    }
  }

  return(max);

}




int SampleDeath(int L){
  double U=as_scalar(arma::randu(1))*L;

  U=floor(U)+1;

  return(U);

}



int SampleBirth(arma::vec s){
  //Vector to store cumulative probility of being in an interval or lower
  arma::vec cumprob(s.n_rows-1);
  int m=0;

  cumprob[0] = s[1]/s[s.n_rows-1];

  for(m=1;m<cumprob.n_rows;m++){
    cumprob[m]=s[m+1]/s[s.n_rows-1];
  }


  //Randomly Draw Uniform that we'll use to determine what interval to add a move to.
  double U=as_scalar(arma::randu(1));

  // Find which interval the random unif is in.
  int Which1 =0;
  if(U<cumprob[0]){
    Which1=0;
  }else{

    for(m=1;m<cumprob.n_rows;m++){
      if( (U>cumprob[m-1]) && (U<cumprob[m])){
        Which1=m;
      }
    }

  }

  return(Which1);

}



//' Returns posterior distribution for key mixture model parameters
//'
//' This function performs MCMC with Metropolis-Hastings-Green steps for the baseline hazard function and is used in the functions Reoptimize, SimPhase123 and SimPhase3.
//' @param Y Patient survival or followup times.
//' @param I Patient event indicators.
//' @param YE Vector of indicators for patient efficacy.
//' @param YT Vector of indicators for patient toxicity.
//' @param Doses Vector of standardized doses given to patients.
//' @param Dose Vector of standardized doses considered in trial.
//' @param B Number of iterations to perform in MCMC.
//' @param prob  length(Doses) X 4 matrix containing the estimated posterior probabilities for each dose and each (Efficacy, Toxicity) outcomes.
//' @param MaxObs length(Doses) X 4 matrix containing the maximum observed survival time we want to evaluate the means to.
//'@return Returns a list containing a matrix of posterior means for each dose, regression coefficients in the cox models, locations of the split points, log hazard heights on each interval, and the number of intervals in the baseline hazard.
//'@examples
//'n=100
//'Y=rexp(n,1)
//'I = rbinom(n,1,.9)
//'YE = rbinom(n,1,.5)
//'YT = rbinom(n,1,.5)
//'Dose = c(1,2,3,3.5,5)
//'Dose=(Dose-mean(Dose))/sd(Dose)
//'Doses = sample(1:5,n,replace=TRUE)
//'Doses=Dose[Doses]
//'B=2000
//'MaxObs = matrix(rep(0,length(Dose)*4),nrow=4)
//'prob=matrix(rep(0,length(Dose)*4),ncol=4)
//'prob=prob+1/4
//'MaxObs=MaxObs+max(Y)
//'G=PieceMCMC(Y,I,YE,YT,Doses,Dose,B,prob,MaxObs)
//'@export
//[[Rcpp::export]]
List PieceMCMC( arma::vec Y, //Survival Times
                arma::vec I, //Censoring Indicators
                arma::vec YE, //Indicators of short term patient Efficacy
                arma::vec YT, //Indicators of short term patient Toxicity
                arma::vec Doses,  //Vector of Standardized Dose values given to patients
                arma::vec Dose, //Set of standardized doses to obtain means for
                int B, // Number of iterations to perform
                arma::mat prob, //Matrix of dose specific probabilities for the for groups
                arma::mat MaxObs //vector containing the maximum

){

  double YEmax = MaxVec(YE);
  double YEmin = MinVec(YE);

  double YTmax = MaxVec(YT);
  double YTmin = MinVec(YT);



  //Prior Params, make user controlled later
  double Poi = 5;

  //For loop quantiites
  int m=0;
  int b=0;
  int k=0;
  int j;
  //Important double quantities used in MCMC
  int B1=B;

  NumericVector z11(5);


  //Innitialize Parameters and storage matrices
  arma::vec beta(4);
  beta.zeros();

  beta[0]=0;
  beta[1]=0;
  beta[2]=0;
  beta[3]=0;


  arma::mat betastore(B,4);
  arma::mat MeanStore(B,Dose.n_rows);
  //Max To Store, really is 1 less here
  int Lmax = 21;
  arma::mat sstore(B,Lmax+2);
  arma::mat lamstore(B,Lmax+1);
  arma::vec Lstore(B);
  arma::vec sigstore(B);
  //Initialize S, lam and J
  int L = 3;
  arma::vec s(Lmax+2);
  arma::vec lam(Lmax+1);
  s.zeros();
  lam.zeros();
  double sig=1;
  double Birth=0;
  int Spot=0;

  arma::vec betaprop=beta;
  double U=0;
  double alpha=0;
  double varbeta=100;
  double varbeta1=100;
  double svar=1;
  double Ints=1;
  double Nums=2;
  //Copy Paste these to top
  arma::vec Means(Dose.n_rows);
  double mean2=0;
  double cum1=0;
  double mean1=0;
  double Con=0;
  //Fill \pi(Y_E,Y_T) with the empirical frequencies in the data




  NumericVector z9(2);


  double m1 = MaxVec(Y);

  arma::vec beta1=beta;


  //Start with split points at every $L/4$
  s[1]=m1/4;
  s[2]=2*m1/4;
  s[3]=3*m1/4;
  s[4]=m1;

  double NewLam1=0;
  double NewLam2=0;
  double U1=0;


  arma::vec lam1=lam;

  arma::vec bprop(4);
  arma::vec Intb(4);
  arma::vec Numb(4);
  arma::vec lamprop=lam;
  arma::vec sprop=s;

  for(j=0;j<4;j++){
    bprop[j]=1;
    Intb[j]=1;
    Numb[j]=2;
  }




  int StoreInx=0;

  for(m=0;m<B;m++){


    if(m<(B/2 + 2)){
      if(m%250==0){

        for(j=0;j<4;j++){
          if((Intb[j]/Numb[j])>.8){
            bprop[j]=bprop[j]*2;
          }

          if((Intb[j]/Numb[j])<.2){
            bprop[j]=bprop[j]/2;
          }


          Intb[j]=1;
          Numb[j]=2;

        }


        //svar
        if(Ints/Nums>.8){
          svar=svar*2;
        }

        if(Ints/Nums<.2){
          svar=svar/2;
        }

        Ints=1;
        Nums=2;



      }
    }







    betaprop = beta;


    //  betaprop[0]=-max1(exp(log(-betaprop[0])+as_scalar(arma::randn(1))*bprop[0]),.01);


    betaprop[0]=betaprop[0]+as_scalar(arma::randn(1))*bprop[0];




    alpha = - .5*(pow(betaprop[0],2) -   pow(beta[0],2))/varbeta  +  Like2(Y,I, YE, YT, Doses, betaprop, s, lam , L) -    Like2(Y,I, YE, YT, Doses, beta, s, lam , L);

    U = log(as_scalar(arma::randu(1)));



    if(U<alpha){
      beta[0]=betaprop[0];
      Intb[0]=Intb[0]+1;
    }
    Numb[0]=Numb[0]+1;







    betaprop = beta;




    betaprop[3]=betaprop[3]+as_scalar(arma::randn(1))*bprop[3];




    alpha = - .5*(pow(betaprop[3],2) -   pow(beta[3],2))/varbeta  +  Like2(Y,I, YE, YT, Doses, betaprop, s, lam , L) -    Like2(Y,I, YE, YT, Doses, beta, s, lam , L);

    U = log(as_scalar(arma::randu(1)));




    if(U<alpha){
      beta[3]=betaprop[3];
      Intb[3]=Intb[3]+1;
    }
    Numb[3]=Numb[3]+1;









    betaprop = beta;




    betaprop[1]=betaprop[1]+as_scalar(arma::randn(1))*bprop[1];




    alpha = - .5*(pow(betaprop[1],2) -   pow(beta[1],2))/varbeta1  +  Like2(Y,I, YE, YT, Doses, betaprop, s, lam , L) -    Like2(Y,I, YE, YT, Doses, beta, s, lam , L);

    U = log(as_scalar(arma::randu(1)));




    if(U<alpha){
      beta[1]=betaprop[1];
      Intb[1]=Intb[1]+1;
    }

    Numb[1]=Numb[1]+1;





    /// \beta_3: Toxicity





    betaprop = beta;


    betaprop[2]=betaprop[2]+as_scalar(arma::randn(1))*bprop[2];



    alpha = - .5*(pow(betaprop[2],2) -   pow(beta[2],2))/varbeta1  +  Like2(Y,I, YE, YT, Doses, betaprop, s, lam , L) -    Like2(Y,I, YE, YT, Doses, beta, s, lam , L);

    U = log(as_scalar(arma::randu(1)));




    if(U<alpha){
      beta[2]=betaprop[2];
      Intb[2]=Intb[2]+1;
    }

    Numb[2]=Numb[2]+1;

    //Sample lam



    if(L>0){

      for(j=0;j<(L+1);j++){

        lamprop = lam;


        lamprop[j]=lam[j]+as_scalar(arma::randu(1))*.5 - .25;


        alpha=Like2(Y,I,YE,YT,  Doses, beta, s, lamprop , L) -    Like2(Y,I,YE,YT,  Doses, beta, s, lam , L);

        if(j==0){
          alpha = alpha - .5*(pow(lamprop[j],2)-pow(lam[j],2))/sig - .5*pow(lamprop[j+1]-lamprop[j],2)/sig + .5*pow(lam[j+1]-lam[j],2)/sig;
        }else{
          alpha = alpha - .5*(pow(lamprop[j]-lamprop[j-1],2) + pow(lam[j]-lam[j-1],2))/sig  - .5*(pow(lamprop[j+1]-lamprop[j],2) + pow(lam[j+1]-lam[j],2))/sig   ;
        }

        U = log(as_scalar(arma::randu(1)));





        if(U<alpha){
          lam=lamprop;
        }


      }

    }else{
      j=0;


      lamprop = lam;


      lamprop[j]=lam[j]+as_scalar(arma::randu(1))*.5 - .25;


      alpha=Like2(Y,I,YE,YT,  Doses, beta, s, lamprop , L) -    Like2(Y,I,YE,YT,  Doses, beta, s, lam , L);

      alpha = alpha - .5*(pow(lamprop[j],2))/25+.5*pow(lam[j],2)/25 ;


      U = log(as_scalar(arma::randu(1)));





      if(U<alpha){
        lam=lamprop;
      }

    }

    //Sample Siglam
    cum1=0;
    if(L>0){
      for(j=1;j<L;j++){
        cum1 =cum1 + pow(lam[j]-lam[j-1],2);
      }

      sig = 1/R::rgamma(L/2 + 1, cum1/2 );

    }else{
      sig=10000;
    }


    //Shuffle Splits

    if(L>0){

      for(j=1; j<(L+1); j++){
        sprop = s;
        //Draw new proposal for s_j, but make sure it's between s_j and s_j+1
        //Draws unif[-c,c] proposal addition
        //    sprop[j]=max1(min1(as_scalar(arma::randu(1))*2*svar-svar+sprop[j],sprop[j+1]-.01),sprop[j-1]+.01);
        sprop[j]=as_scalar(arma::randu(1))*(sprop[j+1]-sprop[j-1]) + sprop[j-1];

        alpha = log(sprop[j+1]-sprop[j])+log(sprop[j]-sprop[j-1])-log(s[j+1]-s[j])+log(s[j]-s[j-1]) +   Like2(Y,I,YE,YT,  Doses, beta, sprop, lam , L) -    Like2(Y,I,YE,YT,  Doses, beta, s, lam , L);

        U = log(as_scalar(arma::randu(1)));





        if(U<alpha){

          s[j]=sprop[j];
          Ints = Ints+1;
        }

        Nums = Nums+1;



      }

    }



    //Add proposal and MHG ratios here

    //Birth move Here
    if(L<(Lmax-1)){
      //What interval is it in?

      lamprop.zeros();

      Spot = SampleBirth(s)+1;
      //     Birth  = R::runif(s[Spot-1],s[Spot+1]);


      Birth = as_scalar(arma::randu(1))*(s[Spot]-s[Spot-1])+s[Spot-1];
      //Random Perturbation
      U1=as_scalar(arma::randu(1));


      NewLam1 =  lam[Spot-1] - log((1-U1)/U1)*(s[Spot+1]-Birth)/(s[Spot+1]-s[Spot-1]);
      NewLam2 =  lam[Spot-1] + log((1-U1)/U1)*(Birth-s[Spot])/(s[Spot+1]-s[Spot-1]);



      //Now we have the interval of the new split in Spot and the actual location in Birth
      sprop.zeros();
      //Let's add Birth to the Spot location and push back the rest of sprop
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
        lamprop[j]=lam[j];
      }
      sprop[Spot]=Birth;
      for(j=(Spot+1);j<(s.n_rows-1);j++){
        sprop[j]=s[j-1];
      }

      lamprop[Spot-1]=NewLam1;
      lamprop[Spot]=NewLam2;

      for(j=(Spot+2);j<(lam.n_rows-1);j++){
        lamprop[j]=lam[j-1];
      }







      //Propose a new lambda corresponding to this split point and a random perturbation


      //Now we have our new proposal vector, evaluate it!
      //Like2 Ratio
      alpha =    Like2(Y,I,YE,YT,  Doses, beta, sprop, lamprop , L+1) -    Like2(Y,I,YE,YT, Doses, beta, s, lam , L);
      //Add proposal ratio
      //Poisson
      alpha= alpha + logINT(Poi) - logINT(L+1);
      // S proposal
      alpha = alpha + logINT(2*L+3)+logINT(2*L+2)+log(Birth-s[Spot-1])+log(s[Spot]-Birth);
      alpha = alpha - 2*log(m1)-log(s[Spot]-s[Spot-1]);
      //Perturbation
      alpha=alpha-log(U1*(1-U1)) ;
      //Add proposal ratio for \lambda

      if(L==0){
        alpha = alpha - .5*pow(lamprop[1]-lamprop[0],2)/sig   ;



      }else{

        if(Spot==(L+1)){
          alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig + pow(lam[Spot-1]-lam[Spot-2],2)/sig;


        }else{


          if(Spot==1){




            alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig ;
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig ;


          }else{





            alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;





          }



        }




      }



      //Add last interval
      alpha = alpha  - .5*pow(lamprop[0],2)/25 + .5*pow(lam[0],2)/25 ;










      U=log(as_scalar(arma::randu(1)));


      if(U<alpha){
        s=sprop;
        lam=lamprop;
        L=L+1;
      }



    }



    //Death move here
    if(L>0){
      //Which one to delete??
      Spot=SampleDeath(L);
      //Setup sprop with deleted value here
      sprop=s;
      sprop.zeros();
      lamprop.zeros();

      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
      }



      if(L>1){

        for(j=0;j<(Spot-1);j++){
          lamprop[j]=lam[j];
        }

        for(j=Spot;j<(lam.n_rows-1);j++){
          lamprop[j]=lam[j+1];
        }

        //New lambda is a weighted average of the old lambdas
        lamprop[Spot-1] = ((s[Spot+1]-s[Spot])*lam[Spot] + (s[Spot]-s[Spot-1])*lam[Spot-1])/(s[Spot+1]-s[Spot-1]);


      }else{
        lamprop[0]=((s[2]-s[1])*lam[1]+(s[1]-s[0])*lam[0])/(s[2]-s[0]);
      }

      for(j=Spot; j<(s.n_rows-1);j++){
        sprop[j]=s[j+1];
      }



      //Now we have our new proposal vector, evaluate it!
      //Like2 Ratio
      alpha =    Like2(Y,I,YE,YT,  Doses, beta, sprop, lamprop , L-1) -    Like2(Y,I,YE,YT,  Doses, beta, s, lam , L);
      //Prior Ratio
      //Poisson
      //Poisson
      alpha = alpha  -logINT(Poi) + logINT(L);
      //S Prior
      alpha = alpha + 2*log(m1) + log(s[Spot+1]-s[Spot-1]) - logINT(2*L+1) - logINT(2L) - log(s[Spot]-s[Spot-1])-log(s[Spot+1]-s[Spot]);
      //Lambda Prior, we DROPPED one

      if(L==1){
        //Removing will drop the sampler to 0 split points

        alpha = alpha + .5*pow(lam[1]-lam[0],2)/sig ;



      }else{



        if(Spot>1){

          if(Spot==(L+1)){
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig  - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;





          }else{

            //Split point is not first or last

            alpha = alpha +  .5*pow(lam[Spot+1]-lam[Spot],2)/sig +  .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;

            alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig ;

            alpha = alpha - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig ;



          }


        }else{

          //First Interval is changing

          alpha = alpha + .5*pow(lam[Spot+1]-lam[Spot],2)/sig + .5*pow(lam[Spot]-lam[Spot-1],2)/sig;

          alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig  ;



        }





      }



      //Add the 0 entry here
      alpha= alpha  - .5*pow(lamprop[0],2)/25 + .5*pow(lam[0],2)/25 ;


      //Random Perturbation
      U1 = as_scalar(arma::randu(1));
      alpha = alpha + log(U1)+log(1-U1);


      U=log(as_scalar(arma::randu(1)));


      if(U<alpha){
        s=sprop;
        L=L-1;
        lam=lamprop;
      }




    }














    //Calculate Means Here
    lam1=exp(lam);

    beta1=beta;

    beta1[1]=-exp(beta1[1]);

    beta1[2]=exp(beta1[2]);



    //Wrap around each dose
    for(j=0;j<Dose.n_rows;j++){

      //Holds Running mean total of the four outcomes
      mean1=0;




      if(L>0){



        Con = exp(beta1[0]*Dose[j]+beta1[3]*pow(Dose[j],2)+beta1[1]*YEmin+beta1[2]*YTmin);

        cum1=0;
        mean2=0;





        for(k=0;k<(L+1);k++){
          cum1 = cum1+ lam1[k]*(s[k+1]-s[k]);

          mean2=mean2+ min1(s[k+1]-s[k],exp(-Con*cum1)*exp(lam1[k]*Con*s[k])*(exp(-Con*lam1[k]*s[k])-exp(-Con*lam1[k]*s[k+1]))/(Con*lam1[k]));









        }


        mean1=mean1 +  prob(j,0)*mean2;







        Con = exp(beta1[0]*Dose[j]+beta1[3]*pow(Dose[j],2)  +beta1[1]*YEmax+beta1[2]*YTmin);

        cum1=0;
        mean2=0;




        for(k=0;k<(L+1);k++){
          cum1 = cum1+ lam1[k]*(s[k+1]-s[k]);

          mean2=mean2+ min1(s[k+1]-s[k],exp(-Con*cum1)*exp(lam1[k]*Con*s[k])*(exp(-Con*lam1[k]*s[k])-exp(-Con*lam1[k]*s[k+1]))/(Con*lam1[k]));









        }




        mean1=mean1 +  prob(j,1)*mean2;



        Con = exp(beta1[0]*Dose[j]+beta1[3]*pow(Dose[j],2)+ beta1[2]*YTmax + beta1[1]*YEmin);

        cum1=0;
        mean2=0;





        for(k=0;k<(L+1);k++){
          cum1 = cum1+ lam1[k]*(s[k+1]-s[k]);

          mean2=mean2+ min1(s[k+1]-s[k],exp(-Con*cum1)*exp(lam1[k]*Con*s[k])*(exp(-Con*lam1[k]*s[k])-exp(-Con*lam1[k]*s[k+1]))/(Con*lam1[k]));









        }



        mean1=mean1 +  prob(j,2)*mean2;





        Con = exp(beta1[0]*Dose[j]+beta1[3]*pow(Dose[j],2) + beta1[1]*YEmax+beta1[2]*YTmax);

        cum1=0;
        mean2=0;




        for(k=0;k<(L+1);k++){
          cum1 = cum1+ lam1[k]*(s[k+1]-s[k]);

          mean2=mean2+ min1(s[k+1]-s[k],exp(-Con*cum1)*exp(lam1[k]*Con*s[k])*(exp(-Con*lam1[k]*s[k])-exp(-Con*lam1[k]*s[k+1]))/(Con*lam1[k]));









        }



        mean1=mean1 +  prob(j,3)*mean2;


      }else{







        k=0;


        Con = exp(beta1[0]*Dose[j]+beta1[3]*pow(Dose[j],2)+YTmin*beta1[3]+YEmin*beta1[2]+lam[0]);


        mean2= min1((1-exp(-Con*MaxObs(0,j)))/Con,s[1]);







        mean1=mean1 +  prob(j,0)*mean2;

        z9[0]=mean2;






        Con = exp(beta1[0]*Dose[j]+beta1[3]*pow(Dose[j],2) + YEmax*beta1[1]+YTmin*beta1[2]+lam[0]);

        cum1=0;
        mean2=0;



        mean2=  min1((1-exp(-Con*MaxObs(1,j)))/Con,s[1]);






        mean1=mean1 +  prob(j,1)*mean2;

        z9[1]=mean2;







        //Outcome 3, toxicity and NO efficacy
        Con = exp(beta1[0]*Dose[j]+beta1[3]*pow(Dose[j],2) + beta1[2]*YTmax+beta1[1]*YEmin+lam[0]);





        mean2=  min1((1-exp(-Con*MaxObs(2,j)))/Con,s[1]);










        mean1=mean1 +  prob(j,2)*mean2;

        z9[0]=mean2;






        Con = exp(beta1[0]*Dose[j]+beta1[3]*pow(Dose[j],2) +beta1[1]*YEmax+beta1[2]*YTmax+lam[0]);



        mean2=  min1((1-exp(-Con*MaxObs(3,j)))/Con,s[1]);






        mean1=mean1 +  prob(j,3)*mean2;

        z9[1]=mean2;







      }


      //Fill with mean of each dose
      Means(j)=mean1;






    }




    //   Rf_PrintValue(wrap(Means));



    //Store Values in Matrix
    StoreInx = m;

    for(j=0;j<sstore.n_cols;j++){
      sstore(StoreInx,j) = s(j);
    }

    for(j=0;j<lamstore.n_cols;j++){
      lamstore(StoreInx,j) = lam(j);
    }

    for(j=0;j<Dose.n_rows;j++){
      MeanStore(StoreInx,j)=Means(j);
    }

    Lstore[StoreInx]=L;
    sigstore[StoreInx]=sig;

    for(j=0;j<4;j++){
      betastore(StoreInx,j)=beta(j);
    }












  }








  List z1 = List::create(MeanStore,betastore,sstore,lamstore,Lstore,sigstore);


  return(z1);


}












//' Obtains estimated posterior probabilities of the four outcomes of (YE,YT) for each dose.
//'
//' This function is used in Reoptimize, SimPhase123 and SimPhase3, here we estimate the mixture probabilities over the four outcomes for efficacy and toxicity.
//' @param YE   Vector containing observed efficacy indicators.
//' @param YT   Vector containing observed toxicity indicators.
//' @param Doses Vector containing Standardized doses of patients in trial.
//' @param Dose Vector containing the standardized doses considered.
//' @param Hypermeans Vector containing prior hypermeans of length 6 for Eff-Tox parameters.
//' @param Hypervars Vector containing prior hypervariances of length 6 for Eff-Tox parameters.
//' @param B Number of iterations to perform in the MCMC.
//' @return The posterior probability matrix for the events (YE,YT) in each row corresponding to a dose level.
//'@examples
//'##Doses, YE,YT
//'Doses= c(1,1,1,2,2,2,1,1,1,3,3,3,1,1,1,2,2,2)
//'YE = c(0,0,1,1,1,0,0,0,0,1,1,1,0,0,1,1,1,0)
//'YT=c(0,0,0,1,1,0,1,0,0,1,1,1,0,0,0,1,0,0)
//'##Vector of Numerical Doses
//'Dose = c(1,2,3,3.5,5)
//'Dose=(Dose-mean(Dose))/sd(Dose)
//'Doses=Dose[Doses]
//'##Hypermeans
//'Hypermeans = c(.022,3.45,0,-4.23,3.1,0)
//'Hypervars = c(2.6761, 2.6852, .2, 3.1304, 3.1165, 1)
//'Hypervars=Hypervars^2
//'##Number of iterations
//'B=2000
//'EFFTOX(YE,YT, Doses, Dose, Hypermeans,  Hypervars, B )
//'@export
//[[Rcpp::export]]
arma::mat EFFTOX(arma::vec YE, //Binary indicators of efficacy
            arma::vec YT, //Binary indicators of Toxicity
            arma::vec Doses, //Vector of standardized Doses given to patients
            arma::vec Dose, //Vector of Doses considered in the trial
            arma::vec Hypermeans, //6 vector of prior means
            arma::vec Hypervars, //6 vector of prior standard deviations
            int B // Number of reps to perform in MCMC
){

  arma::vec MeanUt(Dose.n_rows);
  arma::vec PiEff(Dose.n_rows);
  arma::vec PiTox(Dose.n_rows);

  //What Doses are Acceptable

  arma::vec ACCTOX(Dose.n_rows);
  arma::vec ACCEFF=ACCTOX;

  int B1 = B/2;

  //Used for for loop
  int b=0;


  //Storage Matrices
  arma::mat piEStore(B1,Dose.n_rows);
  arma::mat piTStore(B1,Dose.n_rows);
  arma::mat AcceptTox(B1,Dose.n_rows);
  arma::mat AcceptEff(B1,Dose.n_rows);
  AcceptTox.zeros();
  AcceptEff.zeros();

int conum=YE.n_rows;
  arma::mat ParamStore(B1,6);
  //Setup innitial
  arma::vec beta(5);
  beta.zeros();
  arma::vec betaprop=beta;
  double psi=0;
  double psiprop=0;
  //Quantities for MCMC
  double U=0;
  double alpha=0;
  int m=0;
  int j=0;
  int StoreInx=0;
  double PE=0;
  double PT=0;

  arma::vec bvar(6);
  bvar.zeros();
  bvar=bvar+1;
  arma::vec Intb = bvar;
  arma::vec Numb = bvar + 1;

  NumericVector z9(2);



  //Perform MCMC
  for(b=0;b<B;b++){


    if(b<(B/2 + 2)){
      if(b%250==0){

        for(j=0;j<6;j++){
          if((Intb[j]/Numb[j])>.8){
            bvar[j]=bvar[j]*2;
          }

          if((Intb[j]/Numb[j])<.2){
            bvar[j]=bvar[j]/2;
          }


          Intb[j]=1;
          Numb[j]=2;

        }





      }
    }




    //All but beta_T can be less than 0, sample those
    for(j=0;j<4;j++){

      betaprop=beta;
      betaprop[j]=as_scalar(arma::randn(1))*bvar[j] + betaprop[j];

      alpha =   Like(YE,YT,Doses,betaprop,psi,conum) - Like(YE,YT,Doses, beta, psi,conum);
      alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        beta[j]=betaprop[j];
        Intb[j]=Intb[j]+1;
      }

      Numb[j] = Numb[j]+1;





    }


    //toxicity slope
    j=4;
    betaprop=beta;
    betaprop[j]=exp(as_scalar(arma::randn(1))*bvar[j] + log(betaprop[j]));

    alpha =   Like(YE,YT,Doses,betaprop,psi,conum) - Like(YE,YT,Doses, beta, psi,conum);
    alpha = alpha - .5*pow(betaprop[j]-Hypermeans[j],2)/Hypervars[j] + .5*pow(beta[j]-Hypermeans[j],2)/Hypervars[j];

    alpha = alpha + log(betaprop[j]) - log(beta[j]);


    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      beta[j]=betaprop[j];
      Intb[j]=Intb[j]+1;
    }

    Numb[j] = Numb[j]+1;




    //Get new Psi

    //toxicity slope
    j=5;

    psiprop=as_scalar(arma::randn(1))*bvar[j] + psi;

    alpha =   Like(YE,YT,Doses,beta,psiprop,conum) - Like(YE,YT,Doses, beta, psi,conum);
    alpha = alpha - .5*pow(psiprop-Hypermeans[j],2)/Hypervars[j] + .5*pow(psi-Hypermeans[j],2)/Hypervars[j];

    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      psi=psiprop;
      Intb[j]=Intb[j]+1;
    }

    Numb[j]=Numb[j]+1;









    if(b>(B-B1-1)){

      StoreInx = b-B1;

      for(j=0;j<5;j++){
        ParamStore(StoreInx,j)=beta[j];
      }
      ParamStore(StoreInx,5)=psi;

      //Fill in Toxcity and Efficacy probabilities


      for(j=0;j<Dose.n_rows;j++){
        PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
        PT = beta[3]+beta[4]*Dose[j];

        PE = exp(PE);
        PT=exp(PT);
        PE=PE/(1+PE);
        PT = PT/(1+PT);
        piEStore(StoreInx,j)=PE;
        piTStore(StoreInx,j)=PT;








      }



    }



  }



  //Compute posterior mean

  for(j=0;j<Dose.n_rows;j++){
    PiEff[j]=sum(piEStore.col(j))/piEStore.n_rows;
    PiTox[j]=sum(piTStore.col(j))/piEStore.n_rows;
  }

  //Get Posterior mean of parameters
  for(j=0;j<5;j++){
    beta[j]=sum( ParamStore.col(j))/ParamStore.n_rows;
  }


  for(j=0;j<Dose.n_rows;j++){
    PE = beta[0]+beta[1]*Dose[j]+beta[2]*pow(Dose[j],2);
    PT = beta[3]+beta[4]*Dose[j];

    PE = exp(PE);
    PT=exp(PT);
    PE=PE/(1+PE);
    PT = PT/(1+PT);
    PiEff(j)=PE;
    PiTox(j)=PT;

  }

  psi =   sum(ParamStore.col(5))/ParamStore.n_rows;



  //Get Outcome probs from this
  arma::mat OutCome(Dose.n_rows,4);
  double const1=0;


  for(j=0;j<Dose.n_rows;j++){
    const1 = PiEff[j]*(1-PiEff[j])*PiTox[j]*(1-PiTox[j])*(exp(psi)-1)/(exp(psi)+1);

    OutCome(j,0) = (1-PiEff[j])*(1-PiTox[j]) + const1;

    OutCome(j,1) = PiEff[j]*(1-PiTox[j]) - const1;

    OutCome(j,2) = (1-PiEff[j])*(PiTox[j]) - const1;

    OutCome(j,3) = (PiEff[j])*(PiTox[j]) + const1;



  }












  return(OutCome);
}
