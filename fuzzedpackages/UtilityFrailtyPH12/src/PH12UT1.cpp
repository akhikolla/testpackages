#include <time.h>
#include <Rmath.h>
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <Rcpp.h>
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


//'Returns the desireability value of a dose.
//'
//'Takes estimated posterior mean efficacy and toxicity values and returns the posterior mean desireability score for a given tradeoff contour.
//'@param PE True or estimated probability of efficacy.
//'@param PT True of estimated probability of toxicity.
//'@param Contour Vector containing 4 entries used to make the desireability function. Contour(1) contains a desired toxicity probability given efficacy, Countour(2) contains a desired efficacy probability given toxicity, and (Contour(3),Contour(4)) is an equally desireable pair of efficacy and toxicity probabilities that are non zero or one.
//'@useDynLib UtilityFrailtyPH12
//'@importFrom Rcpp evalCpp
//'@examples
//'PE=.6
//'PT=.2
//'##Contour values
//'Contour=c(.35,.7,.8,.6)
//'GetDesire(PE,PT,Contour)
//'@export
//[[Rcpp::export]]
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
  int Mean = sum(X)/X.n_rows;
  
  return(pow(sum(pow(X-Mean,2))/(X.n_rows-1),.5));
  
  
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











//' Samples from the posterior of the utility based phase12 model.
//' 
//' Takes arguments of data, hypermens and hypervariance vectors and returns a list of posterior samples from the Utility based phase12 model decribed by Chapple and Thall (2019).
//'
//'@param YE Binary indicator vector of efficacy status.
//'@param YT Binary indicator vector of toxicity status.
//'@param Doses Vector of integer Doses given to patients.
//'@param HypermeansEFF Vector of length nDose for dose prior means for efficacy.
//'@param HypermeansTOX Vector of length nDose for dose prior means for toxicity
//'@param Hypervars Length 5 vector of hypervariances. Hypervars(1) and Hypervars(2) contains the Latent parameter variance for normal probability of efficacy and toxicity. Hypervars(3) and Hypervars(4) contains the hypervariance on dose specific mean efficacy and toxicity parameters and Hypervars(5) contains the frailty variance parameter.
//'@param B Number of iterations to run for the MCMC.
//'@return A list of posterior samples after burnin in order: Posterior efficacy dose-vector, Posterior toxicity dose-vector, Posterior correlation.
//'@examples
//'n=100  #Generate Data
//'YE=rbinom(n,1,.6)
//'YT=rbinom(n,1,.2)
//'nDose=5
//'Doses=sample(1:nDose,n,replace=TRUE)
//'##Hyperparameters 
//'HypermeansEFF=c(-1,-.5,0,.5,1,2)
//'HypermeansTOX=HypermeansEFF
//'Hypervars=c(1,1,36,36,1)
//'B=100
//'UTEFFTOX(YE, YT,Doses,HypermeansEFF,HypermeansTOX, Hypervars, B)
//'@export
//[[Rcpp::export]]
List UTEFFTOX(arma::vec YE, //Binary indicators of efficacy
                 arma::vec YT, //Binary indicators of Toxicity
                 arma::vec Doses, //Vector of Numerical Doses given to patients
                 arma::vec HypermeansEFF, //Vector of length nDose + 1 for dose prior means for efficacy
                 arma::vec HypermeansTOX, //Vector of length nDose + 1 for dose prior means for toxicity
                 arma::vec Hypervars, // (\sigma^2_E, \sigma^2_T, \sigma^2_{E,0},\sigma^2_{T,0}, \tau^2, \corr between EFF-TOX )
                 int B // Number of reps to perform in MCMC
){
  
  
  //This should be a bit different in the actual trial conduct
  Doses=Doses-1;  //Formats dose vector for c++
  
  
  
  //Extract hyperparameters
  int nDose = HypermeansEFF.n_rows;

 //Get The hypervars
 double SigmaE = Hypervars[0];  //Dose effect EFF vars
 double SigmaT = Hypervars[1];  //Dose effect TOX vars
 double SigmaE0 = Hypervars[2]; //MEAN Dose effect EFF vars
 double SigmaT0 = Hypervars[3]; //MEAN Dose effect TOX vars
 double tau = Hypervars[4]; //Variance parameter for frailty
 double phi=0;
 double phinew=0;
 double rho = 2*exp(phi)/(1+exp(phi))-1;  //Correlation between EFF-TOX within frailties.
 
 
 
 
 
 //Important Integers
  int B1 = B/2;
  //Used for for loop
  int b=0;
  int i=0; //Patient level index
  int j =0; //For indexing doses
  
  
  int StoreInx = 0;
  
  
  //Ok First we need to set up our storage
  int n = Doses.n_rows; //Sample Size
  
  arma::vec epsilonE(n); //Vector containing latent effiacy
  arma::vec epsilonT(n); //Vector containing latent toxicity
  //Ok let's fill these with starting values depending on patient (Y_E, Y_T)
  epsilonE=YE - .5;
  epsilonT=YT - .5;
  //Starting values are all -.5 or .5

  //Get sums of counts within each dose
arma::vec NumDose(nDose);
int count=0; //count for number of treated at each dose
  //Loop over each dose
  for(j=0;j<nDose;j++){
    //Loop over each individual and count
    count=0;
    for(i=0;i<n;i++){
     if(Doses[i]==j){
       count++;
     } 
      
      
    }
    
    NumDose[j]=count;
    
  }

  //Now NumDose contains the number treated at each dose
  
  
  
  
  //Dose mean vectors
  arma::vec MuE(nDose);
  arma::vec MuT(nDose);
  MuE.zeros();
  MuT.zeros();
  
  //Storage matrix for these:
  arma::mat MuEStore(B1,nDose);
  arma::mat MuTStore(B1,nDose);
  //Storage for epsilons
  arma::mat epsilonEStore(B1,n);
  arma::mat epsilonTStore(B1,n);
  
  
  //MCMC Quantities
  double U = 0; //For metropolis hastings
  double alpha = 0; //For metropolis hastings
 double epprop = 0; //proposal values for epsilon 
double  muprop = 0; //proposal for mu_epsilon
  
  
  
  //Get solved covariance matrix that we need for evaluation of individual 
  //Latent EFF-TOX variables
 arma::mat SIGMA(2,2);
//Fill the matrix
SIGMA(0,0)=SigmaE + tau;
SIGMA(0,1)=rho*tau;
SIGMA(1,0)=rho*tau;
SIGMA(1,1)=SigmaT+tau;
//Ok now let's invert this matrix
arma::mat SIGMAINV=inv(SIGMA);


//Fill in prior variance matrix for dose specific parameters
arma::mat Sigma0(2,2);
Sigma0.zeros();
Sigma0(0,0)=1/SigmaE0;
Sigma0(1,1)=1/SigmaT0;
//Fill in prior mean vector for dose specific parameters
arma::vec MEANS0(2);
MEANS0[0]=0;
MEANS0[1]=0;

//Let's get all our sigma matrices so we don't need to  do them later.

arma::mat SigmaNEW(2,2); //Used in sampling
arma::mat SigmaNEWINV=SigmaNEW;

  
  //Make meanvec for MCMC sampling each latent vector
arma::vec MEANVEC(2);
MEANVEC.zeros();

//For counting sum of epsilons
arma::vec epVEC(2);
  
  //These are the efficacy and toxicity tuning params
  double epEC = 1;
  double epTC = 1;
  

  //Used for sampling condiitonal norm, med
  double munorm=0;
  double signorm=0;
  
  double rhonew=0;
  
  arma::mat SIGMANEW(2,2);
  arma::mat SIGMANEWINV(2,2);
  
  
  arma::vec rhoStore(B1);
  
  arma::vec z9(2);
  double corr=0;
  
  //counters
  double nrho = 1;
  double arho=1;
  
  double crho = 1;
  
  
  //Perform MCMC
  for(b=0;b<B;b++){

    
    if((b%100==0) & (b<(B-B1))){
      if(arho/nrho<.2){
        crho=crho/2;
        
      }
      
      if(arho/nrho>.5){
        crho=crho*2;
      }
      

      
      arho=1;
      nrho=1;
      
    }
    
  
  //Sample individual epsilons
  
  

    
    //Now let's sample the mean parameters
    for(j=0;j<nDose;j++){
      epVEC.zeros();
      //Get sum of epsilon vectors
      for(i=0;i<n;i++){
        if(Doses[i]==j){
          epVEC[0]=epVEC[0]+epsilonE[i];
          epVEC[1]=epVEC[1]+epsilonT[i];
        }
      }
      
      //This is dose specific hyperparameters
      MEANS0[0]=HypermeansEFF[j];
      MEANS0[1]=HypermeansTOX[j];
      
      //First make the mean and covariance matrices
      SigmaNEW = Sigma0 + NumDose[j]*SIGMAINV;
      
      

      SigmaNEWINV = inv(SigmaNEW);
   

      MEANVEC = SigmaNEWINV*(Sigma0*MEANS0 + SIGMAINV*epVEC);
  
  
  
  //Posterior distribution is truncated normal with above mean vector and covariance matrix
  //This comes straight from a conjugate normal setup.
  
  
  
corr=SigmaNEWINV(0,1)/pow(SigmaNEWINV(0,0)*SigmaNEWINV(1,1),.5);
      
      
      //Get the conditional distribution
      signorm = pow((1-corr)*SigmaNEWINV(0,0),.5);
      //Conditional mean
      munorm = MEANVEC[0] + corr*(MuT[j]-MEANVEC[1])*pow(SigmaNEWINV(0,0)/SigmaNEWINV(1,1),.5);
    

      
   
      
      //Generate new variable, where it's truncated
      
      if(j>0){
        
        if(j==(nDose-1)){
          MuE[j] = max1(MuE[j-1]+.001,as_scalar(arma::randn(1))*signorm+munorm);
          //Last dose level
        }else{
          //Now were in the middle doses
          MuE[j] = max1(MuE[j-1]+.001,min1(MuE[j+1]-.001,as_scalar(arma::randn(1))*signorm+munorm));
        }
        
        
        
      }else{
        //First dose mean efficacy parameter
        MuE[j] = min1(MuE[j+1]-.001,as_scalar(arma::randn(1))*signorm+munorm);
        
      }
      
      
      
      
      
      //Do MeanT for this dose next

      
      //Get the conditional distribution
      signorm = pow((1-corr)*SigmaNEWINV(1,1),.5);
      //Conditional mean
      munorm = MEANVEC[1] + corr*(MuE[j]-MEANVEC[0])*pow(SigmaNEWINV(1,1)/SigmaNEWINV(0,0),.5);
      
      
      
      //Generate new variable, where it's truncated
      if(j>0){
        
        if(j==(nDose-1)){
          MuT[j] = max1(MuT[j-1]+.001,as_scalar(arma::randn(1))*signorm+munorm);
          //Last dose level
        }else{
          //Now were in the middle doses
          MuT[j] = max1(MuT[j-1]+.001,min1(MuT[j+1]-.001,as_scalar(arma::randn(1))*signorm+munorm));
        }
        
        
        
      }else{
        //First dose mean efficacy parameter
        MuT[j] = min1(MuT[j+1]-.001,as_scalar(arma::randn(1))*signorm+munorm);
        
      }
      
      
      
      
      
      
      
    }
    
    
    
    
    
    
    
    for(i=0;i<n;i++){
      //Lets get the Meanvector that we need first
      MEANVEC[0]=MuE[Doses[i]];
      MEANVEC[1]=MuT[Doses[i]];
      //Sample the effiacy latent parameter first
      
      // (\epsilon_E,\epsilon_T) comes from 
      // A bivariate normal random vector with mean
      // MuE[j], MuT[j] and covariance matrix
      // \sigma_\epsilon_E^2 + \tau^2    \tau^2 \rho 
      // \tau^2 \rho       \sigma_\epsilon_T^2 +\tau^2
      
      //Correlation of this bivariate normal
      corr=SIGMA(0,1)/pow(SIGMA(0,0)*SIGMA(1,1),.5);
      
      
      //Get the conditional distribution
      signorm = pow((1-corr)*SIGMA(0,0),.5);
      //Conditional mean
      munorm = MEANVEC[0] + corr*(epsilonT[i]-MEANVEC[1])*pow(SIGMA(0,0)/SIGMA(1,1),.5);
      
      
      
      
      //Generate new variable, where it's truncated
      if(YE[i]>0){
        //Truncate this above 0
        epsilonE[i] = max1(.001,as_scalar(arma::randn(1))*signorm+munorm);
        
      }else{
        //Truncate this above 0
        epsilonE[i] = min1(-.001,as_scalar(arma::randn(1))*signorm+munorm);
        
      }
      
      
      
      
      //Sample the toxicity latent parameter next
      //Get the conditional distribution
      
      
      
      
      
      //Get the conditional distribution
      signorm = pow((1-corr)*SIGMA(1,1),.5);
      //Conditional mean
      munorm = MEANVEC[1] + corr*(epsilonE[i]-MEANVEC[0])*pow(SIGMA(1,1)/SIGMA(0,0),.5);
      
      
      //Generate new variable, where it's truncated
      if(YT[i]>0){
        //Truncate this above 0
        epsilonT[i] = max1(.001,as_scalar(arma::randn(1))*signorm+munorm);
        
      }else{
        //Truncate this above 0
        epsilonT[i] = min1(-.001,as_scalar(arma::randn(1))*signorm+munorm);
        
      }
      
      
      
      
      
      
      
      
      
      
    }
    
    
    
    
    
    
    
    
 
    phinew=phi + as_scalar(arma::randn(1))*crho;
    
    z9[0]=phi;
    z9[1]=phinew;
    
// Rf_PrintValue(wrap(z9));
    
    rhonew=2*exp(phinew)/(1+exp(phinew))-1;

    
    
    SIGMANEW(0,0)=SigmaE +tau;
    SIGMANEW(0,1)=rhonew*tau;
    SIGMANEW(1,0)=rhonew*tau;
    SIGMANEW(1,1)=SigmaT +tau;
    SIGMANEWINV=inv(SIGMANEW);
    
    
    //Compute the loglikelihood:
    
    alpha = 0;
    
    
    for(i=0;i<n;i++){
      MEANVEC[0]=MuE[Doses[i]];
      MEANVEC[1]=MuT[Doses[i]];
      
      epVEC[0]=epsilonE[i];
      epVEC[1]=epsilonT[i];
      
      alpha = alpha - as_scalar(.5*(epVEC-MEANVEC).t()*SIGMANEWINV*(epVEC-MEANVEC)) + as_scalar(.5*(epVEC-MEANVEC).t()*SIGMAINV*(epVEC-MEANVEC));
      //Determinant
      alpha = alpha - .5*log(abs1(SIGMANEW(1,1)*SIGMANEW(0,0)-SIGMANEW(0,1)*SIGMANEW(1,0))) + .5*log(abs1(SIGMA(1,1)*SIGMA(0,0)-SIGMA(0,1)*SIGMA(1,0))) ;
      
      alpha=alpha + .5*pow(phi,2) - .5*pow(phinew,2);
      
      
    }
    
    
    
    
    U=log(as_scalar(arma::randu(1)));
    
    
    
    if(U<alpha){
      
      SIGMA=SIGMANEW;
      SIGMAINV=SIGMANEWINV;
      rho=rhonew;
      phi=phinew;
      arho=arho+1;
    }
    nrho=nrho+1;

    
    
    
    
    
    
    
    
    
    
    if(b>(B-B1-1)){
      
      StoreInx = b-B1;
      
      for(j=0;j<nDose;j++){
        MuEStore(StoreInx,j)=MuE[j];
        MuTStore(StoreInx,j)=MuT[j];
      }

 
 for(i=0; i<n; i++){
  epsilonTStore(StoreInx,i) = epsilonT[i];
   epsilonEStore(StoreInx,i) = epsilonE[i];
   
 }
 
 rhoStore(StoreInx)=rho;
 
 
      
    }
    
    
    
  }
  
  

  
  
  //Return the two storage matrices and discard the rest
  List z1 = List::create(MuEStore,MuTStore,rhoStore);
  return(z1);
}
