
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
#include <tgmath.h>
#include <complex.h>
#include <cmath>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]



using namespace Rcpp;



double log1(int Y){
  double Y1=Y;
  return(log(Y1));
}




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


double abs1(double a){
  
  if(a>=0){
    return(a);
  }else{
    return(-a);
  }
  
}


arma::vec GetSlopePLLH(arma::vec s, //split point locations
                       arma::vec lam, //Log-hazard heights
                       int J //Number of split points
){
  
  int m;
  
  //Make Slope
  arma::vec slopes(J+1);
  for(m=0;m<slopes.n_rows;m++){
    slopes(m)=(lam(m+1)-lam(m))/(s(m+1)-s(m));
  }
  
  
  return(slopes);
  
  
}








double LikePLLHTrt( arma::vec Y, //Survival Times
                    arma::vec I1, //Censoring Indicators
                    arma::vec Trt, //Treatment Indicators
                    arma::vec s, //Vector containing the split point locations
                    arma::vec lam, //Vector containing the log hazard heights
                    int J, //Numer of split points
                    double Beta // Beta for hazard multiplication
){
  
  
  
  
  int m=0;
  int l=0;
  double LogL=0;
  double Y1=0;
  //Make Slope
  
  NumericVector z9(4);
  
  
  arma::vec slopes = GetSlopePLLH(s,lam,J);
  //Cycle through each interval and obtain survival estimates and add hazard heights
  for(l=0; l<(J+1);l++){
    for(m=0; m<Y.n_rows;m++){
      //Y1 here contains the max time that a patient was in interval l
      Y1=min1(s(l+1),Y(m));
      if(Y1>s(l)){
        //Means patient has survived past time s_l
        //This is the survival contribution for patient i in interval l
        
        LogL = LogL + exp(Trt[m]*Beta)*exp(lam(l))*(1-exp(slopes(l)*(Y1-s(l))))/slopes(l);
        
        //If this holds below, then the patient died in the interval [s_l, s_{l+1}]
        if((Y1<s(l+1)) && (I1(m)==1)){
          //Patient died in interval, need hazard here
          LogL = LogL + lam(l)+(Y1-s(l))*slopes(l) + Trt(m)*Beta;
        }
        
        
      }
      
      
      
      
    }
    
    
  }
  
  
  
  return(LogL);
  
  
  
}

double LikePLLHCOV( arma::vec Y, //Survival Times
                    arma::vec I1, //Censoring Indicators
                    arma::mat COV, //Treatment Indicators
                    arma::vec s, //Vector containing the split point locations
                    arma::vec lam, //Vector containing the log hazard heights
                    int J, //Numer of split points
                    arma::vec Beta // Beta for hazard multiplication
){
  
  
  
  
  int m=0;
  int l=0;
  double LogL=0;
  double Y1=0;
  //Make Slope
  
  NumericVector z9(4);
  arma::vec eta=COV*Beta;
  
  
  arma::vec slopes = GetSlopePLLH(s,lam,J);
  
  
  //Cycle through each interval and obtain survival estimates and add hazard heights
  for(l=0; l<(J+1);l++){
    for(m=0; m<Y.n_rows;m++){
      
      
      
      //Y1 here contains the max time that a patient was in interval l
      Y1=min1(s(l+1),Y(m));
      if(Y1>s(l)){
        //Means patient has survived past time s_l
        //This is the survival contribution for patient i in interval l
        
        
        
        LogL = LogL + exp(eta[m])*exp(lam(l))*(1-exp(slopes(l)*(Y1-s(l))))/slopes(l);
        
        
        
        
        
        //If this holds below, then the patient died in the interval [s_l, s_{l+1}]
        if((Y1<s(l+1)) && (I1(m)==1)){
          
          //Patient died in interval, need hazard here
          LogL = LogL + lam(l)+(Y1-s(l))*slopes(l) + eta(m);
        }
        
        
      }
      
      
      
      
    }
    
    
  }
  
  
  
  return(LogL);
  
  
  
}














//' Returns the approximate restricted posterior mean survival for the PLLH model.
//'
//' Uses a grid and parameter values to approximate the restricted posterior mean survival for the PLLH model using the integral of the survival function.
//' @param Y Sequence from 0.01 to the maximum observed event time used to compute the approximate restricted mean survival time. Smaller spaced sequences results in better approximation but longer computation time.
//' @param s Vector of split points. The first and last entries must be 0 and max(Y).
//' @param lam Vector of log-hazard values at each split point location. Must be same length as s.
//' @param J Number of split points. 
//' @return Returns the approximate restricted posterior mean survival time for the PLLH model.
//' @importFrom Rcpp evalCpp
//' @useDynLib BayesReversePLLH
//' @examples
//' ##Generate Data
//' Y1=rweibull(100,4,1)
//' ##Create sequence from (0,max(Y1)) for approximation
//' Y=seq(.01,max(Y1),.01)
//' ##Parameters used to approximate the mean
//' s=c(0,1,max(Y1))
//' lam=c(-2,0,-2)
//' J=1
//' ApproxMean( Y, s, lam, J)
//' @export
//[[Rcpp::export]]
double ApproxMean(arma::vec Y , //Time points for approximation
                  arma::vec s, //Vector containing the split point locations
                  arma::vec lam, //Vector containing the log hazard heights
                  int J //Numer of split points
){
  
  NumericVector z9(2);
  
  int m=0;
  int l=0;
  double mean=0;
  double Y1=0;
  //Make Slope
  arma::vec slopes = GetSlopePLLH(s,lam,J);
  double cum1 = 0;
  //Cycle through each interval and obtain survival estimates and add hazard heights
  for(m=0; m<Y.n_rows;m++){
    cum1=0;
    
    for(l=0; l<(J+1);l++){
      //Cycle over the intervals within each observation 
      //Y1 here contains the max time that a patient was in interval l
      Y1=min1(s(l+1),Y(m));
      if(Y1>s(l)){
        //Means patient has survived past time s_l
        //This is the survival contribution for patient i in interval l
        
        //       cum1 = cum1 - exp(lam(l) - slopes(l)*s(l))*(exp(slopes(l)*Y1)-exp(slopes(l)*s(l)))/slopes(l);
        
        cum1 = cum1 + exp(lam(l))*(1-exp(slopes(l)*(Y1-s(l))))/slopes(l);
        
        
      }
      
    }
    
    
    mean = mean + exp(cum1);
    
    
  }
  
  
  
  
  //Returns Approx Mean
  return(mean*(Y(1)-Y(0)));
  
  
  
}





//[[Rcpp::export]]
arma::vec SurvPLLH(arma::vec Y , //Time points for Getting Survival Estimates
                   arma::vec s, //Vector containing the split point locations
                   arma::vec lam, //Vector containing the log hazard heights
                   int J //Numer of split points
){
  
  NumericVector z9(2);
  arma::vec SurvHold=Y;
  
  
  
  int m=0;
  int l=0;
  double mean=0;
  double Y1=0;
  //Make Slope
  arma::vec slopes = GetSlopePLLH(s,lam,J);
  double cum1 = 0;
  //Cycle through each interval and obtain survival estimates and add hazard heights
  for(m=0; m<Y.n_rows;m++){
    cum1=0;
    
    for(l=0; l<(J+1);l++){
      //Cycle over the intervals within each observation 
      //Y1 here contains the max time that a patient was in interval l
      Y1=min1(s(l+1),Y(m));
      if(Y1>s(l)){
        //Means patient has survived past time s_l
        //This is the survival contribution for patient i in interval l
        
        //       cum1 = cum1 - exp(lam(l) - slopes(l)*s(l))*(exp(slopes(l)*Y1)-exp(slopes(l)*s(l)))/slopes(l);
        
        cum1 = cum1 + exp(lam(l))*(1-exp(slopes(l)*(Y1-s(l))))/slopes(l);
        
        
      }
      
    }
    
    
    SurvHold(m)=exp(cum1);
    
    
  }
  
  
  
  
  //Returns Approx Mean
  return(SurvHold);
  
  
  
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


int SampleBirth(arma::vec s,int L){
  //Vector to store cumulative probility of being in an interval or lower
  arma::vec cumprob(L+1);
  int m=0;
  
  cumprob[0] = s[1]/s[L+1];
  
  for(m=1;m<cumprob.n_rows;m++){
    cumprob[m]=s[m+1]/s[L+1];
  }
  
  
  //Randomly Draw Uniform that we'll use to determine what interval to add a move to.
  double U=as_scalar(arma::randu(1));
  
  // Find which interval the random unif is in.
  int Which1 =0;
  if(U<cumprob[0]){
    Which1=0;
  }else{
    
    for(m=1;m<cumprob.n_rows;m++){
      if((U>cumprob[m-1]) & (U<cumprob[m])){
        Which1=m;
      }
    }
    
  }
  
  return(Which1);
  
}


double factorial(int n)
{
  
  return(tgamma(n));
  
}

















//' Samples from the PEH Cox model with a treatment indicator.
//'
//' Samples from the Piecewise Exponential Hazard (PEH) Cox model with a treatment indicator and returns a list containing posterior parameters and posterior restricted mean survival.
//' @param Y Vector of event or censoring times.
//' @param I1 Vector of event indicators.
//' @param Trt Vector containing patient treatment/control assignment.
//' @param Poi Prior mean number of split points.
//' @param B Number of iterations for MCMC.
//' @return Returns a list containing posterior samples of (1) the split point locations, (2) the log-hazards at each split point, (3) the number of split points, (4) the variance parameter for the log-hazard values, (5) the treatment coefficient, (6) the mean restricted survivial time of the control therapy, (7) the mean restricted survival time of the treatment therapy.
//' @examples
//' ##Generate Data
//' Y=rweibull(20,4,1)
//' I=rbinom(20,1,.5)
//' Trt=rbinom(20,1,.5)
//' ##Hyperparameter for number of split points
//' Poi=5
//'##Number of iterations for MCMC
//'B=200
//'BayesPiecewiseLinearLogHazardTrt( Y, I,Trt, Poi,  B)
//'@export
//[[Rcpp::export]]
List BayesPiecewiseLinearLogHazardTrt( arma::vec Y, //Survival Times
                                       arma::vec I1, //Censoring Indicators
                                       arma::vec Trt, //Treatment Indicators
                                       double Poi, //Prior Number of split points
                                       int B // Number of iterations to perform
){
  
  //Upper Limit on Accept
  double High=.4;
  //Lower Limit on Accept
  double Low=.2;
  
  
  //Prior variance on first log-hazard height
  double Var1=25;
  //Prior probability of slab
  double PSlab = .9;
  //Prior probability of spike
  double PSpike=PSlab/(1-PSlab);
  
  
  //Prior Params, make user controlled later
  //For loop quantiites
  int m=0;
  int b=0;
  int k=0;
  int j=0;
  //Important double quantities used in MCMC
  //  int B1=B/10;  //Burnin first 90 %
  //  int B1=B/10;
  int B1=B/2;
  
  
  
  //Storing the restricted mean for control group
  arma::vec MeanStore1(B1);
  //Stores the restricted mean for treatment group.
  arma::vec MeanStore2(B1);
  
  
  MeanStore1.zeros();
  MeanStore2.zeros();
  
  //Max To Store, really is 1 less here
  int Lmax = 30;
  //Storage for S vector
  arma::mat sstore(B1,Lmax+2);
  //Storage for lambda
  arma::mat lamstore(B1,Lmax+2);
  //Storage for L
  arma::vec Lstore(B1);
  //Storage for sigma (dependence)
  arma::vec sigstore(B1);
  //Initialize S, lam and J
  //Make Vector of split points
  arma::vec s(Lmax+2);
  arma::vec lam(Lmax+2);
  s.zeros();
  lam.zeros();
  //Initialize sigma_\lambda and 
  double sig=1;
  //Used for birth moves
  double Birth=0;
  int Spot=0;
  //Used for metropolis hastings moves
  double U=0;
  double alpha=0;
  //This needed for siglam sampling
  double cum1=0;
  //multiplicative perturbation
  double U1=0;
  //used to find what interval something is in
  int which1=0;  
  
  
  
  //Copy Paste these to top
  double mean2=0;
  
  
  
  
  
  double mean1=0;
  double Con=0;
  
  
  
  //Storage for printout
  NumericVector z9(2);
  //Contains the maximum observed event time
  double m1=MaxVec(Y) ;
  //Start with split points at every $L/4$
  s[1]=m1/4;
  s[2]=2*m1/4;
  s[3]=3*m1/4;
  s[4]=m1;
  //Innitialize L with 3 split points
  int L = 3;
  
  
  //Make a sequence x1 that we'll use to get approximate means
  double max12 = ceil(m1*100)/100;
  //Seq size
  int size1 = max12/.01;
  
  
  
  arma::vec x1(size1);
  for(m=0;m<size1;m++){
    x1(m)=.01*(m+1);
  }
  
  
  
  
  
  //Make A sequence for 
  
  
  
  //Now Our initial split point vector is set
  lam(0)=as_scalar(arma::randn(1));
  lam(1)=as_scalar(arma::randn(1));
  lam(2)=as_scalar(arma::randn(1));
  lam(3)=as_scalar(arma::randn(1));
  lam(4)=as_scalar(arma::randn(1));
  
  
  //Double for adding new lambda
  double NewLam1=0;
  double NewLam2=0;
  
  
  
  
  
  //This piece vector contains indicators for if \lambda_l = \lambda_{l-1}
  arma::vec Piece=lam;
  //Sets the piece vector to 0s
  Piece.zeros();
  arma::vec lam1=lam;
  //Sampling counters
  double bprop=1;
  double Intb=1;
  double Numb=2;
  arma::vec lamprop=lam;
  arma::vec sprop=s;
  
  double prob1=0;
  double med2=0;
  
  
  
  //Holds storage index
  int StoreInx=0;
  
  //Lambda Tuning Parameters
  double LamC = .5;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam = 1;  //Number of proposals made
  double ALam = 1; //Number of proposals accepted
  //Here's for lam0, which ALWAYS exists
  double LamC1 = .5;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam1 = 1;  //Number of proposals made
  double ALam1 = 1; //Number of proposals accepted
  //Here's for lam_J+1 which ALWAYS exists
  double LamC2 = .5;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam2 = 1;  //Number of proposals made
  double ALam2 = 1; //Number of proposals accepted
  
  
  //Do these for beta now
  double BetaC = .02; //Beta coefficient for changing value, really starts at .25
  double NBeta = 1;  //Number of Beta proposals made
  double ABeta = 1; //Number of Beta proposals accepted
  
  //Storage
  arma::vec BetaStore(B1);
  
  
  
  
  //Beta Starting Value
  double Beta = as_scalar(arma::randn(1));
  double BetaProp=0;
  
  
  
  for(m=0;m<B;m++){
    
    //Every 250 observations, let's adjust variance until burnin
    if( (m%250==0) & (m<(B-B1))){
      //Tune proposal variance for lambda here
      if((ALam/NLam)>High){
        //Double the variance, acceptance is too high
        LamC = min1(LamC*2,2);
        
      }
      
      
      if((ALam/NLam)<Low){
        //Halve the variance, acceptance is too Low
        LamC = max1(LamC/2,.0625);
        
      }
      
      //Reset counters
      ALam =1;
      NLam=1;
      
      
      
      
      
      //Tune proposal variance for lambda here
      if((ALam1/NLam1)>High){
        //Double the variance, acceptance is too high
        LamC1 = min1(LamC1*2,2);
        
      }
      
      
      if((ALam1/NLam1)<Low){
        //Halve the variance, acceptance is too Low
        LamC1 = max1(LamC1/2,.0625);
        
      }
      
      
      
      
      
      //Reset counters
      ALam1 =1;
      NLam1=1;
      
      
      
      //Tune proposal variance for lambda here
      if((ALam2/NLam2)>High){
        //Double the variance, acceptance is too high
        LamC2 = min1(LamC2*2,2);
        
      }
      
      
      if((ALam2/NLam2)<Low){
        //Halve the variance, acceptance is too Low
        LamC2 = max1(LamC2/2,.0625);
        
      }
      
      
      //Reset counters
      ALam2 =1;
      NLam2=1;
      
      
      
      //Tune Proposals For \beta
      
      if((ABeta/NBeta)>High){
        //Double the variance, acceptance is too high
        BetaC = min1(BetaC*2,.16);
        
      }
      
      
      if((ABeta/NBeta)<Low){
        //Halve the variance, acceptance is too Low
        BetaC = max1(BetaC/2,.005);
        
      }
      
      //Reset counters
      ABeta =1;
      NBeta=1;
      
      
      
      
      
      
    }
    
    
    
    
    //Sample Lambda_{0}
    lamprop = lam;
    //Default for now, could make it adaptive
    lamprop(0)=lam(0)+as_scalar(arma::randu(1))*LamC1*2 - LamC1;
    //Log-Acceptance ratio
    alpha=LikePLLHTrt(Y,I1, Trt, s, lamprop , L,Beta) -    LikePLLHTrt(Y,I1, Trt, s, lam , L,Beta);
    //adjust for proposal ratio
    alpha = alpha - .5*pow(lamprop[0],2)/Var1 + .5*pow(lam[0],2)/Var1 - .5*pow(lamprop[1]-lamprop[0],2)/sig + .5*pow(lam[1]-lam[0],2)/sig;
    
    
    
    U=log(as_scalar(arma::randu(1)));
    if(U<alpha){
      lam(0)=lamprop(0);
      ALam1=ALam1+1;
    }
    NLam1=NLam1+1;
    
    
    
    
    
    //Sample Lambda_{J+1}
    lamprop = lam;
    //Default for now, could make it adaptive
    lamprop(L+1)=lam(L+1)+as_scalar(arma::randu(1))*LamC2*2 - LamC2;
    //Log-Acceptance ratio
    alpha=LikePLLHTrt(Y,I1, Trt, s, lamprop , L,Beta) -    LikePLLHTrt(Y,I1, Trt, s, lam , L,Beta);
    //Adjust for proposal ratio
    alpha = alpha - .5*pow(lamprop[L+1]-lamprop[L],2)/sig + .5*pow(lam[L+1]-lam[L],2)/sig   ;
    
    U=log(as_scalar(arma::randu(1)));
    if(U<alpha){
      lam(L+1)=lamprop(L+1);
      ALam2=ALam2+1;
    }
    NLam2=NLam2+1;
    
    
    
    
    //Sample Lambda | L>0, Here's the middle ones
    if(L>0){
      for(j=1;j<(L+1);j++){
        lamprop = lam;
        //Default for now, could make it adaptive
        lamprop(j)=lam(j)+as_scalar(arma::randu(1))*LamC*2 - LamC;
        //Log-Acceptance ratio
        alpha=LikePLLHTrt(Y,I1, Trt, s, lamprop , L, Beta) -    LikePLLHTrt(Y,I1, Trt, s, lam , L, Beta);
        //Adjust for proposal ratio
        alpha = alpha - .5*pow(lamprop[j]-lamprop[j-1],2)/sig + .5*pow(lam[j]-lam[j-1],2)/sig  - .5*pow(lamprop[j+1]-lamprop[j],2)/sig + .5*pow(lam[j+1]-lam[j],2)/sig   ;
        U=log(as_scalar(arma::randu(1)));
        if(U<alpha){
          lam(j)=lamprop(j);
          ALam=ALam+1;
        }
        NLam=NLam+1;
        
        
      }
    }
    
    
    
    
    
    
    
    
    //Sample Siglam directly
    //Set sum of quadratic lambda differences to 0
    cum1=0;
    //Do a for loop over the intervals for (lambda_j-lambda_{j-1})^2
    for(j=1;j<(L+2);j++){
      cum1 =cum1 + pow(lam[j]-lam[j-1],2);
    }
    //Sample sigma_lambda directly from a Gamma distribution
    sig = 1/R::rgamma((L+2)/2 , cum1/2  );
    
    
    
    //Shuffle Splitpoints
    if(L>0){
      for(j=1; j<(L+1); j++){
        sprop = s;
        //Randomly draw s_j \sim U[s_{j-1},s_{j+1}]
        sprop[j]=as_scalar(arma::randu(1))*(sprop[j+1]-sprop[j-1]) + sprop[j-1];
        //Calculate Acceptance Ratio
        alpha =  log(sprop[j+1]-sprop[j])+log(sprop[j]-sprop[j-1])-log(s[j+1]-s[j])-log(s[j]-s[j-1]) + 
          LikePLLHTrt(Y,I1, Trt,  sprop, lam , L, Beta) -    LikePLLHTrt(Y,I1, Trt,s, lam , L, Beta);
        //Draw Uniform 
        U = log(as_scalar(arma::randu(1)));
        if(U<alpha){
          //Metropolis Hastings
          s[j]=sprop[j];
        }
        
        
        
        
      }
    }
    
    
    //Death move
    if(L>0){
      //Draw proposed split point to delete
      Spot=SampleDeath(L);
      sprop.zeros();
      lamprop.zeros();
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
      }
      
      //We just straight up delete \lamda_spot and s[spot] and make adjustments
      if(L>1){
        for(j=0;j<Spot;j++){
          lamprop[j]=lam[j];
        }
        for(j=Spot;j<(lam.n_rows-1);j++){
          lamprop[j]=lam[j+1];
        }
      }else{
        lamprop[0]=lam[0];
        lamprop[1]=lam[2];
      }
      for(j=Spot; j<(s.n_rows-1);j++){
        sprop[j]=s[j+1];
      }
      //Likelihood ratio
      alpha =    LikePLLHTrt(Y,I1, Trt,  sprop, lamprop , L-1,Beta) -    LikePLLHTrt(Y,I1, Trt, s, lam , L,Beta);
      //Prior Ratio
      //Poisson
      alpha = alpha  -log(Poi) + log1(L);
      //S Prior
      alpha = alpha +2*log(m1) + log(s[Spot+1]-s[Spot-1]) - log1(2*L+1) - log1(2L) - log(s[Spot]-s[Spot-1])-log(s[Spot+1]-s[Spot]);
      //Lambda Prior, we DROPPED one
      if(L==1){
        //Removing will drop the sampler to 0 split points
        alpha = alpha - .5*pow(lam[2]-lam[0],2)/sig + .5*pow(lam[2]-lam[1],2)/sig + .5*pow(lam[1]-lam[0],2)/sig ;
      }else{
        alpha = alpha -.5*pow(lam[Spot+1]-lam[Spot-1],2)/sig + .5*pow(lam[Spot+1]-lam[Spot],2)/sig+.5*pow(lam[Spot]-lam[Spot-1],2)/sig;
      }
      
      
      //Random Perturbation to maintain balance between parameter spaces
      U1 = as_scalar(arma::randu(1));
      ///Add log Jacobian Here
      alpha = alpha +log(U1*(1-U1));
      //Metropolis Hastings Green
      U=log(as_scalar(arma::randu(1)));
      
      if(U<alpha){
        s=sprop;
        L=L-1;
        lam=lamprop;
      }
    }
    
    
    
    //Birth move Here
    if(L<(Lmax-1)){
      //What interval is it in?
      //This generates the interval location that the new proposed split point is located
      Spot = SampleBirth(s,L)+1;
      //Set Lambda to 0
      lamprop.zeros();
      //Now generate the new proposed split point location
      Birth = as_scalar(arma::randu(1))*(s[Spot]-s[Spot-1])+s[Spot-1];
      //Random Perturbation for dimension matching
      U1=as_scalar(arma::randu(1));
      //Slope
      
      //Old Slope for log hazard
      NewLam1= (lam[Spot]-lam[Spot-1])/(s[Spot]-s[Spot-1]);
      //Compute the new slope that we can use to get the new lambda value
      //This quantity comes from the MCMC section on how to obtain new slopes
      NewLam1 = NewLam1 - log((1-U1)/U1)*(s[Spot]-Birth)/(s[Spot]-s[Spot-1]);
      //We can now use this slope to obtain the new hazard height at the time point Birth
      NewLam2 = NewLam1*(Birth-s[Spot-1])+lam[Spot-1];
      //Now NewLam2 contains the hazard at the time point Birth
      //Now we have the interval of the new split in Spot and the actual location in Birth
      sprop.zeros();
      //Here we're going to fill in our new vector
      //Let's add Birth to the Spot location and push back the rest of sprop
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
        lamprop[j]=lam[j];
      }
      sprop[Spot]=Birth;
      for(j=(Spot+1);j<(s.n_rows-1);j++){
        sprop[j]=s[j-1];
      }
      //Log of hazard height at Birth.
      lamprop[Spot]=NewLam2;
      for(j=(Spot+1);j<(lam.n_rows-1);j++){
        lamprop[j]=lam[j-1];
      }
      
      
      //Now we have our new proposal vector, evaluate it!
      //Like Ratio
      alpha =    LikePLLHTrt(Y,I1, Trt, sprop, lam , L+1, Beta) -    LikePLLHTrt(Y,I1, Trt, s, lam , L, Beta);
      //Add proposal ratio
      //Poisson
      alpha= alpha + log(Poi) - log1(L+1) ;
      // S proposal
      alpha = alpha + log1(2*L+3)+log1(2*L+2)+log(Birth-s[Spot-1])+log(s[Spot]-Birth) - 2*log(m1)-log(s[Spot]-s[Spot-1]);
      //Add proposal ratio for \lambda
      if(L==0){
        alpha = alpha  - .5*pow(lamprop[1]-lamprop[0],2)/sig - .5*pow(lamprop[2]-lamprop[1],2)/sig + .5*pow(lam[1]-lam[0],2)/sig ;
      }else{
        alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig + .5*pow(lam[Spot]-lam[Spot-1],2)/sig;
      }
      
      //Add log Jacobian
      alpha = alpha -log(U1*(1-U1));
      //Uniform for Metropolis-Hastings
      U=log(as_scalar(arma::randu(1)));
      if(U<alpha){
        //We have to set the entire vector here
        s=sprop;
        lam=lamprop;
        L=L+1;
      }
    }
    
    
    
    
    //Beta Sampler
    
    BetaProp = Beta  +as_scalar(arma::randu(1))*BetaC*2 - BetaC;
    
    
    
    //Default for now, could make it adaptive
    //Log-Acceptance ratio
    alpha=LikePLLHTrt(Y,I1, Trt, s, lam , L, BetaProp) -    LikePLLHTrt(Y,I1, Trt, s, lam , L, Beta);
    //Adjust for proposal ratio
    alpha = alpha - .5*pow(BetaProp,2)/25 +.5*pow(Beta,2)/25    ;
    U=log(as_scalar(arma::randu(1)));
    if(U<alpha){
      Beta=BetaProp;
      ABeta=ABeta+1;
    }
    NBeta=NBeta+1;
    
    
    
    
    
    
    
    if(m>(B-B1-1)){
      //Store Values in Matrix only after certain burnin
      StoreInx = m-B+B1;
      
      for(j=0;j<sstore.n_cols;j++){
        sstore(StoreInx,j) = s(j);
      }
      for(j=0;j<lamstore.n_cols;j++){
        lamstore(StoreInx,j) = lam(j);
      }
      Lstore[StoreInx]=L;
      //  sigstore[StoreInx]=sig;
      sigstore[StoreInx]=sig;
      
      BetaStore[StoreInx]=Beta;
      
      //Can compute the mean here.
      MeanStore1(StoreInx)=ApproxMean(x1,s,lam,L);
      MeanStore2(StoreInx)=ApproxMean(x1,s,lam+Beta,L);
      
      
      
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
  
  
  
  
  
  
  
  
  List z1 = List::create(sstore,lamstore,Lstore,sigstore,BetaStore,MeanStore1,MeanStore2);
  
  
  return(z1);
  
  
}














//Function for computing the likelihood of Piecewise Exponential Hazard
double LikePEHTrt( arma::vec Y, //Survival Times
                   arma::vec I1, //Censoring Indicators
                   arma::vec Trt,
                   arma::vec s, //Vector containing the split point locations
                   arma::vec lam, //Vector containing the log hazard heights
                   int J,
                   double Beta){
  
  int m=0;
  int l=0;
  double LogL=0;
  
  
  //Cycle through each interval and obtain survival estimates and add hazard heights
  for(l=0; l<(J+1);l++){
    for(m=0; m<Y.n_rows;m++){
      LogL = LogL - max1(0,min1(s(l+1),Y(m))-s(l))*exp(lam[l]+Beta*Trt[m]);
      
      
      if((Y(m)>s[l]) & (Y(m)<=s[l+1]) & (I1[m]==1)){
        LogL = LogL + lam[l] + Beta*Trt[m];
      }
      
      
      
      
    }
    
  }
  
  
  
  
  return(LogL);  
  
  
  
}



//Function for computing the likelihood of Piecewise Exponential Hazard
double LikePEHCOV( arma::vec Y, //Survival Times
                   arma::vec I1, //Censoring Indicators
                   arma::mat COV,
                   arma::vec s, //Vector containing the split point locations
                   arma::vec lam, //Vector containing the log hazard heights
                   int J,
                   arma::vec Beta){
  
  int m=0;
  int l=0;
  double LogL=0;
  
  arma::vec eta=COV*Beta;
  
  //Cycle through each interval and obtain survival estimates and add hazard heights
  for(l=0; l<(J+1);l++){
    for(m=0; m<Y.n_rows;m++){
      LogL = LogL - max1(0,min1(s(l+1),Y(m))-s(l))*exp(lam[l]+eta[m]);
      
      
      if((Y(m)>s[l]) & (Y(m)<=s[l+1]) & (I1[m]==1)){
        LogL = LogL + lam[l] + eta[m];
      }
      
      
      
      
    }
    
  }
  
  
  
  
  return(LogL);  
  
  
  
}






//' Samples from the PEH Cox model with a patient covariate vector.
//'
//' Samples from the Piecewise Linear Log-Hazard (PLLH) Cox model and returns a list containing posterior parameters and posterior restricted mean survival.
//' @param Y Vector of event or censoring times.
//' @param I1 Vector of event indicators.
//' @param Trt Vector containing patient treatment/control assignment.
//' @param Poi Prior mean number of split points.
//' @param B Number of iterations for MCMC.
//' @return Returns a list containing posterior samples of (1) the split point locations, (2) the log-hazards at each split point, (3) the number of split points, (4) the variance parameter for the log-hazard values, (5) the treatment coefficient, (6) the mean restricted survivial time of the control therapy, (7) the mean restricted survival time of the treatment therapy.
//' @examples
//' ##Generate Data
//' Y=rweibull(20,4,1)
//' I=rbinom(20,1,.5)
//' Trt=rbinom(20,1,.5)
//' ##Hyperparameter for number of split points
//' Poi=5
//'##Number of iterations for MCMC
//'B=200
//'BayesPiecewiseHazardTrt( Y, I,Trt, Poi,  B)
//'@export
//[[Rcpp::export]]
List  BayesPiecewiseHazardTrt( arma::vec Y, //Survival Times
                               arma::vec I1, //Censoring Indicators
                               arma::vec Trt, // Treatment Indicators
                               double Poi, //Prior Mean Number of split points
                               int B ){
  
  
  
  //Prior Params, make user controlled later
  double Var1=25;
  
  
  //
  //Upper Limit on Accept
  double High=.4;
  //Lower Limit on Accept
  double Low=.2;
  
  //For loop quantiites
  int m=0;
  int b=0;
  int k=0;
  int j;
  //Important double quantities used in MCMC
  int B1=B/10;
  
  
  //Innitialize Parameters and storage matrices
  double beta = 1  ;
  
  
  
  arma::vec betastore(B1);
  arma::vec MeanStore1(B1);
  arma::vec MeanStore2(B1);
  
  
  
  arma::vec MedianStore=MeanStore1;
  //Max To Store, really is 1 less here
  int Lmax = 30;
  arma::mat sstore(B1,Lmax+2);
  arma::mat lamstore(B1,Lmax+1);
  arma::vec Lstore(B1);
  arma::vec sigstore(B1);
  //Initialize S, lam and J
  arma::vec s(Lmax+2);
  arma::vec lam(Lmax+1);
  s.zeros();
  lam.zeros();
  double sig=1;
  double Birth=0;
  int Spot=0;
  
  double betaprop=0;
  double U=0;
  double alpha=0;
  double svar=1;
  double Ints=1;
  double Nums=2;
  //Copy Paste these to top
  double mean2=0;
  double cum1=0;
  double mean1=0;
  double Con=0;
  
  
  
  
  NumericVector z9(2);
  
  //Innitialize l with 3 split points
  int L = 3;
  
  //Contains the maximum observed censoring time
  double m1 = MaxVec(Y);
  //Start with split points at every $L/4$
  s[1]=m1/4;
  s[2]=2*m1/4;
  s[3]=3*m1/4;
  s[4]=m1;
  
  double NewLam1=0;
  double NewLam2=0;
  double U1=0;
  int which1=0;
  
  
  arma::vec lam1=lam;
  
  double bprop=1;
  double Intb=1;
  double Numb=2;
  arma::vec lamprop=lam;
  arma::vec sprop=s;
  
  double prob1=0;
  double med2=0;
  
  
  
  
  int StoreInx=0;
  
  
  //Lambda Tuning Parameters
  double LamC = .25;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam = 1;  //Number of proposals made
  double ALam = 1; //Number of proposals accepted
  
  
  
  //Do these for beta now
  double BetaC = .02; //Beta coefficient for changing value, really starts at .25
  double NBeta = 1;  //Number of Beta proposals made
  double ABeta = 1; //Number of Beta proposals accepted
  
  //Storage
  arma::vec BetaStore(B1);
  
  
  
  
  //Beta Starting Value
  double Beta = as_scalar(arma::randn(1));
  double BetaProp=0;
  
  
  
  
  
  for(m=0;m<B;m++){
    
    
    if( (m%250==0) & (m<(B-B1))){
      //Tune proposal variance for lambda here
      if((ALam/NLam)>High){
        //Double the variance, acceptance is too high
        LamC = min1(LamC*2,2);
        
      }
      
      
      if((ALam/NLam)<Low){
        //Halve the variance, acceptance is too Low
        LamC = max1(LamC/2,.0625);
        
      }
      
      //Reset counters
      ALam =1;
      NLam=1;
      
      
      
      
      //Tune Proposals For \beta
      
      if((ABeta/NBeta)>High){
        //Double the variance, acceptance is too high
        BetaC = min1(BetaC*2,.16);
        
      }
      
      
      if((ABeta/NBeta)<Low){
        //Halve the variance, acceptance is too Low
        BetaC = max1(BetaC/2,.005);
        
      }
      
      //Reset counters
      ABeta =1;
      NBeta=1;
      
      //Reset counters
      ABeta =1;
      NBeta=1;
      
      
      
      
      
      
    }
    
    
    
    if(L>0){
      
      for(j=0;j<(L+1);j++){
        
        lamprop = lam;
        //Propose new value for lambda
        lamprop(j)=lam[j]+as_scalar(arma::randu(1))*2*LamC - LamC;
        //Likelihood Ratio
        alpha=LikePEHTrt(Y,I1,Trt,s, lamprop , L,Beta) -    LikePEHTrt(Y,I1,Trt, s, lam , L,Beta);
        //Prior Ratio
        if(j==0){
          alpha = alpha - .5*(pow(lamprop[j],2)-pow(lam[j],2))/Var1 - .5*pow(lamprop[j+1]-lamprop[j],2)/sig + .5*pow(lam[j+1]-lam[j],2);
        }else{
          if(j==L){
            alpha = alpha - .5*(pow(lamprop[j]-lamprop[j-1],2) + .5*pow(lam[j]-lam[j-1],2))/sig    ;
            
          }else{
            alpha = alpha - .5*(pow(lamprop[j]-lamprop[j-1],2) + .5*pow(lam[j]-lam[j-1],2))/sig  - .5*(pow(lamprop[j+1]-lamprop[j],2) + .5*pow(lam[j+1]-lam[j],2))/sig   ;
          }
        }
        
        
        
        //Metropolis Hastings Draw
        U = log(as_scalar(arma::randu(1)));
        //Accept/rej
        if(U<alpha){
          //Accept or Reject
          lam(j)=lamprop(j);
          ALam=ALam+1;
        }
        NLam=NLam+1;
        
      }
      
    }else{
      
      
      lamprop = lam;
      lamprop[0]=lam[0]+as_scalar(arma::randu(1))*2*LamC - LamC;
      //Likelihood Ratio
      alpha=LikePEHTrt(Y,I1, Trt,s, lamprop , L,Beta) -    LikePEHTrt(Y,I1,Trt, s, lam , L,Beta);
      //Prior Ratio Here for lambda 0
      alpha = alpha - .5*pow(lamprop[0],2)/Var1+.5*pow(lam[0],2)/Var1 ;
      U = log(as_scalar(arma::randu(1)));
      
      
      
      
      
      if(U<alpha){
        //accept proposal
        lam(0)=lamprop(0);
        ALam=ALam+1;
        
      }
      
      NLam=NLam+1;
      
    }
    
    //Sample Siglam
    cum1=0;
    if(L>0){
      for(j=1;j<L;j++){
        cum1 =cum1 + pow(lam[j]-lam[j-1],2);
      }
      
      sig = 1/R::rgamma(L/2 + 1, cum1/2 );
      
    }else{
      //Why is this here??
      sig=10000;
    }
    
    
    //Shuffle Splits
    
    if(L>0){
      
      for(j=1; j<(L+1); j++){
        sprop = s;
        //Draw new proposal for s_j, but make sure it's between s_j and s_j+1
        sprop[j]=as_scalar(arma::randu(1))*(sprop[j+1]-sprop[j-1]) + sprop[j-1];
        //Acceptance ratio
        alpha = log(sprop[j+1]-sprop[j])+log(sprop[j]-sprop[j-1])-log(s[j+1]-s[j])+
          log(s[j]-s[j-1]) +   LikePEHTrt(Y,I1,Trt, sprop, lam , L,Beta) -    LikePEHTrt(Y,I1,Trt, s, lam , L,Beta);
        
        //Metropolis Hastings
        U = log(as_scalar(arma::randu(1)));
        
        if(U<alpha){
          s[j]=sprop[j];
        }
        
        
        
      }
      
    }
    
    
    
    //Add proposal and MHG ratios here
    
    //Birth move Here
    if(L<(Lmax-1)){
      //What interval is it in?
      
      lamprop.zeros();
      //Find birth location      
      Spot = SampleBirth(s,L)+1;
      
      //Sample birthed split point
      Birth = as_scalar(arma::randu(1))*(s[Spot]-s[Spot-1])+s[Spot-1]; 
      //Now we have the interval of the new split in Spot and the actual location in Birth
      
      //Random Perturbation for detailed balance when dimension matching
      U1=as_scalar(arma::randu(1));
      
      //Find new lambdas
      NewLam1 =  lam[Spot-1] - log((1-U1)/U1)*(s[Spot+1]-Birth)/(s[Spot+1]-s[Spot-1]);
      NewLam2 =  lam[Spot-1] + log((1-U1)/U1)*(Birth-s[Spot])/(s[Spot+1]-s[Spot-1]);
      
      //Set sprop to 0s so we can fill it
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
      
      //Now we have our new proposal vector, evaluate it!
      //Like Ratio
      alpha =    LikePEHTrt(Y,I1,Trt,sprop, lamprop , L+1,Beta) -    LikePEHTrt(Y,I1,Trt, s, lam , L,Beta);
      //Add proposal ratio
      //Poisson
      alpha= alpha + log(Poi) - log1(L+1) ;
      // S proposal
      alpha = alpha + log1(2*L+3)+log1(2*L+2)+log(Birth-s[Spot-1])+log(s[Spot]-Birth);
      alpha = alpha - 2*log(m1)-log(s[Spot]-s[Spot-1]);
      //log of Jacobian for detailed balance
      alpha=alpha-log(U1*(1-U1)) ;
      
      //Add proposal ratio for \lambda
      if(L==0){
        alpha = alpha - .5*pow(lamprop[1]-lamprop[0],2)/sig   ;
      }else{
        
        if(Spot==(L+1)){
          //We're in the last interval
          alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig + pow(lam[Spot-1]-lam[Spot-2],2)/sig;
        }else{
          
          
          if(Spot==1){
            
            
            
            //Birthed is in the first interval      
            alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig ;
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig ;
            //First interval proposal
            alpha = alpha - .5*(pow(lamprop[0],2)-pow(lam[0],2))/Var1;
          }else{
            //The interval chosen is NOT the first or the last one, so we're ok here.
            alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;
          }  
          
          
        }
      }
      
      
      
      
      
      
      
      
      //Metropolis Hastings      
      U=log(as_scalar(arma::randu(1)));
      
      //Accept/Reject
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
      //Fill in Storage
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
      }
      
      //Fill in lambda proposal vector
      
      if(L>1){
        
        for(j=0;j<(Spot-1);j++){
          lamprop[j]=lam[j];
        }
        
        for(j=Spot;j<(lam.n_rows-1);j++){
          lamprop[j]=lam[j+1];
        }
        
        //New lambda is a weighted average of the old lambdas
        lamprop[Spot-1] = ((s[Spot+1]-s[Spot])*lam[Spot] +(s[Spot]-s[Spot-1])*lam[Spot-1])/(s[Spot+1]-s[Spot-1]);
        
        
      }else{
        lamprop[0]=((s[2]-s[1])*lam[1]+(s[1]-s[0])*lam[0])/(s[2]-s[0]);
      }
      //Finish sproposal fill in
      for(j=Spot; j<(s.n_rows-1);j++){
        sprop[j]=s[j+1];
      }
      
      //Now we have our new proposal vectors, evaluate then!
      //Like Ratio
      alpha =    LikePEHTrt(Y,I1,Trt, sprop, lamprop , L-1,Beta) -    LikePEHTrt(Y,I1,Trt, s, lam , L,Beta);
      //Prior Ratio
      //Poisson
      alpha = alpha  -log(Poi) + log1(L);
      //S Prior
      alpha = alpha + 2*log(m1) + log(s[Spot+1]-s[Spot-1]) - log1(2*L+1) - log1(2L) - log(s[Spot]-s[Spot-1])-log(s[Spot+1]-s[Spot]);
      //Lambda Prior, we DROPPED one
      if(L==1){
        //Removing will drop the sampler to 0 split points
        alpha = alpha + .5*pow(lam[1]-lam[0],2)/sig ;
      }else{
        
        
        
        if(Spot>1){
          
          if(Spot==(L+1)){
            //We dropped the last interval, so we can't evaluate intervals past this
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig  - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;
          }else{
            
            //Split point is not first or last
            alpha = alpha +  .5*pow(lam[Spot+1]-lam[Spot],2)/sig +  .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;
            alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig ;
            
          }
          
          
        }else{
          
          //First Interval is changing so we must evaluate it
          alpha = alpha + .5*pow(lam[Spot+1]-lam[Spot],2)/sig + .5*pow(lam[Spot]-lam[Spot-1],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig  ;
          //Add the 0 entry here
          alpha= alpha  - .5*pow(lamprop[0],2)/Var1 + .5*pow(lam[0],2)/Var1 ;
          
          
        }
        
        
        
        
        
      }
      
      
      
      
      
      
      //Random Perturbation for dimension matching
      U1 = as_scalar(arma::randu(1));
      alpha = alpha + log(U1)+log(1-U1);
      
      //Metropolis Draw
      U=log(as_scalar(arma::randu(1)));
      //Accept/Reject
      if(U<alpha){
        s=sprop;
        L=L-1;
        lam=lamprop;
      }
      
      
      
      
    }
    
    
    
    
    //Beta Sampler
    
    BetaProp = Beta  +as_scalar(arma::randu(1))*BetaC*2 - BetaC;
    //Default for now, could make it adaptive
    //Log-Acceptance ratio
    alpha=LikePEHTrt(Y,I1, Trt, s, lam , L, BetaProp) -    LikePEHTrt(Y,I1, Trt, s, lam , L, Beta);
    //Adjust for proposal ratio
    alpha = alpha - .5*pow(BetaProp,2)/25 +.5*pow(Beta,2)/25    ;
    U=log(as_scalar(arma::randu(1)));
    if(U<alpha){
      Beta=BetaProp;
      ABeta=ABeta+1;
    }
    NBeta=NBeta+1;
    
    
    
    
    
    
    if(m>(B-B1-1)){
      //Store Values in Matrix
      StoreInx = m-B+B1;
      
      for(j=0;j<sstore.n_cols;j++){
        sstore(StoreInx,j) = s(j);
      }
      
      for(j=0;j<lamstore.n_cols;j++){
        lamstore(StoreInx,j) = lam(j);
      }
      
      Lstore[StoreInx]=L;
      sigstore[StoreInx]=sig;
      BetaStore[StoreInx]=Beta;
      //Calculate Means Here
      lam1=exp(lam);
      mean2=0;
      cum1=0;
      
      if(L>0){
        for(k=0;k<(L+1);k++){
          
          mean2= mean2+ exp(-cum1+lam1[k]*s[k])*(exp(-lam1[k]*s[k])-exp(-lam1[k]*s[k+1]))/(lam1[k]);
          
          cum1 = cum1+ lam1[k]*(s[k+1]-s[k]);
          
          
          
        }
      }else{
        k=0;
        
        mean2 = (1-exp(-lam1[0]*s[1]))/lam1[0];
        
      }
      
      
      
      MeanStore1(StoreInx)=mean2;
      
      
      
      
      //Calculate Means Here
      lam1=exp(lam + Beta);
      mean2=0;
      cum1=0;
      
      if(L>0){
        for(k=0;k<(L+1);k++){
          
          mean2= mean2+ exp(-cum1+lam1[k]*s[k])*(exp(-lam1[k]*s[k])-exp(-lam1[k]*s[k+1]))/(lam1[k]);
          
          cum1 = cum1+ lam1[k]*(s[k+1]-s[k]);
          
          
          
        }
      }else{
        k=0;
        
        mean2 = (1-exp(-lam1[0]*s[1]))/lam1[0];
        
      }
      
      
      
      MeanStore2(StoreInx)=mean2;
      
      
      
    }
    
    
    
    
  }
  
  
  List z1 = List::create(sstore,lamstore,Lstore,sigstore,BetaStore,MeanStore1,MeanStore2);
  
  
  return(z1);
  
  
}








//' Samples from the PEH Cox model with a patient covariate vector.
//'
//' Samples from the Piecewise Exponential Hazard (PEH) Cox model with a patient covariate vector and returns a list containing posterior parameters and posterior restricted mean survival.
//' @param Y Vector of event or censoring times.
//' @param I1 Vector of event indicators.
//' @param COV Matrix of size nxp containing p patient covariates.
//' @param Poi Prior mean number of split points.
//' @param B Number of iterations for MCMC.
//' @return Returns a list containing posterior samples of (1) the split point locations, (2) the log-hazards at each split point, (3) the number of split points, (4) the variance parameter for the log-hazard values, (5) the coefficients in the Cox model.
//' @examples
//' ##Generate Data
//' Y=rweibull(20,4,1)
//' I=rbinom(20,1,.5)
//' COV = matrix(rnorm(40,0,1),ncol=2)
//' ##Hyperparameter for number of split points
//' Poi=5
//'##Number of iterations for MCMC
//'B=200
//'BayesPiecewiseHazardCOV( Y, I,COV, Poi,  B)
//'@export
//[[Rcpp::export]]
List  BayesPiecewiseHazardCOV( arma::vec Y, //Survival Times
                               arma::vec I1, //Censoring Indicators
                               arma::mat COV, // Covariate Matrix Indicators
                               double Poi, //Prior Mean Number of split points
                               int B ){
  
  
  
  //Prior Params, make user controlled later
  double Var1=25;
  
  
  //
  //Upper Limit on Accept
  double High=.4;
  //Lower Limit on Accept
  double Low=.2;
  
  //For loop quantiites
  int m=0;
  int b=0;
  int k=0;
  int j;
  //Important double quantities used in MCMC
  int B1=B/2;
  
  
  //Innitialize Parameters and storage matrices
  double beta = 1  ;
  
  
  
  arma::vec betastore(B1);
  arma::vec MeanStore1(B1);
  arma::vec MeanStore2(B1);
  
  
  
  arma::vec MedianStore=MeanStore1;
  //Max To Store, really is 1 less here
  int Lmax = 30;
  arma::mat sstore(B1,Lmax+2);
  arma::mat lamstore(B1,Lmax+1);
  arma::vec Lstore(B1);
  arma::vec sigstore(B1);
  //Initialize S, lam and J
  arma::vec s(Lmax+2);
  arma::vec lam(Lmax+1);
  s.zeros();
  lam.zeros();
  double sig=1;
  double Birth=0;
  int Spot=0;
  
  double betaprop=0;
  double U=0;
  double alpha=0;
  double svar=1;
  double Ints=1;
  double Nums=2;
  //Copy Paste these to top
  double mean2=0;
  double cum1=0;
  double mean1=0;
  double Con=0;
  
  
  
  
  NumericVector z9(2);
  
  //Innitialize l with 3 split points
  int L = 3;
  
  //Contains the maximum observed censoring time
  double m1 = MaxVec(Y);
  //Start with split points at every $L/4$
  s[1]=m1/4;
  s[2]=2*m1/4;
  s[3]=3*m1/4;
  s[4]=m1;
  
  double NewLam1=0;
  double NewLam2=0;
  double U1=0;
  int which1=0;
  
  
  arma::vec lam1=lam;
  
  double bprop=1;
  double Intb=1;
  double Numb=2;
  arma::vec lamprop=lam;
  arma::vec sprop=s;
  
  double prob1=0;
  double med2=0;
  
  
  
  
  int StoreInx=0;
  
  
  //Lambda Tuning Parameters
  double LamC = .25;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam = 1;  //Number of proposals made
  double ALam = 1; //Number of proposals accepted
  
  
  
  //Do these for beta now
  arma::vec BetaC(COV.n_cols); //Beta Proposal Vector
  arma::vec NBeta = BetaC; //Counts number of Beta Proposals
  arma::vec ABeta=BetaC; //Counts Number of Beta Acceptances
  
  BetaC.zeros();
  NBeta.zeros();
  ABeta.zeros();
  
  BetaC = BetaC+.02;
  NBeta=NBeta+1;
  ABeta=ABeta+1;
  
  //Storage
  arma::mat BetaStore(B1,COV.n_cols);
  
  
  
  
  //Beta Starting Value
  arma::vec Beta(COV.n_cols);
  arma::vec BetaProp(COV.n_cols);
  Beta.zeros();
  
  for(j=0;j<COV.n_cols;j++){
    Beta[j]=as_scalar(arma::randn(1));
  }
  
  
  BetaProp.zeros();
  
  
  arma::vec z1(2);
  
  
  for(m=0;m<B;m++){
    
    
    if( (m%250==0) & (m<(B-B1))){
      //Tune proposal variance for lambda here
      if((ALam/NLam)>High){
        //Double the variance, acceptance is too high
        LamC = min1(LamC*2,2);
        
      }
      
      
      if((ALam/NLam)<Low){
        //Halve the variance, acceptance is too Low
        LamC = max1(LamC/2,.0625);
        
      }
      
      //Reset counters
      ALam =1;
      NLam=1;
      
      
      
      
      //Tune Proposals For \beta
      for(j=0;j<COV.n_cols;j++){
        if((ABeta[j]/NBeta[j])>High){
          //Double the variance, acceptance is too high
          BetaC[j] = min1(BetaC[j]*2,.16);
          
        }
        
        
        if((ABeta[j]/NBeta[j])<Low){
          //Halve the variance, acceptance is too Low
          BetaC[j] = max1(BetaC[j]/2,.005);
          
        }
        
        //Reset counters
        ABeta[j] =1;
        NBeta[j]=1;
        
      }
      
      
      
      
      
      
    }
    
    
    
    if(L>0){
      
      for(j=0;j<(L+1);j++){
        
        lamprop = lam;
        //Propose new value for lambda
        lamprop(j)=lam[j]+as_scalar(arma::randu(1))*2*LamC - LamC;
        //Likelihood Ratio
        alpha=LikePEHCOV(Y,I1,COV,s, lamprop , L,Beta) -    LikePEHCOV(Y,I1,COV, s, lam , L,Beta);
        //Prior Ratio
        if(j==0){
          alpha = alpha - .5*(pow(lamprop[j],2)-pow(lam[j],2))/Var1 - .5*pow(lamprop[j+1]-lamprop[j],2)/sig + .5*pow(lam[j+1]-lam[j],2);
        }else{
          if(j==L){
            alpha = alpha - .5*(pow(lamprop[j]-lamprop[j-1],2) + .5*pow(lam[j]-lam[j-1],2))/sig    ;
            
          }else{
            alpha = alpha - .5*(pow(lamprop[j]-lamprop[j-1],2) + .5*pow(lam[j]-lam[j-1],2))/sig  - .5*(pow(lamprop[j+1]-lamprop[j],2) + .5*pow(lam[j+1]-lam[j],2))/sig   ;
          }
        }
        
        
        
        //Metropolis Hastings Draw
        U = log(as_scalar(arma::randu(1)));
        //Accept/rej
        if(U<alpha){
          //Accept or Reject
          lam(j)=lamprop(j);
          ALam=ALam+1;
        }
        NLam=NLam+1;
        
      }
      
    }else{
      
      
      lamprop = lam;
      lamprop[0]=lam[0]+as_scalar(arma::randu(1))*2*LamC - LamC;
      //Likelihood Ratio
      alpha=LikePEHCOV(Y,I1, COV,s, lamprop , L,Beta) -    LikePEHCOV(Y,I1,COV, s, lam , L,Beta);
      //Prior Ratio Here for lambda 0
      alpha = alpha - .5*pow(lamprop[0],2)/Var1+.5*pow(lam[0],2)/Var1 ;
      U = log(as_scalar(arma::randu(1)));
      
      
      
      
      
      if(U<alpha){
        //accept proposal
        lam(0)=lamprop(0);
        ALam=ALam+1;
        
      }
      
      NLam=NLam+1;
      
    }
    
    //Sample Siglam
    cum1=0;
    if(L>0){
      for(j=1;j<L;j++){
        cum1 =cum1 + pow(lam[j]-lam[j-1],2);
      }
      
      sig = 1/R::rgamma(L/2 + 1, cum1/2 );
      
    }else{
      //Why is this here??
      sig=10000;
    }
    
    
    //Shuffle Splits
    
    if(L>0){
      
      for(j=1; j<(L+1); j++){
        sprop = s;
        //Draw new proposal for s_j, but make sure it's between s_j and s_j+1
        sprop[j]=as_scalar(arma::randu(1))*(sprop[j+1]-sprop[j-1]) + sprop[j-1];
        //Acceptance ratio
        alpha = log(sprop[j+1]-sprop[j])+log(sprop[j]-sprop[j-1])-log(s[j+1]-s[j])+
          log(s[j]-s[j-1]) +   LikePEHCOV(Y,I1,COV, sprop, lam , L,Beta) -    LikePEHCOV(Y,I1,COV, s, lam , L,Beta);
        
        //Metropolis Hastings
        U = log(as_scalar(arma::randu(1)));
        
        if(U<alpha){
          s[j]=sprop[j];
        }
        
        
        
      }
      
    }
    
    
    
    //Add proposal and MHG ratios here
    
    //Birth move Here
    if(L<(Lmax-1)){
      //What interval is it in?
      
      lamprop.zeros();
      //Find birth location      
      Spot = SampleBirth(s,L)+1;
      
      //Sample birthed split point
      Birth = as_scalar(arma::randu(1))*(s[Spot]-s[Spot-1])+s[Spot-1]; 
      //Now we have the interval of the new split in Spot and the actual location in Birth
      
      //Random Perturbation for detailed balance when dimension matching
      U1=as_scalar(arma::randu(1));
      
      //Find new lambdas
      NewLam1 =  lam[Spot-1] - log((1-U1)/U1)*(s[Spot+1]-Birth)/(s[Spot+1]-s[Spot-1]);
      NewLam2 =  lam[Spot-1] + log((1-U1)/U1)*(Birth-s[Spot])/(s[Spot+1]-s[Spot-1]);
      
      //Set sprop to 0s so we can fill it
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
      
      //Now we have our new proposal vector, evaluate it!
      //Like Ratio
      alpha =    LikePEHCOV(Y,I1,COV,sprop, lamprop , L+1,Beta) -    LikePEHCOV(Y,I1,COV, s, lam , L,Beta);
      //Add proposal ratio
      //Poisson
      alpha= alpha + log(Poi) - log1(L+1) ;
      // S proposal
      alpha = alpha + log1(2*L+3)+log1(2*L+2)+log(Birth-s[Spot-1])+log(s[Spot]-Birth);
      alpha = alpha - 2*log(m1)-log(s[Spot]-s[Spot-1]);
      //log of Jacobian for detailed balance
      alpha=alpha-log(U1*(1-U1)) ;
      
      //Add proposal ratio for \lambda
      if(L==0){
        alpha = alpha - .5*pow(lamprop[1]-lamprop[0],2)/sig   ;
      }else{
        
        if(Spot==(L+1)){
          //We're in the last interval
          alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig + pow(lam[Spot-1]-lam[Spot-2],2)/sig;
        }else{
          
          
          if(Spot==1){
            
            
            
            //Birthed is in the first interval      
            alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig ;
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig ;
            //First interval proposal
            alpha = alpha - .5*(pow(lamprop[0],2)-pow(lam[0],2))/Var1;
          }else{
            //The interval chosen is NOT the first or the last one, so we're ok here.
            alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;
          }  
          
          
        }
      }
      
      
      
      
      
      
      
      
      //Metropolis Hastings      
      U=log(as_scalar(arma::randu(1)));
      
      //Accept/Reject
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
      //Fill in Storage
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
      }
      
      //Fill in lambda proposal vector
      
      if(L>1){
        
        for(j=0;j<(Spot-1);j++){
          lamprop[j]=lam[j];
        }
        
        for(j=Spot;j<(lam.n_rows-1);j++){
          lamprop[j]=lam[j+1];
        }
        
        //New lambda is a weighted average of the old lambdas
        lamprop[Spot-1] = ((s[Spot+1]-s[Spot])*lam[Spot] +(s[Spot]-s[Spot-1])*lam[Spot-1])/(s[Spot+1]-s[Spot-1]);
        
        
      }else{
        lamprop[0]=((s[2]-s[1])*lam[1]+(s[1]-s[0])*lam[0])/(s[2]-s[0]);
      }
      //Finish sproposal fill in
      for(j=Spot; j<(s.n_rows-1);j++){
        sprop[j]=s[j+1];
      }
      
      //Now we have our new proposal vectors, evaluate then!
      //Like Ratio
      alpha =    LikePEHCOV(Y,I1,COV, sprop, lamprop , L-1,Beta) -    LikePEHCOV(Y,I1,COV, s, lam , L,Beta);
      //Prior Ratio
      //Poisson
      alpha = alpha  -log(Poi) + log1(L);
      //S Prior
      alpha = alpha + 2*log(m1) + log(s[Spot+1]-s[Spot-1]) - log1(2*L+1) - log1(2L) - log(s[Spot]-s[Spot-1])-log(s[Spot+1]-s[Spot]);
      //Lambda Prior, we DROPPED one
      if(L==1){
        //Removing will drop the sampler to 0 split points
        alpha = alpha + .5*pow(lam[1]-lam[0],2)/sig ;
      }else{
        
        
        
        if(Spot>1){
          
          if(Spot==(L+1)){
            //We dropped the last interval, so we can't evaluate intervals past this
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig  - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;
          }else{
            
            //Split point is not first or last
            alpha = alpha +  .5*pow(lam[Spot+1]-lam[Spot],2)/sig +  .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;
            alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig ;
            
          }
          
          
        }else{
          
          //First Interval is changing so we must evaluate it
          alpha = alpha + .5*pow(lam[Spot+1]-lam[Spot],2)/sig + .5*pow(lam[Spot]-lam[Spot-1],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig  ;
          //Add the 0 entry here
          alpha= alpha  - .5*pow(lamprop[0],2)/Var1 + .5*pow(lam[0],2)/Var1 ;
          
          
        }
        
        
        
        
        
      }
      
      
      
      
      
      
      //Random Perturbation for dimension matching
      U1 = as_scalar(arma::randu(1));
      alpha = alpha + log(U1)+log(1-U1);
      
      //Metropolis Draw
      U=log(as_scalar(arma::randu(1)));
      //Accept/Reject
      if(U<alpha){
        s=sprop;
        L=L-1;
        lam=lamprop;
      }
      
      
      
      
    }
    
    
    
    
    
    //Beta Sampler
    for(j=0;j<COV.n_cols;j++){
      BetaProp=Beta;
      BetaProp[j] = Beta[j]  +as_scalar(arma::randu(1))*BetaC[j]*2 - BetaC[j];
      //Default for now, could make it adaptive
      //Log-Acceptance ratio
      alpha=LikePEHCOV(Y,I1, COV, s, lam , L, BetaProp) -    LikePEHCOV(Y,I1, COV, s, lam , L, Beta);
      //Adjust for proposal ratio
      alpha = alpha - .5*pow(BetaProp[j],2)/25 +.5*pow(Beta[j],2)/25    ;
      
      
      
      U=log(as_scalar(arma::randu(1)));
      if(U<alpha){
        Beta[j]=BetaProp[j];
        ABeta[j]=ABeta[j]+1;
      }
      NBeta[j]=NBeta[j]+1;
    }
    
    
    
    
    
    if(m>(B-B1-1)){
      //Store Values in Matrix
      StoreInx = m-B+B1;
      
      for(j=0;j<sstore.n_cols;j++){
        sstore(StoreInx,j) = s(j);
      }
      
      for(j=0;j<lamstore.n_cols;j++){
        lamstore(StoreInx,j) = lam(j);
      }
      
      Lstore[StoreInx]=L;
      sigstore[StoreInx]=sig;
      for(j=0;j<COV.n_cols;j++){
        BetaStore(StoreInx,j)=Beta[j];
      }
      
    }
  }
  
  
  
  List z101 = List::create(sstore,lamstore,Lstore,sigstore,BetaStore);
  
  
  return(z101);
  
  
}














//' Samples from the PLLH Cox model with a patient covariate vector.
//'
//' Samples from the Piecewise Linear Log-Hazard (PLLH) Cox model with a patient covariate vector and returns a list containing posterior parameters and posterior restricted mean survival.
//' @param Y Vector of event or censoring times.
//' @param I1 Vector of event indicators.
//' @param COV Matrix of size nxp containing p patient covariates.
//' @param Poi Prior mean number of split points.
//' @param B Number of iterations for MCMC.
//' @return Returns a list containing posterior samples of (1) the split point locations, (2) the log-hazards at each split point, (3) the number of split points, (4) the variance parameter for the log-hazard values, (5) the coefficients in the Cox model.
//' @examples
//' ##Generate Data
//' Y=rweibull(20,4,1)
//' I=rbinom(20,1,.5)
//' COV = matrix(rnorm(40,0,1),ncol=2)
//' ##Hyperparameter for number of split points
//' Poi=5
//'##Number of iterations for MCMC
//'B=200
//'BayesPiecewiseLinearLogHazardCOV( Y, I,COV, Poi,  B)
//'@export
//[[Rcpp::export]]
List  BayesPiecewiseLinearLogHazardCOV( arma::vec Y, //Survival Times
                                        arma::vec I1, //Censoring Indicators
                                        arma::mat COV, // Covariate Matrix Indicators
                                        double Poi, //Prior Mean Number of split points
                                        int B ){
  
  
  
  //Prior Params, make user controlled later
  double Var1=25;
  
  
  //
  //Upper Limit on Accept
  double High=.4;
  //Lower Limit on Accept
  double Low=.2;
  
  //For loop quantiites
  int m=0;
  int b=0;
  int k=0;
  int j;
  //Important double quantities used in MCMC
  int B1=B/2;
  
  
  //Innitialize Parameters and storage matrices
  double beta = 1  ;
  
  
  
  arma::vec betastore(B1);
  arma::vec MeanStore1(B1);
  arma::vec MeanStore2(B1);
  
  
  
  arma::vec MedianStore=MeanStore1;
  //Max To Store, really is 1 less here
  int Lmax = 30;
  arma::mat sstore(B1,Lmax+2);
  arma::mat lamstore(B1,Lmax+2);
  arma::vec Lstore(B1);
  arma::vec sigstore(B1);
  //Initialize S, lam and J
  arma::vec s(Lmax+2);
  arma::vec lam(Lmax+2);
  s.zeros();
  lam.zeros();
  double sig=1;
  double Birth=0;
  int Spot=0;
  
  double betaprop=0;
  double U=0;
  double alpha=0;
  double svar=1;
  double Ints=1;
  double Nums=2;
  //Copy Paste these to top
  double mean2=0;
  double cum1=0;
  double mean1=0;
  double Con=0;
  
  
  
  
  NumericVector z9(2);
  
  //Innitialize l with 3 split points
  int L = 3;
  
  //Contains the maximum observed censoring time
  double m1 = MaxVec(Y);
  //Start with split points at every $L/4$
  s[1]=m1/4;
  s[2]=2*m1/4;
  s[3]=3*m1/4;
  s[4]=m1;
  
  double NewLam1=0;
  double NewLam2=0;
  double U1=0;
  int which1=0;
  
  
  arma::vec lam1=lam;
  
  double bprop=1;
  double Intb=1;
  double Numb=2;
  arma::vec lamprop=lam;
  arma::vec sprop=s;
  
  double prob1=0;
  double med2=0;
  
  
  
  
  int StoreInx=0;
  
  
  //Lambda Tuning Parameters
  double LamC = .25;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam = 1;  //Number of proposals made
  double ALam = 1; //Number of proposals accepted
  
  
  
  //Do these for beta now
  arma::vec BetaC(COV.n_cols); //Beta Proposal Vector
  arma::vec NBeta = BetaC; //Counts number of Beta Proposals
  arma::vec ABeta=BetaC; //Counts Number of Beta Acceptances
  
  BetaC.zeros();
  NBeta.zeros();
  ABeta.zeros();
  
  BetaC = BetaC+.02;
  NBeta=NBeta+1;
  ABeta=ABeta+1;
  
  //Storage
  arma::mat BetaStore(B1,COV.n_cols);
  
  //Now Our initial split point vector is set
  lam(0)=as_scalar(arma::randn(1));
  lam(1)=as_scalar(arma::randn(1));
  lam(2)=as_scalar(arma::randn(1));
  lam(3)=as_scalar(arma::randn(1));
  lam(4)=as_scalar(arma::randn(1));
  
  
  
  //Beta Starting Value
  arma::vec Beta(COV.n_cols);
  arma::vec BetaProp(COV.n_cols);
  Beta.zeros();
  
  for(j=0;j<COV.n_cols;j++){
    Beta[j]=as_scalar(arma::randn(1));
  }
  
  
  BetaProp.zeros();
  
  
  arma::vec z1(2);
  
  
  
  
  for(m=0;m<B;m++){
    
    
    if( (m%250==0) & (m<(B-B1))){
      //Tune proposal variance for lambda here
      if((ALam/NLam)>High){
        //Double the variance, acceptance is too high
        LamC = min1(LamC*2,2);
        
      }
      
      
      if((ALam/NLam)<Low){
        //Halve the variance, acceptance is too Low
        LamC = max1(LamC/2,.0625);
        
      }
      
      //Reset counters
      ALam =1;
      NLam=1;
      
      
      
      
      //Tune Proposals For \beta
      for(j=0;j<COV.n_cols;j++){
        if((ABeta[j]/NBeta[j])>High){
          //Double the variance, acceptance is too high
          BetaC[j] = min1(BetaC[j]*2,.16);
          
        }
        
        
        if((ABeta[j]/NBeta[j])<Low){
          //Halve the variance, acceptance is too Low
          BetaC[j] = max1(BetaC[j]/2,.005);
          
        }
        
        //Reset counters
        ABeta[j] =1;
        NBeta[j]=1;
        
      }
      
      
      
      
      
      
    }
    
    
    
    
    
    if(L>0){
      
      for(j=0;j<(L+2);j++){
        
        lamprop = lam;
        //Propose new value for lambda
        lamprop(j)=lam[j]+as_scalar(arma::randu(1))*2*LamC - LamC;
        //Likelihood Ratio
        
        alpha=LikePLLHCOV(Y,I1,COV,s, lamprop , L,Beta) -    LikePLLHCOV(Y,I1,COV, s, lam , L,Beta);
        
        
        //Prior Ratio
        if(j==0){
          alpha = alpha - .5*(pow(lamprop[j],2)-pow(lam[j],2))/Var1 - .5*pow(lamprop[j+1]-lamprop[j],2)/sig + .5*pow(lam[j+1]-lam[j],2);
        }else{
          if(j==L){
            alpha = alpha - .5*(pow(lamprop[j]-lamprop[j-1],2) + .5*pow(lam[j]-lam[j-1],2))/sig    ;
            
          }else{
            alpha = alpha - .5*(pow(lamprop[j]-lamprop[j-1],2) + .5*pow(lam[j]-lam[j-1],2))/sig  - .5*(pow(lamprop[j+1]-lamprop[j],2) + .5*pow(lam[j+1]-lam[j],2))/sig   ;
          }
        }
        
        
        
        //Metropolis Hastings Draw
        U = log(as_scalar(arma::randu(1)));
        //Accept/rej
        if(U<alpha){
          //Accept or Reject
          lam(j)=lamprop(j);
          ALam=ALam+1;
        }
        NLam=NLam+1;
        
      }
      
    }else{
      
      
      lamprop = lam;
      lamprop[0]=lam[0]+as_scalar(arma::randu(1))*2*LamC - LamC;
      //Likelihood Ratio
      alpha=LikePLLHCOV(Y,I1, COV,s, lamprop , L,Beta) -    LikePLLHCOV(Y,I1,COV, s, lam , L,Beta);
      //Prior Ratio Here for lambda 0
      alpha = alpha - .5*pow(lamprop[0],2)/Var1+.5*pow(lam[0],2)/Var1 ;
      
      alpha = alpha - .5*pow(lamprop[1]-lamprop[0],2)/sig + .5*pow(lam[1]-lam[0],2);
      
      U = log(as_scalar(arma::randu(1)));
      
      
      
      
      
      if(U<alpha){
        //accept proposal
        lam(0)=lamprop(0);
        ALam=ALam+1;
        
      }
      
      NLam=NLam+1;
      
      
      
      lamprop = lam;
      lamprop[1]=lam[1]+as_scalar(arma::randu(1))*2*LamC - LamC;
      //Likelihood Ratio
      alpha=LikePLLHCOV(Y,I1, COV,s, lamprop , L,Beta) -    LikePLLHCOV(Y,I1,COV, s, lam , L,Beta);
      //Prior Ratio Here for lambda 0
      alpha = alpha - .5*pow(lamprop[0],2)/Var1+.5*pow(lam[0],2)/Var1 ;
      
      alpha = alpha - .5*pow(lamprop[1]-lamprop[0],2)/sig + .5*pow(lam[1]-lam[0],2);
      
      U = log(as_scalar(arma::randu(1)));
      
      
      
      
      
      if(U<alpha){
        //accept proposal
        lam(1)=lamprop(1);
        ALam=ALam+1;
        
      }
      
      NLam=NLam+1;
      
      
      
    }
    
    
    
    //Sample Siglam
    cum1=0;
    for(j=1;j<(L+1);j++){
      cum1 =cum1 + pow(lam[j]-lam[j-1],2);
    }
    
    sig = 1/R::rgamma(L/2 + 1, cum1/2 );
    
    
    
    
    //Shuffle Splits
    
    if(L>0){
      
      for(j=1; j<(L+1); j++){
        sprop = s;
        //Draw new proposal for s_j, but make sure it's between s_j and s_j+1
        sprop[j]=as_scalar(arma::randu(1))*(sprop[j+1]-sprop[j-1]) + sprop[j-1];
        //Acceptance ratio
        alpha = log(sprop[j+1]-sprop[j])+log(sprop[j]-sprop[j-1])-log(s[j+1]-s[j])+
          log(s[j]-s[j-1]) +   LikePLLHCOV(Y,I1,COV, sprop, lam , L,Beta) -    LikePLLHCOV(Y,I1,COV, s, lam , L,Beta);
        
        //Metropolis Hastings
        U = log(as_scalar(arma::randu(1)));
        
        if(U<alpha){
          s[j]=sprop[j];
        }
        
        
        
      }
      
    }
    
    
    
    //Add proposal and MHG ratios here
    
    //Death move
    if(L>0){
      //Draw proposed split point to delete
      Spot=SampleDeath(L);
      sprop.zeros();
      lamprop.zeros();
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
      }
      
      //We just straight up delete \lamda_spot and s[spot] and make adjustments
      if(L>1){
        for(j=0;j<Spot;j++){
          lamprop[j]=lam[j];
        }
        for(j=Spot;j<(lam.n_rows-1);j++){
          lamprop[j]=lam[j+1];
        }
      }else{
        lamprop[0]=lam[0];
        lamprop[1]=lam[2];
      }
      for(j=Spot; j<(s.n_rows-1);j++){
        sprop[j]=s[j+1];
      }
      //Likelihood ratio
      alpha =    LikePLLHCOV(Y,I1, COV,  sprop, lamprop , L-1,Beta) -    LikePLLHCOV(Y,I1, COV, s, lam , L,Beta);
      //Prior Ratio
      //Poisson
      alpha = alpha  -log(Poi) + log1(L);
      //S Prior
      alpha = alpha +2*log(m1) + log(s[Spot+1]-s[Spot-1]) - log1(2*L+1) - log1(2L) - log(s[Spot]-s[Spot-1])-log(s[Spot+1]-s[Spot]);
      //Lambda Prior, we DROPPED one
      if(L==1){
        //Removing will drop the sampler to 0 split points
        alpha = alpha - .5*pow(lam[2]-lam[0],2)/sig + .5*pow(lam[2]-lam[1],2)/sig + .5*pow(lam[1]-lam[0],2)/sig ;
      }else{
        alpha = alpha -.5*pow(lam[Spot+1]-lam[Spot-1],2)/sig + .5*pow(lam[Spot+1]-lam[Spot],2)/sig+.5*pow(lam[Spot]-lam[Spot-1],2)/sig;
      }
      
      
      //Random Perturbation to maintain balance between parameter spaces
      U1 = as_scalar(arma::randu(1));
      ///Add log Jacobian Here
      alpha = alpha +log(U1*(1-U1));
      //Metropolis Hastings Green
      U=log(as_scalar(arma::randu(1)));
      
      if(U<alpha){
        s=sprop;
        L=L-1;
        lam=lamprop;
      }
    }
    
    
    
    //Birth move Here
    if(L<(Lmax-1)){
      //What interval is it in?
      //This generates the interval location that the new proposed split point is located
      Spot = SampleBirth(s,L)+1;
      //Set Lambda to 0
      lamprop.zeros();
      //Now generate the new proposed split point location
      Birth = as_scalar(arma::randu(1))*(s[Spot]-s[Spot-1])+s[Spot-1];
      //Random Perturbation for dimension matching
      U1=as_scalar(arma::randu(1));
      //Slope
      
      //Old Slope for log hazard
      NewLam1= (lam[Spot]-lam[Spot-1])/(s[Spot]-s[Spot-1]);
      //Compute the new slope that we can use to get the new lambda value
      //This quantity comes from the MCMC section on how to obtain new slopes
      NewLam1 = NewLam1 - log((1-U1)/U1)*(s[Spot]-Birth)/(s[Spot]-s[Spot-1]);
      //We can now use this slope to obtain the new hazard height at the time point Birth
      NewLam2 = NewLam1*(Birth-s[Spot-1])+lam[Spot-1];
      //Now NewLam2 contains the hazard at the time point Birth
      //Now we have the interval of the new split in Spot and the actual location in Birth
      sprop.zeros();
      //Here we're going to fill in our new vector
      //Let's add Birth to the Spot location and push back the rest of sprop
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
        lamprop[j]=lam[j];
      }
      sprop[Spot]=Birth;
      for(j=(Spot+1);j<(s.n_rows-1);j++){
        sprop[j]=s[j-1];
      }
      //Log of hazard height at Birth.
      lamprop[Spot]=NewLam2;
      for(j=(Spot+1);j<(lam.n_rows-1);j++){
        lamprop[j]=lam[j-1];
      }
      
      
      //Now we have our new proposal vector, evaluate it!
      //Like Ratio
      alpha =    LikePLLHCOV(Y,I1, COV, sprop, lam , L+1, Beta) -    LikePLLHCOV(Y,I1, COV, s, lam , L, Beta);
      //Add proposal ratio
      //Poisson
      alpha= alpha + log(Poi) - log1(L+1) ;
      // S proposal
      alpha = alpha + log1(2*L+3)+log1(2*L+2)+log(Birth-s[Spot-1])+log(s[Spot]-Birth) - 2*log(m1)-log(s[Spot]-s[Spot-1]);
      //Add proposal ratio for \lambda
      if(L==0){
        alpha = alpha  - .5*pow(lamprop[1]-lamprop[0],2)/sig - .5*pow(lamprop[2]-lamprop[1],2)/sig + .5*pow(lam[1]-lam[0],2)/sig ;
      }else{
        alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig + .5*pow(lam[Spot]-lam[Spot-1],2)/sig;
      }
      
      //Add log Jacobian
      alpha = alpha -log(U1*(1-U1));
      //Uniform for Metropolis-Hastings
      U=log(as_scalar(arma::randu(1)));
      if(U<alpha){
        //We have to set the entire vector here
        s=sprop;
        lam=lamprop;
        L=L+1;
      }
    }
    
    
    
    
    //Beta Sampler
    for(j=0;j<COV.n_cols;j++){
      BetaProp=Beta;
      BetaProp[j] = Beta[j]  +as_scalar(arma::randu(1))*BetaC[j]*2 - BetaC[j];
      //Default for now, could make it adaptive
      //Log-Acceptance ratio
      alpha=LikePLLHCOV(Y,I1, COV, s, lam , L, BetaProp) -    LikePLLHCOV(Y,I1, COV, s, lam , L, Beta);
      //Adjust for proposal ratio
      alpha = alpha - .5*pow(BetaProp[j],2)/25 +.5*pow(Beta[j],2)/25    ;
      
      
      
      U=log(as_scalar(arma::randu(1)));
      if(U<alpha){
        Beta[j]=BetaProp[j];
        ABeta[j]=ABeta[j]+1;
      }
      NBeta[j]=NBeta[j]+1;
    }
    
    
    
    
    
    if(m>(B-B1-1)){
      //Store Values in Matrix
      
      StoreInx = m-B+B1;
      
      for(j=0;j<sstore.n_cols;j++){
        sstore(StoreInx,j) = s(j);
      }
      
      for(j=0;j<lamstore.n_cols;j++){
        lamstore(StoreInx,j) = lam(j);
      }
      
      Lstore[StoreInx]=L;
      sigstore[StoreInx]=sig;
      for(j=0;j<COV.n_cols;j++){
        BetaStore(StoreInx,j)=Beta[j];
      }
      
    }
  }
  
  
  
  List z101 = List::create(sstore,lamstore,Lstore,sigstore,BetaStore);
  
  
  return(z101);
  
  
}









//Function for computing the likelihood of Piecewise Exponential Hazard
double LikePEH( arma::vec Y, //Survival Times
                arma::vec I1, //Censoring Indicators
                arma::vec s, //Vector containing the split point locations
                arma::vec lam, //Vector containing the log hazard heights
                int J){
  
  int m=0;
  int l=0;
  double LogL=0;
  
  
  //Cycle through each interval and obtain survival estimates and add hazard heights
  for(l=0; l<(J+1);l++){
    for(m=0; m<Y.n_rows;m++){
      LogL = LogL - max1(0,min1(s(l+1),Y(m))-s(l))*exp(lam[l]);
      
      
      if((Y(m)>s[l]) & (Y(m)<=s[l+1]) & (I1[m]==1)){
        LogL = LogL + lam[l];
      }
      
      
      
      
    }
    
  }
  
  
  
  
  return(LogL);  
  
  
  
}






//' Samples from the PEH model without covariates.
//'
//' Samples from the Piecewise Exponential Hazard (PEH) model and returns a list containing posterior parameters and posterior restricted mean survival.
//' @param Y Vector of event or censoring times.
//' @param I1 Vector of event indicators.
//' @param Poi Prior mean number of split points.
//' @param B Number of iterations for MCMC.
//' @return Returns a list containing posterior samples of (1) the split point locations, (2) the log-hazards at each split point, (3) the number of split points, (4) the variance parameter for the log-hazard values, (5) the posterior mean restricted survivial time.
//' @examples
//' ##Generate Data
//' Y=rweibull(20,4,1)
//' I=rbinom(20,1,.5)
//' ##Hyperparameter for number of split points
//' Poi=5
//'##Number of iterations for MCMC
//'B=200
//'BayesPiecewiseHazard( Y, I, Poi,  B)
//'@export
//[[Rcpp::export]]
List  BayesPiecewiseHazard( arma::vec Y, //Survival Times
                            arma::vec I1, //Censoring Indicators
                            double Poi, //Prior Mean Number of split points
                            int B ){
  
  
  
  //Prior Params, make user controlled later
  double Var1=25;
  
  
  //
  //Upper Limit on Accept
  double High=.4;
  //Lower Limit on Accept
  double Low=.2;
  
  //For loop quantiites
  int m=0;
  int b=0;
  int k=0;
  int j;
  //Important double quantities used in MCMC
  int B1=B/10;
  
  
  //Innitialize Parameters and storage matrices
  double beta = 1  ;
  
  
  
  arma::vec betastore(B1);
  arma::vec MeanStore(B1);
  arma::vec MedianStore=MeanStore;
  //Max To Store, really is 1 less here
  int Lmax = 30;
  arma::mat sstore(B1,Lmax+2);
  arma::mat lamstore(B1,Lmax+1);
  arma::vec Lstore(B1);
  arma::vec sigstore(B1);
  //Initialize S, lam and J
  arma::vec s(Lmax+2);
  arma::vec lam(Lmax+1);
  s.zeros();
  lam.zeros();
  double sig=1;
  double Birth=0;
  int Spot=0;
  
  double betaprop=0;
  double U=0;
  double alpha=0;
  double svar=1;
  double Ints=1;
  double Nums=2;
  //Copy Paste these to top
  double mean2=0;
  double cum1=0;
  double mean1=0;
  double Con=0;
  
  
  
  
  NumericVector z9(2);
  
  //Innitialize l with 3 split points
  int L = 3;
  
  //Contains the maximum observed censoring time
  double m1 = MaxVec(Y);
  //Start with split points at every $L/4$
  s[1]=m1/4;
  s[2]=2*m1/4;
  s[3]=3*m1/4;
  s[4]=m1;
  
  double NewLam1=0;
  double NewLam2=0;
  double U1=0;
  int which1=0;
  
  
  arma::vec lam1=lam;
  
  double bprop=1;
  double Intb=1;
  double Numb=2;
  arma::vec lamprop=lam;
  arma::vec sprop=s;
  
  double prob1=0;
  double med2=0;
  
  
  
  
  int StoreInx=0;
  
  
  //Lambda Tuning Parameters
  double LamC = .25;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam = 1;  //Number of proposals made
  double ALam = 1; //Number of proposals accepted
  
  
  for(m=0;m<B;m++){
    
    
    if( (m%250==0) & (m<(B-B1))){
      //Tune proposal variance for lambda here
      if((ALam/NLam)>High){
        //Double the variance, acceptance is too high
        LamC = min1(LamC*2,2);
        
      }
      
      
      if((ALam/NLam)<Low){
        //Halve the variance, acceptance is too Low
        LamC = max1(LamC/2,.0625);
        
      }
      
      //Reset counters
      ALam =1;
      NLam=1;
      
    }
    
    
    
    if(L>0){
      
      for(j=0;j<(L+1);j++){
        
        lamprop = lam;
        //Propose new value for lambda
        lamprop(j)=lam[j]+as_scalar(arma::randu(1))*2*LamC - LamC;
        //Likelihood Ratio
        alpha=LikePEH(Y,I1,s, lamprop , L) -    LikePEH(Y,I1, s, lam , L);
        //Prior Ratio
        if(j==0){
          alpha = alpha - .5*(pow(lamprop[j],2)-pow(lam[j],2))/Var1 - .5*pow(lamprop[j+1]-lamprop[j],2)/sig + .5*pow(lam[j+1]-lam[j],2);
        }else{
          if(j==L){
            alpha = alpha - .5*(pow(lamprop[j]-lamprop[j-1],2) + .5*pow(lam[j]-lam[j-1],2))/sig    ;
            
          }else{
            alpha = alpha - .5*(pow(lamprop[j]-lamprop[j-1],2) + .5*pow(lam[j]-lam[j-1],2))/sig  - .5*(pow(lamprop[j+1]-lamprop[j],2) + .5*pow(lam[j+1]-lam[j],2))/sig   ;
          }
        }
        
        
        
        //Metropolis Hastings Draw
        U = log(as_scalar(arma::randu(1)));
        //Accept/rej
        if(U<alpha){
          //Accept or Reject
          lam(j)=lamprop(j);
          ALam=ALam+1;
        }
        NLam=NLam+1;
        
      }
      
    }else{
      
      
      lamprop = lam;
      lamprop[0]=lam[0]+as_scalar(arma::randu(1))*2*LamC - LamC;
      //Likelihood Ratio
      alpha=LikePEH(Y,I1,s, lamprop , L) -    LikePEH(Y,I1, s, lam , L);
      //Prior Ratio Here for lambda 0
      alpha = alpha - .5*pow(lamprop[0],2)/Var1+.5*pow(lam[0],2)/Var1 ;
      U = log(as_scalar(arma::randu(1)));
      
      
      
      
      
      if(U<alpha){
        //accept proposal
        lam(0)=lamprop(0);
        ALam=ALam+1;
        
      }
      
      NLam=NLam+1;
      
    }
    
    //Sample Siglam
    cum1=0;
    if(L>0){
      for(j=1;j<L;j++){
        cum1 =cum1 + pow(lam[j]-lam[j-1],2);
      }
      
      sig = 1/R::rgamma(L/2 + 1, cum1/2 );
      
    }else{
      //Why is this here??
      sig=10000;
    }
    
    
    //Shuffle Splits
    
    if(L>0){
      
      for(j=1; j<(L+1); j++){
        sprop = s;
        //Draw new proposal for s_j, but make sure it's between s_j and s_j+1
        sprop[j]=as_scalar(arma::randu(1))*(sprop[j+1]-sprop[j-1]) + sprop[j-1];
        //Acceptance ratio
        alpha = log(sprop[j+1]-sprop[j])+log(sprop[j]-sprop[j-1])-log(s[j+1]-s[j])+log(s[j]-s[j-1]) +   LikePEH(Y,I1, sprop, lam , L) -    LikePEH(Y,I1, s, lam , L);
        
        //Metropolis Hastings
        U = log(as_scalar(arma::randu(1)));
        
        if(U<alpha){
          s[j]=sprop[j];
        }
        
        
        
      }
      
    }
    
    
    
    //Add proposal and MHG ratios here
    
    //Birth move Here
    if(L<(Lmax-1)){
      //What interval is it in?
      
      lamprop.zeros();
      //Find birth location      
      Spot = SampleBirth(s,L)+1;
      
      //Sample birthed split point
      Birth = as_scalar(arma::randu(1))*(s[Spot]-s[Spot-1])+s[Spot-1]; 
      //Now we have the interval of the new split in Spot and the actual location in Birth
      
      //Random Perturbation for detailed balance when dimension matching
      U1=as_scalar(arma::randu(1));
      
      //Find new lambdas
      NewLam1 =  lam[Spot-1] - log((1-U1)/U1)*(s[Spot+1]-Birth)/(s[Spot+1]-s[Spot-1]);
      NewLam2 =  lam[Spot-1] + log((1-U1)/U1)*(Birth-s[Spot])/(s[Spot+1]-s[Spot-1]);
      
      //Set sprop to 0s so we can fill it
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
      
      //Now we have our new proposal vector, evaluate it!
      //Like Ratio
      alpha =    LikePEH(Y,I1,sprop, lamprop , L+1) -    LikePEH(Y,I1, s, lam , L);
      //Add proposal ratio
      //Poisson
      alpha= alpha + log(Poi) - log1(L+1) ;
      // S proposal
      alpha = alpha + log1(2*L+3)+log1(2*L+2)+log(Birth-s[Spot-1])+log(s[Spot]-Birth);
      alpha = alpha - 2*log(m1)-log(s[Spot]-s[Spot-1]);
      //log of Jacobian for detailed balance
      alpha=alpha-log(U1*(1-U1)) ;
      
      //Add proposal ratio for \lambda
      if(L==0){
        alpha = alpha - .5*pow(lamprop[1]-lamprop[0],2)/sig   ;
      }else{
        
        if(Spot==(L+1)){
          //We're in the last interval
          alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig + pow(lam[Spot-1]-lam[Spot-2],2)/sig;
        }else{
          
          
          if(Spot==1){
            
            
            
            //Birthed is in the first interval      
            alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig ;
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig ;
            //First interval proposal
            alpha = alpha - .5*(pow(lamprop[0],2)-pow(lam[0],2))/Var1;
          }else{
            //The interval chosen is NOT the first or the last one, so we're ok here.
            alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;
          }  
          
          
        }
      }
      
      
      
      
      
      
      
      
      //Metropolis Hastings      
      U=log(as_scalar(arma::randu(1)));
      
      //Accept/Reject
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
      //Fill in Storage
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
      }
      
      //Fill in lambda proposal vector
      
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
      //Finish sproposal fill in
      for(j=Spot; j<(s.n_rows-1);j++){
        sprop[j]=s[j+1];
      }
      
      //Now we have our new proposal vectors, evaluate then!
      //Like Ratio
      alpha =    LikePEH(Y,I1, sprop, lamprop , L-1) -    LikePEH(Y,I1, s, lam , L);
      //Prior Ratio
      //Poisson
      alpha = alpha  -log(Poi) + log1(L);
      //S Prior
      alpha = alpha + 2*log(m1) + log(s[Spot+1]-s[Spot-1]) - log1(2*L+1) - log1(2L) - log(s[Spot]-s[Spot-1])-log(s[Spot+1]-s[Spot]);
      //Lambda Prior, we DROPPED one
      if(L==1){
        //Removing will drop the sampler to 0 split points
        alpha = alpha + .5*pow(lam[1]-lam[0],2)/sig ;
      }else{
        
        
        
        if(Spot>1){
          
          if(Spot==(L+1)){
            //We dropped the last interval, so we can't evaluate intervals past this
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig  - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;
          }else{
            
            //Split point is not first or last
            alpha = alpha +  .5*pow(lam[Spot+1]-lam[Spot],2)/sig +  .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;
            alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig ;
            
          }
          
          
        }else{
          
          //First Interval is changing so we must evaluate it
          alpha = alpha + .5*pow(lam[Spot+1]-lam[Spot],2)/sig + .5*pow(lam[Spot]-lam[Spot-1],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig  ;
          //Add the 0 entry here
          alpha= alpha  - .5*pow(lamprop[0],2)/Var1 + .5*pow(lam[0],2)/Var1 ;
          
          
        }
        
        
        
        
        
      }
      
      
      
      
      
      
      //Random Perturbation for dimension matching
      U1 = as_scalar(arma::randu(1));
      alpha = alpha + log(U1)+log(1-U1);
      
      //Metropolis Draw
      U=log(as_scalar(arma::randu(1)));
      //Accept/Reject
      if(U<alpha){
        s=sprop;
        L=L-1;
        lam=lamprop;
      }
      
      
      
      
    }
    
    
    
    if(m>(B-B1-1)){
      //Store Values in Matrix
      StoreInx = m-B+B1;
      
      for(j=0;j<sstore.n_cols;j++){
        sstore(StoreInx,j) = s(j);
      }
      
      for(j=0;j<lamstore.n_cols;j++){
        lamstore(StoreInx,j) = lam(j);
      }
      
      Lstore[StoreInx]=L;
      sigstore[StoreInx]=sig;
      
      //Calculate Means Here
      lam1=exp(lam);
      mean2=0;
      cum1=0;
      
      if(L>0){
        for(k=0;k<(L+1);k++){
          
          mean2= mean2+ exp(-cum1+lam1[k]*s[k])*(exp(-lam1[k]*s[k])-exp(-lam1[k]*s[k+1]))/(lam1[k]);
          
          cum1 = cum1+ lam1[k]*(s[k+1]-s[k]);
          
          
          
        }
      }else{
        k=0;
        
        mean2 = (1-exp(-lam1[0]*s[1]))/lam1[0];
        
      }
      
      
      
      MeanStore(StoreInx)=mean2;
      
      
      
      
    }
    
    
    
    
  }
  
  
  List z1 = List::create(sstore,lamstore,Lstore,sigstore,MeanStore);
  
  
  return(z1);
  
  
}







double LikePLLH( arma::vec Y, //Survival Times
                 arma::vec I1, //Censoring Indicators
                 arma::vec s, //Vector containing the split point locations
                 arma::vec lam, //Vector containing the log hazard heights
                 int J //Numer of split points
){
  
  
  
  
  int m=0;
  int l=0;
  double LogL=0;
  double Y1=0;
  //Make Slope
  
  NumericVector z9(4);
  
  
  arma::vec slopes = GetSlopePLLH(s,lam,J);
  //Cycle through each interval and obtain survival estimates and add hazard heights
  for(l=0; l<(J+1);l++){
    for(m=0; m<Y.n_rows;m++){
      //Y1 here contains the max time that a patient was in interval l
      Y1=min1(s(l+1),Y(m));
      if(Y1>s(l)){
        //Means patient has survived past time s_l
        //This is the survival contribution for patient i in interval l
        
        LogL = LogL + exp(lam(l))*(1-exp(slopes(l)*(Y1-s(l))))/slopes(l);
        
        //If this holds below, then the patient died in the interval [s_l, s_{l+1}]
        if( (Y1<s(l+1)) && (I1(m)==1)){
          //Patient died in interval, need hazard here
          LogL = LogL + lam(l)+(Y1-s(l))*slopes(l);
        }
        
        
      }
      
      
      
      
    }
    
    
  }
  
  
  
  return(LogL);
  
  
  
}


//' Samples from the PLLH model without covariates.
//'
//' Samples from the Piecewise Linear Log-Hazard (PLLH) model and returns a list containing posterior parameters and posterior restricted mean survival.
//' @param Y Vector of event or censoring times.
//' @param I1 Vector of event indicators.
//' @param Poi Prior mean number of split points.
//' @param B Number of iterations for MCMC.
//' @return Returns a list containing posterior samples of (1) the split point locations, (2) the log-hazards at each split point, (3) the number of split points, (4) the variance parameter for the log-hazard values, (5) the posterior mean restricted survivial time.
//' @examples
//' ##Generate Data
//' Y=rweibull(20,4,1)
//' I=rbinom(20,1,.5)
//' ##Hyperparameter for number of split points
//' Poi=5
//'##Number of iterations for MCMC
//'B=200
//'BayesPiecewiseLinearLogHazard( Y, I, Poi,  B)
//'@export
//[[Rcpp::export]]
List BayesPiecewiseLinearLogHazard( arma::vec Y, //Survival Times
                                    arma::vec I1, //Censoring Indicators
                                    double Poi, //Prior Number of split points
                                    int B // Number of iterations to perform
){
  
  //Upper Limit on Accept
  double High=.4;
  //Lower Limit on Accept
  double Low=.2;
  
  
  //Prior variance on first log-hazard height
  double Var1=25;
  //Prior probability of slab
  double PSlab = .9;
  //Prior probability of spike
  double PSpike=PSlab/(1-PSlab);
  
  
  //Prior Params, make user controlled later
  //For loop quantiites
  int m=0;
  int b=0;
  int k=0;
  int j=0;
  //Important double quantities used in MCMC
  //  int B1=B/10;  //Burnin first 90 %
  //  int B1=B/10;
  int B1=B/10;
  
  
  //Storing the restricted mean
  arma::vec MeanStore(B1);
  MeanStore.zeros();
  //Max To Store, really is 1 less here
  int Lmax = 30;
  //Storage for S vector
  arma::mat sstore(B1,Lmax+2);
  //Storage for lambda
  arma::mat lamstore(B1,Lmax+2);
  //Storage for L
  arma::vec Lstore(B1);
  //Storage for sigma (dependence)
  arma::vec sigstore(B1);
  //Initialize S, lam and J
  //Make Vector of split points
  arma::vec s(Lmax+2);
  arma::vec lam(Lmax+2);
  s.zeros();
  lam.zeros();
  //Initialize sigma_\lambda and 
  double sig=1;
  //Used for birth moves
  double Birth=0;
  int Spot=0;
  //Used for metropolis hastings moves
  double U=0;
  double alpha=0;
  //This needed for siglam sampling
  double cum1=0;
  //multiplicative perturbation
  double U1=0;
  //used to find what interval something is in
  int which1=0;  
  
  
  
  //Copy Paste these to top
  double mean2=0;
  
  
  
  
  
  double mean1=0;
  double Con=0;
  
  
  
  //Storage for printout
  NumericVector z9(2);
  //Contains the maximum observed event time
  double m1=MaxVec(Y) ;
  //Start with split points at every $L/4$
  s[1]=m1/4;
  s[2]=2*m1/4;
  s[3]=3*m1/4;
  s[4]=m1;
  //Innitialize L with 3 split points
  int L = 3;
  
  
  //Make a sequence x1 that we'll use to get approximate means
  double max12 = ceil(m1*100)/100;
  //Seq size
  int size1 = max12/.01;
  
  
  
  arma::vec x1(size1);
  for(m=0;m<size1;m++){
    x1(m)=.01*(m+1);
  }
  
  
  
  
  
  //Make A sequence for 
  
  
  
  //Now Our initial split point vector is set
  lam(0)=as_scalar(arma::randn(1));
  lam(1)=as_scalar(arma::randn(1));
  lam(2)=as_scalar(arma::randn(1));
  lam(3)=as_scalar(arma::randn(1));
  lam(4)=as_scalar(arma::randn(1));
  
  
  //Double for adding new lambda
  double NewLam1=0;
  double NewLam2=0;
  
  
  
  
  
  //This piece vector contains indicators for if \lambda_l = \lambda_{l-1}
  arma::vec Piece=lam;
  //Sets the piece vector to 0s
  Piece.zeros();
  arma::vec lam1=lam;
  //Sampling counters
  double bprop=1;
  double Intb=1;
  double Numb=2;
  arma::vec lamprop=lam;
  arma::vec sprop=s;
  
  double prob1=0;
  double med2=0;
  
  
  
  //Holds storage index
  int StoreInx=0;
  
  //Lambda Tuning Parameters
  double LamC = .5;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam = 1;  //Number of proposals made
  double ALam = 1; //Number of proposals accepted
  //Here's for lam0, which ALWAYS exists
  double LamC1 = .5;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam1 = 1;  //Number of proposals made
  double ALam1 = 1; //Number of proposals accepted
  //Here's for lam_J+1 which ALWAYS exists
  double LamC2 = .5;
  //Note, it really starts at .25 since it will be halved immediately
  double NLam2 = 1;  //Number of proposals made
  double ALam2 = 1; //Number of proposals accepted
  
  
  
  
  for(m=0;m<B;m++){
    
    //Every 250 observations, let's adjust variance until burnin
    if( (m%250==0) & (m<(B-B1))){
      //Tune proposal variance for lambda here
      if((ALam/NLam)>High){
        //Double the variance, acceptance is too high
        LamC = min1(LamC*2,2);
        
      }
      
      
      if((ALam/NLam)<Low){
        //Halve the variance, acceptance is too Low
        LamC = max1(LamC/2,.0625);
        
      }
      
      //Reset counters
      ALam =1;
      NLam=1;
      
      
      
      
      
      //Tune proposal variance for lambda here
      if((ALam1/NLam1)>High){
        //Double the variance, acceptance is too high
        LamC1 = min1(LamC1*2,2);
        
      }
      
      
      if((ALam1/NLam1)<Low){
        //Halve the variance, acceptance is too Low
        LamC1 = max1(LamC1/2,.0625);
        
      }
      
      
      
      
      
      //Reset counters
      ALam1 =1;
      NLam1=1;
      
      
      
      //Tune proposal variance for lambda here
      if((ALam2/NLam2)>High){
        //Double the variance, acceptance is too high
        LamC2 = min1(LamC2*2,2);
        
      }
      
      
      if((ALam2/NLam2)<Low){
        //Halve the variance, acceptance is too Low
        LamC2 = max1(LamC2/2,.0625);
        
      }
      
      
      //Reset counters
      ALam2 =1;
      NLam2=1;
      
      
      
      
      
    }
    
    
    
    
    //Sample Lambda_{0}
    lamprop = lam;
    //Default for now, could make it adaptive
    lamprop(0)=lam(0)+as_scalar(arma::randu(1))*LamC1*2 - LamC1;
    //Log-Acceptance ratio
    alpha=LikePLLH(Y,I1, s, lamprop , L) -    LikePLLH(Y,I1, s, lam , L);
    //adjust for proposal ratio
    alpha = alpha - .5*pow(lamprop[0],2)/Var1 + .5*pow(lam[0],2)/Var1 - .5*pow(lamprop[1]-lamprop[0],2)/sig + .5*pow(lam[1]-lam[0],2)/sig;
    
    
    
    U=log(as_scalar(arma::randu(1)));
    if(U<alpha){
      lam(0)=lamprop(0);
      ALam1=ALam1+1;
    }
    NLam1=NLam1+1;
    
    
    
    
    
    //Sample Lambda_{J+1}
    lamprop = lam;
    //Default for now, could make it adaptive
    lamprop(L+1)=lam(L+1)+as_scalar(arma::randu(1))*LamC2*2 - LamC2;
    //Log-Acceptance ratio
    alpha=LikePLLH(Y,I1, s, lamprop , L) -    LikePLLH(Y,I1, s, lam , L);
    //Adjust for proposal ratio
    alpha = alpha - .5*pow(lamprop[L+1]-lamprop[L],2)/sig + .5*pow(lam[L+1]-lam[L],2)/sig   ;
    
    U=log(as_scalar(arma::randu(1)));
    
    z9[0]=alpha;
    Rf_PrintValue(wrap(z9));
    
    if(U<alpha){
      lam(L+1)=lamprop(L+1);
      ALam2=ALam2+1;
    }
    NLam2=NLam2+1;
    
    
    
    
    //Sample Lambda | L>0, Here's the middle ones
    if(L>0){
      for(j=1;j<(L+1);j++){
        lamprop = lam;
        //Default for now, could make it adaptive
        lamprop(j)=lam(j)+as_scalar(arma::randu(1))*LamC*2 - LamC;
        //Log-Acceptance ratio
        alpha=LikePLLH(Y,I1, s, lamprop , L) -    LikePLLH(Y,I1, s, lam , L);
        //Adjust for proposal ratio
        alpha = alpha - .5*pow(lamprop[j]-lamprop[j-1],2)/sig + .5*pow(lam[j]-lam[j-1],2)/sig  - .5*pow(lamprop[j+1]-lamprop[j],2)/sig + .5*pow(lam[j+1]-lam[j],2)/sig   ;
        U=log(as_scalar(arma::randu(1)));
        if(U<alpha){
          lam(j)=lamprop(j);
          ALam=ALam+1;
        }
        NLam=NLam+1;
        
        
      }
    }
    
    
    
    
    
    
    
    
    //Sample Siglam directly
    //Set sum of quadratic lambda differences to 0
    cum1=0;
    //Do a for loop over the intervals for (lambda_j-lambda_{j-1})^2
    for(j=1;j<(L+2);j++){
      cum1 =cum1 + pow(lam[j]-lam[j-1],2);
    }
    //Sample sigma_lambda directly from a Gamma distribution
    sig = 1/R::rgamma((L+2)/2 , cum1/2  );
    
    
    
    //Shuffle Splitpoints
    if(L>0){
      for(j=1; j<(L+1); j++){
        sprop = s;
        //Randomly draw s_j \sim U[s_{j-1},s_{j+1}]
        sprop[j]=as_scalar(arma::randu(1))*(sprop[j+1]-sprop[j-1]) + sprop[j-1];
        //Calculate Acceptance Ratio
        alpha =  log(sprop[j+1]-sprop[j])+log(sprop[j]-sprop[j-1])-log(s[j+1]-s[j])-log(s[j]-s[j-1]) +   LikePLLH(Y,I1,  sprop, lam , L) -    LikePLLH(Y,I1,   s, lam , L);
        //Draw Uniform 
        U = log(as_scalar(arma::randu(1)));
        if(U<alpha){
          //Metropolis Hastings
          s[j]=sprop[j];
        }
        
        
        
        
      }
    }
    
    
    //Death move
    if(L>0){
      //Draw proposed split point to delete
      Spot=SampleDeath(L);
      sprop.zeros();
      lamprop.zeros();
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
      }
      
      //We just straight up delete \lamda_spot and s[spot] and make adjustments
      if(L>1){
        for(j=0;j<Spot;j++){
          lamprop[j]=lam[j];
        }
        for(j=Spot;j<(lam.n_rows-1);j++){
          lamprop[j]=lam[j+1];
        }
      }else{
        lamprop[0]=lam[0];
        lamprop[1]=lam[2];
      }
      for(j=Spot; j<(s.n_rows-1);j++){
        sprop[j]=s[j+1];
      }
      //Likelihood ratio
      alpha =    LikePLLH(Y,I1,  sprop, lamprop , L-1) -    LikePLLH(Y,I1, s, lam , L);
      //Prior Ratio
      //Poisson
      alpha = alpha  -log(Poi) + log1(L);
      //S Prior
      alpha = alpha +2*log(m1) + log(s[Spot+1]-s[Spot-1]) - log1(2*L+1) - log1(2L) - log(s[Spot]-s[Spot-1])-log(s[Spot+1]-s[Spot]);
      //Lambda Prior, we DROPPED one
      if(L==1){
        //Removing will drop the sampler to 0 split points
        alpha = alpha - .5*pow(lam[2]-lam[0],2)/sig + .5*pow(lam[2]-lam[1],2)/sig + .5*pow(lam[1]-lam[0],2)/sig ;
      }else{
        alpha = alpha -.5*pow(lam[Spot+1]-lam[Spot-1],2)/sig + .5*pow(lam[Spot+1]-lam[Spot],2)/sig+.5*pow(lam[Spot]-lam[Spot-1],2)/sig;
      }
      
      
      //Random Perturbation to maintain balance between parameter spaces
      U1 = as_scalar(arma::randu(1));
      ///Add log Jacobian Here
      alpha = alpha +log(U1*(1-U1));
      //Metropolis Hastings Green
      U=log(as_scalar(arma::randu(1)));
      
      if(U<alpha){
        s=sprop;
        L=L-1;
        lam=lamprop;
      }
    }
    
    
    
    //Birth move Here
    if(L<(Lmax-1)){
      //What interval is it in?
      //This generates the interval location that the new proposed split point is located
      Spot = SampleBirth(s,L)+1;
      //Set Lambda to 0
      lamprop.zeros();
      //Now generate the new proposed split point location
      Birth = as_scalar(arma::randu(1))*(s[Spot]-s[Spot-1])+s[Spot-1];
      //Random Perturbation for dimension matching
      U1=as_scalar(arma::randu(1));
      //Slope
      
      //Old Slope for log hazard
      NewLam1= (lam[Spot]-lam[Spot-1])/(s[Spot]-s[Spot-1]);
      //Compute the new slope that we can use to get the new lambda value
      //This quantity comes from the MCMC section on how to obtain new slopes
      NewLam1 = NewLam1 - log((1-U1)/U1)*(s[Spot]-Birth)/(s[Spot]-s[Spot-1]);
      //We can now use this slope to obtain the new hazard height at the time point Birth
      NewLam2 = NewLam1*(Birth-s[Spot-1])+lam[Spot-1];
      //Now NewLam2 contains the hazard at the time point Birth
      //Now we have the interval of the new split in Spot and the actual location in Birth
      sprop.zeros();
      //Here we're going to fill in our new vector
      //Let's add Birth to the Spot location and push back the rest of sprop
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
        lamprop[j]=lam[j];
      }
      sprop[Spot]=Birth;
      for(j=(Spot+1);j<(s.n_rows-1);j++){
        sprop[j]=s[j-1];
      }
      //Log of hazard height at Birth.
      lamprop[Spot]=NewLam2;
      for(j=(Spot+1);j<(lam.n_rows-1);j++){
        lamprop[j]=lam[j-1];
      }
      
      
      //Now we have our new proposal vector, evaluate it!
      //Like Ratio
      alpha =    LikePLLH(Y,I1, sprop, lam , L+1) -    LikePLLH(Y,I1, s, lam , L);
      //Add proposal ratio
      //Poisson
      alpha= alpha + log(Poi) - log1(L+1) ;
      // S proposal
      alpha = alpha + log1(2*L+3)+log1(2*L+2)+log(Birth-s[Spot-1])+log(s[Spot]-Birth) - 2*log(m1)-log(s[Spot]-s[Spot-1]);
      //Add proposal ratio for \lambda
      if(L==0){
        alpha = alpha  - .5*pow(lamprop[1]-lamprop[0],2)/sig - .5*pow(lamprop[2]-lamprop[1],2)/sig + .5*pow(lam[1]-lam[0],2)/sig ;
      }else{
        alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig + .5*pow(lam[Spot]-lam[Spot-1],2)/sig;
      }
      
      //Add log Jacobian
      alpha = alpha -log(U1*(1-U1));
      //Uniform for Metropolis-Hastings
      U=log(as_scalar(arma::randu(1)));
      if(U<alpha){
        //We have to set the entire vector here
        s=sprop;
        lam=lamprop;
        L=L+1;
      }
    }
    
    
    
    
    
    
    
    
    
    if(m>(B-B1-1)){
      //Store Values in Matrix only after certain burnin
      StoreInx = m-B+B1;
      
      for(j=0;j<sstore.n_cols;j++){
        sstore(StoreInx,j) = s(j);
      }
      for(j=0;j<lamstore.n_cols;j++){
        lamstore(StoreInx,j) = lam(j);
      }
      Lstore[StoreInx]=L;
      //  sigstore[StoreInx]=sig;
      sigstore[StoreInx]=sig;
      
      //Can compute the mean here.
      MeanStore(StoreInx)=ApproxMean(x1,s,lam,L);
      
      
      
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
  
  
  
  
  
  
  
  
  List z1 = List::create(sstore,lamstore,Lstore,sigstore,MeanStore);
  
  
  return(z1);
  
  
}
