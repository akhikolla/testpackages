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






//Returns the maximum from a vector
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


double cumsumlog(int Rate){
  double Prod =0;
  int m=1;
  double m1=1;

  for(m=2;m<(Rate+1);m++){
    m1=m;
    Prod = Prod + log(m1);
  }

  return(Prod);

}




//Function for computing the likelihood of a Lognormal GLM
double Like( arma::vec Y, //Survival Times
             arma::vec Rates, //Rates, Number of severity
             arma::vec s, //Vector containing the split point locations
             arma::vec lam, //Vector containing the log hazard heights
             arma::vec poi, //Vector of log poisson rates for the # of casualties for a given event time
             int J //Numer of split points
){

  int m=0;
  int l=0;
  double LogL=0;

arma::vec sum1(J+1);
sum1.zeros();

double m1=0;

    //Cycle through each interval and obtain survival estimates and add hazard heights
  for(l=0; l<(J+1);l++){
    for(m=0; m<Y.n_rows;m++){
      LogL = LogL - max1(0,min1(s(l+1),Y(m))-s(l))*exp(lam[l]);


      if( (Y(m)>=s[l]) && (Y(m)<s[l+1]) ){

        m1=Rates[m];
      LogL = LogL + lam[l] + Rates[m]*poi[l]-exp(poi[l]) - cumsumlog(m1);

        //If interval HAS NONE BIG PENALTY!!!

        sum1(l)=sum1(l)+1;



      }




    }

  }

if(5<3){
  for(m=0;m<(J+1);m++){
    if(sum1(m)<4){
      LogL = -100000;
    }
  }
}

  return(LogL);



}




arma::vec GetBounds(arma::vec Y, double s1, double s2){

  arma::vec Bounds(2);

  double c=s1;
  //Get first bound
  double sum1 = 0;


  while(sum1<5 ){
    c=c+.1;
    sum1 = sum(Y<c)-sum(Y<=s1);
    if(c>s2){
      c=-10;
      break;

    }

  }

  Bounds(0)=c;

  sum1=0;
  c=s2;

  while(sum1<5){
    c=c-.1;
    sum1 = sum(c<=Y) - sum(Y>s2);
    if(c<s1){
      c=-10;
      break;

    }
  }

  Bounds(1)=c;


  return(Bounds);




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
      if( (U>cumprob[m-1]) && (U<cumprob[m]) ){
        Which1=m;
      }
    }

  }

  return(Which1);

}


int factorial(int N){
  int m=0;
  int cum1=1;
  for(m=2;m<(N+1);m++){
    cum1 = cum1*m;
  }
  return(cum1);
}

//'  C++ Sampling Function for MCMC
//'
//' C++ Sampling Function used in the PieceExpIntensity function.
//' @param Y   Vector containing observed event times.
//' @param Rates Vector containing poisson count intensities.
//' @param B Number of iterations to run the MCMC with half burned in.
//' @param Poi Prior mean number of split points,
//' @return A list of all posterior quantities.
//' @examples
//' B=1000
//' n=100
//' Y=rexp(n,1)
//' Rates=Y
//' Rates[Y<.5]=rpois(sum(Y<.5),20)
//' Rates[Y>.5]=rpois(sum(Y>.5),3)
//' Poi=10
//' PieceExpIntensity2(Y,Rates,B,Poi)
//' @export
//[[Rcpp::export]]
List PieceExpIntensity2( arma::vec Y, //Survival Times
                arma::vec Rates, //Number of casualties/etc at each event time
                int B, // Number of iterations to perform
                double Poi
){


  //Prior Params, make user controlled later

  //For loop quantiites
  int m=0;
  int b=0;
  int k=0;
  int j=0;
  //Important double quantities used in MCMC
  int B1=B/2;


  //Innitialize Parameters and storage matrices
  double beta = 1  ;



  arma::vec betastore(B1);
  arma::vec MeanStore(B1);
  arma::vec MedianStore=MeanStore;
  //Max To Store, really is 1 less here
  int Lmax = 20;
  arma::mat sstore(B1,Lmax+2);
  arma::mat lamstore(B1,Lmax+1);
  arma::mat poistore(B1,Lmax+1);
  arma::vec Lstore(B1);
  arma::vec sigstore(B1);
  arma::vec sig1store(B1);
  //Initialize S, lam and J
  int L = 3;
  arma::vec s(Lmax+2);
  arma::vec lam(Lmax+1);
  s.zeros();
  lam.zeros();
  arma::vec poi = lam;
  arma::vec poiprop=poi;
//Now set the first few entries of the poisson and lambda vector to random values
lam[0]=as_scalar(arma::randn(1));
lam[1]=as_scalar(arma::randn(1));
lam[2]=as_scalar(arma::randn(1));
lam[3]=as_scalar(arma::randn(1));
poi[0]=as_scalar(arma::randn(1));
poi[1]=as_scalar(arma::randn(1));
poi[2]=as_scalar(arma::randn(1));
poi[3]=as_scalar(arma::randn(1));







  arma::vec I = Rates;


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


  //Contains the maximum observed event time
  double m1 = MaxVec(Y);
  //Start with split points at every $L/4$
  s[1]=m1/4;
  s[2]=2*m1/4;
  s[3]=3*m1/4;
  s[4]=m1;
  double U2=0;

  double NewLam1=0;
  double NewLam2=0;
  double U1=0;
  int which1=0;


  arma::vec lam1=lam;

  double bprop=1;
  double Intb=1;
  double Numb=2;
  double sig1=1;
  arma::vec lamprop=lam;
  arma::vec sprop=s;

  double prob1=0;
  double med2=0;
arma::vec Bounds(2);

double L1=1;


  int StoreInx=0;


  for(m=0;m<B;m++){

    if(m%1000==0){
      Bounds[0]=m;
      Bounds[1]=m;
      Rf_PrintValue(wrap(Bounds));
      Rprintf("Iterations");
    }






    for(j=0;j<(L+1);j++){

      lamprop = lam;


      lamprop(j)=lam(j)+as_scalar(arma::randu(1))*.5 - .25;





      alpha=Like(Y,Rates, s, lamprop , poi, L) -    Like(Y,Rates, s, lam ,poi, L);

      if(L>0){
      if(j==0){
        alpha = alpha - .5*(pow(lamprop[j],2)-pow(lam[j],2))/25 - .5*pow(lamprop[j+1]-lamprop[j],2)/sig + .5*pow(lam[j+1]-lam[j],2)/sig;
      }else{
        if(j==L){
          alpha = alpha - .5*pow(lamprop[j]-lamprop[j-1],2)/sig + pow(lam[j]-lam[j-1],2)/sig    ;

        }else{
        alpha = alpha - .5*pow(lamprop[j]-lamprop[j-1],2)/sig + .5*pow(lam[j]-lam[j-1],2)/sig  - .5*pow(lamprop[j+1]-lamprop[j],2)/sig + .5*pow(lam[j+1]-lam[j],2)/sig   ;
      }
      }

      }else{

        alpha = alpha - .5*(pow(lamprop[j],2)-pow(lam[j],2))/25 ;



      }


      U = log(as_scalar(arma::randu(1)));





      if(U<alpha){
        lam=lamprop;
      }


    }



    if(L>0){
      //Sample Siglam
      cum1=0;
      for(j=1;j<(L+1);j++){
        cum1 =cum1 + pow(lam(j)-lam(j-1),2);
      }



      sig = 1/R::rgamma(L/2 + 1, cum1/2 );

    }










    for(j=0;j<(L+1);j++){

      poiprop = poi;


      poiprop[j]=poi[j]+as_scalar(arma::randu(1))*.5 - .25;



      alpha=Like(Y,Rates, s, lam , poiprop, L) -    Like(Y,Rates, s, lam ,poi, L);



      if(L>0){
      if(j==0){
        alpha = alpha - .5*(pow(poiprop[j],2)-pow(poi[j],2))/25 - .5*pow(poiprop[j+1]-poiprop[j],2)/sig1 + .5*pow(poi[j+1]-poi[j],2)/sig1;
      }else{
        if(j==L){
          alpha = alpha - .5*pow(poiprop[j]-poiprop[j-1],2)/sig1 + pow(poi[j]-poi[j-1],2)/sig1    ;

        }else{
          alpha = alpha - .5*pow(poiprop[j]-poiprop[j-1],2)/sig1 + .5*pow(poi[j]-poi[j-1],2)/sig1  - .5*pow(poiprop[j+1]-poiprop[j],2)/sig1 + .5*pow(poi[j+1]-poi[j],2)/sig1   ;
        }
      }

      }else{

        alpha = alpha - .5*(pow(poiprop[j],2)-pow(poi[j],2))/25 ;
      }




      U = log(as_scalar(arma::randu(1)));





      if(U<alpha){
        poi=poiprop;

      }


    }

    //Sample Siglam
    if(L>0){
      cum1=0;

      for(j=1;j<(L+1);j++){
        cum1 =cum1 + pow(poi[j]-poi[j-1],2);
      }

      sig1 = 1/R::rgamma(L/2 + 1, cum1/2 );

    }






    //Shuffle Splits
if(L>0){
    for(j=1; j<(L+1); j++){
      sprop = s;
      //Draw new proposal for s_j, but make sure it's between s_j and s_j+1
       sprop[j]=as_scalar(arma::randu(1))*(sprop[j+1]-sprop[j-1])+sprop[j-1];




      alpha =  log(sprop[j+1]-sprop[j])+log(sprop[j]-sprop[j-1])-log(s[j+1]-s[j])-log(s[j]-s[j-1]) +   Like(Y,Rates,  sprop, lam ,poi, L) -    Like(Y,Rates,   s, lam ,poi, L);

      U = log(as_scalar(arma::randu(1)));







      if(U<alpha){

        s=sprop;
      }




    }
}

    //Add proposal and MHG ratios here
    //Birth move Here


    if(L<(Lmax-1)){
      //What interval is it in?

      lamprop.zeros();
      poiprop.zeros();

      Spot = SampleBirth(s)+1;


    Birth = as_scalar(arma::randu(1))*(s[Spot]-s[Spot-1])+s[Spot-1];



      //Random Perturbation
      U1=as_scalar(arma::randu(1));





      NewLam1 =  lam[Spot-1] - log((1-U1)/U1)*(s[Spot]-Birth)/(s[Spot]-s[Spot-1]);
      NewLam2 =  lam[Spot-1] + log((1-U1)/U1)*(Birth-s[Spot-1])/(s[Spot]-s[Spot-1]);



      //Now we have the interval of the new split in Spot and the actual location in Birth
      sprop.zeros();
      //Let's add Birth to the Spot location and push back the rest of sprop
      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
        lamprop[j]=lam[j];
        poiprop(j)=poi(j);
      }
      sprop[Spot]=Birth;
      for(j=(Spot+1);j<(s.n_rows-1);j++){
        sprop[j]=s[j-1];
      }

      lamprop[Spot-1]=NewLam1;
      lamprop[Spot]=NewLam2;

      for(j=(Spot+1);j<(lam.n_rows-1);j++){
        lamprop[j]=lam[j-1];
        poiprop(j)=poi[j-1];
      }


      U2=as_scalar(arma::randu(1));

   NewLam1 =  poi(Spot-1) - log((1-U2)/U2)*(s(Spot)-Birth)/(s(Spot)-s(Spot-1));
      NewLam2 =  poi(Spot-1) + log((1-U2)/U2)*(Birth-s(Spot-1))/(s(Spot)-s(Spot-1));


    poiprop(Spot-1)=NewLam1;
      poiprop(Spot)=NewLam2;

      //Propose a new lambda corresponding to this split point and a random perturbation


      //Now we have our new proposal vector, evaluate it!
      //Like Ratio
      alpha =    Like(Y,Rates, sprop, lamprop ,poiprop, L+1) -    Like(Y,Rates, s, lam ,poi, L);
      //Add proposal ratio
      //Poisson
      L1=L;

  //    alpha= alpha + log(Poi) - log(L1+1);
      // S proposal
      alpha = alpha + log(2*L1+3)+log(2*L1+2)+log(Birth-s[Spot-1])+log(s[Spot]-Birth);
      alpha = alpha - 2*log(m1)-log(s[Spot]-s[Spot-1]);
      //Perturbation
      alpha=alpha-log(U1*(1-U1)) - log(U2*(1-U2));
      //Add proposal ratio for \lambda

      if(L==0){
        alpha = alpha - .5*pow(lamprop[1]-lamprop[0],2)/sig   ;

        alpha = alpha - .5*pow(poiprop[1]-poiprop[0],2)/sig1   ;


      }else{

        if(Spot==(L+1)){
          alpha = alpha - .5*pow(poiprop[Spot]-poiprop[Spot-1],2)/sig1 - .5*pow(poiprop[Spot-1]-poiprop[Spot-2],2)/sig1 + pow(poi[Spot-1]-poi[Spot-2],2)/sig1;
          alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig + pow(lam[Spot-1]-lam[Spot-2],2)/sig;


        }else{


if(Spot==1){


  alpha = alpha - .5*pow(poiprop[Spot+1]-poiprop[Spot],2)/sig1 - .5*pow(poiprop[Spot]-poiprop[Spot-1],2)/sig1;
  alpha = alpha + .5*pow(poi[Spot]-poi[Spot-1],2)/sig1;


  alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig ;
  alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig ;


}else{




  alpha = alpha - .5*pow(poiprop[Spot+1]-poiprop[Spot],2)/sig1 - .5*pow(poiprop[Spot]-poiprop[Spot-1],2)/sig1 - pow(poiprop[Spot-1]-poiprop[Spot-2],2)/sig1;
  alpha = alpha + .5*pow(poi[Spot]-poi[Spot-1],2)/sig1 + .5*pow(poi[Spot-1]-poi[Spot-2],2)/sig1;


  alpha = alpha - .5*pow(lamprop[Spot+1]-lamprop[Spot],2)/sig - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig - pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;
  alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;





}



        }




      }



  //Add last interval
  alpha = alpha - .5*pow(poiprop[0],2)/25 + .5*pow(poi[0],2)/25 ;
      alpha = alpha  - .5*pow(lamprop[0],2)/25 + .5*pow(lam[0],2)/25 ;





      U=log(as_scalar(arma::randu(1)));


      if(U<alpha){
        s=sprop;
        lam=lamprop;
        L=L+1;
        poi=poiprop;
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
      poiprop.zeros();

      for(j=0;j<Spot;j++){
        sprop[j]=s[j];
      }


      if(L>1){

        if(Spot>1){
        for(j=0;j<(Spot-1);j++){
          lamprop[j]=lam[j];
          poiprop[j]=poi[j];
        }
        }

        for(j=Spot;j<(lam.n_rows-1);j++){
          lamprop[j]=lam[j+1];
          poiprop[j]=poi[j+1];
        }

        //New lambda is a weighted average of the old lambdas
        lamprop[Spot-1] = ((s[Spot+1]-s[Spot])*lam[Spot] + (s[Spot]-s[Spot-1])*lam[Spot-1])/(s[Spot+1]-s[Spot-1]);

        poiprop[Spot-1] = ((s[Spot+1]-s[Spot])*poi[Spot] + (s[Spot]-s[Spot-1])*poi[Spot-1])/(s[Spot+1]-s[Spot-1]);


      }else{
        lamprop[0]=((s[2]-s[1])*lam[1]+(s[1]-s[0])*lam[0])/(s[2]-s[0]);
        poiprop[0]=((s[2]-s[1])*poi[1]+(s[1]-s[0])*poi[0])/(s[2]-s[0]);

      }

      for(j=Spot; j<(s.n_rows-1);j++){
        sprop[j]=s[j+1];
      }



      //Now we have our new proposal vector, evaluate it!
      //Like Ratio
      //Prior Ratio
      alpha =    Like(Y,Rates,  sprop, lamprop ,poiprop, L-1) -    Like(Y,Rates, s, lam , poi, L);


      L1=L;
      //Poisson
    //  alpha = alpha  -log(Poi) + log(L1);
      //S Prior
      alpha = alpha + 2*log(m1) + log(s[Spot+1]-s[Spot-1]) - log(2*L1+1) - log(2*L1) - log(s[Spot]-s[Spot-1])-log(s[Spot+1]-s[Spot]);
      //Lambda Prior, we DROPPED one

      if(L==1){
        //Removing will drop the sampler to 0 split points

        alpha = alpha + .5*pow(lam[1]-lam[0],2)/sig ;

        alpha = alpha + .5*pow(poi[1]-poi[0],2)/sig1 ;


      }else{



          if(Spot>1){

          if(Spot==(L+1)){
            alpha = alpha  + .5*pow(poi[Spot]-poi[Spot-1],2)/sig1 + .5*pow(poi[Spot-1]-poi[Spot-2],2)/sig1 - .5*pow(poiprop[Spot-1]-poiprop[Spot-2],2)/sig1;
            alpha = alpha + .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig  - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig;





          }else{

            //Split point is not first or last

            alpha = alpha + .5*pow(poi[Spot+1]-poi[Spot],2)/sig1 + .5*pow(poi[Spot]-poi[Spot-1],2)/sig1 + .5*pow(poi[Spot-1]-poi[Spot-2],2)/sig1;
            alpha = alpha +  .5*pow(lam[Spot+1]-lam[Spot],2)/sig +  .5*pow(lam[Spot]-lam[Spot-1],2)/sig + .5*pow(lam[Spot-1]-lam[Spot-2],2)/sig;

          alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig  - .5*pow(poiprop[Spot]-poiprop[Spot-1],2)/sig1;

            alpha = alpha - .5*pow(lamprop[Spot-1]-lamprop[Spot-2],2)/sig  - .5*pow(poiprop[Spot-1]-poiprop[Spot-2],2)/sig1;



          }


          }else{

            //First Interval is changing

            alpha = alpha + .5*pow(poi[Spot+1]-poi[Spot],2)/sig1+ .5*pow(poi[Spot]-poi[Spot-1],2)/sig1 ;
            alpha = alpha + .5*pow(lam[Spot+1]-lam[Spot],2)/sig + .5*pow(lam[Spot]-lam[Spot-1],2)/sig;

            alpha = alpha - .5*pow(lamprop[Spot]-lamprop[Spot-1],2)/sig  - .5*pow(poiprop[Spot]-poiprop[Spot-1],2)/sig1;



          }





      }



      //Add the 0 entry here
    alpha= alpha  - .5*pow(lamprop[0],2)/25 + .5*pow(lam[0],2)/25 ;
      alpha = alpha  - .5*pow(poiprop[0],2)/25 + .5*pow(poi[0],2)/25 ;




      //Random Perturbation
      U1 = as_scalar(arma::randu(1));
      U2 = as_scalar(arma::randu(1));

      alpha = alpha + log(U1)+log(1-U1) + log(U2)+log(1-U2);





      U=log(as_scalar(arma::randu(1)));


      if(U<alpha){
        s=sprop;
        L=L-1;
        lam=lamprop;
        poi=poiprop;
      }




    }













    if(m>(B1-1)){
      //Store Values in Matrix
      StoreInx = m-B1;

      for(j=0;j<sstore.n_cols;j++){
        sstore(StoreInx,j) = s(j);
      }

      for(j=0;j<lamstore.n_cols;j++){
        lamstore(StoreInx,j) = lam(j);
        poistore(StoreInx,j) = poi(j);
      }


      Lstore[StoreInx]=L;
      sigstore[StoreInx]=sig;
sig1store[StoreInx]=sig1;



    }









  }









  List z1 = List::create(sstore,lamstore,poistore,Lstore,sigstore, sig1store);


  return(z1);


}














