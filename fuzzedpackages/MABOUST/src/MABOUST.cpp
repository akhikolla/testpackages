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

arma::mat MATMULT(arma::mat X, arma::mat Y){

  arma::mat PROD(X.n_rows,Y.n_cols);

  int k;
  int j;
  int m=0;

  for(j = 0 ;j<X.n_rows;j++){
    for(k=0;k<PROD.n_cols;k++){
      PROD(j,k)=0;
      for(m=0;m<X.n_cols;m++){
        PROD(j,k)=PROD(j,k)+X(j,m)*Y(m,k);
      }



    }
  }


  return(PROD);


}



int IsAdmissable(arma::mat gamma1,arma::mat gamma2){
  int ADM=1;
  arma::vec z9(2);
  int m1=1;
  int k=1;
  for(k=0;k<gamma1.n_rows;k++){
    for(m1=0;m1<(gamma1.n_cols-1);m1++){


      if((gamma1(k,m1)>gamma1(k,m1+1))|| ((gamma1(k,m1)+gamma2(k,m1))>(gamma1(k,m1+1)+gamma2(k,m1+1)))){
        //Cumulative hazard is violated.
        ADM = 0;
        break;
      }

    }}

  return(ADM);

}


int Sample1(int G){
  arma::vec groupprob(G);
  int m=0;
  double G1=G;
  groupprob.zeros();
  groupprob = groupprob + 1/G1;


  return(Sample2(groupprob));

}


int SampleSpike(arma::vec SPIKEHOLD,int k){

  //How many are unclustered right now?
  int count = 0;
  int m1=0;
  for(m1 =0 ; m1<SPIKEHOLD.n_rows;m1++){
    if(SPIKEHOLD(m1)==m1){
      count++;
    }
  }

  count=count-1;

  if(count==0){
    return(0); //Cluster on 0...
  }else{

    arma::vec UNCLUST(count);

    count = 0;
    for(m1=0;m1<SPIKEHOLD.n_rows;m1++){
      if(SPIKEHOLD(m1)==m1 && m1!=k){
        UNCLUST(count)=m1;
        count++;
      }
    }

    return(UNCLUST(Sample1(UNCLUST.n_rows)));


  }
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



arma::vec GetBoundariesALPHA(int m,  arma::vec beta){

  arma::vec HOLD(2);
  double bound1=0;
  double bound2=0;
  if(m==0){
    bound1 = -10;
    bound2 = beta(m+1);
    HOLD(0)=bound1;
    HOLD(1)=bound2;

  }else{
    if(m==(beta.n_rows-1)){
      bound1 = beta(m-1);
      bound2 = 10;
      HOLD(0)=bound1;
      HOLD(1)=bound2;

    }else{
      //In the middle
      bound1 =beta(m-1);
      bound2 = beta(m+1);
      HOLD(0)=min1(bound1,bound2);
      HOLD(1)=max1(bound1,bound2);

    }


  }

  return(HOLD);

}

double TruncNormALPHA(int m,  arma::vec beta, double c1){
  //Truncation boundaries
  arma::vec BOUNDARIES = GetBoundariesALPHA(m,beta);


  //Generate Uniform variable
  double U = as_scalar(arma::randu(1));

  //Get value to plug into qnorm
  double X = U*R::pnorm5(BOUNDARIES(1),beta(m),c1,1,0)+(1-U)*R::pnorm5(BOUNDARIES(0),beta(m),c1,1,0);

  X=R::qnorm5(X,beta(m),c1,1,0);

  return(X);

}




arma::vec GetBoundariesBETA(int m, arma::vec beta){
  //
  arma::vec alpha=beta;
  alpha.zeros();
  arma::vec HOLD(2);
  double bound1=0;
  double bound2=0;
  if(m==0){
    bound1 = -10;
    bound2 = beta(m+1);
    HOLD(0)=bound1;
    HOLD(1)=bound2;

  }else{
    if(m==(alpha.n_rows-1)){
      bound1 = beta(m-1);
      bound2 = 10;
      HOLD(0)=bound1;
      HOLD(1)=bound2;

    }else{
      //In the middle
      bound1 = beta(m-1);
      bound2 = beta(m+1);
      HOLD(0)=bound1;
      HOLD(1)=bound2;

    }


  }


  return(HOLD);

}

double TruncNormBETA( int m, arma::vec beta, double c1){


  //Truncation boundaries
  arma::vec BOUNDARIES = GetBoundariesBETA(m,beta);

  //Generate Uniform variable
  double U = as_scalar(arma::randu(1));

  //Get value to plug into qnorm
  double X = U*R::pnorm5(BOUNDARIES(1),beta(m),c1,1,0)+(1-U)*R::pnorm5(BOUNDARIES(0),beta(m),c1,1,0);

  X=R::qnorm5(X,beta(m),c1,1,0);

  return(X);

}



















//Clusters on a random group. TEST THIS
int GetRandGroup(arma::vec INC){ //Vector containing stopped indicators)


  //How many are UNCLUSTERED AND not stopped
  double NumGroup = sum(INC);
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
    sum1 = sum1+INC(k);
  }





  return(k);

}







double LIKECOV(arma::vec Y,  //Vector of observed ordinal outcomes
               arma::vec T,  // Vector of observed treatment indicators, = 0,1 ,2, ..., nTreat
               arma::mat X, //MAtrix of covariates
               arma::mat thetavec, //Matrix of  NPO treatment effects
               arma::vec Beta //Vector of PO covariate effects
){


  //vector containing individual log likelihood contributions
  arma::vec eta(Y.n_rows);
  eta.zeros();
  //  eta=eta+1;

  //Needed for For loops
  int m=0;
  int k=1;
  int j=0;
  int J = thetavec.n_cols;


  arma::vec VAL = MATMULT(X,Beta);


  for(m=0;m<eta.n_rows;m++){
    if(Y(m)==0){
      //First outcome
      eta(m) = exp(VAL(m)+ thetavec(T(m),Y(m)))/(1+exp(VAL(m)+thetavec(T(m),Y(m))));

    }else{
      if(Y(m)==J){
        eta(m) = 1-exp(VAL(m)+ thetavec(T(m),Y(m)-1))/(1+exp( VAL(m)+ thetavec(T(m),Y(m)-1)));

      }else{
        //Somewhere in between
        eta(m) = exp(VAL(m)+ thetavec(T(m),Y(m)))/(1+exp(VAL(m)+ thetavec(T(m),Y(m))))-exp(VAL(m)+thetavec(T(m),Y(m)-1))/(1+exp(VAL(m)+thetavec(T(m),Y(m)-1)));
      }

    }

  }


  return(sum(log(eta)));

}









//' Performs posterior sampling for the MABOUST design.
//' @param Y Ordinal Outcome Vector, labeled 1,...,J
//' @param T1 Treatment Indicator, labeled 1,...,K.
//' @param X Matrix of patient covariates.
//' @param NTreat Number of treatments in consideration, i.e. K.
//' @param NOUT Number of ordinal outcome categories, i.e. J.
//' @param PSPIKE Prior probability of a pairwise null effect.
//' @param B Number of MCMC iterations to perform.
//'@importFrom Rcpp evalCpp
//'@useDynLib MABOUST
//'@return List containing posterior distributions of \eqn{\bf{\theta}} and \eqn{\bf{\beta}}.
//'@export
//[[Rcpp::export]]
List MCMC_MABOUST( arma::vec Y,  //Vector of observed ordinal outcomes
               arma::vec T1,  // Vector of observed treatment indicators
               arma::mat X, //Matrix of covariates
               double B, //Number of iterations for MCMC
               double NTreat, //Number of treatments
               double NOUT, //Number of outcomes
               double PSPIKE
){



  double PNPO = .5; //Probability of npo

  int g=0;
  int m1=0;
  int StoreInx=0;



  int MEANS=0;

  int IntIN=0;
  int IntOUT=0;





  //Important For loop integer
  int m =0; //For the inner MCMC looprep
  int i=0; //For the patient index
  int rep=0; //For the simulation repetition.




  double B1=B;

  //Make List objects we'll use



  double NewSlope=0;
  int Which1=0;


  //Important quantities we will need in the MCMC
  double alpha=0;
  double U=0;
  double signew=0;
  double slopenew=0;
  double Munew=0;

  //For loop integers
  int j=0;
  int k=0;




  arma::mat thetavec(NTreat,NOUT-1);
  thetavec.zeros();
  arma::mat thetavecprop=thetavec;
  int count=0;


  NumericVector z9(2);
  NumericVector zprop(5);



  for(j =0;j<thetavec.n_cols;j++){
    thetavec.col(j)=thetavec.col(j)+j;
  }

  arma::mat ThetaStore(B1,(NOUT-1)*NTreat);



  ThetaStore.zeros();

  arma::mat etamat(NTreat,NOUT-1);
  etamat.zeros();
  etamat=etamat+1; //All unclustered to start
  arma::mat etamatprop = etamat;

  int m3=0;


  arma::mat ATheta(NTreat,NOUT-1);
  arma::mat NTheta(NTreat,NOUT-1);
  arma::mat varTheta(NTreat,NOUT-1);
  ATheta.zeros();
  NTheta.zeros();
  varTheta.zeros();
  ATheta=ATheta+1;
  NTheta=NTheta+2;
  varTheta=varTheta+.25;



  ///Things for beta, covariate effects on outcome
  arma::vec Beta(X.n_cols);
  Beta.zeros();
  arma::vec Betaprop=Beta;
  arma::vec ABeta=Beta+1;
  arma::vec NBeta=Beta+2;
  arma::vec VarBeta = Beta+.25;

  arma::vec SPIKEHOLD(thetavec.n_rows);
  arma::vec SPIKEHOLDPROP(thetavec.n_rows);


  for(m=0;m<SPIKEHOLD.n_rows;m++){
    SPIKEHOLD[m]=m;
    SPIKEHOLDPROP[m]=m;
  }

  arma::mat BetaStore(B1,Beta.n_rows);
  arma::mat SPIKESTORE(B1,SPIKEHOLD.n_rows);
  SPIKESTORE.zeros();
  BetaStore.zeros();

  arma::vec STORERAND(B);
  arma::vec STORERAND1(B);
  arma::vec STOREHOLD(B);

  Beta=Beta-1;

  int m2=0;



  for(m=0;m<B;m++){

    if(m<(B/2 + 2)){
      if(m%100==0){



        for(k=0;k<thetavec.n_rows;k++){
          for(m1=0;m1<thetavec.n_cols;m1++){
            if((ATheta(k,m1)/NTheta(k,m1))>.6){
              varTheta(k,m1)=min1(varTheta(k,m1)+.05,2);
            }
            if((ATheta(k,m1)/NTheta(k,m1))<.25){
              varTheta(k,m1)=min1(varTheta(k,m1)-.05,.05);
            }


          }
        }



        ATheta=ATheta.zeros()+1;
        NTheta=NTheta.zeros()+2;


        for(k=0;k<Beta.n_rows;k++){
          if((ABeta[k]/NBeta[k])>.6){
            VarBeta[k]=VarBeta[k]*2;
          }

          if((ABeta[k]/NBeta[k])<.25){
            VarBeta[k]=VarBeta[k]/2;

          }

        }

        ABeta=ABeta.zeros()+1;
        NBeta=NBeta.zeros()+1;




      }
    }

    for(k=0;k<thetavec.n_rows;k++){
      //Treatment Effects
      //Can't move Unless they are unclustered!!
      if(SPIKEHOLD(k)==k){
        //Is this currently unclustered?
        for(m1=0;m1<thetavec.n_cols;m1++){
          //Get our Alphas
          thetavecprop=thetavec;
          //   thetavecprop(k,m1)=TruncNormALPHA(m1, thetavecprop.row(0).t(), thetavecprop.row(k).t(), varTheta(k,m1));
          thetavecprop(k,m1)=TruncNormALPHA(m1, thetavecprop.row(k).t(), varTheta(k,m1));

          //Move anything clustered with this, MOVE IT!!!
          for(j=0;j<SPIKEHOLD.n_rows;j++){
            if(SPIKEHOLD(j)==k){
              thetavecprop.row(j)=thetavecprop.row(k);
            }
          }





          // TruncNormBETA(m1,thetavecprop.row(k).t()+thetavecprop.row(0).t(),varTheta(k,m1))-thetavecprop(0,m1);
          //Prior Ratio
          //  alpha =   -.5*(gammadf+1)*log(1+pow((thetavecprop(k,m1)-0)/gammascale,2)/gammadf)+.5*(gammadf+1)*log(1+pow((thetavec(k,m1)-0)/gammascale,2)/gammadf);
          //Generalized t distribution
          alpha=0;
          //Loglikelihood ratio
          alpha=alpha+LIKECOV( Y,  T1,X, thetavecprop,Beta)-LIKECOV( Y,  T1, X, thetavec,Beta);



          //Metropolis Hastings
          U=log(as_scalar(arma::randu(1)));



          if(U<alpha){
            ATheta(k,m1)=ATheta(k,m1)+1;

            thetavec=thetavecprop;


          }
          NTheta(k,m1)=NTheta(k,m1)+1;





        }
      }




    }



    //Adjusting for all clustering indications
    k=Sample1(thetavec.n_rows);


    STORERAND(m)=k;
    STOREHOLD(m)=SPIKEHOLD(k);

    if(SPIKEHOLD(k)==k){
      //Try SPIKING it at 0

      thetavecprop=thetavec;
      //Sample another random row to cluster it with
      m1 = SampleSpike(SPIKEHOLD,k);

      SPIKEHOLDPROP=SPIKEHOLD;
      SPIKEHOLDPROP(k)=m1;

      STORERAND1(m)=m1;



      //Not clustered on treatment 0
      for(j=0;j<thetavecprop.n_cols;j++){
        thetavecprop(k,j)=thetavecprop(m1,j);
      }






      //Move everything with it...
      //Move anything clustered with this previously...
      for(j=0;j<SPIKEHOLD.n_rows;j++){
        if(SPIKEHOLDPROP(j)==k){
          SPIKEHOLDPROP(j)=m1;

          for(m2=0;m2<thetavecprop.n_cols;m2++){
            thetavecprop(j,m2)=thetavecprop(k,m2);
          }
        }
      }




      alpha=log(PSPIKE)-log(1-PSPIKE); //Prior for spike and slab
      //Loglikelihood ratio
      alpha=alpha+LIKECOV( Y,  T1,X, thetavecprop,Beta)-LIKECOV( Y,  T1, X, thetavec,Beta);

      //Metropolis Hastings
      U=log(as_scalar(arma::randu(1)));





      if(U<alpha){
        thetavec=thetavecprop;
        SPIKEHOLD=SPIKEHOLDPROP;


      }



    }else{

      STORERAND1(m)=-7;

      //Propose unclustering a treatment vector
      thetavecprop=thetavec;
      //uncluster one entry of this vector
      j=Sample1(thetavecprop.n_cols);
      //Uncluster one value
      thetavecprop(k,j)=TruncNormALPHA(j, thetavecprop.row(k).t(),  1);




      //Nothing should be clustered on group k, so it doesn't need to move anything else.
      //Group k is just clustered ON something else.



      alpha=log(1-PSPIKE)-log(PSPIKE); //Prior for spike and slab
      //Loglikelihood ratio
      alpha=alpha+LIKECOV( Y,  T1,X, thetavecprop,Beta)-LIKECOV( Y,  T1, X, thetavec,Beta);



      //Metropolis Hastings
      U=log(as_scalar(arma::randu(1)));



      if(U<alpha){
        thetavec=thetavecprop;
        SPIKEHOLD(k)=k;

      }




    }




    ///Sample Beta
    for(k=0;k<Beta.n_rows;k++){

      Betaprop=Beta;
      Betaprop(k)=-exp(as_scalar(arma::randn(1))*VarBeta(k)+log(-Beta(k)));


      //  alpha =   -.5*(gammadf+1)*log(1+pow((Betaprop(k)-0)/gammascale,2)/gammadf)+.5*(gammadf+1)*log(1+pow((Beta(k)-0)/gammascale,2)/gammadf);
      //Generalized t distribution
      alpha=log(-Betaprop(k))-log(-Beta(k));
      //Loglikelihood ratio
      alpha=alpha+LIKECOV( Y,  T1,X, thetavec,Betaprop)-LIKECOV( Y,  T1, X, thetavec,Beta);



      //Metropolis Hastings
      U=log(as_scalar(arma::randu(1)));



      if(U<alpha){
        ABeta(k)=ABeta(k)+1;

        Beta=Betaprop;

      }
      NBeta(k)=NBeta(k)+1;



    }



    //Storage

    StoreInx=m;

    thetavecprop=thetavec;


    for(j=0;j<thetavec.n_rows;j++){
      for(m1=0;m1<thetavec.n_cols;m1++){
        ThetaStore(StoreInx,(thetavec.n_cols*j+m1))=thetavecprop(j,m1);
      }
    }


    for(k=0;k<Beta.n_rows;k++){
      BetaStore(StoreInx,k)=Beta(k);
    }

    for(k=0;k<SPIKEHOLD.n_rows;k++){
      SPIKESTORE(StoreInx,k)=SPIKEHOLD(k);
    }












  }





  //End of MCMC














  List z1 = List::create(ThetaStore,BetaStore,SPIKESTORE,STORERAND,STORERAND1,STOREHOLD);
  //Entry 8 has NPO/PO
  //Entry 5,6 has




  return(z1);

}



