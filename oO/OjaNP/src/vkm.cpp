// VKM.cpp : Defines the entry point for the console application.
//



#define COMPILE_AS_DLL
//#define COMPILE_AS_EXE

#include "stdafx.h"
#include <math.h>
#include <iostream>
// #include <conio.h> //Nicht ANSI!
#include <time.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "DebugTags.h"

#include "Matrix2D.h"
#include "SubsetGenerator.h" 
#include "MersenneTwister.h"
#include "pipes.h"

using namespace std;

static int DIM = 15; //dimension
static int NUMBER  = 1000; //number of points given
static int INTERVALL = 20; //intervall for adapting sigma, =1 for USE_GRADIENT, =10..20 for USE_RANDOM 
static double SCALE_FACTOR = 0.5;
static int NUM_OF_ITERATIONS = 10000000;// maximal number of iterations allowed
static bool USE_ALL_SUBSETS = false; //= false: take random subsets
static int NUMBER_OF_SUBSETS_USED = 500; // Number of random subsets to use
static int SIGMA_LOG10_DECREASE=8; // Convergence is reached if log10(sigmaStart)-log10(sigma) <= SIGMA_LOG10_DECREASE
static bool STORE_SUBDETERMINANTS=true; // Store Sub-Dets for faster computation, be aware of size of CPU-Cache
static bool USE_GRADIENT=false;
static bool USE_RANDOM=true;
static bool SET_SIGMA_MANUALLY=false;

double normRand(double mean, double std)
{
  /*double v = 2;
  double u1 = 0;
  double u2 = 0;
  double u3 = 0;
  double normRand = 0;
  int count = 0;
  while (v >= 1)
  {
    u1 = rand() / (double)RAND_MAX;
    u2 = rand() / (double)RAND_MAX; 
    v = (2*u1 - 1)*(2*u1 - 1) + (2*u2 - 1)*(2*u2 - 1);
  }
  u3 = ((2*u1 - 1) * sqrt(-2*log(v)/v));
  normRand = std * u3 + mean;
  return normRand;*/
  double normRand = std*norm_rand()+mean;
  return normRand;
}

void genPoints(Matrix2D * points)
{
 double *x=new double[DIM];
 int i;

 for(i=0;i<NUMBER*1/3.0;i++)
 {
	 for(int j=0;j<DIM;j++)
       x[j]=normRand(0,1);

	 points->setColumn(i,x);
 }

 for(;i<NUMBER;i++)
 {
	 for(int j=0;j<DIM;j++)
       x[j]=normRand(0,0.1);

     points->setColumn(i,x);
 }

 delete[] x;
} 



int calc_ojaMedian(Matrix2D & points, 
	double & result_bestValue, double & result_sigma, Vector & result_bestMu, int &usedIter)
{
  #ifdef PRINT_TRACE
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei," \n Anfang der Methode \n");
  fclose(debugDatei); 
  #endif
  
  /** Mersenne Twister **/
  
  unsigned int seed= unif_rand()*32767;

  SubsetGenerator *g1,*g2,*g3,*g4,*g5;

  if (USE_ALL_SUBSETS)
  {
	  g1=new AllSubsets(NUMBER, DIM); 
	  g2=new AllSubsets(NUMBER, DIM); 
	  g3=new AllSubsets(NUMBER, DIM); 
	  g4=new AllSubsets(NUMBER, DIM); 
	  g5=new AllSubsets(NUMBER, DIM); 
  }
  else
  {
	  g1=new RandomSubsets(NUMBER, DIM, seed, NUMBER_OF_SUBSETS_USED);
	  g2=new RandomSubsets(NUMBER, DIM, seed, NUMBER_OF_SUBSETS_USED);
	  g3=new RandomSubsets(NUMBER, DIM, seed, NUMBER_OF_SUBSETS_USED);
	  g4=new RandomSubsets(NUMBER, DIM, seed, NUMBER_OF_SUBSETS_USED);
	  g5=new RandomSubsets(NUMBER, DIM, seed, NUMBER_OF_SUBSETS_USED);
  }

  #ifdef PRINT_TRACE
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei,"\n Subsets wurden initialisiert \n");
  fclose(debugDatei);
  #endif
  	

  ComputeObjectiveFunction objFunc(g5,&points,STORE_SUBDETERMINANTS);
  
  #ifdef PRINT_TRACE
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei,"\n objFunc wurde erstellt \n");
  fclose(debugDatei);
  #endif
  
  ComputeNabla compNabla(g5,&points,STORE_SUBDETERMINANTS);

  #ifdef PRINT_TRACE
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei,"\n  compNabla wurde erstellt \n");
  fclose(debugDatei);
  #endif
  
  /*** Calculate variance and initialize sigma ***/
  
  if (!SET_SIGMA_MANUALLY) result_sigma = points.getVariance();

  /*** Choose initial point for mu ***/
  
   Vector mu(DIM);         

  // try some random points for initial mu
  
  int best_r=0;
  Vector *v;
  v=points.getColumn(0);
  double best=objFunc.compute(v); //calculateObjectiveFunctionValue(g4, points, *v);
  delete v;
  
  
  #ifdef PRINT_TRACE
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei,"\n  Mu Initialisierungsschleife: NUMBER: %i, SEED: %i \n",NUMBER, seed);
  fclose(debugDatei);
  #endif
  
  for(int i=0;i<10;i++)
  {
	  #ifdef PRINT_TRACE
	  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  	  PRINT_LOCATION;
  	  fprintf(debugDatei,"\n    Mu Initialisierungsschleife: Durchlauf %i gestartet\n",i);
  	  fclose(debugDatei);
  	  #endif
  
	  int r=NUMBER+1;
	  while(r>=NUMBER)
 //         r=(int)rr.rand(NUMBER);    //cast um warning zu umgehen
            r=(int) unif_rand()*NUMBER;
	  v=points.getColumn(r);
	  double f=objFunc.compute(v); //calculateObjectiveFunctionValue(g4, points, *v);
	  delete v;
	  if (f<best)
	  {
		  best=f;
		  best_r=r;
	  }
	  
	  #ifdef PRINT_TRACE
	  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  	  PRINT_LOCATION;
  	  fprintf(debugDatei,"\n    Mu Initialisierungsschleife: Durchlauf %i beendet \n",i); 
  	  fclose(debugDatei);
  	  #endif 
  	  
  }
  

  for(int i=0;i<DIM;i++)
	  mu.setValue(i, points.getValue(i,best_r));
	  
	  
  #ifdef PRINT_TRACE
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei,"\n  mu wurde mit zuf�lligen werten initialisiert \n");
  fclose(debugDatei);
  #endif
  	  
//  mu.print();
  

  /** store initial value as current best known value **/

  /** calculate initial objective function value **/
  
  result_bestValue = objFunc.compute(&mu);
  
  result_bestMu.setValues(mu);

  #ifdef COMPILE_AS_EXE
  /*************************************/
  #ifdef PRINT_SUB_RESULTS
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei,"\n Start value = %f \n",result_bestValue);
  fclose(debugDatei); 
  #endif
  /*************************************/
  #endif

  int numOfGoodMutations = 0;
  int numOfIterations = 0;

  double sigmaStart=result_sigma;
 
  bool nablaOutOfDate=true;
  Vector* nablaF_org = new Vector(DIM);
  Vector* nablaF = new Vector(DIM);

  
  #ifdef PRINT_TRACE
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei,"\n Die Initialisierung vor der Hauptschleife wurde erfolgreich abgeschlossen \n");
  fclose(debugDatei);
  #endif
  
  #ifdef PRINT_TRACE 
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei,"\n Nun beginnt die Hauptschleife \n ");
  fclose(debugDatei);
  #endif
  for(int i = 0; i < NUM_OF_ITERATIONS && log10(sigmaStart)-log10(result_sigma)<=SIGMA_LOG10_DECREASE; i++)
  {
	#ifdef PRINT_TRACE 
  	debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  	PRINT_LOCATION;
  	fprintf(debugDatei,"\n Anfang des %i-ten Hauptschleifendurchlaufs",i);
  	fclose(debugDatei); 	
  	#endif
  	
  	#ifdef PRINT_SUB_RESULTS
  	debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  	PRINT_LOCATION;
  	fprintf(debugDatei,"\n Best values = %f \n",result_bestValue);
  	fclose(debugDatei);
  	#endif
	
	usedIter=i;
    numOfIterations++;

    /** calculate nablaF and adapt length **/

	if (USE_GRADIENT && !USE_RANDOM)
	{
		if (nablaOutOfDate)
		{
			nablaOutOfDate=false;

			*nablaF_org=compNabla.compute(&mu);
		}

		*nablaF=*nablaF_org;

		nablaF->setLength(fabs(result_sigma));
	}

	if (USE_RANDOM)
	{
			for(int i=0;i<DIM;i++)
				nablaF->setValue(i,normRand(0, result_sigma));

			double randomValue = normRand(0, result_sigma); 
			nablaF->setLength(fabs(randomValue));
	}

    /** add value to mu **/
    mu += *nablaF;


    /** evaluate objective function **/
    double f;
	
  
	f = objFunc.compute(&mu);
	
    if(f < result_bestValue)
    {
      result_bestValue = f;
      result_bestMu.setValues(mu);
      numOfGoodMutations++;

	  nablaOutOfDate=true;
    }
    else
    {
      mu.setValues(result_bestMu);
	}

    /*************************************/
    
    #ifdef PRINT_SUB_RESULTS
    debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  	PRINT_LOCATION;
  	fprintf(debugDatei,"\n Current value = %f \n",f);
    fclose(debugDatei);
    #endif

    //adapt sigma if necessary
    if(numOfIterations == INTERVALL)
    {
      
      if((double) numOfGoodMutations / (double) numOfIterations >= 0.2)
        result_sigma /= SCALE_FACTOR; // lots of good mutations --> increase
      else if((double) numOfGoodMutations / (double) numOfIterations < 0.2)
        result_sigma *= SCALE_FACTOR; // few good mutations --> decrease
      
	  numOfIterations = 0;
      numOfGoodMutations = 0;
    }
  } 
  
  #ifdef PRINT_TRACE
  debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  PRINT_LOCATION;
  fprintf(debugDatei,"\n  Die Hauptschleife wurde erfolgreich abgeschlossen \n");
  fclose(debugDatei);
  #endif 

  delete nablaF_org;
  delete nablaF;

  delete g1;
  delete g2;
  delete g3;
  delete g4;
  delete g5;

  return 0;
}

#ifdef COMPILE_AS_EXE

int main(int argc, char **argv) 
{

  Matrix2D points(DIM, NUMBER);
  genPoints(&points); 

  double bestValue;
  double sigma;
  Vector bestMu(DIM);
  int usedIter = 0;
  
  calc_ojaMedian(points,bestValue,sigma,bestMu,usedIter);						// Result


  //window stays open
  //cout << "Convergence reached after " << usedIter << " iterations \n";

  //cout << "Best value = ";
  printf("%20f ",bestValue);;
  //cout << "Best mu = ";
  bestMu.print();
  //cout << "Press any key to close program";
  getch();

	return 0;
}

#endif

#ifdef COMPILE_AS_DLL



SEXP answer;

SEXP mkDblAns(double x) 
{
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,1));
	REAL(ans)[0]=x;
	UNPROTECT(1);
	return ans;
}

SEXP mkIntAns(int x) 
{
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=x;
	UNPROTECT(1);
	return ans;
}



extern "C" {
	SEXP ojaEvo(SEXP data,SEXP initialSigma,SEXP sigmaAdaptation,SEXP adaptationFactor,SEXP iterations,SEXP useAllSubsets,SEXP numberOfSubsetsUsed,SEXP sigmaLog10Decrease,SEXP storeSubdeterminants)
	{	
		#ifdef PRINT_TRACE 
		debugDatei=fopen(DEBUG_FILE_NAME,"w"); // erstellt eine leere Debugdatei 	
		fclose(debugDatei);
	    #endif
		
		INTERVALL = (int)REAL(sigmaAdaptation)[0];
		SCALE_FACTOR = REAL(adaptationFactor)[0];
		NUM_OF_ITERATIONS = (int)REAL(iterations)[0];
		USE_ALL_SUBSETS = (bool)LOGICAL(useAllSubsets)[0];  
		NUMBER_OF_SUBSETS_USED = (int)REAL(numberOfSubsetsUsed)[0];
		SIGMA_LOG10_DECREASE= (int)REAL(sigmaLog10Decrease)[0];
		STORE_SUBDETERMINANTS=(bool)LOGICAL(storeSubdeterminants)[0];
		//USE_GRADIENT=(bool)LOGICAL(useGradient)[0];
		//USE_RANDOM=(bool)LOGICAL(useRandom)[0];
        USE_GRADIENT=(bool)false;
		USE_RANDOM=(bool)true;

		
    GetRNGstate();
	  	SEXP names,className,dataClass,column,solution;
	  	
	  	//PROTECT(column);    //warning R
        PROTECT(column=allocVector(REALSXP,DIM));	  
		//PROTECT(solution);  //warning R
        PROTECT(solution=allocVector(REALSXP,DIM));
		PROTECT(answer=allocVector(VECSXP,5));
		PROTECT(names=allocVector(STRSXP,5));
		PROTECT(className = allocVector(STRSXP,1));
		SET_STRING_ELT(className,0,mkChar("Oja Median"));
		SET_STRING_ELT(names,0,mkChar("best"));
		SET_STRING_ELT(names,1,mkChar("crit"));		
		SET_STRING_ELT(names,2,mkChar("Sigma"));	
		SET_STRING_ELT(names,3,mkChar("readData"));
		SET_STRING_ELT(names,4,mkChar("Sigma"));						
		classgets(answer,className);
		namesgets(answer,names);

// --------Ab R-2.4.0 wird EnsureString nicht mehr unterst�tzt (daf�r funktioniert jetzt isFrame)---------	  

	 	PROTECT(dataClass=getAttrib(data,R_ClassSymbol));

							
//	  	if ((strcmp(CHAR(EnsureString(dataClass)),"data.frame")==0)||(isMatrix(data)==1)) {
	  	if ((isFrame(data)==1)||(isMatrix(data)==1)) {	
			#ifdef PRINT_TRACE 
			debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  			PRINT_LOCATION;
  			fprintf(debugDatei,"\n 1. Vorbedingung erf�llt \n");	 
	    	fclose(debugDatei);
	    	#endif
	  		
			bool dataFrame=false;
//			if (strcmp(CHAR(EnsureString(dataClass)),"data.frame")==0) {
			if (isFrame(data)==1) {				
				dataFrame=true;
				DIM=length(data);
				NUMBER=length(VECTOR_ELT(data,0));
				
				#ifdef PRINT_TRACE 
				debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  				PRINT_LOCATION;
				fprintf(debugDatei,"\ndataClass ist ein Data Frame \n");
				fclose(debugDatei);
	    		#endif
			} else {
				dataClass=getAttrib(data,R_DimSymbol);
				DIM= INTEGER(dataClass)[1];
				NUMBER= INTEGER(dataClass)[0];
			}
	    	    
	    // POINTS EINLESEN 
	    
  	    Matrix2D points(DIM, NUMBER);
		
		if (dataFrame) {			
			for(int j=0;j<DIM;j++) {
				column=coerceVector(VECTOR_ELT(data,j),REALSXP);
		    	for(int i=0;i<NUMBER;i++) {			
					points.setValue(j,i,REAL(column)[i]); 
				}
		    }     		    
		} else {		    	
	    	for(int j=0;j<DIM;j++) {							
				for(int i=0;i<NUMBER;i++) {
					if (isReal(data)==1) points.setValue(j,i,REAL(data)[j*NUMBER + i]); else points.setValue(j,i,(double)INTEGER(data)[j*NUMBER + i]);	
				}
		   	}     		    
		}
		
	    #ifdef PRINT_TRACE
	    debugDatei=fopen(DEBUG_FILE_NAME,"a");	
		PRINT_LOCATION;
	    fprintf(debugDatei,"\n Punkte wurden eingelesen \n");	 
	    fclose(debugDatei);
	    #endif


				
		double bestValue;
  		double sigma=REAL(initialSigma)[0];
        SET_SIGMA_MANUALLY = (REAL(initialSigma)[0]!=0.0);
  		Vector bestMu(DIM);
  		int usedIter = 0;
		
	    #ifdef PRINT_TRACE
	    debugDatei=fopen(DEBUG_FILE_NAME,"a");	
		PRINT_LOCATION;
	    fprintf(debugDatei,"\nBerechnung startet jetzt: \n");	 
	    fclose(debugDatei);
		#endif
  
  		calc_ojaMedian(points,bestValue,sigma,bestMu,usedIter);						// Result

	    
	    
	    #ifdef PRINT_TRACE
	    debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  		PRINT_LOCATION;
  		fprintf(debugDatei,"\nBerechnung wurde erfolgreich beendet! \n");	 
	    fclose(debugDatei);
		#endif
  		
		SET_VECTOR_ELT(answer,1,mkDblAns(bestValue));		   		  
		SET_VECTOR_ELT(answer,2,mkDblAns(sigma));		   		  		  
		//SET_VECTOR_ELT(answer,3,readData);
		SET_VECTOR_ELT(answer,3,data);
		SET_VECTOR_ELT(answer,4,mkDblAns(sigma));	  		  
	    solution = allocVector(REALSXP,DIM);
 		for(int i=0;i<DIM;i++) {
 			REAL(solution)[i]=(double)bestMu.getValue(i);
 		}
 		SET_VECTOR_ELT(answer,0,solution);
	  	
	  	#ifdef PRINT_TRACE
	  	PRINT_LOCATION;
	  	debugDatei=fopen(DEBUG_FILE_NAME,"a");	
  		PRINT_LOCATION;
  		fprintf(debugDatei,"\nDie L�sung wurde erfolgreich zur�ck geschrieben \n");	 
	    fclose(debugDatei);
		#endif
	  		  	
	  	}
	    UNPROTECT(6);
		return (answer);
	  PutRNGstate();
	}
	
}

#endif
