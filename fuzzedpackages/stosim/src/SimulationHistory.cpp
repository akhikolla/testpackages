  // SimulatioinHistory.cpp file
 /*
 * Author: David Silkworth
 *         (c) 2011-2018 OpenReliability.org
 */


#include <R.h>
#include <Rcpp.h>
#include "stosim.h"


 SEXP SimulationHistory(SEXP arg1, SEXP arg2, SEXP arg3,
      SEXP arg4,  SEXP arg5, SEXP arg6,
        SEXP arg7,  SEXP arg8, SEXP arg9,
          SEXP arg10,  SEXP arg11, SEXP arg12)

{     using namespace Rcpp;

//src <- '
		// set vector arguments from R into Rcpp vector classes for use in C++
Rcpp::NumericVector OpLine_Vec(arg1);
Rcpp::NumericVector Event_ID_Vec(arg2);
Rcpp::NumericVector FD_Vec(arg3);
Rcpp::NumericVector FP1_Vec(arg4);
Rcpp::NumericVector FP2_Vec(arg5);
Rcpp::NumericVector FP3_Vec(arg6);
Rcpp::NumericVector RD_Vec(arg7);
Rcpp::NumericVector RP1_Vec(arg8);
Rcpp::NumericVector RP2_Vec(arg9);
Rcpp::NumericVector RP3_Vec(arg10);
Rcpp::IntegerVector Seed_Vec(arg11);
Rcpp::NumericVector SimulationYears(arg12);


int ModelSize = Seed_Vec.size();
Rcpp::RNGScope Scope;



		//Determine the max amount of random events to be sampled
		//  per event element, considering all elements in the model.
		//  The Weibull FP1 should still be pretty valid here.
		//  I have removed consideration of lognormals as failure distributions
		//  The revised calculation will take effect on next package rebuild *(1.0+0.25/ModelSize) 
		//  This attempts to over-correct for risks of single element models
int NumRands=0;
 for (int x= 0; x < ModelSize; x++)  {
double num = (8760/(FP1_Vec[x]+FP3_Vec[x])+0.5)*(1.0+0.25/ModelSize)*SimulationYears[0];
if((int)num>NumRands) NumRands=(int)num;
}


int i, j, t;
Rcpp::NumericVector RandCol(NumRands);
		// get the random values for intervals and durations
Rcpp::NumericMatrix RandVals(NumRands, ModelSize*2);

Environment base("package:base");
Function SetSeed = base["set.seed"];
// rlnorm has been fixed in Rcpp as of release 0.95
//Environment stats("package:stats");
//Function rlnorm=stats["rlnorm"];


int mySeed ;
j=0;
for( t=0; t< ModelSize; t++)   {
mySeed =Seed_Vec[t];
SetSeed(Named("seed",mySeed) );

switch((int)FD_Vec[t])  {
case 1:   {
RandCol= rexp(NumRands,1/FP1_Vec[t]);
for(i=0; i<NumRands; i++) {
RandVals(i,j) = RandCol(i);
}
j++;
break; }
case 2:   {
RandCol= rnorm(NumRands,FP1_Vec[t],FP2_Vec[t]);
for(i=0; i<NumRands; i++) {
RandVals(i,j) = RandCol(i);
}
j++;
break; }
case 3:   {
RandCol=rweibull(NumRands,FP2_Vec[t],FP1_Vec[t]);
for(i=0; i<NumRands; i++) {
RandVals(i,j) = RandCol(i);
}
j++;
break; }

}

switch((int)RD_Vec[t])  {
case 2:   {
RandCol= rnorm(NumRands,RP1_Vec[t],RP2_Vec[t]);
for(i=0; i<NumRands; i++) {
RandVals(i,j) = RandCol(i);
}
j++;
break; }
case 3:   {
RandCol=rweibull(NumRands,RP2_Vec[t],RP1_Vec[t]);
for(i=0; i<NumRands; i++) {
RandVals(i,j) = RandCol(i)+RP3_Vec(t);
}
j++;
break; }
case 4:   {
RandCol=rlnorm(NumRands,RP1_Vec[t],RP2_Vec[t]);
for(i=0; i<NumRands; i++) {
RandVals(i,j) = RandCol(i)+RP3_Vec(t);
}
j++;
break; }

}
}		// NumRands Matrix has been built for reference in this scope

		// build an initial event que

		// I really want to build the individual vectors then bind them
		//  into a DataFrame only for export back to R

Rcpp::NumericVector  timeQ;
Rcpp::NumericVector durationQ;
Rcpp::NumericVector oplineQ;
Rcpp::NumericVector eidQ;
Rcpp::NumericVector eeQ;


Rcpp::IntegerVector rPtr(ModelSize);
int ee=0;
int qPtr=0;
int iPtr=0;

		// this initial event que will be built in a sorted fashion

for (ee=0; ee<ModelSize; ee++)    {
if(ee == 0)   {		// just place the first entry
timeQ.push_back(RandVals(rPtr(ee),ee*2));
durationQ.push_back(RandVals(rPtr(ee),ee*2+1));
oplineQ.push_back(OpLine_Vec(ee));
eidQ.push_back(Event_ID_Vec(ee));
eeQ.push_back(ee);
rPtr(ee)=rPtr(ee)+1;
}else{
		// this could be larger than any we've seen yet
if(RandVals(rPtr(ee),ee*2)>timeQ(timeQ.size()-1))  {
timeQ.push_back(RandVals(rPtr(ee),ee*2));
durationQ.push_back(RandVals(rPtr(ee),ee*2+1));
oplineQ.push_back(OpLine_Vec(ee));
eidQ.push_back(Event_ID_Vec(ee));
eeQ.push_back(ee);
rPtr(ee)=rPtr(ee)+1;
}else{

iPtr=0;
		// now I have to find the point to insert
for(int x=qPtr; x<timeQ.size(); x++)   {
if(timeQ(x)>RandVals(rPtr(ee),ee*2))   {
iPtr=x;
break;  }
}
timeQ.insert(iPtr,RandVals(rPtr(ee),ee*2));
durationQ.insert(iPtr,RandVals(rPtr(ee),ee*2+1));
oplineQ.insert(iPtr,OpLine_Vec(ee));
eidQ.insert(iPtr,Event_ID_Vec(ee));
eeQ.insert(iPtr,ee);
rPtr(ee)=rPtr(ee)+1;
}
}
}
		// Now we can enter the main engine loop
		// Always use _Ptr's and timeQ.size() to reflect on where we are or what we have

double SimulationLimit=8760*SimulationYears[0];
double ProposedTime=0;



do   {
ee=eeQ(qPtr);
	// add event downtime to rest of events in  the que
	// only if the Model is larger than single event element
if(ModelSize>1)  {
for(int y=qPtr+1; y<timeQ.size(); y++)  {
timeQ(y)=timeQ(y)+durationQ(qPtr);
} }

qPtr++;		// increment the qPtr (thereby establishing the last entry to SimHistory)
		//  establish a ProposedTime for the next event on the element that has just been secured in history
ProposedTime=timeQ(qPtr-1)+durationQ(qPtr-1)+RandVals(rPtr(ee),ee*2);
		// this could be larger than any we've seen yet
if(ProposedTime>timeQ(timeQ.size()-1))  {
timeQ.push_back(ProposedTime);
durationQ.push_back(RandVals(rPtr(ee),ee*2+1));
oplineQ.push_back(OpLine_Vec(ee));
eidQ.push_back(Event_ID_Vec(ee));
eeQ.push_back(ee);
rPtr(ee)=rPtr(ee)+1;
}else{

iPtr=0;
		// now I have to find the point to insert
for(int x=qPtr; x<timeQ.size(); x++)   {
if(timeQ(x)>ProposedTime)   {
iPtr=x;
break;  }
}
timeQ.insert(iPtr,ProposedTime);
durationQ.insert(iPtr,RandVals(rPtr(ee),ee*2+1));
oplineQ.insert(iPtr,OpLine_Vec(ee));
eidQ.insert(iPtr,Event_ID_Vec(ee));
eeQ.insert(iPtr,ee);
rPtr(ee)=rPtr(ee)+1;
}

		// We should be done now the while test will
		// be on timeQ(qPtr)   Seems that EventTime was never needed.
}
while(timeQ(qPtr)<SimulationLimit);

		// will use pop_back to erase undesired events at end
if(ModelSize>1)  {
for(i=0; i<ModelSize-1; i++)  {
t=timeQ.size()-1;
timeQ.erase(t);
durationQ.erase(t);
oplineQ.erase(t);
eidQ.erase(t);
eeQ.erase(t);
}  }

t=timeQ.size()-1;
timeQ(t)=SimulationLimit;
durationQ(t)= 0.0;
oplineQ(t)=-1;
eidQ(t)=-1;

if(timeQ(t-1)+durationQ(t-1) > SimulationLimit)  {
durationQ(t-1)=SimulationLimit-timeQ(t-1);
}


Rcpp::DataFrame SimHistory=
  Rcpp::DataFrame::create( Rcpp::Named("Time")=timeQ,
     Rcpp::Named("Duration")=durationQ,
     Rcpp::Named("OpLine")=oplineQ,
     Rcpp::Named("EventID")=eidQ);

 return(SimHistory);
 //'
}
