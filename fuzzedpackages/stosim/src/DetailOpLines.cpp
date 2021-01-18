// DetailOpLines.cpp file
/*
* Author: David Silkworth
* (c) 2011-2018 OpenReliability.org
*/
#include <R.h>
#include <Rcpp.h>
#include "stosim.h"
SEXP DetailOpLinesCPP(SEXP arg1, SEXP arg2, SEXP arg3,
SEXP arg4)
{ using namespace Rcpp;
Rcpp::NumericVector TimesVec(arg1);
Rcpp::NumericVector DurationsVec(arg2);
Rcpp::NumericVector LengthsVec(arg3);
Rcpp::CharacterVector NamesVec(arg4);
int num_opl=NamesVec.size();
int s = TimesVec.size()*2;
Rcpp::NumericVector Times2(s+1);
Rcpp::IntegerVector ud(s+1);
Rcpp::NumericVector newTimes(s);
Rcpp::NumericVector newDurations(s);
Rcpp::IntegerMatrix opd(s,num_opl);
Rcpp::IntegerVector oplPtrs(num_opl);
Rcpp::NumericVector TimeQ(num_opl);
int u=1;
Times2[0]=0.0;
ud[0]=u;
int x=0;
do {
u--;
Times2[x*2+1]=TimesVec[x];
ud[x*2+1]=u;
u++;
Times2[x*2+2]=TimesVec[x]+DurationsVec[x];
ud[x*2+2]=u;
x++;
}
while(x*2+2<ud.size());
Rcpp::NumericVector SimLimit(1);
int row=0;
int col=0;
int cur_col=0;
double LastEventTime=0.0;
SimLimit[0]=TimesVec[LengthsVec[0]-1];
int TempCount=0;
oplPtrs[col]=1;
TimeQ[col]=Times2[oplPtrs[col]];
opd(row,col)=1;
for(col=1; col<num_opl; col++) {
oplPtrs[col]=oplPtrs[col-1]+LengthsVec[col-1]*2;
TimeQ[col]=Times2[oplPtrs[col]];
opd(row,col)=1;
}
newTimes[row]=0.0;
do {
cur_col = std::min_element (TimeQ.begin(), TimeQ.end())-TimeQ.begin();
if(TimeQ[cur_col]>LastEventTime) {
row++;
for(col=0; col<num_opl; col++) {
opd(row,col)=opd(row-1,col);
}
newTimes[row]=Times2[oplPtrs[cur_col]];
LastEventTime=newTimes[row];
newDurations[row-1]=Times2[oplPtrs[cur_col]]-newTimes[row-1];
}
opd(row,cur_col)=ud[oplPtrs[cur_col]];
oplPtrs[cur_col]++;
TimeQ[cur_col]=Times2[oplPtrs[cur_col]];
TempCount++;
}
while ( LastEventTime<SimLimit[0]) ;
for(int e=s-1; e>row-1; e--) {newTimes.erase(e);}
for(int e=s-1; e>row-1; e--) {newDurations.erase (e);}
SubMatrix<INTSXP> yy = opd( Range(0,row-1), Range(0,num_opl-1) ) ;
opd = yy;
Rcpp::DataFrame DF=
Rcpp::DataFrame::create(Rcpp::Named("Time")=newTimes,
Named("Duration")=newDurations,opd);
return DF;
}
