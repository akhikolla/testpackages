#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

#include "Functions_mean.h"

#include <vector>

#include <string>
#include "user_interupt.h"
#include "capa.exception.h"

using namespace anomalymv;

std::vector<int> MeanAnomalyMV(SEXP Rx, SEXP Rn, SEXP Rp, SEXP Rl, SEXP Rminlength, SEXP Rbetachange, SEXP Rbetaanomaly, SEXP Rmaxlength, SEXP Ronline)
// SEXP MeanAnomalyMV(SEXP Rx, SEXP Rn, SEXP Rp, SEXP Rl, SEXP Rminlength, SEXP Rbetachange, SEXP Rbetaanomaly, SEXP Rmaxlength, SEXP Ronline)
{
  
  /* 
  Rx    : Data
  Rn    : Length of data
  */

  /*
 	PROTECT(Rx) ; 
 	PROTECT(Rn) ;
 	PROTECT(Rp) ;
	PROTECT(Rl) ;
	PROTECT(Rminlength) ;
	PROTECT(Rbetachange) ;
	PROTECT(Rbetaanomaly) ;
	PROTECT(Rmaxlength) ;
	PROTECT(Ronline) ;
  */	
  	int n = 0, p = 0, l = 0, minlength = 0, ii = 0, error = 0, maxlength = 0, online = 0;
  	double betaanomaly = 0.0;
  	double* x = NULL;
  	double* betachange_DUMMY = NULL;
  	double* betachange = NULL;
	
	maxlength        = *(INTEGER(Rmaxlength));
 	minlength        = *(INTEGER(Rminlength));
	n                = *(INTEGER(Rn));
	p                = *(INTEGER(Rp));
	l  		 = *(INTEGER(Rl));
	online           = *(INTEGER(Ronline));
  	x          	 = REAL(Rx);
  	betachange_DUMMY = REAL(Rbetachange);
  	betaanomaly      = *REAL(Rbetaanomaly);

	std::vector<int> vout;	
	std::string reason;
	// int *out;
	struct orderedobservationlist_mean* mylist;


	try
	{
		betachange = new double[p];
	}
	catch(std::bad_alloc& e)
	{
		reason = "Not enough memory";
		error = 1;
		goto clearup;
	}

  	for (ii = 0; ii < p; ii++)
  	{
  		betachange[ii] = betachange_DUMMY[ii];
  	}

	try
	{
		populate_mean(&mylist, x, n, p, l); 
	}
	catch(std::bad_alloc& e)
	{
		reason = "Not enough memory";
		error = 1;
		goto clearup;
	}

	try
	{
		solveorderedobservationlist_mean(mylist, n, p, l, betachange, betaanomaly, minlength, maxlength);
	}
	catch(user_interupt& a)
	{
		reason = "user interrupt";
		error = 1;
		goto clearup;
	}	

	if (online)
	{

		vout.resize(n*(2 + 3*p));
		
		changepointreturn_mean_online(mylist, n, p, vout);
		
	} 
	else
	{

		int numberofchanges = 0, *changes = NULL, *components = NULL, *startlag = NULL, *endlag = NULL;
		
		try
		{
			changepointreturn_mean(mylist, n, p, &numberofchanges, &changes, &components, &startlag, &endlag);
		}
		catch(std::bad_alloc& e)
		{
			if(components){delete[] components;}
			if(startlag){delete[] startlag;}
			if(endlag){delete[] endlag;}
			if(changes){delete[] changes;}
			reason = "Not enough memory";
			error = 1;
			goto clearup;
		}		


		vout.resize(numberofchanges*(3 + 3*p));
	
		for (ii = 0; ii < 3*numberofchanges; ii++)
		{
			vout[ii] = changes[ii];
		}

		for (ii = 0; ii < numberofchanges*p; ii++)
		{
			vout[ii + 3*numberofchanges] = components[ii];
		}

		for (ii = 0; ii < numberofchanges*p; ii++)
		{
			vout[ii + numberofchanges*(3 + p)] = startlag[ii];
		}

		for (ii = 0; ii < numberofchanges*p; ii++)
		{
			vout[ii + numberofchanges*(3 + 2*p)] = endlag[ii];
		}

		if(components){delete[] components;}
		if(startlag){delete[] startlag;}
		if(endlag){delete[] endlag;}
		if(changes){delete[] changes;}

	}

clearup:	

	if(mylist)
	{
		for (ii = 0; ii < n + l + 2; ii++)
		{

			if(mylist[ii].observation){ delete[] mylist[ii].observation;}
			if(mylist[ii].mean_of_xs){ delete[] mylist[ii].mean_of_xs;}
			if(mylist[ii].segmentcosts){ delete[] mylist[ii].segmentcosts;}
			if(mylist[ii].best_end_costs){  delete[] mylist[ii].best_end_costs;}
			if(mylist[ii].affectedcomponents){  delete[] mylist[ii].affectedcomponents;}
			if(mylist[ii].startlag){ delete[] mylist[ii].startlag;}
			if(mylist[ii].endlag){ delete[] mylist[ii].endlag;}

		}

		delete[] mylist;
	}

	if(betachange){ delete[] betachange;}

	//  	UNPROTECT(9);

	if (error != 0)
	{
	  throw_capa_exception(reason);
	}
	
	return(vout) ;
	
}










