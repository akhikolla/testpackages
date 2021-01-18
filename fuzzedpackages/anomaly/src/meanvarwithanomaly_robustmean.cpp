
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

#include "Functions.h"

#include <vector>

#include <string>
#include "capa.exception.h"
#include "user_interupt.h"

using namespace anomaly;

std::vector<int> RobustMeanAnomaly(SEXP Rx, SEXP Rn, SEXP Rminlength, SEXP Rmaxlength, SEXP Rbetachange, SEXP Rbetaanomaly, SEXP Ronline)
// SEXP MeanAnomaly(SEXP Rx, SEXP Rn, SEXP Rminlength, SEXP Rmaxlength, SEXP Rbetachange, SEXP Rbetaanomaly, SEXP Ronline)
{

  /*
 	PROTECT(Rx) ; 
 	PROTECT(Rn) ;
	PROTECT(Rminlength) ;
	PROTECT(Rmaxlength) ;
	PROTECT(Rbetachange) ;
	PROTECT(Rbetaanomaly) ;
	PROTECT(Ronline) ;
  */	
  	int n = 0, minlength = 0, maxlength = 0, error = 0, online = -1, ii = 0;
  	double betaanomaly = 0.0;
  	double* x = NULL, *betachange = NULL, *betavector = NULL ;
	std::string reason;

  
 	minlength        = *(INTEGER(Rminlength));
	maxlength        = *(INTEGER(Rmaxlength));
	n                = *(INTEGER(Rn));
  	x          	 =   REAL(Rx);
  	betachange       =   REAL(Rbetachange);
  	betaanomaly      = *REAL(Rbetaanomaly);
	online           = *INTEGER(Ronline);

	struct orderedobservationlist_robustmean* mylist;

	int numberofchanges = 0, *changes = NULL;

	std::vector<int> Rout;


	try
	{
		betavector = new double[maxlength];
	}
	catch(std::bad_alloc& e)
	{
		reason = "Not enough memory";
		error = 1;
		goto clearup;
	}
	
	for (ii = 0; ii < minlength-1; ii++){betavector[ii] = 0;}
	for (ii = minlength-1; ii < maxlength; ii++){betavector[ii] = betachange[ii+1-minlength];}

	try
	{
		populateorderedobservationlist_robustmean(&mylist, x, n);
	}
	catch(std::bad_alloc& e)
	{
		reason = "Not enough memory";
		error = 1;
		goto clearup;
	}

	try
	{	
		solveorderedobservationlist_robustmean(mylist, n, betavector, betaanomaly, minlength, maxlength);
	}
	catch(user_interupt& a)
	{
		reason = "user interrupt";
		error = 1;
		goto clearup;
	}
	catch(std::bad_alloc& e)
	{
		reason = "Not enough memory";
		error = 1;
		goto clearup;
	}

	if (online == 0)
	{

		try
		{
			changepointreturn_robustmean(mylist, n, &numberofchanges, &changes);
		}
		catch(std::bad_alloc& e)
		{
			reason = "Not enough memory";
			error = 1;
			goto clearup;
		}
	
		Rout.resize(3*numberofchanges);
  		
		for (ii = 0; ii < 3*numberofchanges; ii++)
		{
			Rout[ii] = changes[ii];
		}

	}
	else
	{

		try
		{
		changepointreturn_online_robustmean(mylist, n, &changes);
		}
		catch(std::bad_alloc& e)
		{
			reason = "Not enough memory";
			error = 1;
			goto clearup;
		}

		Rout.resize(2*n);

		for (ii = 0; ii < 2*n; ii++)
		{
			Rout[ii] = changes[ii];
		}

	}

clearup:

	for (ii = 0; ii < n+2; ii++)
	{	
		if(mylist[ii].Tukey_Stuff){delete mylist[ii].Tukey_Stuff;}
	}

	if(changes){delete[] changes;}
	if(betavector){delete[] betavector;}
	if(mylist){delete[] mylist;}

	// UNPROTECT(7);

	if (error != 0)
	{
	  throw_capa_exception(reason);
	}
	
	return(Rout) ;

}










