#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <vector>


std::vector<int> recursive_mvanomalies(SEXP RPrevious, SEXP ROptions, SEXP Rcomponents, SEXP Rstartlags, SEXP Rendlags, SEXP Rn, SEXP Rp)
// SEXP recursive_mvanomalies(SEXP RPrevious, SEXP ROptions, SEXP Rcomponents, SEXP Rstartlags, SEXP Rendlags, SEXP Rn, SEXP Rp)
{
	 
 	PROTECT(RPrevious) ; 
 	PROTECT(ROptions) ;
 	PROTECT(Rcomponents) ; 
 	PROTECT(Rstartlags) ;
 	PROTECT(Rendlags) ; 
	PROTECT(Rn) ;
	PROTECT(Rp) ;
	
  	int n = 0, ii = 0, p = 0, *Previous = NULL, *Components = NULL, *Startlags = NULL, *Endlags = NULL, *Options = NULL;

	n              = *(INTEGER(Rn));
	p              = *(INTEGER(Rp));
  	Previous       = INTEGER(RPrevious);
  	Options        = INTEGER(ROptions);
	Startlags      = INTEGER(Rstartlags);
	Endlags        = INTEGER(Rendlags);
	Components     = INTEGER(Rcomponents);

	int number_of_anomalies = 1, current_position = n-1;

	while (current_position > -1)
	{

		if (Options[current_position] > 0)
		{
			number_of_anomalies++;
		}

		current_position = Previous[current_position] - 1;


	}

	SEXP Rout;

	PROTECT(Rout = allocVector(INTSXP, (3 + 3*p)*number_of_anomalies));
	
	int *out;
  	out  = INTEGER(Rout);

	for (ii = 0; ii < 3; ii++)
	{
		out[ii] = -1;
	}

	for (ii = 0; ii < p; ii++)
	{
		out[(3 + 0*p)*number_of_anomalies + ii] = -1;
		out[(3 + 1*p)*number_of_anomalies + ii] = -1;
		out[(3 + 2*p)*number_of_anomalies + ii] = -1;
	}

	current_position = n-1;

	int anomaly_counter = 1;

	while (current_position > -1)
	{

		if (Options[current_position] > 0)
		{

			out[3*anomaly_counter+0] = current_position + 1;
			out[3*anomaly_counter+1] = Previous[current_position] + 1;
			out[3*anomaly_counter+2] = Options[current_position];
			
			for (ii = 0; ii < p; ii++)
			{

				out[(3 + 0*p)*number_of_anomalies + anomaly_counter*p + ii] = Components[current_position*p + ii];
				out[(3 + 1*p)*number_of_anomalies + anomaly_counter*p + ii] = Startlags[current_position*p + ii];
				out[(3 + 2*p)*number_of_anomalies + anomaly_counter*p + ii] = Endlags[current_position*p + ii];

			}

			anomaly_counter++;
		}

		current_position = Previous[current_position] - 1;

	}

	std::vector<int> vout((3 + 3*p)*number_of_anomalies);
	for(unsigned int cursor = 0; cursor < vout.size(); cursor++)
	  {
	    vout[cursor] = out[cursor];
	  }


	
  	UNPROTECT(8);

  	// return(Rout);
	return(vout);

}










