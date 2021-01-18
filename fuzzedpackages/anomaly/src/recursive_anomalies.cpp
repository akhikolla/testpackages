#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <vector>




std::vector<int> recursive_anomalies(SEXP RPrevious, SEXP ROptions, SEXP Rn)
// SEXP recursive_anomalies(SEXP RPrevious, SEXP ROptions, SEXP Rn)
{
	 
 	PROTECT(RPrevious) ; 
 	PROTECT(ROptions) ;
	PROTECT(Rn) ;
	
  	int n = 0, ii = 0, *Previous = NULL, *Options = NULL, *Previous_dummy = NULL, *Options_dummy = NULL;

	n              = *(INTEGER(Rn));
  	Previous_dummy = INTEGER(RPrevious);
  	Options_dummy  = INTEGER(ROptions);

	Previous = (int*) calloc(n, sizeof(int));
	Options  = (int*) calloc(n, sizeof(int));
	
	for (ii = 0; ii < n; ii++)
	{

		Previous[ii] = Previous_dummy[ii];
		Options[ii]  = Options_dummy[ii];

	}

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
	PROTECT(Rout = allocVector(INTSXP, 3*number_of_anomalies));


	
	int *out;
  	out  = INTEGER(Rout);

	out[0] = -1;
	out[1] = -1;
	out[2] = -1;

	current_position = n-1;

	int anomaly_counter = 1;

	while (current_position > -1)
	{

		if (Options[current_position] > 0)
		{

			out[3*anomaly_counter+0] = Options[current_position];
			out[3*anomaly_counter+1] = Previous[current_position] + 1;
			out[3*anomaly_counter+2] = current_position + 1;

			anomaly_counter++;
		}

		current_position = Previous[current_position] - 1;

	}

	std::vector<int> vout(3*number_of_anomalies);
	for(unsigned int cursor = 0; cursor < vout.size(); cursor++)
	  {
	    vout[cursor] = out[cursor];
	  }
	
	free(Previous);
	free(Options);

  	UNPROTECT(4);

  	// return(Rout); 
	return(vout);
}










