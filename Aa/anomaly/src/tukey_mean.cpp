#include <math.h>
#include <stdlib.h>
#include <list>
#include <vector>
#include "tukey.h"
#include "tukey.h"
#include "Online_tukey.h"


// [[Rcpp::export]]
double tukey_mean(std::vector<double> x, double th){

	int n = x.size();

	double th_squared = th*th;
	
	Online_tukey tmp;

	for (int ii = 0; ii < n; ii++)
	{

		tmp.Add_observation(x[ii], x[ii]*x[ii], th, th_squared);		

	}	

	return(tmp.Find_mean());

};
