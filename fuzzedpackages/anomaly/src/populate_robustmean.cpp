#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "Functions_robustmean.h"

namespace anomalymv
{

void populate_robustmean(struct orderedobservationlist_robustmean **list, double* x , int n, int p, int l) 
{

	int ii = 0, jj = 0;

	*list = new orderedobservationlist_robustmean[n+l+2];

	struct orderedobservationlist_robustmean *mylist = *list;

	for (ii = 0; ii < n+l+2; ii ++){

		mylist[ii].numberofobservation = ii;

		mylist[ii].observation        = NULL;
		mylist[ii].observationsquared = NULL;
		mylist[ii].Tukey_Stuff        = NULL;

		mylist[ii].segmentcosts          = NULL;
		mylist[ii].best_end_costs        = NULL;

		mylist[ii].optimalcostofprevious =    0;
		mylist[ii].optimalcost           =    0;
		mylist[ii].costofstartingsegment =    0;
	
		mylist[ii].affectedcomponents = NULL;
		mylist[ii].startlag           = NULL;  
		mylist[ii].endlag             = NULL;

		mylist[ii].optimalcut         = NULL;
		mylist[ii].option             =   -1;

		mylist[ii].destruction = n + 100;
		mylist[ii].next        =    NULL;
		mylist[ii].previous    =    NULL;

	}


	mylist[0].next       = &(mylist[1]);
	mylist[n+l+1].previous = &(mylist[n+l]);


	for (ii = 1; ii < n+l+1; ii++)
	{

		mylist[ii].observation        = new double[p];
		mylist[ii].observationsquared = new double[p];
		mylist[ii].Tukey_Stuff        = new Online_tukey[p];

		mylist[ii].segmentcosts   = new double[p * (l+1)]; 
		mylist[ii].best_end_costs = new double[p];
	
		mylist[ii].affectedcomponents = new int[p];
		mylist[ii].startlag           = new int[p];
		mylist[ii].endlag             = new int[p];

		for (jj = 0; jj < p; jj ++)
		{

			mylist[ii].best_end_costs[jj]     = 100;

			mylist[ii].affectedcomponents[jj] = 0;
			mylist[ii].startlag[jj]           = 0;
			mylist[ii].endlag[jj]             = 0;

		}

		for (jj = 0; jj < p* (l+1); jj ++)
		{

			mylist[ii].segmentcosts[jj] = 100;

		}

		mylist[ii].next        = &(mylist[ii+1]);
		mylist[ii].previous    = &(mylist[ii-1]);

	}

	for (ii = 1; ii < n+1; ii++)
	{

		for (jj = 0; jj < p; jj ++)
		{

			mylist[ii].observation[jj]               = x[n*jj+ii-1];
			mylist[ii].observationsquared[jj]        = x[n*jj+ii-1] * x[n*jj+ii-1];

		}

	}


}

} // namespace anomalymv
