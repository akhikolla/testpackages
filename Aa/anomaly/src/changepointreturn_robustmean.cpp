#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "Functions_robustmean.h"

namespace anomalymv
{

void changepointreturn_robustmean(struct orderedobservationlist_robustmean *list, int n, int p, int* numberofchanges, int** changepoints, int** components, int** startlag, int** endlag)
{

	*numberofchanges = 1;
	int  ii = 1, jj = 0;
	struct orderedobservationlist_robustmean* current;


	current = list[n+1].previous;


		
	while (current->numberofobservation > 0)
	{	
		if (current->option > 0){*numberofchanges = *numberofchanges + 1;}
		current = current->optimalcut;
	}

	
	*changepoints = new int[(*numberofchanges)*3];
	*components   = new int[(*numberofchanges)*p];
	*startlag     = new int[(*numberofchanges)*p];
	*endlag       = new int[(*numberofchanges)*p];


	(*changepoints)[0] = -1; 
	(*changepoints)[1] = -1;
	(*changepoints)[2] = -1;

	for (jj = 0; jj < p; jj++)
	{

		(*components)[jj] = -1;
		(*startlag)[jj]   = -1;
		(*endlag)[jj]     = -1;

	}

	
 
	current = list[n+1].previous;
	
	ii = 1;
	
	while (current->numberofobservation > 0)
	{	

		if (current->option > 0)
		{

			(*changepoints)[3*ii]   = current->numberofobservation;
			(*changepoints)[3*ii+1] = current->optimalcut->numberofobservation + 1;		
			(*changepoints)[3*ii+2] = current->option;	

			for (jj = 0; jj < p; jj++)
			{
				
				(*components)[ii*p+jj] = current->affectedcomponents[jj];
				(*startlag)[ii*p+jj]   = current->startlag[jj];
				(*endlag)[ii*p+jj]     = current->endlag[jj]; 

			}

			ii++;

		}

		current = current->optimalcut;

	}


}


} // namespace anomalymv











