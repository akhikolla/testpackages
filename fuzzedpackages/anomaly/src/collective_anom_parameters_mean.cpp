#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "Functions_mean.h"

namespace anomalymv
{

void collective_anom_parameters_mean(struct orderedobservationlist_mean *list, int ii, int p, int l, int minseglength, double *penaltycomponent, struct position_saving *savingvector)
{

	int jj, kk, bestobservation, numaffected, relevant_component, examined, tmp; 
	struct orderedobservationlist_mean *current = NULL, *examine = NULL, *tobefilled = NULL;
	double mincost = 100, currentcost = 0;
	double *costs;

	bestobservation = list[ii].optimalcut->numberofobservation + 1;
	current    = &(list[bestobservation]);
	tobefilled = &(list[ii]);

	examine = current;
	costs   = examine->best_end_costs;

	for (jj = 0; jj < p; jj++)
	{

	 	savingvector[jj].saving   = costs[jj];
	 	savingvector[jj].position = jj;
	 	tobefilled->startlag[jj] = 0;

	}

	for (kk = 0; kk < l; kk++)
	{

		examine = examine->next;
		costs   = examine->best_end_costs;

		for (jj = 0; jj < p; jj++)
		{

			if(savingvector[jj].saving > costs[jj])
			{

				savingvector[jj].saving   = costs[jj];
				tobefilled->startlag[jj]  = kk+1;

			}

		}

	}

	
		
	qsort(savingvector, p, sizeof(struct position_saving), cmpfunc_sorting);

	numaffected = 1;

	for (jj = 0; jj < p; jj++)
	{

		currentcost = currentcost + savingvector[jj].saving + penaltycomponent[jj];

		if (currentcost < mincost)
		{

			mincost = currentcost;
			numaffected = jj + 1;

		}

	}

	tmp = (ii -1)% (l+1);

	for (jj = 0; jj < numaffected; jj++)
	{

		relevant_component = savingvector[jj].position;
		examined           = relevant_component;
		
		tobefilled->affectedcomponents[relevant_component] = 1;

		

		examine = &(list[bestobservation + tobefilled->startlag[relevant_component]]);

		tobefilled->endlag[relevant_component] = -1;

		mincost = 100;

		for (kk = 0; kk < tmp + 1; kk++)
		{

			if (mincost > examine->segmentcosts[examined])
			{

				mincost = examine->segmentcosts[examined];
				tobefilled->endlag[relevant_component] = tmp - kk;		

			}

			examined = examined + p  ;

		}

		for (kk = tmp + 1; kk < (l + 1); kk++)
		{

			if (mincost > examine->segmentcosts[examined])
			{

				mincost = examine->segmentcosts[examined];
				tobefilled->endlag[relevant_component] = tmp - kk + l + 1;	

			}

			examined = examined + p  ;

		}

	}

	
}


} // namespace anomalymv
