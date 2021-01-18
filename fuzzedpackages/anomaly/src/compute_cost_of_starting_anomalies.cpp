#include "Functions.h"
#include <math.h>
#include <stdlib.h>

namespace anomalymv
{


void compute_cost_of_starting_anomalies(struct orderedobservationlist *list , int ii, int n, int p, int l, int minseglength, double *penaltycomponent, double *componentcost)
{

	int jj = 0, kk = 0;
	struct orderedobservationlist *current = NULL, *examine = NULL;
	double *costs, min, total;

	current = list[0].next;

	while ((current->numberofobservation) < (ii - minseglength + 2) )
	{

		examine = current;
		costs   = examine->best_end_costs;

		for (jj = 0; jj < p; jj++)
		{

	 		componentcost[jj] = costs[jj];

		}

		for (kk = 0; kk < l; kk++)
		{

			examine = examine->next;
			costs   = examine->best_end_costs;

			for (jj = 0; jj < p; jj++)
			{

				if(componentcost[jj] > costs[jj])
				{

					componentcost[jj] = costs[jj];

				}

			}

		}
		
		qsort(componentcost, p, sizeof(double), cmpfunc_nosorting);

		min = 100;
		total = 0;

		for (jj = 0; jj < p; jj++)
		{

			total = total + componentcost[jj] + penaltycomponent[jj];

			if (total < min)
			{
				min = total;
			}


		}

		current->costofstartingsegment = min + current->optimalcostofprevious;

		current = current->next;

	}

}

} // namespace anomalymv
