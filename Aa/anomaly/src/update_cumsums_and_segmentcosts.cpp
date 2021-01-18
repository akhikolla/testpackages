#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "Functions.h"


namespace anomalymv
{

void update_cumsums_and_segmentcosts(struct orderedobservationlist *list, int ii, int n, int p, int l, int minseglength)
{
	
	int jj = 0, start;
	struct orderedobservationlist* current = NULL;
	double factor, varianceestimate;
	double *means, *means_squared, *segmentcosts, *new_observation, *new_observationsquared;

	current = list[0].next;

	new_observation        = list[ii].observation;
	new_observationsquared = list[ii].observationsquared;

	start = ((ii-1) % (l+1))*p;

	while ((current->numberofobservation) < (ii - minseglength + 2) )
	{

		factor = (ii - current->numberofobservation + 1);

		means          = current->mean_of_xs;
		means_squared  = current->mean_of_xs_squared;
		segmentcosts   = current->segmentcosts;

		for (jj = 0; jj < p; jj++)
		{

			means[jj] = means[jj] + (new_observation[jj] - means[jj])/factor;
			means_squared[jj] = means_squared[jj] + (new_observationsquared[jj] - means_squared[jj])/factor;

			varianceestimate = means_squared[jj] - (means[jj])*(means[jj]);
		
			if (varianceestimate <= DBL_MIN)
			{
				varianceestimate = DBL_MIN;
			}

			segmentcosts[start + jj] = (log(varianceestimate) + 1 - means_squared[jj]) * factor;

			current->best_end_costs[jj] = find_lowest_end_cost(segmentcosts,jj,p,l);

		}

		current = current->next;

	}

	while ((current->numberofobservation) < (ii + 1) )
	{

		factor = (ii - current->numberofobservation + 1);

		means          = current->mean_of_xs;
		means_squared  = current->mean_of_xs_squared;

		for (jj = 0; jj < p; jj++)
		{

			means[jj] = means[jj] + (new_observation[jj] - means[jj])/factor;
			means_squared[jj] = means_squared[jj] + (new_observationsquared[jj] - means_squared[jj])/factor;

		}

		current = current->next;

	}

}

} // namespace anomalymv
