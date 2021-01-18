#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "Functions_mean.h"

namespace anomalymv
{

void update_cumsums_and_segmentcosts_mean(struct orderedobservationlist_mean *list, int ii, int n, int p, int l, int minseglength)
{
	
	int jj = 0, start;
	struct orderedobservationlist_mean* current = NULL;
	double factor;
	double *means, *segmentcosts, *new_observation;

	current = list[0].next;

	new_observation        = list[ii].observation;

	start = ((ii-1) % (l+1))*p;

	while ((current->numberofobservation) < (ii - minseglength + 2) )
	{

		factor = (ii - current->numberofobservation + 1);

		means          = current->mean_of_xs;
		segmentcosts   = current->segmentcosts;

		for (jj = 0; jj < p; jj++)
		{

			means[jj] = means[jj] + (new_observation[jj] - means[jj])/factor;

			segmentcosts[start + jj] = -(means[jj] * means[jj] ) * factor;

			current->best_end_costs[jj] = find_lowest_end_cost(segmentcosts,jj,p,l);

		}

		current = current->next;

	}

	while ((current->numberofobservation) < (ii + 1) )
	{

		factor = (ii - current->numberofobservation + 1);

		means          = current->mean_of_xs;

		for (jj = 0; jj < p; jj++)
		{

			means[jj] = means[jj] + (new_observation[jj] - means[jj])/factor;

		}

		current = current->next;

	}

}

} // namespace anomalymv
