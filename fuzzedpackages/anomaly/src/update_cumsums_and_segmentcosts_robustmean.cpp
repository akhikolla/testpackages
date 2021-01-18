#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "Functions_robustmean.h"

namespace anomalymv
{

void update_cumsums_and_segmentcosts_robustmean(struct orderedobservationlist_robustmean *list, int ii, int n, int p, int l, int minseglength, const double threshold, const double threshold_squared)
{
	
	int jj = 0, start;
	struct orderedobservationlist_robustmean* current = NULL;
	double *new_observation_squared, *segmentcosts, *new_observation;

	current = list[0].next;

	new_observation         = list[ii].observation;
	new_observation_squared = list[ii].observationsquared;

	start = ((ii-1) % (l+1))*p;

	while ((current->numberofobservation) < (ii - minseglength + 2) )
	{

		segmentcosts   = current->segmentcosts;

		for (jj = 0; jj < p; jj++)
		{

			current->Tukey_Stuff[jj].Add_observation(new_observation[jj],new_observation_squared[jj],threshold,threshold_squared);

			segmentcosts[start + jj] = -current->Tukey_Stuff[jj].Find_minimum();

			current->best_end_costs[jj] = find_lowest_end_cost(segmentcosts,jj,p,l);

		}

		current = current->next;

	}

	while ((current->numberofobservation) < (ii + 1) )
	{

		for (jj = 0; jj < p; jj++)
		{

			current->Tukey_Stuff[jj].Add_observation(new_observation[jj],new_observation_squared[jj],threshold,threshold_squared);

		}

		current = current->next;

	}

}

} // namespace anomalymv
