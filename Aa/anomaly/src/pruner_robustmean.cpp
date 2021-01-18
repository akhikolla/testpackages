#include "Functions_robustmean.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace anomalymv
{

void pruner_robustmean(struct orderedobservationlist_robustmean *list, int ii, int p, int l, int minseglength, int maxseglength, double totalpenalty)
{

	double threshold;
	threshold = totalpenalty + list[ii].optimalcost;

	struct orderedobservationlist_robustmean* current = NULL;
	current = list[0].next;

	int destroy = 1;

	if (maxseglength < ii - current->numberofobservation + 2) 
	{

		current->previous->next = current->next;
		current->next->previous = current->previous;
		current = current->next;

	}

	while (current->numberofobservation < ii - minseglength - l + 2)
	{


		if (current->costofstartingsegment > threshold)
		{

			if (current->destruction > ii + minseglength + l)
			{

				current->destruction = ii + minseglength + l;

			}

		}

		if (destroy == 1)
		{

			destroy = 0;

			if (current->destruction < ii + 1 )
			{

				delete[] current->Tukey_Stuff;
				current->Tukey_Stuff = NULL;
				current->previous->next = current->next;
				current->next->previous = current->previous;
				destroy = 1;

			}

		}


		current = current->next;

	}

	
}


} // namespace anomalymv 
