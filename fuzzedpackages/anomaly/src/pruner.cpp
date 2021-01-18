#include "Functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace anomalymv
{

void pruner(struct orderedobservationlist *list, int ii, int p, int l, int minseglength, int maxseglength, double totalpenalty)
{

	double threshold;
	threshold = totalpenalty + list[ii].optimalcost;

	struct orderedobservationlist* current = NULL;
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

				current->previous->next = current->next;
				current->next->previous = current->previous;
				destroy = 1;

			}

		}


		current = current->next;

	}

	
}

} // namespace anomalymv
