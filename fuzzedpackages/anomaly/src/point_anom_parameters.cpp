#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "Functions.h"

namespace anomalymv
{

void point_anom_parameters(struct orderedobservationlist *list, int ii, int p, double penaltyanomaly)
{

	int jj;
	double obs, extra;

	for (jj = 0; jj < p; jj++)
	{

		obs = list[ii].observationsquared[jj];

		if (obs <= DBL_MIN)
		{
			obs = DBL_MIN;
		}

		extra = penaltyanomaly + log(obs) + 1 - obs;

		if (extra < 0)
		{
			list[ii].affectedcomponents[jj] = 1;
		}

	}

}

} // namespace anomalymv
