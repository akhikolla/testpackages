#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "Functions_mean.h"

namespace anomalymv
{

void point_anom_parameters_mean(struct orderedobservationlist_mean *list, int ii, int p, double penaltyanomaly)
{

	int jj;
	double obs, extra;

	for (jj = 0; jj < p; jj++)
	{

		obs = list[ii].observation[jj];

		extra = penaltyanomaly - obs*obs;

		if (extra < 0)
		{
			list[ii].affectedcomponents[jj] = 1;
		}

	}

}

} // namespace anomalymv
