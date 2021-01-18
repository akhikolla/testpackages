#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "Functions.h"
#include <vector>
namespace anomalymv
{

void changepointreturn_online(struct orderedobservationlist *mylist, int n, int p, std::vector<int> &out )
{

	struct orderedobservationlist current; 

	int ii = 1, jj = 0;

	for (ii = 1; ii < n+1; ii++)
	{

		current = mylist[ii];

		out[(ii-1)*(3*p+2) + 0] = current.option;
		out[(ii-1)*(3*p+2) + 1] = current.optimalcut->numberofobservation;

		for (jj = 0; jj<p; jj++)
		{

			out[(ii-1)*(3*p+2) + 2 + 0*p + jj] = current.affectedcomponents[jj];
			out[(ii-1)*(3*p+2) + 2 + 1*p + jj] = current.startlag[jj];
			out[(ii-1)*(3*p+2) + 2 + 2*p + jj] = current.endlag[jj];

		}

		
	}
}

} // namespace anomalymv
