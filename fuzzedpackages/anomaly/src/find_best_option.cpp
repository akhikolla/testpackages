#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "Functions.h"


namespace anomalymv
{

void find_best_option(struct orderedobservationlist *list, int ii, int n, int p, int l, int minseglength, double *penaltycomponent, double penaltyanomaly, struct position_saving *savingvector)
{

	int jj = 0, option = 0, bestcut = 0;
	struct orderedobservationlist *current = NULL;
	double cost_pt_anom, cost_coll_anom, mincost, extra, obs;

	mincost = list[ii].optimalcostofprevious;

	current = list[0].next;

	while ((current->numberofobservation) < (ii - minseglength + 2) )
	{

		cost_coll_anom = current->costofstartingsegment;

		if (cost_coll_anom < mincost)
		{
			mincost = cost_coll_anom;
			option  = 2;
			bestcut = current->numberofobservation - 1;
		}

		current = current->next;

	}

	cost_pt_anom = list[ii].optimalcostofprevious;

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
			cost_pt_anom = cost_pt_anom + extra;
		}

	}

	if (cost_pt_anom < mincost)
	{
		mincost = cost_pt_anom;
		option  = 1;
	}


	
	list[ii].option                  = option;
	list[ii].optimalcost             = mincost;
	list[ii+1].optimalcostofprevious = mincost;


	list[ii].optimalcut = &(list[ii-1]);
	
	if (option == 2)
	{

		list[ii].optimalcut = &(list[bestcut]);
		collective_anom_parameters(list,ii,p,l,minseglength,penaltycomponent,savingvector);
	}
	

	if (option == 1)
	{

		point_anom_parameters(list,ii,p,penaltyanomaly);
		
	}


}



} // namespace anomalymv
