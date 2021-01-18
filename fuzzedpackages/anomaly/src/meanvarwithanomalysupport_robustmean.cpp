#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> 
#include <math.h> 
#include <stdlib.h>
#include <float.h>
#include <stdbool.h>

#include "Functions.h"

#include "user_interupt.h"
#include "check_user_interrupt.h"

namespace anomaly
{

void populateorderedobservationlist_robustmean(struct orderedobservationlist_robustmean **list, double* x , int n) 
{

	int ii = 0;

	*list = new orderedobservationlist_robustmean[n+2];
	
	(*list)[0].numberofobservation = 0;
	(*list)[0].observation = 0;
	(*list)[0].observationsquared = 0;
	(*list)[0].Tukey_Stuff = new Online_tukey;

	(*list)[0].optimalcostofprevious = 0;
	(*list)[0].segmentcost = 0;

	(*list)[0].optimalcost = 0;
	(*list)[0].optimalcut = NULL;
	(*list)[0].option = -99;

	(*list)[0].destruction = n+100;
	(*list)[0].next = (orderedobservationlist_robustmean*)&((*list)[1]);
	(*list)[0].previous = NULL;

	for (ii = 1; ii < n+1; ii++)
	{

		(*list)[ii].numberofobservation = ii;
		(*list)[ii].observation = x[ii-1];
		(*list)[ii].observationsquared = x[ii-1]*x[ii-1];
		(*list)[ii].Tukey_Stuff = new Online_tukey;

		(*list)[ii].optimalcostofprevious = 0;
		(*list)[ii].segmentcost = 0;

		(*list)[ii].optimalcost = 0;
		(*list)[ii].optimalcut = NULL;
		(*list)[ii].option = -99;

		(*list)[ii].destruction = n+100;
		(*list)[ii].next = &(*list)[ii+1];
		(*list)[ii].previous = &(*list)[ii-1];

	}

	(*list)[n+1].numberofobservation = n+1;
	(*list)[n+1].observation = 0;
	(*list)[n+1].observationsquared = 0;
	(*list)[n+1].Tukey_Stuff = new Online_tukey;

	(*list)[n+1].optimalcostofprevious = 0;
	(*list)[n+1].segmentcost = 0;

	(*list)[n+1].optimalcost = 0;
	(*list)[n+1].optimalcut = NULL;
	(*list)[n+1].option = -99;

	(*list)[n+1].destruction = n+100;
	(*list)[n+1].next = NULL;
	(*list)[n+1].previous = &(*list)[n];

}

void updatewithobservation_robustmean(int ii, struct orderedobservationlist_robustmean *list, double* penaltychange, const double threshold, const double threshold_squared)
{
	
	double x, x_squared, saving;
	int factor;

	x         = list[ii].observation;
	x_squared = list[ii].observationsquared;

     	struct orderedobservationlist_robustmean* current = NULL;
	current = list[0].next;

	while (current->numberofobservation < ii+1)
	{
		
		current->Tukey_Stuff->Add_observation(x,x_squared,threshold,threshold_squared);

		saving  = current->Tukey_Stuff->Find_minimum();
		factor  = (ii - current->numberofobservation + 1);

		current->segmentcost = current->optimalcostofprevious - saving + penaltychange[factor - 1];
		current = current->next;
	}

}

void findoptimaloption_robustmean(int ii, struct orderedobservationlist_robustmean *list, int minseglength, double penaltyoutlier)
{
	
	int option = 0;
	struct orderedobservationlist_robustmean *bestcut = NULL;
	double optimalscore = 0, scoreanomalous = 0, squareestimate = 0;
	
	optimalscore = list[ii].optimalcostofprevious;
	bestcut= &(list[ii-1]);
	option = 0;

	squareestimate  = list[ii].observation * list[ii].observation;

	scoreanomalous = list[ii].optimalcostofprevious - squareestimate + penaltyoutlier;
	
	if (scoreanomalous < optimalscore)
	{
		optimalscore = scoreanomalous;
		option = 1;
	}

	struct orderedobservationlist_robustmean* currentcheck = NULL;
	currentcheck = list[0].next;

	while (currentcheck->numberofobservation < ii - minseglength + 2)
	{
		if (currentcheck->segmentcost < optimalscore)
		{
			bestcut   = &(list[currentcheck->numberofobservation-1]);
			option = 2;
			optimalscore = currentcheck->segmentcost;
		}

		currentcheck = currentcheck->next;
	}	
	
	list[ii].optimalcut              = bestcut;
	list[ii].optimalcost             = optimalscore;
	list[ii].option                  = option;
	list[ii+1].optimalcostofprevious = optimalscore;
}

void pruner_robustmean(struct orderedobservationlist_robustmean *list, int ii, double penaltychange_max, int minseglength, int maxseglength)
{

	double threshold;
	threshold = penaltychange_max + list[ii].optimalcost;

     	struct orderedobservationlist_robustmean* current = NULL;
	current = list[0].next;

	if (maxseglength < ii - current->numberofobservation + 2) 
	{

		current->previous->next = current->next;
		current->next->previous = current->previous;
		current = current->next;

	}

	while (current->numberofobservation < ii - minseglength + 2)
	{

		if (current->segmentcost > threshold)
		{

			if (current->destruction > ii + minseglength){current->destruction = ii + minseglength;}

		}

		if (current->destruction < ii + 1 )
		{
			delete current->Tukey_Stuff;
			current->Tukey_Stuff = NULL;
			current->previous->next = current->next;
			current->next->previous = current->previous;
		}

		current = current->next;

	}

	
}


void solveorderedobservationlist_robustmean(struct orderedobservationlist_robustmean *list, int n, double *penaltychange, double penaltyoutlier, int minseglength, int maxseglength)
{

	int ii = 1;

	double penaltychange_max = 0.0, sqrt_penaltyoutlier = sqrt(penaltyoutlier);
	
	for (ii = 0; ii < maxseglength; ii++)
	{
		if (penaltychange_max < penaltychange[ii]){ penaltychange_max = penaltychange[ii];}
	}

	for (ii = 1; ii < n+1; ii++)
	{
	  
		updatewithobservation_robustmean(ii,list,penaltychange,sqrt_penaltyoutlier,penaltyoutlier);
		findoptimaloption_robustmean(ii,list,minseglength,penaltyoutlier);
		pruner_robustmean(list,ii,penaltychange_max,minseglength,maxseglength);
		
		if (ii % 128 == 0)
		{
		  	if(check_user_interrupt())
		  	{
		     		user_interupt e;
		     		throw(e);
		  	}
		}
		
	}

}


void changepointreturn_robustmean(struct orderedobservationlist_robustmean *list, int n, int* numberofchanges, int** changepoints)
{

	*numberofchanges = 1;
	int  ii = 1;
	struct orderedobservationlist_robustmean* current;

	current = list[n+1].previous;
	
	while (current->numberofobservation > 0)
	{	
		if (current->option > 0){*numberofchanges = *numberofchanges + 1;}
		current = current->optimalcut;
	}

	
        *changepoints = new int[3*(*numberofchanges)];

	(*changepoints)[0] = -1; 
	(*changepoints)[1] = -1;
	(*changepoints)[2] = -1;
 
	current = list[n+1].previous;
	
	ii = 1;
	
	while (current->numberofobservation > 0)
	{	

		if (current->option > 0)
		{
			(*changepoints)[3*ii  ] = current->numberofobservation;
			(*changepoints)[3*ii+1] = current->optimalcut->numberofobservation + 1;
			(*changepoints)[3*ii+2] = current->option;
			ii++;
		}

		current = current->optimalcut;

	}
	

}

void changepointreturn_online_robustmean(struct orderedobservationlist_robustmean *list, int n, int** changepoints)
{
	
  	*changepoints = new int[2*n];
	
	int ii = 0;

	for (ii = 1; ii < n+1; ii++)
	{
		(*changepoints)[2*ii-2] = list[ii].option;
		(*changepoints)[2*ii-1] = list[ii].optimalcut->numberofobservation;
	}
	

}


} // namespace anomaly











