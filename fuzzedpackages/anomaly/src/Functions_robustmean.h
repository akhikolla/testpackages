#include <stdbool.h>
#include "Online_tukey.h"

#include "check_user_interrupt.h"
#include "user_interupt.h"
#include <vector>

namespace anomalymv
{

typedef struct position_saving
{

	double saving;
	int position;

} position_saving;

typedef struct orderedobservationlist_robustmean 
{

	int numberofobservation;
	
	double* observation;
	double* observationsquared;
	Online_tukey* Tukey_Stuff;

	double* segmentcosts;
	double* best_end_costs;

	double optimalcostofprevious;
	double costofstartingsegment;
	double optimalcost;
	
	int* affectedcomponents;
	int* startlag;
	int* endlag;

	struct orderedobservationlist_robustmean* optimalcut;
	int option;

	int    destruction;
  	struct orderedobservationlist_robustmean* next;
  	struct orderedobservationlist_robustmean* previous;

} orderedobservationlist_robustmean;


void populate_robustmean(struct orderedobservationlist_robustmean **list, double* x , int n, int p, int l );

int cmpfunc_nosorting (const void * a, const void * b);

int cmpfunc_sorting (const void * a, const void * b);

void solveorderedobservationlist_robustmean(struct orderedobservationlist_robustmean *list, int n, int p, int l, double* penaltycomponent, double penaltyanomaly, int minseglength, int maxseglength);

void compute_cost_of_starting_anomalies_robustmean(struct orderedobservationlist_robustmean *list , int ii, int n, int p, int l, int minseglength, double *penaltycomponent, double *componentcost);

void update_cumsums_and_segmentcosts_robustmean(struct orderedobservationlist_robustmean *list, int ii, int n, int p, int l, int minseglength, const double threshold, const double threshold_squared);

double find_lowest_end_cost(double* segmentcosts, int jj, int p, int l);

void find_best_option_robustmean(struct orderedobservationlist_robustmean *list, int ii, int n, int p, int l, int minseglength, double *penaltycomponent, double penaltyanomaly, struct position_saving *savingvector);

void collective_anom_parameters_robustmean(struct orderedobservationlist_robustmean *list, int ii, int p, int l, int minseglength, double *penaltycomponent, struct position_saving *savingvector);

void point_anom_parameters_robustmean(struct orderedobservationlist_robustmean *list, int ii, int p, double penaltyanomaly);

void changepointreturn_robustmean(struct orderedobservationlist_robustmean *list, int n, int p, int* numberofchanges, int** changepoints, int** components, int** startlag, int** endlag);

void pruner_robustmean(struct orderedobservationlist_robustmean *list, int ii, int p, int l, int minseglength, int maxseglength, double totalpenalty);

void changepointreturn_robustmean_online(struct orderedobservationlist_robustmean *mylist, int n, int p, std::vector<int> &out);



} // namespace anomalymv
