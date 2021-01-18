#include <stdbool.h>
#include "user_interupt.h"
#include "check_user_interrupt.h"
#include <vector>

namespace anomalymv
{

typedef struct position_saving
{

	double saving;
	int position;

} position_saving;

typedef struct orderedobservationlist_mean 
{

	int numberofobservation;
	
	double* observation;
	double* mean_of_xs;

	double* segmentcosts;
	double* best_end_costs;

	double optimalcostofprevious;
	double costofstartingsegment;
	double optimalcost;
	
	int* affectedcomponents;
	int* startlag;
	int* endlag;

	struct orderedobservationlist_mean* optimalcut;
	int option;

	int    destruction;
  	struct orderedobservationlist_mean* next;
  	struct orderedobservationlist_mean* previous;

} orderedobservationlist;


void populate_mean(struct orderedobservationlist_mean **list, double* x , int n, int p, int l );

int cmpfunc_nosorting (const void * a, const void * b);

int cmpfunc_sorting (const void * a, const void * b);

void solveorderedobservationlist_mean(struct orderedobservationlist_mean *list, int n, int p, int l, double* penaltycomponent, double penaltyanomaly, int minseglength, int maxseglength);

void compute_cost_of_starting_anomalies_mean(struct orderedobservationlist_mean *list , int ii, int n, int p, int l, int minseglength, double *penaltycomponent, double *componentcost);

void update_cumsums_and_segmentcosts_mean(struct orderedobservationlist_mean *list, int ii, int n, int p, int l, int minseglength);

double find_lowest_end_cost(double* segmentcosts, int jj, int p, int l);

void find_best_option_mean(struct orderedobservationlist_mean *list, int ii, int n, int p, int l, int minseglength, double *penaltycomponent, double penaltyanomaly, struct position_saving *savingvector);

void collective_anom_parameters_mean(struct orderedobservationlist_mean *list, int ii, int p, int l, int minseglength, double *penaltycomponent, struct position_saving *savingvector);

void point_anom_parameters_mean(struct orderedobservationlist_mean *list, int ii, int p, double penaltyanomaly);

void changepointreturn_mean(struct orderedobservationlist_mean *list, int n, int p, int* numberofchanges, int** changepoints, int** components, int** startlag, int** endlag);

void pruner_mean(struct orderedobservationlist_mean *list, int ii, int p, int l, int minseglength, int maxseglength, double totalpenalty);

void changepointreturn_mean_online(struct orderedobservationlist_mean *mylist, int n, int p, std::vector<int> &out);



} // namespace anomalymv
