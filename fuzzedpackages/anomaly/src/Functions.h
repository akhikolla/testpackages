#include <stdbool.h>
#include "tukey.h"
#include "Online_tukey.h"
#include "user_interupt.h"
#include <vector>

namespace anomaly
{

typedef struct orderedobservationlist 
{
	int    numberofobservation;
	double observation;
  	double observationsquared;

	double cumulativesum;
	double cumulativesumofsquares;
	double optimalcostofprevious;
	double segmentcost;
	
	double optimalcost;
	struct orderedobservationlist* optimalcut;
	int option;

	int    destruction;
  	struct orderedobservationlist* next;
  	struct orderedobservationlist* previous;
} orderedobservationlist;

void populateorderedobservationlist(struct orderedobservationlist **list, double* x , int n);

void updatewithobservation(int ii, struct orderedobservationlist *list, double* penaltychange);

void findoptimaloption(int ii, struct orderedobservationlist *list, int minseglength, double penaltyoutlier);

void solveorderedobservationlist(struct orderedobservationlist *list, int n, double* penaltychange, double penaltyoutlier, int minseglength, int maxseglength);

void changepointreturn(struct orderedobservationlist *list, int n, int* numberofchanges, int** changepoints);

void changepointreturn_online(struct orderedobservationlist *list, int n, int** changepoints);

void pruner(struct orderedobservationlist *list, int ii, double penaltychange_max, int minseglength, int maxseglength);

typedef struct orderedobservationlist_mean 
{
	int    numberofobservation;
	double observation;

	double cumulativesum;
	double optimalcostofprevious;
	double segmentcost;
	
	double optimalcost;
	struct orderedobservationlist_mean* optimalcut;
	int option;

	int    destruction;
  	struct orderedobservationlist_mean* next;
  	struct orderedobservationlist_mean* previous;
} orderedobservationlist_mean;

void populateorderedobservationlist_mean(struct orderedobservationlist_mean **list, double* x , int n);

void updatewithobservation_mean(int ii, struct orderedobservationlist_mean *list, double* penaltychange);

void findoptimaloption_mean(int ii, struct orderedobservationlist_mean *list, int minseglength, double penaltyoutlier);

void solveorderedobservationlist_mean(struct orderedobservationlist_mean *list, int n, double* penaltychange, double penaltyoutlier, int minseglength, int maxseglength);

void changepointreturn_mean(struct orderedobservationlist_mean *list, int n, int* numberofchanges, int** changepoints);

void changepointreturn_online_mean(struct orderedobservationlist_mean *list, int n, int** changepoints);

void pruner_mean(struct orderedobservationlist_mean *list, int ii, double penaltychange_max, int minseglength, int maxseglength);

typedef struct orderedobservationlist_robustmean 
{
	int numberofobservation;
	double observation;
	double observationsquared;
	Online_tukey *Tukey_Stuff;

	double optimalcostofprevious;
	double segmentcost;
	
	double optimalcost;
	struct orderedobservationlist_robustmean* optimalcut;
	int option;

	int    destruction;
  	struct orderedobservationlist_robustmean* next;
  	struct orderedobservationlist_robustmean* previous;
} orderedobservationlist_robustmean;

void populateorderedobservationlist_robustmean(struct orderedobservationlist_robustmean **list, double* x , int n);

void updatewithobservation_robustmean(int ii, struct orderedobservationlist_robustmean *list, double* penaltychange, const double threshold, const double threshold_squared);

void findoptimaloption_robustmean(int ii, struct orderedobservationlist_robustmean *list, int minseglength, double penaltyoutlier);

void solveorderedobservationlist_robustmean(struct orderedobservationlist_robustmean *list, int n, double* penaltychange, double penaltyoutlier, int minseglength, int maxseglength);

void changepointreturn_robustmean(struct orderedobservationlist_robustmean *list, int n, int* numberofchanges, int** changepoints);

void changepointreturn_online_robustmean(struct orderedobservationlist_robustmean *list, int n, int** changepoints);

void pruner_robustmean(struct orderedobservationlist_robustmean *list, int ii, double penaltychange_max, int minseglength, int maxseglength);


} // namespace anomaly


namespace anomalymv
{

typedef struct position_saving
{

	double saving;
	int position;

} position_saving;

typedef struct orderedobservationlist 
{

	int numberofobservation;
	
	double* observation;
	double* observationsquared;
	double* mean_of_xs;
	double* mean_of_xs_squared;

	double* segmentcosts;
	double* best_end_costs;

	double optimalcostofprevious;
	double costofstartingsegment;
	double optimalcost;
	
	int* affectedcomponents;
	int* startlag;
	int* endlag;

	struct orderedobservationlist* optimalcut;
	int option;

	int    destruction;
  	struct orderedobservationlist* next;
  	struct orderedobservationlist* previous;

} orderedobservationlist;


void populate(struct orderedobservationlist **list, double* x , int n, int p, int l );

int cmpfunc_nosorting (const void * a, const void * b);

int cmpfunc_sorting (const void * a, const void * b);

void solveorderedobservationlist(struct orderedobservationlist *list, int n, int p, int l, double* penaltycomponent, double penaltyanomaly, int minseglength, int maxseglength);

void compute_cost_of_starting_anomalies(struct orderedobservationlist *list , int ii, int n, int p, int l, int minseglength, double *penaltycomponent, double *componentcost);

void update_cumsums_and_segmentcosts(struct orderedobservationlist *list, int ii, int n, int p, int l, int minseglength);

double find_lowest_end_cost(double* segmentcosts, int jj, int p, int l);

void find_best_option(struct orderedobservationlist *list, int ii, int n, int p, int l, int minseglength, double *penaltycomponent, double penaltyanomaly, struct position_saving *savingvector);

void collective_anom_parameters(struct orderedobservationlist *list, int ii, int p, int l, int minseglength, double *penaltycomponent, struct position_saving *savingvector);

void point_anom_parameters(struct orderedobservationlist *list, int ii, int p, double penaltyanomaly);

void changepointreturn(struct orderedobservationlist *list, int n, int p, int* numberofchanges, int** changepoints, int** components, int** startlag, int** endlag);

void pruner(struct orderedobservationlist *list, int ii, int p, int l, int minseglength, int maxseglength, double totalpenalty);

void changepointreturn_online(struct orderedobservationlist *mylist, int n, int p, std::vector<int> &out);



} // namespace anomalymv
