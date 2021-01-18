#include <math.h>
#include <stdlib.h>

#ifndef TUKEY_OBJECT
#define TUKEY_OBJECT

class tukey_object
{

	public:
	double mean_of_x_squared,mean_of_x,penalty_contribution,start,end;
	double min_cost;
	int number_of_observations;
	tukey_object();
	tukey_object(tukey_object,double);
	tukey_object(double,tukey_object);
	void add_an_x_and_x_squared(double,double);
	void add_constant(double);

};

#endif
