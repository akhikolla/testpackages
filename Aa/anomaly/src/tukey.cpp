#include <math.h>
#include <stdlib.h>
#include <list>
#include "tukey.h"

void tukey_object::add_an_x_and_x_squared(double x, double x_squared)
{

	number_of_observations++; 
	mean_of_x_squared += (x_squared-mean_of_x_squared)/number_of_observations;
	mean_of_x         += (x-mean_of_x)/number_of_observations;	

	min_cost  = number_of_observations*(mean_of_x_squared - mean_of_x*mean_of_x) + penalty_contribution;

	if (mean_of_x > end)
	{

		min_cost +=  number_of_observations*(mean_of_x-end)*(mean_of_x-end);

	}

	if (mean_of_x < start)
	{

		min_cost += number_of_observations*(start-mean_of_x)*(start - mean_of_x);

	}

};

void tukey_object::add_constant(double h)
{

	min_cost += h;
	penalty_contribution += h;	

};

tukey_object::tukey_object()
{

	mean_of_x_squared      = 0;
	mean_of_x              = 0;
	penalty_contribution   = 0;
	number_of_observations = 0;
	start                  = -999999999999;
	end                    = +999999999999;
	penalty_contribution   = 0;
	min_cost               = 0;
	
};

tukey_object::tukey_object(double new_start, tukey_object old)
{

	mean_of_x_squared      = old.mean_of_x_squared;
	mean_of_x              = old.mean_of_x;
	penalty_contribution   = old.penalty_contribution;
	number_of_observations = old.number_of_observations;
	start                  = new_start;
	end                    = old.end;
	min_cost               = old.min_cost;
	if (mean_of_x < start)
	{

		min_cost += number_of_observations*(start-mean_of_x)*(start - mean_of_x);

		if (mean_of_x < old.start)
		{
			min_cost -= number_of_observations*(mean_of_x-old.start)*(mean_of_x-old.start);
		}

	} 
	
};

tukey_object::tukey_object(tukey_object old ,double new_end)
{

	mean_of_x_squared      = old.mean_of_x_squared;
	mean_of_x              = old.mean_of_x;
	penalty_contribution   = old.penalty_contribution;
	number_of_observations = old.number_of_observations;
	start                  = old.start;
	end                    = new_end;
	min_cost               = old.min_cost;
	if (mean_of_x > end)
	{

		min_cost += number_of_observations*(mean_of_x-end)*(mean_of_x-end);

		if (mean_of_x > old.end)
		{
			min_cost -= number_of_observations*(mean_of_x-old.end)*(mean_of_x-old.end);
		}

	} 
	
};
