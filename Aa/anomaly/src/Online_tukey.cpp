#include <math.h>
#include <stdlib.h>
#include <list>
#include "tukey.h"
#include "Online_tukey.h"

Online_tukey::Online_tukey()
{

	cumsumofsquares = 0;
	tukey_object New;
	object_list.push_back(New);

};

void Online_tukey::Add_observation(double observation, double observation_squared, const double threshold, const double threshold_squared)
{

	double lower = observation-threshold;
	double upper = observation+threshold;

	cumsumofsquares += observation_squared;

	std::list<tukey_object>::iterator it = object_list.begin();

	while (it->end < lower)
	{
		it->add_constant(threshold_squared);
		it++;
	}

	tukey_object New_left(*it,lower);
	tukey_object New_right(lower,*it);

	it = object_list.erase(it);
	object_list.insert(it,New_left);
	object_list.insert(it,New_right);

	it--;
	it--;

	it->add_constant(threshold_squared);
	it++;

	while (it->end < upper)
	{
		it->add_an_x_and_x_squared(observation,observation_squared);
		it++;
	}

	tukey_object New_New_left(*it,upper);
	tukey_object New_New_right(upper,*it);

	it = object_list.erase(it);
	object_list.insert(it,New_New_left);
	object_list.insert(it,New_New_right);

	it--;
	it--;

	it->add_an_x_and_x_squared(observation,observation_squared);
	it++;

	while (it != object_list.end())
	{
		it->add_constant(threshold_squared);
		it++;
	}

}

double Online_tukey::Find_mean()
{

	double Min = object_list.begin()->min_cost;
	double Mean = object_list.begin()->mean_of_x;

	if (Mean > object_list.begin()->end){Mean = object_list.begin()->end;}
	if (Mean < object_list.begin()->start){Mean = object_list.begin()->start;}

	for (std::list<tukey_object>::iterator it = object_list.begin(); it != object_list.end(); it ++)
	{
		
		if (it->min_cost < Min)
		{

			Min  = it->min_cost;
			Mean = it->mean_of_x;

			if (Mean > it->end){Mean = it->end;}
			if (Mean < it->start){Mean = it->start;}
		}

	}

	return(Mean);

};

double Online_tukey::Find_minimum()
{

	double Min = object_list.begin()->min_cost;

	for (std::list<tukey_object>::iterator it = object_list.begin(); it != object_list.end(); it ++)
	{
		
		if (it->min_cost < Min){Min = it->min_cost;}

	}

	return(cumsumofsquares - Min);

};
