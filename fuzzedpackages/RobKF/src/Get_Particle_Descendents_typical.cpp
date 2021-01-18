#include "Particle.h"
#include <RcppEigen.h>

std::list < struct Particle > Get_Particle_Descendents_typical(const struct Particle & Ancestor, double likelihood)
{

	std::list < struct Particle > out;

	out.push_back(initiate_Particle(Ancestor.log_prob + likelihood, 1, &Ancestor, 0, -1, -1));

	return(out);
	
};
