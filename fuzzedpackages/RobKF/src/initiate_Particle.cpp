#include "Particle.h"

struct Particle initiate_Particle(const double & log_prob, const int & horizon, const struct Particle * ancestor, const int & anomaly_type, const int & anomaly_comp, const double & strength)
{
	
	struct Particle out;

	out.log_prob = log_prob;

	out.horizon  = horizon;
	out.ancestor = ancestor;

	out.anomaly_type = anomaly_type;
	out.anomaly_comp = anomaly_comp;

	out.anomaly_strength = strength;

	return(out);	

};
