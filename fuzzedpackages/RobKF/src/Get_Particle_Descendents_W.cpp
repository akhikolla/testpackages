#include "Particle.h"
#include <RcppEigen.h>
#include <R.h>

std::list < struct Particle > Get_Particle_Descendents_W(const struct Particle & Ancestor, int ii, int Number, double likelihood, const double & sigma_hat, double Sigma_Inn_comp, double Numerator, double Denominator, 
double s, double General_weight, int horizon)
{

	std::list < struct Particle > out;

	double sampled_scale, Common_weight, specific_weight, parameter, tmp, factor;

	tmp = Sigma_Inn_comp*Denominator;

	factor = Numerator/Denominator;

	parameter = 1/(s + 0.5*sigma_hat*factor/tmp);

	Common_weight = (s+0.5)*std::log(parameter);
	
	for (int jj = 0; jj < Number; jj++)
	{
		
		sampled_scale  = sigma_hat*(R::rgamma(s + 0.5, parameter));

		specific_weight = -0.5*(std::log(tmp+sampled_scale) - (1 +  sampled_scale*sampled_scale/tmp/(tmp+sampled_scale))  * factor );

		out.push_back(initiate_Particle(Ancestor.log_prob + likelihood + General_weight + Common_weight + specific_weight, horizon, &Ancestor, 1, ii, sampled_scale));

	}

	return(out);
	
};
