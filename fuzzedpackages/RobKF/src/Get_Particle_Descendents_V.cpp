#include "Particle.h"
#include <RcppEigen.h>
#include <R.h>
#include <iostream>

std::list < struct Particle > Get_Particle_Descendents_V(const struct Particle & Ancestor, int ii, int Number, double likelihood, const double & sigma_tilde, double Sigma_Add_comp, double Numerator, double Denominator, 
double s, double General_weight)
{

	std::list < struct Particle > out;

	double sampled_scale, Common_weight, specific_weight, parameter, tmp, factor;

	tmp = Sigma_Add_comp*Denominator;

	factor = Numerator/Denominator;

	parameter = 1/(s + 0.5*sigma_tilde*factor/tmp);

	Common_weight = (s+0.5)*std::log(parameter);



	for (int jj = 0; jj < Number; jj++)
	{
		
		sampled_scale  = sigma_tilde*(R::rgamma(s + 0.5, parameter));

		specific_weight = -0.5*(  std::log(tmp+sampled_scale) - (1 +  sampled_scale*sampled_scale/tmp/(tmp+sampled_scale)) * factor );

		/*std::cout << std::endl;
		std::cout << "sampled_scale_stuff";
		std::cout << std::endl;
		std::cout << Ancestor.log_prob;
		std::cout << std::endl;
		std::cout << s + 0.5;		
		std::cout << std::endl;
		std::cout << 1/parameter;
		std::cout << std::endl;
		std::cout << likelihood;
		std::cout << std::endl;
		std::cout << General_weight;
		std::cout << std::endl;
		std::cout << Common_weight;
		std::cout << std::endl;
		std::cout << specific_weight;
		std::cout << std::endl;
		std::cout << 0.5*(1 +  sampled_scale*sampled_scale/tmp/tmp)*tmp/(tmp+sampled_scale) * factor;
		std::cout << std::endl;
		std::cout << -0.5*(  log(tmp+sampled_scale) ) ;
		std::cout << std::endl;
		std::cout << specific_weight+likelihood;
		std::cout << std::endl; */

		out.push_back(initiate_Particle(Ancestor.log_prob + likelihood + General_weight + Common_weight + specific_weight, 1, &Ancestor, 2, ii, sampled_scale));

	}

	return(out);
	
};
