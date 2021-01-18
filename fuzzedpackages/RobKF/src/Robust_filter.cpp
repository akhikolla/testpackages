#include "Particle.h"
#include <RcppEigen.h>

#include <list>
#include <vector>

#include "exception.h"
#include "user_interupt.h"
#include "check_user_interrupt.h"

// [[Rcpp::export]]
std::list< std::list <  std::list < Eigen::MatrixXd > > > Robust_filter( const std::list<std::list<Eigen::MatrixXd> > & Y_expanded, const std::list<Eigen::MatrixXd> & C_list, const std::list<Eigen::MatrixXd> & Sigma_Add_list, 
const std::list<Eigen::MatrixXd> & Sigma_Inn_Contribution, const Eigen::MatrixXd & A, const Eigen::MatrixXd & Sigma_Inn, const Eigen::MatrixXd & Sigma_Add, const double & s, const int & Num_Descendents, const int & Num_Particles,
const std::list<std::vector<int> > & to_sample, const std::vector<int> & Number_of_resamples, const std::vector<double> & sigma_tilde, const std::vector<double> & sigma_hat, const Eigen::MatrixXd & mu_0, const Eigen::MatrixXd & Sigma_0,
const int & horizon, const std::vector <double> & prob_inn, const std::vector <double> & prob_add, int Particle_Number, const std::list<Eigen::MatrixXd> & Y_Full_list)
{

	std::list< std::list < struct Particle > > pre_out = Kalman_filter(Y_expanded, C_list, Sigma_Add_list, Sigma_Inn_Contribution, A, Sigma_Inn, Sigma_Add, s, Num_Descendents, Num_Particles, to_sample, Number_of_resamples,
                                                                           sigma_tilde, sigma_hat, mu_0, Sigma_0, horizon, prob_inn, prob_add, Particle_Number, Y_Full_list);



	std::list< std::list < std::list < Eigen::MatrixXd > > > Out;

	std::list< struct Particle >::iterator it_List;

	std::list< std::list < struct Particle > >::iterator it_List_of_List;

	it_List_of_List = pre_out.begin();

	std::list < std::list < Eigen::MatrixXd > > State;

	for (it_List = it_List_of_List->begin(); it_List != it_List_of_List->end(); it_List++)
		{

      			if(check_user_interrupt())
      			{
	  			throw_exception("User interrupt");
      			}

			std::list < Eigen::MatrixXd > Particle;

			Particle.push_back(it_List->mu);

			Particle.push_back(it_List->Sigma);

			Eigen::MatrixXd Type(1,1);

			Type << it_List->anomaly_type;

			Particle.push_back(Type);

			Eigen::MatrixXd Comp(1,1);

			Comp << it_List->anomaly_comp;

			Particle.push_back(Comp);

			Eigen::MatrixXd id(1,1);

			id << it_List->id;

			Particle.push_back(id);	

			Eigen::MatrixXd position(1,1);

			position << -1;

			Particle.push_back(position);

			Eigen::MatrixXd horizon(1,1);

			horizon << it_List->horizon;

			Particle.push_back(horizon);		

			Eigen::MatrixXd Strength(1,1);

			Strength << it_List->anomaly_strength;

			Particle.push_back(Strength);	

			State.push_back(Particle);

		}

	Out.push_back(State);

	it_List_of_List++;

	for (; it_List_of_List != pre_out.end(); it_List_of_List++)
	{

		std::list < std::list < Eigen::MatrixXd > > State;

		for (it_List = it_List_of_List->begin(); it_List != it_List_of_List->end(); it_List++)
		{

			std::list < Eigen::MatrixXd > Particle;

			Particle.push_back(it_List->mu);

			Particle.push_back(it_List->Sigma);

			Eigen::MatrixXd Type(1,1);

			Type << it_List->anomaly_type;

			Particle.push_back(Type);

			Eigen::MatrixXd Comp(1,1);

			Comp << it_List->anomaly_comp;

			Particle.push_back(Comp);

			Eigen::MatrixXd id(1,1);

			id << it_List->id;

			Particle.push_back(id);	

			Eigen::MatrixXd position(1,1);

			position << it_List->ancestor->position;

			Particle.push_back(position);

			Eigen::MatrixXd horizon(1,1);

			horizon << it_List->horizon;

			Particle.push_back(horizon);		

			Eigen::MatrixXd Strength(1,1);

			Strength << it_List->anomaly_strength;

			Particle.push_back(Strength);	

			State.push_back(Particle);

		}

		Out.push_back(State);


	}

	return(Out);

};
