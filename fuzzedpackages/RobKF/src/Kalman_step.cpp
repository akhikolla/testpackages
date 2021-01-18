#include "Particle.h"
#include <RcppEigen.h>

std::list < struct Particle > Kalman_step(const std::list<Eigen::MatrixXd> & considered_Y, int horizon, int counter, std::list<std::list<struct Particle> > & Particle_history,  const Eigen::MatrixXd & A, 
const std::list<Eigen::MatrixXd> & C_list, const std::list<Eigen::MatrixXd> & Sigma_Add_list, const std::list<Eigen::MatrixXd> & Sigma_Inn_Contribution, const std::vector < double > & sigma_hat, 
const std::vector < double > & sigma_tilde, const std::list<std::vector <int> > & Sample_From_List,const Eigen::MatrixXd & Sigma_Inn, const Eigen::MatrixXd & Sigma_Add, const int & Num_Descendents, const int & s,
const std::vector <double> & General_Weight_Inn, const std::vector <double> & General_Weight_Add, int Particle_Number,  std::list<Eigen::MatrixXd> Y_list) 
{

	std::list < struct Particle > Candidates;
	std::list < struct Particle > Output;
	std::list < struct Particle > Candidates_addition;

	std::list < int > Candidate_data;
	std::list < int > Candidate_data_addition;

	std::list< std::list<struct Particle> >::iterator History_iterator       = Particle_history.end();
	std::list<struct Particle>::iterator              Particle_iterator;

	std::list< std::vector<int> >::const_iterator         Sampling_iterator  = Sample_From_List.begin();

	std::list<Eigen::MatrixXd>::const_iterator            Y_iterator         = considered_Y.begin();
	std::list<Eigen::MatrixXd>::const_iterator            C_iterator         = C_list.begin();
	std::list<Eigen::MatrixXd>::const_iterator            Sigma_Inn_iterator = Sigma_Inn_Contribution.begin();
	std::list<Eigen::MatrixXd>::const_iterator            Sigma_Add_iterator = Sigma_Add_list.begin();

	// Get W updates first

	for (int ii = 0; ii < horizon; ii++)
	{

		History_iterator--;

		prepare_Particles(*History_iterator,A,*C_iterator,*Sigma_Inn_iterator,*Sigma_Add_iterator);

		for (Particle_iterator = History_iterator->begin(); Particle_iterator != History_iterator->end(); Particle_iterator++)
		{

			Candidates_addition = Get_Particle_Innovative_Descendents(*Particle_iterator, *Y_iterator, *C_iterator, sigma_hat, Sigma_Inn, Num_Descendents, s, General_Weight_Inn,*Sampling_iterator, ii+1);

			Candidates.splice(Candidates.end(),Candidates_addition);

		}

		Y_iterator++;
		C_iterator++;
		Sampling_iterator++;
		Sigma_Inn_iterator++;
		Sigma_Add_iterator++;
		
	}

	for (Particle_iterator = Particle_history.back().begin(); Particle_iterator != Particle_history.back().end(); Particle_iterator++)
	{

		Candidates_addition = Get_Particle_Descendents(*Particle_iterator, considered_Y.front(), sigma_tilde, Sigma_Add, Num_Descendents, s, General_Weight_Add);

		Candidates.splice(Candidates.end(),Candidates_addition);

	}

	Output = Subsample_Particles(Candidates,Particle_Number);

	int position = 0;

	for (Particle_iterator = Output.begin(); Particle_iterator != Output.end(); Particle_iterator++)
	{

		update_Particle(*Particle_iterator, A, C_list.front(), Sigma_Inn, Sigma_Add, Y_list);
		Particle_iterator->id       = counter;
		Particle_iterator->position = position;

		position++;

	}

	return Output;

};

