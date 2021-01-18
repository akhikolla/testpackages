#include <RcppEigen.h>

typedef struct Particle 
{

	double log_prob;

	Eigen::MatrixXd mu;
	Eigen::MatrixXd Sigma;

	int horizon;
	int id;
	int position;
	const struct Particle* ancestor;

	int anomaly_type;
	int anomaly_comp;
	
	double anomaly_strength; 

	Eigen::MatrixXd obs_Pred;
	Eigen::MatrixXd obs_Prec;

} Particle;

std::list < struct Particle > Get_Particle_Innovative_Descendents(const struct Particle & Ancestor, const Eigen::MatrixXd & relevant_Y, const Eigen::MatrixXd & C_Augmented, const std::vector < double > & sigma_hat, 
const Eigen::MatrixXd & Sigma_Inn, int Number, double s, const std::vector <double> & General_Weights, const std::vector <int> & Sample_From, const int & horizon, const std::vector<int> & Number_of_resamples);

std::list < struct Particle > Get_Particle_Additive_Descendents(const struct Particle & Ancestor, const double & log_likelihood, const int & Number, const double & s, const std::vector<double> sigma_tilde, const Eigen::MatrixXd & Sigma_Add, 
const Eigen::MatrixXd & Pre_Numerator, const Eigen::MatrixXd & Pre_Denominator, const std::vector <double> General_Weights_Add);

std::list < struct Particle > Get_Particle_Descendents(const struct Particle & Ancestor, const Eigen::MatrixXd  & relevant_Y, const std::vector < double > & sigma_tilde, const Eigen::MatrixXd & Sigma_Add,
 int Number, double s, const std::vector <double> & General_Weights_Add);

std::list < struct Particle > Get_Particle_Descendents_typical(const struct Particle & Ancestor, double likelihood);

std::list < struct Particle > Get_Particle_Descendents_V(const struct Particle & Ancestor, int ii, int Number, double likelihood, const double & sigma_tilde, double Sigma_Add_comp, double Numerator, double Denominator, 
double s, double General_weight);

std::list < struct Particle > Get_Particle_Descendents_W(const struct Particle & Ancestor, int ii, int Number, double likelihood, const double & sigma_hat, double Sigma_Inn_comp, double Numerator, double Denominator, 
double s, double General_weight, int horizon);

std::list < struct Particle > Get_Particle_Innovative_Descendents(const struct Particle & Ancestor, const Eigen::MatrixXd & relevant_Y, const Eigen::MatrixXd & C_Augmented, const std::vector < double > & sigma_hat, 
const Eigen::MatrixXd & Sigma_Inn, int Number, double s, const std::vector <double> & prob_inn, const std::vector <int> & Sample_From, const int & horizon);

std::list < struct Particle > Initial_list(const Eigen::MatrixXd & mu_0, const Eigen::MatrixXd & Sigma_0, int Number_of_Particles);

struct Particle initiate_Particle(const double & log_prob, const int & horizon, const struct Particle * ancestor, const int & anomaly_type, const int & anomaly_comp, const double & strength);

std::list< std::list < struct Particle > > Kalman_filter( const std::list<std::list<Eigen::MatrixXd> > & Y_expanded, const std::list<Eigen::MatrixXd> & C_list, const std::list<Eigen::MatrixXd> & Sigma_Add_list, 
const std::list<Eigen::MatrixXd> & Sigma_Inn_Contribution, const Eigen::MatrixXd & A, const Eigen::MatrixXd & Sigma_Inn, const Eigen::MatrixXd & Sigma_Add, const double & s, const int & Num_Descendents, const int & Num_Particles,
const std::list<std::vector<int> > & to_sample, const std::vector<int> & Number_of_resamples, const std::vector<double> & sigma_tilde, const std::vector<double> & sigma_hat, const Eigen::MatrixXd & mu_0, const Eigen::MatrixXd & Sigma_0,
const int & horizon,  const std::vector <double> & prob_inn, const std::vector <double> & prob_add, int Particle_Number, const std::list<Eigen::MatrixXd> & Y_Full_list);

std::list < struct Particle > Kalman_step(const std::list<Eigen::MatrixXd> & considered_Y, int horizon, int counter, std::list<std::list<struct Particle> > & Particle_history,  const Eigen::MatrixXd & A, 
const std::list<Eigen::MatrixXd> & C_list, const std::list<Eigen::MatrixXd> & Sigma_Add_list, const std::list<Eigen::MatrixXd> & Sigma_Inn_Contribution, const std::vector < double > & sigma_hat, 
const std::vector < double > & sigma_tilde, const std::list<std::vector <int> > & Sample_From_List,const Eigen::MatrixXd & Sigma_Inn, const Eigen::MatrixXd & Sigma_Add, const int & Num_Descendents, const int & s,
const std::vector <double> & General_Weight_Inn, const std::vector <double> & General_Weight_Add, int Particle_Number,  std::list<Eigen::MatrixXd> Y_list);

void prepare_Particles(std::list <struct Particle> & Particle_List, const Eigen::MatrixXd & A, const Eigen::MatrixXd & C, const Eigen::MatrixXd & Sigma_Inn_contribution, const Eigen::MatrixXd & Sigma_Add_contribution);

std::list< std::list < std::list < Eigen::MatrixXd > > > Robust_filter( const std::list<std::list<Eigen::MatrixXd> > & Y_expanded, const std::list<Eigen::MatrixXd> & C_list, const std::list<Eigen::MatrixXd> & Sigma_Add_list, 
const std::list<Eigen::MatrixXd> & Sigma_Inn_Contribution, const Eigen::MatrixXd & A, const Eigen::MatrixXd & Sigma_Inn, const Eigen::MatrixXd & Sigma_Add, const double & s, const int & Num_Descendents, const int & Num_Particles,
const std::list<std::vector<int> > & to_sample, const std::vector<int> & Number_of_resamples, const std::vector<double> & sigma_tilde, const std::vector<double> & sigma_hat, const Eigen::MatrixXd & mu_0, const Eigen::MatrixXd & Sigma_0,
const int & horizon, const std::vector <double> & prob_inn, const std::vector <double> & prob_add, int Particle_Number, const std::list<Eigen::MatrixXd> & Y_Full_list);

std::list < struct Particle > Subsample_Particles(std::list < struct Particle > & candidates, const int &N);

void update_Particle(struct Particle & Sampled_Particle, const Eigen::MatrixXd & A, const Eigen::MatrixXd & C, const Eigen::MatrixXd & Sigma_Inn, const Eigen::MatrixXd & Sigma_Add, std::list<Eigen::MatrixXd> Y);

