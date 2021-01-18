#ifndef _CPARAM_H
#define _CPARAM_H

#include "CHeader.h"  // this code uses c++ library "newmat11" and "Newran03" by R B Davies

class CHyperParam {
  public:
    CHyperParam();
    ~CHyperParam();
    void init(double h, double aphi, double bphi, double aalpha, double balpha, double fsigma);
    
    double h_Mu;                        				// (compare 1.0, 5.0, 10.0) Small h_mu makes the drawn value spread
  	double a_Phi; double b_Phi;          // Phi_ll ~ Gamma(a_Phi,b_Phi)
    
  	double a_alpha; double b_alpha;       //  nu_j ~ Beta(1,alpha), alpha ~ Gamma(a_alpha,b_alpha)
    
    double f_Sigma;        // Sigma ~ IW(f_Sigma,Phi)
    
};

class CParam {
public:
  CParam();
  
	virtual ~CParam(); //Destructor
  // constants copied over for convienence 
	int n_var;
  int n_sample;
  int n_var_independent;
  int K; //number of components
  
  //output parameters
  double alpha;
  SymmetricMatrix Phi;
  
  Matrix Mu;
  vector<SymmetricMatrix> SIGMA; //Sigma
  ColumnVector pi; 
  
  Matrix Y_in;
  ColumnVector logUnif_y_tilde;
  
  ColumnVector z_in;
  ColumnVector n_z; 
  ColumnVector n_z_in; // Ver 1.0.3
  
  Matrix S_Mat;
  
  RowVector accept_rate;
  RowVector is_accept;
  double Prob_A;
  
  void iterate(int iter, CData &Data, CFeasibilityMap &FM, CHyperParam &hyper, Uniform &randUnif, int n_simul);
  void initialize(CData &Data, int nComp,CFeasibilityMap &FM,Uniform &randUnif, int n_simul);
  
  ColumnVector GetComponentCounts();
  int msg_level; //0 errors only; 1: errors and warnings; 2: errors, warnings and information
  
private:
  
  ColumnVector logpi;
  ColumnVector z_aug;
  
  ColumnVector nu_short;
  ColumnVector mu_bar;
  
  vector<LowerTriangularMatrix> LSIGMA;
  vector<LowerTriangularMatrix> LSIGMA_i; //inverse of LSIGMA
  
  // initial log conditional densities  
  Matrix Sigma_k_inv_ll;
  Matrix Y_in_compact;
  Matrix Y_aug_compact;
  
  void CalculateInitProbA(CData &Data);
  void init_pi(); //init nu_short  pi_short pi logpi
  void init_sigmamu(CData &Data); //init mu mu_bar LSIGMA
  void init_Y_in(CData &Data);
  void init_z_in();
  void init_logUnif_y_tilde(CData &Data, CFeasibilityMap &FM, Uniform &randUnif, int n_simul);
  double calculate_log_cond_norm(CData &Data, int i_original, ColumnVector &item_by_rnorm, 
  ColumnVector &tilde_y_i, ColumnVector &y_q, bool is_q, LowerTriangularMatrix &LSigma_1_i, ColumnVector &s_q); // MODIFIED 2015/02/16;
  
  //samplers
  //input: z_in, Mu, LSIGMA,  logUnif_y_tilde
  //output: logUnif_y_tilde, Y_in, S_Mat
  void S1(int iter, Uniform &randUnif, CData &Data, CFeasibilityMap &FM, int n_simul);
          
  //input: Y_in, pi, Mu, LSIGMA
  //output: Y_in_compact,  z_in
  void S2_add(Uniform &randUnif,CData &Data);
  
  void S3_Z_in(CData &Data);
  
  //input: Y_in, pi, Mu, LSIGMA,z_in
  //output: Y_in_compact,  Prob_A,Y_aug_compact,X_bar,z_aug,n_z
  //return ProbA
  void S4_Z_out(CData &Data);
  
  //input: n_z,X_bar,Phi,Y_aug_compact,z_aug
  //output: Mu,  LSIGMA,Sigma_k_inv_ll
  void S5_MuSigma(CData &Data, double f_Sigma,double h_Mu);
  
  //input: n_z, alpha
  //output: nu_short, logpi,pi_short,pi
  void S6_pi();
  
  //input: logpi
  //output: alpha
  void S7_alpha(double a_alpha, double b_alpha);
  
  //input: Sigma_k_inv_ll
  //output: Phi
  void S8_Phi(double f_Sigma, double a_Phi, double b_Phi);
  
  Matrix X_bar; //holds cluster mean in columns
  ColumnVector logdet_and_more;
  
  //tuning parameter, make it public later on;
  double min_Prob_A;
  int toolarge_nout;
};

#endif
