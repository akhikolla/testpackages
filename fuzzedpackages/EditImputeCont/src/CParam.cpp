#include <R.h>
#include "CData.h"
#include "CFeasibilityMap.h"
#include "CParam.h"
#define LOG_2_PI 1.83787706640935

CHyperParam::CHyperParam()  {
  h_Mu=5.0;
  a_Phi=0.25;
  b_Phi=0.25 ;
  a_alpha=0.25;
  b_alpha=0.25;
  f_Sigma = -1;  //need to be inilized by data
}

CHyperParam::~CHyperParam()  {
}

void CHyperParam::init(double h, double aphi, double bphi, double aalpha, double balpha, double fsigma)  {
  h_Mu=h;
  a_Phi=aphi;
  b_Phi=bphi;
  a_alpha=aalpha;
  b_alpha=balpha;
  f_Sigma = fsigma;  //need to be inilized by data
}

CParam::CParam()  {
  min_Prob_A = 0.01;
  accept_rate = RowVector(2) ; accept_rate = 0 ;
  is_accept = RowVector(2) ; is_accept = 0;
  Prob_A = 0;
  K = -1; //to indicate no number of components has been set
  msg_level = 0; //0 errors only; 1: errors and warnings; 2: errors, warnings and information
}

CParam::~CParam() {
}

void CParam::iterate(int iter, CData &Data, CFeasibilityMap &FM, CHyperParam &hyper, Uniform &randUnif, int n_simul) {
  if (hyper.f_Sigma < 0) { //if we use default hyper parameters
  hyper.f_Sigma  = Data.n_var-Data.n_balance_edit+1;
  }
  if (K == -1) {
	Rprintf( "The number of components has not been set\n");
    return;
  }

  if (K < 0) { //the K has been changed and delayed re-initilization is needed
    //cout << "Initilizing parameters" << endl;
    initialize(Data,-K,FM,randUnif,n_simul);
  }
  //clock_t start = clock(), diff;

  // if ( msg_level >= 2 ) cout << "iter=" << iter << endl ; // Ver 1.1.0
    
  // S1. // vector_r
  // if ( msg_level >= 2 ) cout << "S1" << endl ; // Ver 1.1.0
  S1(iter, randUnif, Data, FM, n_simul);
  /*
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  
  << "s1:" << msec << endl;
  start = clock();
  */

  // S2-add. // Y_in | S_i for sum(s_i) >= 2
  // if ( msg_level >= 2 ) cout << "S2_add" << endl ; // Ver 1.1.0
  S2_add(randUnif,Data);

  // if ( msg_level >= 2 ) cout << "S3_Z_in" << endl ; // Ver 1.1.0
  S3_Z_in(Data);
  
  // if ( msg_level >= 2 ) cout << "S4_Z_out" << endl ; // Ver 1.1.0
  S4_Z_out(Data);

  // if ( msg_level >= 2 ) cout << "S5_MuSigma" << endl ; // Ver 1.1.0
  S5_MuSigma(Data, hyper.f_Sigma,hyper.h_Mu);
  
  // if ( msg_level >= 2 ) cout << "S6_pi" << endl ; // Ver 1.1.0
  S6_pi();
  
  // if ( msg_level >= 2 ) cout << "S7_alpha" << endl ; // Ver 1.1.0
  S7_alpha(hyper.a_alpha, hyper.b_alpha);
  
  // if ( msg_level >= 2 ) cout << "S8_Phi" << endl ; // Ver 1.1.0
  S8_Phi(hyper.f_Sigma, hyper.a_Phi, hyper.b_Phi);
  
}

void CParam::initialize(CData &Data, int nComp, CFeasibilityMap &FM, Uniform &randUnif, int n_simul){
  // Set the minimum value of (no. of correct records)/(no. of all records)
  toolarge_nout = (1-min_Prob_A)/min_Prob_A * Data.n_sample ;

  n_var = Data.n_var;
  n_var_independent = n_var - Data.n_balance_edit;
  K = nComp;
  n_sample = Data.n_sample;
  
  n_z_in = ColumnVector(K); n_z_in = 0 ; // Ver 1.1.0

  alpha = 1.0 ;
  Phi = IdentityMatrix(n_var_independent) * 5.0;
  init_pi();
  init_sigmamu(Data);
  init_Y_in(Data);
  init_z_in();

  S_Mat = Data.initial_S_Mat.t();
  init_logUnif_y_tilde(Data,FM,randUnif,n_simul);

  n_z = ColumnVector(K); n_z = 0 ;
  Sigma_k_inv_ll = Matrix(K,n_var_independent);  Sigma_k_inv_ll=0.0;
  X_bar= Matrix(n_var_independent,K); X_bar = 0.0 ;
  CalculateInitProbA(Data);
  
} // CParam::initialize

void CParam::init_pi() {
  nu_short = ColumnVector(K - 1);
  nu_short = 0.1 ;    // WAS 0.1 starting value

  logpi = ColumnVector(K);
  pi= ColumnVector(K);

  double Sum_logOneMinusVg = 0.0 ;
  for (int i= 1; i<= (K-1); i++) {
    double nu_k = nu_short(i) ;
    logpi(i) = log(nu_k) + Sum_logOneMinusVg ;
    pi(i) = exp(logpi(i)) ;

    double one_minus_V = ( 1.0 - nu_short(i) ) ;
    Sum_logOneMinusVg = Sum_logOneMinusVg + log(one_minus_V) ;
  }
  logpi(K) = Sum_logOneMinusVg ;
  pi(K) = exp(Sum_logOneMinusVg) ;
}

void CParam::init_sigmamu(CData &Data) {
  mu_bar = ColumnVector(n_var_independent);
  ColumnVector log_obs_mean(n_var), log_obs_sd(n_var) ;
  Matrix logD_editpassing(n_sample-Data.n_faulty,n_var) ;
  for (int i_sample=1, count = 1; i_sample<=n_sample; i_sample++){
    if (Data.is_case(i_sample,0)) {
      logD_editpassing.row(count++) = Data.log_D_Observed.row(i_sample);
    }
  }

  for (int i_var=1; i_var<=n_var; i_var++){
    ColumnVector temp_col = logD_editpassing.column(i_var) ;
    log_obs_mean(i_var) = 1.0/(n_sample-Data.n_faulty)*(temp_col.sum()) ;
    ColumnVector temp_mean(n_sample-Data.n_faulty) ; temp_mean = log_obs_mean(i_var) ;
    Matrix SumSq = (temp_col-temp_mean).t()*(temp_col-temp_mean) ;
    log_obs_sd(i_var) = sqrt( 1.0/(n_sample-Data.n_faulty-1)*SumSq(1,1) ) ;
  }

  Data.UpdateCompactVector(mu_bar,log_obs_mean);
  DiagonalMatrix LSigma_temp(n_var_independent);
  for (int i=1; i<= n_var_independent; i++){
    LSigma_temp(i) = log_obs_sd(i);
  }
  DiagonalMatrix LSigma_temp_i = LSigma_temp.i();

  LSIGMA = vector<LowerTriangularMatrix>(K);
  LSIGMA_i = vector<LowerTriangularMatrix>(K);
  SIGMA = vector<SymmetricMatrix>(K);
  for (int k=0 ; k<K; k++) {
    LSIGMA[k] = LSigma_temp;
    LSIGMA_i[k] = LSigma_temp_i;
    SIGMA[k] << LSigma_temp * LSigma_temp; //lazy
  }
  logdet_and_more = ColumnVector(K); logdet_and_more = -0.5*n_var* LOG_2_PI + logdet(LSIGMA_i[0]);

  Mu = Matrix((n_var_independent),K);      // Note that mu_k is Mu.column(k)
  for (int k=1; k<=K; k++){
    Mu.column(k) = mu_bar ;  // starting value
  }
}

void CParam::init_Y_in(CData &Data) {
  Y_in = Matrix(n_sample,n_var) ;
  Real* a = Y_in.data();
  Real* b = Data.D_initial.data();
  for (int i =0; i <n_sample * n_var; i++) { *a++ = log(*b++);}
  Y_in_compact= Matrix(n_sample,n_var_independent);
  Data.UpdateCompactMatrix(Y_in_compact,Y_in);
}

// Ver 1.1.0 moves "void CParam::S6_pi()" below // // Ver 1.1.0

void CParam::init_z_in() {
  
  z_in = ColumnVector(n_sample) ;
  n_z_in = 0 ; // Ver 1.1.0
  
  for (int i_sample=1; i_sample<=n_sample; i_sample++) {
    // calculate pi_tilde_in
    ColumnVector logN_unnorm(K);
    ColumnVector y_compact = Y_in_compact.row(i_sample).t() ;
    for (int k=1; k<=K; k++) {
      ColumnVector mu_k = Mu.column(k);
      logN_unnorm(k) = log_MVN_fn( y_compact, mu_k, LSIGMA_i[k-1],logdet_and_more(k));
    }
    double max_logN_unnorm = logN_unnorm.maximum();

    ColumnVector pi_tilde_in_unnorm(K); pi_tilde_in_unnorm = 0.0 ;
    for (int k=1; k<=K; k++){
      pi_tilde_in_unnorm(k) = pi(k) * exp(logN_unnorm(k)-max_logN_unnorm);
    }
    ColumnVector pi_tilde_in = (1.0/pi_tilde_in_unnorm.sum()) * pi_tilde_in_unnorm ;

    z_in(i_sample) = rdiscrete_fn(pi_tilde_in);
    n_z_in(z_in(i_sample)) = n_z_in(z_in(i_sample)) + 1 ; // Ver 1.1.0
  }
  
} // void CParam::init_z_in()

void CParam::init_logUnif_y_tilde(CData &Data,CFeasibilityMap &FM, Uniform &randUnif, int n_simul) {
  logUnif_y_tilde = ColumnVector(Data.n_faulty);
  for (int i = 1; i <= Data.n_faulty; i++) {
    ColumnVector s_i = S_Mat.column(i);
    int i_original = Data.Faulty2Original[i-1];
    // double log_f_y_tilde_q = 0;  Commented by HANG on 5/16/2015
    if (FM.useMap && Data.is_case(i_original,2)) {
      logUnif_y_tilde(i) = FM.Simulate_logUnif_case2(CFeasibilityMap::s_to_tau_fn(s_i),i_original,n_simul,randUnif);
    } else {
      logUnif_y_tilde(i) = Data.get_logUnif_y_tilde(s_i, CFeasibilityMap::s_to_tau_fn(s_i),i_original); 
    } 
    //logUnif_y_tilde(i) = Data.get_logUnif_y_tilde(s_i, CFeasibilityMap::s_to_tau_fn(s_i),i_original);
  }
}

double CParam::calculate_log_cond_norm(CData &Data, int i_original, ColumnVector &item_by_rnorm, ColumnVector &tilde_y_i, ColumnVector &y_q, bool is_q, LowerTriangularMatrix &LSigma_1_i, ColumnVector &s_q) { // MODIFIED 2015/02/16

  double log_cond_norm;

  if ( item_by_rnorm.sum() >= 1 ) {

    ColumnVector mu_z_i = Mu.column(z_in(i_original));
    ColumnVector s_1_compact = Data.get_compact_vector(item_by_rnorm); 
    ColumnVector Mu_1 = subvector(mu_z_i,s_1_compact);
    Matrix Sigma_1 = Submatrix_elem_2(SIGMA[z_in(i_original)-1],s_1_compact,s_1_compact);

    // ADDED 2015/01/27

    ColumnVector s_q_compact = Data.get_compact_vector(s_q) ; // MODIFIED 2015/02/16 
    ColumnVector VectorOne = s_q_compact ; VectorOne = 1 ; // MODIFIED 2015/02/16
    ColumnVector s_0_compact = VectorOne - s_q_compact ; // MODIFIED 2015/02/16
    int sum_s_0_comp = s_0_compact.sum() ;
    LowerTriangularMatrix LSigma_cond ;
    ColumnVector Mu_cond ;

    if ( sum_s_0_comp>0 ){

      ColumnVector Mu_0 = subvector(mu_z_i,s_0_compact); // (s_1_compact.sum()) vector
      Matrix Sigma_0 = Submatrix_elem_2(SIGMA[z_in(i_original)-1],s_0_compact,s_0_compact);
      Matrix Sigma_10 = Submatrix_elem_2(SIGMA[z_in(i_original)-1],s_1_compact,s_0_compact);
      ColumnVector y_tilde_compact = Data.get_compact_vector(tilde_y_i) ;
      ColumnVector y_tilde_0 = subvector(y_tilde_compact,s_0_compact) ;
      SymmetricMatrix Sigma_0_symm ; Sigma_0_symm << Sigma_0 ;

      LowerTriangularMatrix LSigma_0 = Cholesky(Sigma_0_symm) ;
      Mu_cond = Mu_1 + Sigma_10 * (LSigma_0.i()).t()*LSigma_0.i() * ( y_tilde_0-Mu_0 ) ;
      Matrix Sigma_cond = Sigma_1 - Sigma_10 * (LSigma_0.i()).t()*LSigma_0.i() * Sigma_10.t() ;
      SymmetricMatrix Sigma_cond_symm ; Sigma_cond_symm << Sigma_cond ;
      int sum_s_1_comp = s_1_compact.sum() ;
      DiagonalMatrix D(sum_s_1_comp) ; Matrix V(sum_s_1_comp,sum_s_1_comp) ;
      Jacobi(Sigma_cond_symm,D,V) ;
      int is_zero_exist = 0 ;
      for (int i_var=1; i_var<=sum_s_1_comp; i_var++){
    	if ( D(i_var) < 1e-9 ){
    	  D(i_var) = 1e-9 ;
    	  is_zero_exist = 1 ;
    	}
      } // for (int i_var=1; i_var<=sum_s_1_comp; i_var++)
      if ( is_zero_exist == 1 ){
    	Sigma_cond_symm << V * D * V.t() ;
    	if ( msg_level >= 1 ) {
		  Rprintf( "   Warning: When generating y_j from conditional normal(Mu_-j,Sigma_-j), Sigma_-j is non-positive definite because of computation precision. The eigenvalues D(j,j) smaller than 1e-9 is replaced with 1e-9, and let Sigma_-j = V D V.t().\n");
    	}
      } //
      LSigma_cond = Cholesky(Sigma_cond_symm);
      // y_part = rMVN_fn(Mu_cond,LSigma_cond);
      // log_cond_norm = log_MVN_fn(y_part,Mu_cond,LSigma_cond) ;

    } else {

      Mu_cond = Mu_1 ;
      SymmetricMatrix Sigma_1_symm = Submatrix_elem(SIGMA[z_in(i_original)-1],s_1_compact);
      LSigma_cond = Cholesky(Sigma_1_symm) ;
      // SymmetricMatrix Sigma_1_symm ; Sigma_1_symm << Sigma_1 ;
      // LowerTriangularMatrix LSigma_1 = Cholesky_Sigma_star_symm(Sigma_1_symm);
      // y_part = rMVN_fn(Mu_1,LSigma_1);
      // log_cond_norm = log_MVN_fn(y_part,Mu_1,LSigma_1) ;

    } // if ( sum_s_0_comp>0 ) else ...

    // ADDED 2015/01/26

    LowerTriangularMatrix LSigma_cond_i = LSigma_cond.i() ;
    // LowerTriangularMatrix LSigma_1 = Cholesky(Sigma_1);
    // LSigma_1_i = LSigma_1.i();

    ColumnVector y_part;
    if (is_q) {
      y_part = rMVN_fn(Mu_cond,LSigma_cond);
    } else {
      ColumnVector y_i = (Y_in.row(i_original)).t();
      y_part = subvector(y_i,item_by_rnorm);
    }
    
    log_cond_norm = log_MVN_fn(y_part,Mu_cond,LSigma_cond_i);

    if (is_q) {
      y_q = tilde_y_i;
      for ( int temp_j = 1,temp_count1 = 0; temp_j<=n_var; temp_j++ ){
        if ( item_by_rnorm(temp_j)==1 ){ 
          y_q(temp_j) = y_part(++temp_count1);
	}
      }
    } // if (is_q)

  } else {

    log_cond_norm = 0;
    if (is_q) { y_q = tilde_y_i;}

  } // if ( item_by_rnorm.sum() > = 1 ) else ..

  return log_cond_norm;
}

void CParam::S1(int iter, Uniform &randUnif, CData &Data, CFeasibilityMap &FM, int n_simul) {
  is_accept(1) = 0;
  for (int i_faulty=1; i_faulty<=Data.n_faulty; i_faulty++) {
    if (FM.useMap) {FM.pmm->iter = iter;}
    int i_original = Data.Faulty2Original[i_faulty-1];
    double g_option_q, g_mode_q;
    int what_type_move;
    ColumnVector s_i = S_Mat.column(i_faulty);
    int tau_q = FM.EvaluateMove(i_original, Data, s_i, iter, what_type_move, g_option_q, g_mode_q, true);
    // calculate g_mode_i & g_option_i
    ColumnVector s_q = CFeasibilityMap::tau_to_s_fn(tau_q, n_var);
    double g_option_i, g_mode_i;
    FM.EvaluateMove(i_original, Data, s_q, iter, what_type_move, g_option_i, g_mode_i, false);
    // calculate g_s_q & g_s_i
    double g_s_q = g_mode_q * g_option_q;
    double g_s_i = g_mode_i * g_option_i;

    // whether drawn by item_by_rnorm and balance edits
    ColumnVector item_by_bal;
    ColumnVector item_by_rnorm;
    if (Data.is_case(i_original,1)) {
      item_by_rnorm = Data.get_item_by_norm_indicator(s_q,item_by_bal);
    } else {
      item_by_rnorm = Data.copy_non_balance_edit(s_q);
    }
    
    //Generate from normal distribution
    ColumnVector mu_z_i = Mu.column(z_in(i_original)) ;
    ColumnVector tilde_y_i = Data.log_D_Observed.row(i_original).t();
    ColumnVector y_q(n_var);
    LowerTriangularMatrix LSigma_i;
    double log_cond_norm_q = calculate_log_cond_norm(Data, i_original, item_by_rnorm, tilde_y_i, y_q, true, LSigma_i, s_q); // MODIFIED 2015/02/16
    
    // Put values from balance edits
    ColumnVector x_q = exp_ColumnVector(y_q);
    if (Data.is_case(i_original,1)) {
      Data.set_balance_edit_values_for_x_q(s_q, x_q, item_by_bal); // CHANGED by Hang, 2014/12/29
    } else {
      Data.update_full_x_for_balance_edit(x_q);
    }

    // Acceptance/Rejection
    if (Data.PassEdits(x_q)) {  // Check constraints
      ColumnVector item_i_by_rnorm;
      if (Data.is_case(i_original,1)) {
        item_i_by_rnorm = Data.get_item_by_norm_indicator(s_i,item_by_bal);
      } else {
        item_i_by_rnorm = Data.copy_non_balance_edit(s_i);
      }
    
    double log_cond_norm_i = calculate_log_cond_norm(Data, i_original, item_i_by_rnorm, tilde_y_i, y_q, false, LSigma_i, s_i); // MODIFIED 2015/02/16

      y_q = log_ColumnVector(x_q) ;
      ColumnVector y_i = (Y_in.row(i_original)).t();
  
      // Calculate acceptance ratio
      ColumnVector y_compact_q = Data.get_compact_vector(y_q);
      ColumnVector y_compact_i = Data.get_compact_vector(y_i);
      double log_full_norm_q = log_MVN_fn(y_compact_q,mu_z_i,LSIGMA_i[z_in(i_original)-1],logdet_and_more(z_in(i_original)));
      double log_full_norm_i = log_MVN_fn(y_compact_i,mu_z_i,LSIGMA_i[z_in(i_original)-1],logdet_and_more(z_in(i_original)));
      
    
      double log_f_y_tilde_q = 0;
      if (FM.useMap && Data.is_case(i_original,2)) {
        log_f_y_tilde_q = FM.Simulate_logUnif_case2(CFeasibilityMap::s_to_tau_fn(s_q),i_original,n_simul,randUnif);
      } else {
        log_f_y_tilde_q = Data.get_logUnif_y_tilde(s_q, CFeasibilityMap::s_to_tau_fn(s_q),i_original); 
      }
      double logNum = log_full_norm_q + log_f_y_tilde_q - log_cond_norm_q - log(g_s_q);
      double logDen = log_full_norm_i + logUnif_y_tilde(i_faulty) - log_cond_norm_i - log(g_s_i);
    
      accept_rate(1) = exp( logNum - logDen ) ;
      if ( randUnif.Next() < accept_rate(1) ){
        Y_in.row(i_original) = y_q.t() ;
        S_Mat.column(i_faulty) = s_q;
        logUnif_y_tilde(i_faulty) = log_f_y_tilde_q;
        is_accept(1)++;
      }
    }
    
  }
  is_accept(1) = is_accept(1) / Data.n_faulty;
}

void CParam::S2_add(Uniform &randUnif,CData &Data) {
  int n_needtoupdate = 0;
  for (int i_faulty=1; i_faulty<=Data.n_faulty; i_faulty++){
    int i_original = Data.Faulty2Original[i_faulty-1];
    ColumnVector item_by_bal;
    ColumnVector s_i = S_Mat.column(i_faulty);
    ColumnVector item_by_rnorm = Data.get_item_by_norm_indicator(s_i,item_by_bal);

    //Generate from normal distribution
    if ( item_by_rnorm.sum() >= 1 ) { // if no random number, other values by balanc edits remain same
    n_needtoupdate++;
    ColumnVector mu_z_i = Mu.column(z_in(i_original)) ;
    ColumnVector tilde_y_i = Data.log_D_Observed.row(i_original).t();
    ColumnVector s_1_compact = Data.get_compact_vector(item_by_rnorm);
    ColumnVector Mu_1i = subvector(mu_z_i,s_1_compact);

    LowerTriangularMatrix LSigma_1i_i;
    ColumnVector y_q(n_var);
    double log_cond_norm_q = calculate_log_cond_norm(Data, i_original, item_by_rnorm, tilde_y_i, y_q, true, LSigma_1i_i, s_i); // MODIFIED 2015/02/16
    ColumnVector y_i = (Y_in.row(i_original)).t() ;
    // ColumnVector y_part_i = subvector(y_i,item_by_rnorm);

    // Put values from balance edits
    ColumnVector x_q = exp_ColumnVector(y_q) ;
    Data.set_balance_edit_values_for_x_q(s_i, x_q, item_by_bal); // CHANGED by Hang, 2014/12/29

    // double log_cond_norm_i = log_MVN_fn(y_part_i,Mu_1i,LSigma_1i_i);
    double log_cond_norm_i = calculate_log_cond_norm(Data, i_original, item_by_rnorm, tilde_y_i, y_q, false, LSigma_1i_i, s_i); // CHANGED 2015/01/27 , // MODIFIED 2015/02/16

    // Acceptance/Rejection
    if (Data.PassEdits(x_q)) {  // Check constraints
    y_q = log_ColumnVector(x_q) ;
    ColumnVector y_compact_q = Data.get_compact_vector(y_q);
    ColumnVector y_compact_i = Data.get_compact_vector(y_i);
    double log_full_norm_q = log_MVN_fn(y_compact_q,mu_z_i,LSIGMA_i[z_in(i_original)-1],logdet_and_more(z_in(i_original)));
    double log_full_norm_i = log_MVN_fn(y_compact_i,mu_z_i,LSIGMA_i[z_in(i_original)-1],logdet_and_more(z_in(i_original)));

    // Calculate acceptance ratio
    double logNum = log_full_norm_q - log_cond_norm_q;
    double logDen = log_full_norm_i - log_cond_norm_i;
    accept_rate(2) = exp( logNum - logDen );

    if (randUnif.Next() < accept_rate(2)){
      Y_in.row(i_original) = y_q.t();
      is_accept(2)++;
    }
    }
    }
  }
  is_accept(2) = is_accept(2) / n_needtoupdate;
}

void CParam::S3_Z_in(CData &Data) {
  
  Data.UpdateCompactMatrix(Y_in_compact,Y_in);
  
  n_z_in = 0 ; // Ver 1.1.0
  
  for (int i_sample=1; i_sample<=n_sample; i_sample++) {
    // calculate pi_tilde_in
    ColumnVector logN_unnorm(K);
    ColumnVector y_compact = (Y_in_compact.row(i_sample)).t() ;

    for (int k=1; k<=K; k++) {
      ColumnVector mu_k = Mu.column(k);
      logN_unnorm(k) = log_MVN_fn(y_compact, mu_k, LSIGMA_i[k-1],logdet_and_more(k));
    }
    double max_logN_unnorm = logN_unnorm.maximum();

    ColumnVector pi_tilde_in_unnorm(K); pi_tilde_in_unnorm = 0.0 ;
    for (int k=1; k<=K; k++){
      pi_tilde_in_unnorm(k) = pi(k) * exp(logN_unnorm(k)-max_logN_unnorm);
    }
    ColumnVector pi_tilde_in = (1.0/pi_tilde_in_unnorm.sum()) * pi_tilde_in_unnorm ;

    z_in(i_sample) = rdiscrete_fn( pi_tilde_in );
    n_z_in(z_in(i_sample)) = n_z_in(z_in(i_sample)) + 1 ; // Ver 1.1.0
  } // for
  
} // void CParam::S3_Z_in

void CParam::S4_Z_out(CData &Data) {
  int count_out = 0,   count_in = 0;

  ColumnVector z_out_large(toolarge_nout) ;
  ColumnVector n_z_out_large(K) ; n_z_out_large = 0 ; // Ver 1.1.0
  
  Matrix Y_out_compact_large(toolarge_nout,n_var_independent);
  Matrix X_aux(n_sample,n_var);

  while ( (count_in < n_sample) && (count_out < toolarge_nout) ) {

    int k = rdiscrete_fn(pi);
    ColumnVector mu_k = Mu.column(k);			// Note that mu_k is Mu.column(k)
    ColumnVector y_compact_i = rMVN_fn( mu_k, LSIGMA[k-1] );

    // ADDED by HANG to check infinity value of x_full_i
    int check_infinity = 0 ;

    if (y_compact_i.maximum()>700){
      check_infinity = 1 ;
      if ( msg_level >= 1 ) {
		Rprintf( "   Warning: x_out from N(Mu_k,Sigma_k) > exp(700). There is no harm for convergence and inference, but the computation may get slower if you see this warning too often, e.g. every iteration\n");
      }
    }

    if (check_infinity==0){

      ColumnVector x_compact_i = exp_ColumnVector(y_compact_i);
      ColumnVector x_full_i(n_var) ;
      Data.UpdateFullVector(x_compact_i,x_full_i);
      Data.update_full_x_for_balance_edit(x_full_i);

      // do not need to consdier min(x_full_i)>0 since x_compact_i=exp(y_compact_i)
      if (Data.PassEdits(x_full_i)) {
        X_aux.row(++count_in) = x_full_i.t();
      } else {
        if ( n_z_out_large(k) < n_z_in(k) ){ // Ver 1.1.0
        Y_out_compact_large.row(++count_out) = y_compact_i.t();
        z_out_large(count_out)=k;
        n_z_out_large(k) = n_z_out_large(k) + 1 ; 
        }
      } // if ... else ...

    } // if (check_infinity==0) : ADDED by HANG

  } // while ( (count_in < n_sample) && (count_out < toolarge_nout) )

  Matrix Y_out_compact = Y_out_compact_large.rows(1,count_out) ;    // cut extra space
  ColumnVector z_out = z_out_large.rows(1,count_out) ;    // cut extra space
  Prob_A = 1.0 * n_sample / (n_sample+count_out) ;

  // calculate n_z and Sum_groupX
  Data.UpdateCompactMatrix(Y_in_compact,Y_in);
  Y_aug_compact = Y_in_compact & Y_out_compact;
  z_aug = z_in & z_out;
  int n_aug = z_aug.nrows() ;

  n_z = 0 ;
  X_bar = Matrix(n_var_independent,K); X_bar= 0.0 ;
  for (int i_aug=1; i_aug<=n_aug; i_aug++) {
    int k = z_aug(i_aug);
    n_z(k) = n_z(k) + 1 ;
    RowVector y_compact_t = Y_aug_compact.row(i_aug);
    ColumnVector y_i_compact = y_compact_t.t() ;
    X_bar.column(k) += y_i_compact;
  }
  for (int k = 1; k <= K; k++) {
    if (n_z(k) > 0) {
      X_bar.column(k) *= (1.0/n_z(k));
    }
  }

}


void CParam::S5_MuSigma(CData &Data, double f_Sigma,double h_Mu) {
  vector<Matrix> X = vector<Matrix>(K);
  int *Counts = new int[K];
  for (int k =0; k < K; k++) {
    if (n_z(k+1) > 0) {
      X[k] = Matrix(n_z(k+1),n_var_independent); X[k] = 0;
      Counts[k] = 0;
    }
  }
  for (int i=1; i<=Y_aug_compact.nrows(); i++){
    int k = z_aug(i);
    X[k-1].row(++Counts[k-1]) = Y_aug_compact.row(i) -  X_bar.column(k).t();
  }
  SymmetricMatrix SqMatrix;
  for (int k=1; k<=K ; k++) {
    // propose Sigma_k_q
    double f_Sigma_tilde_k = f_Sigma + n_z(k);
    double h_k = h_Mu + n_z(k);
    SymmetricMatrix Phi_tilde_k = Phi;
    ColumnVector mu_tilde_k = mu_bar;
    if ( n_z(k) > 0) {
      mu_tilde_k = (h_Mu * mu_bar +  X_bar.column(k) * n_z(k)) / h_k;
      SqMatrix << X[k-1].t() * X[k-1];
      Phi_tilde_k += SqMatrix; //can be further optimized
      ColumnVector xbar_mubar = X_bar.column(k) - mu_bar;
      SqMatrix << (h_Mu*n_z(k)/h_k) * ( xbar_mubar * xbar_mubar.t());
      Phi_tilde_k += SqMatrix ;
    }
    LowerTriangularMatrix LPhi_tilde_k = Cholesky(Phi_tilde_k);
    LowerTriangularMatrix LSigma_k_q = rIW_w_pd_check_fn( f_Sigma_tilde_k, LPhi_tilde_k );

    // propose mu_k_q
    LowerTriangularMatrix LSigma_k_tilde = (1.0/sqrt(h_k)) * LSigma_k_q ;
    ColumnVector mu_k_q = rMVN_fn( mu_tilde_k, LSigma_k_tilde ); // Modified

    // Gibbs update
    Mu.column(k) =	mu_k_q ;
    LSIGMA[k-1] = LSigma_k_q;
    LSIGMA_i[k-1] = LSigma_k_q.i();
    SIGMA[k-1] << LSigma_k_q * LSigma_k_q.t();
    logdet_and_more(k) = -0.5*n_var* LOG_2_PI + logdet(LSIGMA_i[k-1]);

    // S = L * L.t() ;   S.i() = (L.i()).t() * L.i() ;
    Matrix Sigma_k_inv = LSIGMA_i[k-1].t() * LSIGMA_i[k-1];
    for (int i_var=1; i_var<= n_var_independent; i_var++) {
      Sigma_k_inv_ll(k,i_var) = Sigma_k_inv(i_var,i_var);
    }
  }
  delete [] Counts;
}

void CParam::S6_pi() {
  double Sum_n_m = n_z.sum();
  for (int k=1; k<=(K-1); k++) {
    ColumnVector nu_short_q = nu_short;
    double n_z_k = n_z(k) ;
    double one_tilde = 1.0 + n_z_k ;
    Sum_n_m = Sum_n_m - n_z_k ;         // start from Sum_n_m - n_z_1 when k=1
    // ->  Sum_n_m - sum(n_z_1 + n_z_2) when k=2 -> ...
    // Sum_n_m - sum(n_z_1 + ... + n_z_k)) i.e. sum_{g=k+1}^K n_z_g
    double alpha_tilde = alpha + Sum_n_m ;
    
    nu_short_q(k) = rbeta_fn( one_tilde, alpha_tilde );
    
    double Sum_logOneMinusVg = 0.0 ;
    for (int m=1; m<=(K-1); m++) {
      double nu_k_q = nu_short_q(m) ;
      logpi(m) = log(nu_k_q) + Sum_logOneMinusVg ;
      pi(m) = exp( logpi(m) ) ;
      
      double one_minus_V = ( 1.0 - nu_short_q(m) ) ;
      Sum_logOneMinusVg = Sum_logOneMinusVg + log(one_minus_V) ;
    }
    logpi(K) = Sum_logOneMinusVg ;
    pi(K) = exp(Sum_logOneMinusVg) ;
    
    nu_short = nu_short_q;
  }
}

void CParam::S7_alpha(double a_alpha, double b_alpha) {
  double a_alpha_tilde = a_alpha + K - 1.0 ;
  double b_alpha_tilde = b_alpha - logpi(K) ;
  double alpha_q = rgamma_fn( a_alpha_tilde, b_alpha_tilde );
  alpha = alpha_q;
}

void CParam::S8_Phi(double f_Sigma, double a_Phi, double b_Phi) {
  for (int i_var=1; i_var<=n_var_independent; i_var++) {
    ColumnVector Sigma_inv_temp = Sigma_k_inv_ll.column(i_var) ;
    double a_Phi_tilde = a_Phi + 0.5 * K * f_Sigma ;			// (n_var-#. of balance edit+1)
    double b_Phi_tilde = b_Phi + 0.5 * Sigma_inv_temp.sum();
    Phi(i_var,i_var) = rgamma_fn( a_Phi_tilde, b_Phi_tilde );
  }       // end: for (int i_var)
}

void CParam::CalculateInitProbA(CData &Data) {
  // initial value of Y_out_compact and z_out
  int count_out = 0;  int count_in = 0;
  ColumnVector z_out_large(toolarge_nout);
  Matrix Y_out_compact_large(toolarge_nout,n_var_independent);
  while ( (count_in < Data.n_sample) && (count_out < toolarge_nout) ) {
    int k = rdiscrete_fn(pi);
    ColumnVector mu_k = Mu.column(k);             // Note that mu_k is Mu.column(k)
    ColumnVector y_compact_i = rMVN_fn( mu_k,LSIGMA[k-1] );
    ColumnVector x_compact_i = exp_ColumnVector(y_compact_i);

    ColumnVector x_i(Data.n_var) ;
    Data.UpdateFullVector(x_compact_i, x_i);
    Data.update_full_x_for_balance_edit(x_i);

    if (Data.PassEdits(x_i)) {
      count_in++;
    } else {
      count_out++;
      Y_out_compact_large.row(count_out) = y_compact_i.t();
      z_out_large(count_out)=k;
    }
  }
  Matrix Y_out_compact = Y_out_compact_large.rows(1,count_out) ;    // cut extra space
  ColumnVector z_out = z_out_large.rows(1,count_out) ;    // cut extra space
  // calculate n_z and Sum_groupX
  Matrix Y_aug_compact = Y_in_compact & Y_out_compact;
  ColumnVector z_aug = z_in & z_out;
  int n_out = z_out.nrows();
  Prob_A =  (1.0 * Data.n_sample / (Data.n_sample+n_out));
}

ColumnVector CParam::GetComponentCounts() {
  // ColumnVector n_z_in(K) ; // Ver 1.1.0
  n_z_in = 0 ; // Ver 1.1.0
  for (int i=1; i<=n_sample ; i++ ){
    n_z_in(z_in(i))++;
  }
  return n_z_in;
}
