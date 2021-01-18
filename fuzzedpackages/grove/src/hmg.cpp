#include <RcppArmadillo.h>
#include "helpers.h"   
#include "hmg.h"  

using namespace Rcpp;
using namespace arma;
using namespace std;


HMG::HMG(   mat D, 
            mat X, 
            vec p, 
            vec tau_par, 
            vec eta_par, 
            vec gamma_par,  
            vec init_state,
            double nu, 
            double sigma0, 
            double alpha, 
            double beta,  
            int transition_mode ):
  X(X), p(p), tau_par(tau_par), eta_par(eta_par), gamma_par(gamma_par), 
  alpha(alpha), beta(beta), transition_mode(transition_mode)
{
    n_factors = p.n_elem;
    tot_states = (int)pow(2, n_factors);
    initial_state = init_init_state(init_state);
    n_tot = D.n_rows;
    J = log2(D.n_cols+1);
    a_0 = nu + 1;
    b_0 = pow(sigma0, 2) * nu;
    w = init_data(D);
    FullLambda = FullPrecision();
    Phi = init_marg();
    Psi = update_marg();
    PostTrans = post_trans();
    PostStates = post_state();  
    
    PriorTrans = prior_trans();
    PriorStates = prior_state();                                  
}


// organize initial state of the product space
vec HMG::init_init_state(vec input_init_state)
{
  vec output(tot_states); output.fill(1);
  for( int s=0; s<tot_states; s++)
  {
    for(int t=0; t<n_factors; t++)
    {
      if( ((s >> t) & 1) == 1 )
        output(s) *= input_init_state(t);
      else if( ((s >> t) & 1) == 0 )
        output(s) *= (1 - input_init_state(t));
    }
    
  }
  return output;
}  

// organize data on a tree
std::vector<mat> HMG::init_data(mat D)
{
  std::vector<mat> data(J);
  int j, k, it=0;
  for(j=0; j<J; j++)
  {
    mat M(n_tot, (int)pow(2,j));
    for(k=0; k < (int)pow(2,j); k++)
    {
      M.col(k) = D.col(it);
      it++;
    }
    data[j] = M;      
  }
  return data;    
}  

// compute marginal likelihood of every node for every state
std::vector<mat> HMG::init_marg()
{
  int j, s;
  std::vector<mat> data(J);

  for(j=0; j<J; j++)
  {
    mat M(tot_states, (int)pow(2,j));
    for(int k=0; k<(int)pow(2,j); k++)
    {
      for(s=0; s < tot_states; s++)
        M(s,k) = MargLike(j, k, s);
    }
    data[j] = M;      
  }
  return data;
}

// vector of active factors from integer
uvec HMG::active_columns(int s)       
{
  vec temp = zeros<vec>(sum(p));
  vec par_p = cumsum(p);
  
  if( (s & 1) == 1 )
    temp(0) = 1;
    
  for(int i=1; i<n_factors; i++)
  {
    if( (s>>i & 1) == 1 )
    {
      if(p(i)==1)
      {
        temp(par_p(i)-1) = 1;
      }
      else
      {
        temp.subvec(par_p(i-1), (par_p(i)-1)  ) = ones<vec>(p(i));
      }
      
    }
      
  }
  uvec output = find(temp>0);
  return output;
}

// extract submatrix according to active factors
mat HMG::DesignMatrix(int s)
{
  // if(s == 0)
  //   cout << "error - DesignMatrix s==0" << endl;    
  uvec ac = active_columns(s);  
  return(X.cols(ac));    
}

// compute the full precision matrix
mat HMG::FullPrecision()
{
  mat output(J, sum(p));
  vec par_p = cumsum(p);
  for(int j=0; j<J; j++)
  {
    output(j,0) = 1.0 / ( tau_par(0)*pow(2,-alpha*j) );
    
    for(int i=1; i<n_factors; i++)
    {
      rowvec temp(p(i)); 
      temp.fill( 1/( tau_par(i)*pow(2,-alpha*j) )); 
      output.row(j).cols(par_p(i-1), par_p(i)-1) = temp;
    }
  }
  return output;
}


// extract submatrix for computing marginal likelihood 
mat HMG::Lambda(int j, int s)
{  
  uvec ac = active_columns(s);  
  uvec scale(1); scale(0)=(unsigned int)j;
  return diagmat( FullLambda.submat(scale,ac) );
}

// evaluate marginal likelihood at given node and given state
double HMG::MargLike(int j, int k, int s)     
{

  vec y = w.at(j).col(k);
  double a_n = a_0 + n_tot/2.0;
  double b_n = b_0 + 0.5 * dot(y, y)  ;

  double output = -n_tot/(2.0) * log2pi  +  a_0 * log(b_0) + lgamma(a_n) - lgamma(a_0);
  if(s == 0)
  {
    output -= a_n * log(b_n);
  }
  else
  {
    mat Lambda_0 = Lambda(j, s);
    mat X = DesignMatrix(s);
    mat Lambda_n = X.t() * X + Lambda_0;
    vec mu_n = inv_sympd( Lambda_n ) * ( X.t() * y );
    b_n -= as_scalar(mu_n.t() * Lambda_n * mu_n)/2.0 ;
    double logDetLambda_n;
    double sign;
    log_det(logDetLambda_n, sign, Lambda_n);
    double logDetLambda_0 = sum( log(diagvec(Lambda_0) ) );
    
    output += 0.5*( logDetLambda_0 - logDetLambda_n ) - a_n * log(b_n);

  }
  return( output );
}

// compute prior transition probability
double HMG::prior_trans_elem(int j,  int s, int t)   
{
  double out = 1.0;  
  for(int i=0; (i<n_factors && out>0) ; i++)
  {
    vec temp1; temp1 << 1 <<  ( eta_par(i) * pow(2,-beta*j) );
    vec temp2; temp2 << (1-gamma_par(i)) << ( eta_par(i) * pow(2,-beta*j) );
  
    if (transition_mode == 1) { // Markov transition
      if( ((s >> i) & 1)==0 && ((t >> i) & 1)==0)
        out *= ( 1 -  min(temp1) );   
      else if(((s >> i) & 1)==0 && ((t >> i) & 1)==1)
        out *=  min(temp1);  
      else if(((s >> i) & 1)==1 && ((t >> i) & 1)==0)
        out *= ( 1 - gamma_par(i) - min(temp2) ); 
      else if(((s >> i) & 1)==1 && ((t >> i) & 1)==1)
        out *= ( gamma_par(i) + min(temp2) );
    }
    
    if (transition_mode == 0) { // Independent transition
      
      if( ((t >> i) & 1)==0)
        out *= ( 1 -  min(temp1) );   
      else if(((t >> i) & 1)==1)
        out *=  min(temp1);  
    }
  }      
  return out;
  
}


// compute marginal likelihood recursively
std::vector<mat> HMG::update_marg() 
{
  std::vector<mat> output(J);    
  for (int j=0; j<J; j++) 
  { 
    mat M(tot_states, (int)pow(2,j));    
    output[j] = M;
  }
      
  for (int j=(J-1); j>=0; j--) 
  { //start from the leaves of the tree        
    for (int k = 0; k < (int)pow(2,j) ; k++) 
    {
      for(int s=0; s<tot_states; s++)
      {
        bool first_it = true;
        for(int t = 0; t<tot_states ; t++)
        {
          double rho  = prior_trans_elem(j, s, t);               
          if(rho > 0.0)
          {
            if(j == (J-1))
            {   
              if(first_it)
              {
                output.at(j).at(s,k) = log(rho) + Phi.at(j).at(t,k); 
                first_it=false;
              }
              else
                output.at(j).at(s,k) = log_exp_x_plus_exp_y( output.at(j).at(s,k), log(rho) + Phi.at(j).at(t,k) );
            }
            else
            {                          
              if(first_it)
              {
                output.at(j).at(s,k) = log(rho)  + Phi.at(j).at(t,k) 
                  + output.at(j+1).at(t,2*k) + output.at(j+1).at(t,2*k+1);
                first_it =false;
              }
              else
                output.at(j).at(s,k) = log_exp_x_plus_exp_y( output.at(j).at(s,k), log(rho)
                  + Phi.at(j).at(t,k) + output.at(j+1).at(t,2*k) + output.at(j+1).at(t,2*k+1) );
            }
          }                         
        }          
      }
    }
  }
  return output;
}
  
// compute posterior transition probabilities  
std::vector<cube> HMG::post_trans()
{
  std::vector<cube> output(J);    
  for (int j=0; j<J; j++) 
  { 
    cube M(tot_states, tot_states, (int)pow(2,j));    
    output[j] = M;
  }
      
  for (int j=(J-1); j>=0; j--) 
  { //start from the leaves of the tree      
    for (int k = 0; k < (int)pow(2,j) ; k++) 
    {
          for(int s=0; s<tot_states; s++)
          {
              for(int t = 0; t<tot_states ; t++)
              {
                output.at(j).at(s, t, k) = post_trans_elem(j, k, s, t);
              }                
          }     
    }
  }
  return output;
}

// compute prior transition probabilities  
std::vector<mat> HMG::prior_trans()
{
  std::vector<mat> output(J);    
  for (int j=0; j<J; j++) 
  { 
    mat M(tot_states, tot_states);    
    output[j] = M;
  }
      
  for (int j=(J-1); j>=0; j--) 
  { //start from the leaves of the tree      

          for(int s=0; s<tot_states; s++)
          {
              for(int t = 0; t<tot_states ; t++)
              {
                output.at(j).at(s, t) = prior_trans_elem(j, s, t);
              }                
          }     

  }
  return output;
}

// compute posterior transition probability for a given node
double HMG::post_trans_elem(int j, int k, int s, int t)
{
  double rho  =  prior_trans_elem(j, s, t);
  double log_numerator;

  if(rho>0.0 && rho<1.0)
  {
      
      if(j==(J-1))
        log_numerator = log(rho)  + Phi.at(j).at(t,k); 
      else
        log_numerator = log(rho) + Phi.at(j).at(t,k) + Psi.at(j+1).at(t,2*k) + Psi.at(j+1).at(t,2*k+1);    
  
      return exp( log_numerator - Psi.at(j).at(s,k) );
  }
  else
      return rho;

}

// compute posterior state probabilities
std::vector<mat> HMG::post_state()
{
  std::vector<mat> output(J);    
  
  mat M(tot_states,1);
  M.col(0) = PostTrans.at(0).slice(0).t() * initial_state;
  output[0] = M;
  
  for (int j=1; j<J; j++) 
  { 
    mat M(tot_states, (int)pow(2,j));    
    for (int k = 0; k < (int)pow(2,j) ; k++) 
    {
      M.col(k) = PostTrans.at(j).slice(k).t() * output.at(j-1).col((int)k/2);
    }
    output[j] = M;
  }
  
  return output;
  
}    

// compute prior marginal alternative probabilities
std::vector<mat> HMG::prior_state()
{
  std::vector<mat> output(J);    
  
  mat M(tot_states,1);
  M.col(0) = PriorTrans.at(0).t() * initial_state;
  output[0] = M;
  
  for (int j=1; j<J; j++) 
  { 
    mat M(tot_states, (int)pow(2,j));    
    for (int k = 0; k < (int)pow(2,j) ; k++) 
    {
      M.col(k) = PriorTrans.at(j).t() * output.at(j-1).col((int)k/2);
    }
    output[j] = M;
  }
  
  return output;
  
}    

// get marginal likelihood
double HMG::get_marginal_likelihood()
{
  mat like = Psi.at(0).col(0).t() * initial_state;
  return like(0,0); 
}

// get posterior null for each factor
vec HMG::get_posterior_null()
{
  vec post_null(n_factors); post_null.fill(0.0);    
  for(int g=0; g<n_factors; g++)
  {
    for(int i=0; i<tot_states; i++)
    {
      if(initial_state(i)>0)
        post_null(g) += posterior_null(0,0,i,g)*initial_state(i);
    }
  }
  return post_null;
}

// get posterior null for each node and each factor
double HMG::posterior_null(int j, int k, int s, int g)
{
  if( j == (J-1) )        //    if( j == (J-1) )
      return 1.0;
  else
  {
    double temp = 0;
    for(int t=0; t<tot_states; t++)
    {
      if( (PostTrans.at(j).at(s,t,k)>0) && ( (t>>g & 1 )==0 )  )
        temp += PostTrans.at(j).at(s,t,k)  * posterior_null(j+1, 2*k, t, g) * posterior_null(j+1, 2*k+1, t, g);
    }
    return temp;
  }           
}

// get prior null probability for each factor
vec HMG::get_prior_null(vec init_state)
{
  vec output(n_factors);
  for(int s=0; s<n_factors; s++)
  {
    if(eta_par(s) >= 1.0)
      output(s) = 0;         
    else
    {
      double temp=0; vec temp2(2);
      for(int j=1; j<J; j++)
      {
        temp += pow(2,j)*log( 1 - eta_par(s)*pow(2,-beta*j) );
      }
      temp2 << 1.0 - gamma_par(s) << eta_par(s) ;
      output(s) = ( init_state(s)*( 1- gamma_par(s) - min(temp2)  )   + (1-init_state(s))*( 1 - eta_par(s)) ) * exp(temp);

    }      

  }
  return output;
}

// get marginal posterior probability of hidden states
std::vector<NumericMatrix> HMG::get_post_states()
{
  std::vector<NumericMatrix> output(J);
  for(int j=0; j<J; j++)
  {
    NumericMatrix M_new = as<NumericMatrix>(wrap(PostStates.at(j)));
    output[j] = M_new;
  }
  return output;
}

// get marginal prior probability of hidden states
std::vector<NumericMatrix> HMG::get_prior_states()
{
  std::vector<NumericMatrix> output(J);
  for(int j=0; j<J; j++)
  {
    NumericMatrix M_new = as<NumericMatrix>(wrap(PriorStates.at(j)));
    output[j] = M_new;
  }
  return output;
}


std::vector<NumericMatrix> HMG::Sample_States(int n_samples) {
  std::vector<NumericMatrix> output(J);    
  NumericMatrix v(n_samples, 1);
  vec probs = PostTrans.at(0).slice(0).t() * initial_state;    
  for (int i=0; i<n_samples; i++) {
    v(i,0) = sampling(tot_states, probs);
  }
  output[0] = v;
  int it = 0;
  for(int j=1; j<J; j++) {
    NumericMatrix v(n_samples, pow(2,j)); 
    for(int k=0; k<pow(2,j); k++) {
      NumericMatrix temp = output.at(j-1);
      for(int i=0; i<n_samples; i++) {
        probs = vectorise(PostTrans.at(j).slice(k).row(temp(i,k/2)));
        v(i,k) = sampling(tot_states, probs);
      }
    } 
    it++;
    output[it] = v;
  }    
  return output;    
}  


std::vector<mat> HMG::Count_Sample_States(std::vector<NumericMatrix> StatesSample)
{
  std::vector<mat> output(J);
  int n_samples = StatesSample.at(0).nrow();
  for(unsigned int j=0; j < (unsigned int)J; j++)
  {
    mat v(tot_states, pow(2,j)); v.zeros();
    NumericMatrix temp = StatesSample.at(j);
    for(int i = 0; i < n_samples; i++)
    {
      for(int k=0; k<pow(2,j); k++)
        v(temp(i,k),k)++;
    }
    output[j] = v;
  }
  return output;
}


mat HMG::post_sample_single_coeff(int j, int k, int s, int n)  
{
  mat mean_samples = zeros<mat>(n, sum(p));
  vec y = w.at(j).col(k);
  double a_n = a_0 + n_tot/2.0;
  double b_n = b_0 + 0.5 * dot(y, y)  ;    
  if(s==0)
    return mean_samples;
  else
  {
    uvec ac = active_columns(s); 
    vec z(ac.n_rows); z.fill(0);  
    mat Lambda_0 = Lambda(j, s);
    mat X = DesignMatrix(s);
    mat Lambda_n = X.t() * X + Lambda_0;
    mat Sigma_n = inv_sympd( Lambda_n );
    vec mu_n = Sigma_n * ( X.t() * y );
    vec mean_n(sum(p)); mean_n.fill(0); mean_n(ac) = mu_n;
    b_n -= as_scalar(mu_n.t() * Lambda_n * mu_n)/2.0 ;
    vec sigma = pow( rgamma(n, a_n, 1.0/b_n)  , -0.5);
    mean_samples.cols(ac) = mvrnormArma(n, z, Sigma_n);   
    for(int i=0; i<n; i++)
      mean_samples.row(i) = mean_n.t() + mean_samples.row(i)  * as_scalar(sigma(i)); 
    return mean_samples;
  }
    
}

std::vector<cube>  HMG::post_sample_coeff(int n_samp)
{
  std::vector<NumericMatrix> StatesSample = Sample_States(n_samp);
  std::vector<mat> Counts = Count_Sample_States(StatesSample);
  std::vector<cube> output(J);
  for(int j=0; j<J; j++)
  {
    cube NodeSamples(n_samp, sum(p), pow(2,j));
    for(int k=0; k<pow(2,j); k++)
    {
      for(int s=0; s<tot_states; s++)
      {
        int n_local = as_scalar( Counts.at(j).at(s,k) );
        if( n_local > 0)
        {
          mat temp = Rcpp::as<arma::mat>(StatesSample.at(j)); 
          uvec idx = find( temp.col(k) == s );
          NodeSamples.slice(k).rows(idx) = post_sample_single_coeff(j, k, s, n_local);
        }
      }
    }
    output[j] = NodeSamples;
  }
  return output;
}


cube HMG::get_post_sample_coeff(std::vector<cube> input)
{
  cube tt = input.at(0);  mat tt2 = tt.slice(0);
  int n_samp = tt2.n_rows;
  cube output(sum(p), pow(2,J)-1, n_samp);
  for(int i=0; i<n_samp; i++)
  {
    cube Q = input.at(0);
    mat aaa = Q.tube( i, 0, size(1, sum(p)) );
    output.slice(i).col(0) = aaa.t();
    for(int j=1; j<J; j++)
    {
      Q = input.at(j);
      mat aaa = Q.tube( i, 0, size(1, sum(p)) );
      output.slice(i).cols(pow(2,j)-1, pow(2,j+1)-2) = aaa;
    }
  }
  return output;
}
