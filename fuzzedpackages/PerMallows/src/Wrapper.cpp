//c++ naming conventions
/*
 they all must return void
 parameters by reference
 include the R.h
 */
#include "Cayley.h"
#include "Kendall.h"
#include "Hamming.h"
#include "Ulam_disk.h"
#include "Ulam.h"
#include "Newton_raphson.h"
#include "Exponential_model.h"
#include "Lap.h"

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

/* dist:id
 0: cayley
 1: kendall
 2: hamming
 3: ulam
 4: ULAM_DISK_DISTANCE
 */
extern "C" {
    
    void count_permus_with_at_least_k_unfixed_points ( int*n, int *k, double * res){
        GetRNGstate();
        Generic gen;
        *res = (double)gen.count_permus_with_at_least_k_unfixed_points(*n, *k);
        PutRNGstate();
    }
    
    void random_permutation(int *n, int *sigma){
      GetRNGstate();
        Generic gen;
        gen.generate_random_permutation(*n, FIRST_ITEM, sigma);
        PutRNGstate();
    }
    
    void compute_distance (int*dist_id, int*n,int*s1, int*s2, int*d) {
      GetRNGstate();
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(*dist_id, *n);
        *d=exp_mod->distance(s1,s2);
        delete exp_mod;
        PutRNGstate();
    }
    
    void count_permus_at_dist (int*dist_id, int*n, int*d, double*res) {
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(*dist_id, *n);
        if (*d < 0 || *d > exp_mod->maximum_distance()) *res = 0;
        else *res= (double) exp_mod->num_permus_at_distance(*d);
        delete exp_mod;
        PutRNGstate();
    }
    
    void probability(int*dist_id, int*n, int*sigma, int*sigma_0, double*theta, double*prob){
      GetRNGstate();
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(*dist_id, *n);
        *prob = exp_mod->probability(sigma, sigma_0, theta);
        delete exp_mod;
        PutRNGstate();
    }
    
    void get_altern_repre_for_permu(int*dist_id, int*n, int*sigma, int*vec){
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(*dist_id, *n);
        exp_mod->perm2dist_decomp_vector(sigma, vec);
        delete exp_mod;
    }
    
    void get_permu_given_altern_repre(int*dist_id, int*n, int*vec, int*sigma){
      GetRNGstate();
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(*dist_id, *n);
        exp_mod->dist_decomp_vector2perm(vec, sigma);
        delete exp_mod;
        PutRNGstate();
    }
    
    void expectation(int*dist_id, int* model, int *n, double * theta, double* h_avg){
      GetRNGstate();
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(*dist_id, *n);
        if ( *model == MALLOWS_MODEL )
          h_avg[0] = exp_mod->expectation(theta[ 0 ]);
        else exp_mod->expectation(theta, h_avg);
        delete exp_mod;
        PutRNGstate();
    }
    
    void save_counts_to_files ( int*n){
      GetRNGstate();
        Ulam_disk * ul = new Ulam_disk(*n);
        ul->save_counts_to_file_bin();
        delete ul;
        PutRNGstate();
    }
    
    void marginals ( int*n, int *dist_id, int *h, double * theta, double *res){
      GetRNGstate();
        Hamming *ham = new Hamming(*n);
        *res = ham->compute_marginal(h, theta);
        delete ham;
        PutRNGstate();
    }
    
    
    /***************       SAMPLING (GENERATING)      *******************/
    
    SEXP get_random_sample_at_dist_d(SEXP dist_id_var, SEXP n_var, SEXP m_var, SEXP d_var){
      GetRNGstate();
        
        SEXP Rval;
        int n = INTEGER_VALUE(n_var);
        int m  = INTEGER_VALUE(m_var);
        int d = INTEGER_VALUE(d_var);
        int dist_id = INTEGER_VALUE(dist_id_var);
        
        int**sample=new int*[m];
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(dist_id, n);
        exp_mod->random_sample_at_dist(d, m, sample);
        
        PROTECT(Rval = allocMatrix(REALSXP, m, n));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                REAL(Rval)[i + m * j] = sample[i][j];
        UNPROTECT(1);
        for (int i = 0 ; i < m ;  i ++ ) delete [] sample[ i ];
        delete [] sample;
        delete exp_mod;
        PutRNGstate();
        return Rval;
    }
    
    SEXP distances_sampling(SEXP dist_id_var, SEXP n_var, SEXP m_var, SEXP theta_var){
      GetRNGstate();
        int     m  =    INTEGER_VALUE(m_var);
        int     n  =    INTEGER_VALUE(n_var);
        int dist_id  =  INTEGER_VALUE(dist_id_var);
        double theta =  NUMERIC_VALUE(theta_var);
        int  **sample=  new int*[m];
        SEXP Rval;
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(dist_id, n);
        exp_mod->distances_sampling(m,theta,sample);

        PROTECT(Rval = allocMatrix(REALSXP, m, n));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                REAL(Rval)[i + m * j] = sample[i][j];
        UNPROTECT(1);
        for (int i = 0 ; i < m ;  i ++ ) delete [] sample[ i ];
        delete [] sample;
        delete exp_mod;
        PutRNGstate();
        return Rval;
    }
    
    
    
    //sam <- .Call("sampling_multi_gibbs_cayley",dist_id,  permu.length, num.permus, theta, 0, algorithm_id)
    SEXP sampling_multi_gibbs_cayley(SEXP dist_id_var,SEXP n_var, SEXP m_var,
                                     SEXP theta_var, SEXP model_var,SEXP method_var){
         GetRNGstate();
        //method_var == 1 multistage
        //method_var == 2 gibbs
        int m  = INTEGER_VALUE(m_var);
        int n  = INTEGER_VALUE(n_var);
        int model    = INTEGER_VALUE(model_var);//GMM MM
        int method ;
        method = INTEGER_VALUE(method_var);//gibbs...
        int dist_id  = INTEGER_VALUE(dist_id_var);
        PROTECT(theta_var = AS_NUMERIC(theta_var));
        double*theta=NUMERIC_POINTER(theta_var);
        UNPROTECT(1);
        //return(R_NilValue);
        
        int**sample=new int*[m];
        SEXP Rval;
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(dist_id, n);
        if(method == 1)
            exp_mod->multistage_sampling(m,theta, sample);
        else
            exp_mod->gibbs_sampling(m, theta, model,sample);
            
        PROTECT(Rval = allocMatrix(REALSXP, m, n));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                REAL(Rval)[i + m * j] = sample[i][j];
        UNPROTECT(1);
        for (int i = 0 ; i < m ;  i ++ ) delete [] sample[ i ];
        delete [] sample;
        delete exp_mod;
        PutRNGstate();
        return Rval;
    }
    
    /***************       LEARNING (FITTING)     *******************/
    
    SEXP consensus (SEXP dist_id_var, SEXP samples, SEXP mm_var, SEXP estim_var, SEXP sigma_0_ini_var){
        //mm_var==0 mm
        //mm_var==1 gmm
        //estim_var == 0 exact
        
        GetRNGstate();
        
        SEXP Rdim = getAttrib(samples, R_DimSymbol);
        
        //PROTECT(sigma_0_ini_var = AS_INTEGER(sigma_0_ini_var));
        //int*sigma_0_ini=INTEGER_POINTER(sigma_0_ini_var);
        //UNPROTECT(1);
        PROTECT(samples = AS_INTEGER(samples));
        int m = INTEGER(Rdim)[0];
        int n = INTEGER(Rdim)[1];
        int model = INTEGER_VALUE(mm_var);
        int estim = INTEGER_VALUE(estim_var); // 0 exact;
        int dist_id = INTEGER_VALUE(dist_id_var);
        int**c_samples=new int*[m];
        int*sigma_0=new int[n];
        
        for(int i=0;i<m;i++){
            c_samples[i]=new int[n];
            for(int j=0;j<n;j++){
                c_samples[i][j]=INTEGER(samples)[i + m * j];
            }
        }
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(dist_id, n);
        if (estim == 0)
            exp_mod->estimate_consensus_exact( model, m,c_samples, NULL, sigma_0);
        else
            exp_mod->estimate_consensus_approx(model, m,c_samples, sigma_0);

        SEXP myint;
        int *p_myint;
        int len = n;
        // Allocating storage space:
        PROTECT(myint = NEW_INTEGER(len));
        p_myint = INTEGER_POINTER(myint);
        for(int i=0;i<n;i++) p_myint[i] = sigma_0[i];
        UNPROTECT(2);
        delete exp_mod;
        for(int i=0;i<m;i++) delete [] c_samples[ i ];
        delete [] c_samples;
        delete [] sigma_0;
        PutRNGstate();
        return myint;
    }
    
    SEXP estimate_theta(SEXP dist_id_var, SEXP n_var, SEXP m_var, SEXP sigma_0_var, SEXP samples_var, SEXP model_var){
        //model (model_var)  MM=0 ;;;; 1 =GMM
        
        GetRNGstate();
        SEXP Rval;
        int i;
        PROTECT(Rval = allocVector(INTSXP, 1));
        for (i = 0; i < 1; i++)
            INTEGER(Rval)[i ] = i;
        UNPROTECT(1);
        int n = INTEGER_VALUE(n_var);
        int m = INTEGER_VALUE(m_var);
        int dist_id = INTEGER_VALUE(dist_id_var);
        int model = INTEGER_VALUE(model_var);
        int**samples=new int*[m];
        PROTECT(samples_var = AS_INTEGER(samples_var));
        for(int i=0;i<m;i++){
            samples[i]=new int[n];
            for(int j=0;j<n;j++)
                samples[i][j]=INTEGER(samples_var)[i + m * j];
        }
        double*theta= new double[n];
        PROTECT(sigma_0_var = AS_INTEGER(sigma_0_var));
        int*sigma_0 = INTEGER_POINTER(sigma_0_var);
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(dist_id, n);
        exp_mod->estimate_theta(m, sigma_0, samples, model, theta);

        UNPROTECT(2);
        PROTECT(Rval = allocVector(REALSXP, n));
        for (int i = 0; i < n ; i++) REAL(Rval)[i] = theta[i];
        UNPROTECT(1);
        delete [] theta;
        for (int i = 0 ; i < m ;  i ++ ) delete [] samples[ i ];
        delete [] samples;
        delete exp_mod;
        PutRNGstate();
        return Rval;
    }
    //Call("get_log_likelihood",dist_id,    permu.length, num.permus, sigma_0, theta, samples, model)
    
    SEXP get_log_likelihood(SEXP dist_id_var, SEXP n_var, SEXP m_var, SEXP sigma_0_var, SEXP theta_var,
                            SEXP samples_var, SEXP model_var){
        
        GetRNGstate();
        //model (model_var)  MM=0 ;;;; 1 =GMM
        SEXP Rval;
        int i;
        double likeli = -1.0;
        PROTECT(Rval = allocVector(INTSXP, 1));
        for (i = 0; i < 1; i++)
            INTEGER(Rval)[i ] = i;
        UNPROTECT(1);
        int n         = INTEGER_VALUE( n_var );
        int m         = INTEGER_VALUE( m_var );
        int dist_id   = INTEGER_VALUE( dist_id_var );
        int model     = INTEGER_VALUE( model_var );
        
        int**samples  = new int*[ m ];
        PROTECT(samples_var = AS_INTEGER( samples_var ));
        for(int i=0;i<m;i++){
            samples[i]=new int[n];
            for(int j=0;j<n;j++)
                samples[ i ][ j ] = INTEGER(samples_var)[i + m * j];
        }
        int last_theta = n;
        if( dist_id == CAYLEY_DISTANCE || dist_id == KENDALL_DISTANCE ) last_theta = n-1;
        
        /*double * theta = new double[ n ];
        PROTECT(theta_var = AS_NUMERIC(samples_var));
        for(int i = 0 ; i < last_theta ; i ++)
            theta[ i ] = REAL(theta_var)[ i ];
        */
        int * sigma_0 = new int[ n ];
        PROTECT(sigma_0_var = AS_INTEGER(sigma_0_var));
        for(int i = 0 ; i < n ; i ++)
            sigma_0[ i ] = INTEGER(sigma_0_var)[ i ];
        
        Generic gen;
        Exponential_model * exp_mod = gen.new_instance(dist_id, n);
        likeli = exp_mod->get_likelihood(m, samples, model, sigma_0);

        UNPROTECT(2);
        //cout<<"liekli .. "<<likeli<<endl;
        PROTECT(Rval = allocVector(REALSXP, 1));//vector len 1
        REAL(Rval)[0] = likeli;
        //for (int i = 0; i < n - 1; i++) REAL(Rval)[i] = theta[i];
        UNPROTECT(1);
        delete exp_mod;
        for (int i = 0 ; i < m ;  i ++ ) delete [] samples[ i ];
        delete [] samples;
        delete [] sigma_0;
        PutRNGstate();
        return Rval;
    }
    
} // extern "C"
