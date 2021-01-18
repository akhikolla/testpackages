#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

#include <string>
#include <iostream>
#include <fstream>
using namespace Rcpp;



#define array_size 100000
#define K_num 3



int GetN(double t)
{

    return (int)(t/0.5+1);

}

double Min(const double t1, const double t2)
{
    if(t1<t2) return t1;
    else return t2;
}


void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta)
{
    int p = XX->size1;
    int q = XX->size2;

    int i=0,j=0;
    double temp=0.0;

    for(i=0;i<p;i++)
    {
        temp=0;
        for(j=0;j<q;j++)  temp+=gsl_matrix_get(XX,i,j)*gsl_vector_get(X,j);
        gsl_vector_set(beta,i,temp);
    }

}



//' Simulation of ordinal longitudinal outcome and competing risks data
//' Currently, only the simulation in Li et al(2010) is implemented. 
//' @title Data simulation of ordinal outcomes and competing risks
//' @param k_val The number of subjects in study.
//' @param p1_val  The dimension of fixed effects in longitudinal measurements. 
//' @param p1a_val  The dimension of random effects in longitudinal measurements. 
//' @param p2_val  The dimension of fixed effects in competing risks failure time data.
//' @param g_val  The number of type of failure in competing risks data.
//' @param truebeta True values for beta, the longitudinal coefficients. 
//' @param truetheta True values for theta, subset of the non-proportional odds longitudinal coefficients.
//' @param truegamma True values for gamma, the survival coefficients. 
//' @param randeffect True values for random effects in longitudinal  and competing risks parts,namely in the order of \eqn{\sigma_b},\eqn{\nu_2},\eqn{\sigma_u}.
//' @param yfn  Filename of genereated Y matrix for longitudinal measurements in long format.
//' @param cfn  Filename of genereated C matrix for competing risks failure time data.
//' @param mfn  Filename of genereated M vector to indicate the number of longitudinal measurements per subject.
//' @return Files with names yfn, cfn and mfn.
//' \tabular{ll}{
//' \code{censoring_rate} \tab Censoring rate of the survival data. \cr
//' \code{rate1} \tab Censoring rate of competing risk 1. \cr
//' \code{rate2} \tab Censoring rate of competing risk 2. \cr
//' \code{yfn} \tab Filename of genereated Y matrix for longitudinal measurements. \cr
//' \code{cfn} \tab Filename of genereated C matrix for competing risks failure time data. \cr
//' \code{mfn} \tab  Filename of genereated M vector to indicate the number of longitudinal measurements per subject.
//'  } 
//' @examples
//' # A toy example testint data generations
//' require(JMcmprsk)
//' set.seed(123)
//' yfn=tempfile(pattern = "", fileext = ".txt")
//' cfn=tempfile(pattern = "", fileext = ".txt")
//' mfn=tempfile(pattern = "", fileext = ".txt")
//'  k_val=50;p1_val=3;p1a_val=1; p2_val=2;g_val=2;
//' truebeta=c(-1,1.5,0.8);truetheta=c(-0.5,1);truegamma=c(0.8,-1,0.5,-1); randeffect=c(1,0.5,0.5);
//' #writing files
//' SimDataO(k_val, p1_val, p1a_val, p2_val, g_val,
//'       truebeta, truetheta, truegamma, randeffect, yfn,  cfn,  mfn)
//' \donttest{
//' jmo_0(p=3,s=1, yfn,cfn,mfn,point=10,do.trace = TRUE)
//'}
//' @references
//' \itemize{
//' \item Ning Li,Robert M. Elashoff,Gang Li and Jeffrey Saver. "Joint modeling of longitudinal ordinal data and competing risks survival times and analysis of the NINDS rt-PA stroke trial." Statistics in medicine 29.5 (2010): 546-557.
//' }
//' @seealso \code{\link{SimDataC}}
//' @export
// [[Rcpp::export]]
Rcpp::List   SimDataO(SEXP k_val,SEXP p1_val,SEXP p1a_val,SEXP p2_val, SEXP g_val, SEXP truebeta, SEXP truetheta,SEXP truegamma,SEXP randeffect,SEXP yfn, SEXP cfn, SEXP mfn)
{

    
    double degree_v=5.0;
	
	size_t k=as<int> (k_val);
	size_t p1=as<int> (p1_val);
	size_t p2=as<int> (p2_val);
	size_t g=as<int> (g_val);
	size_t p1a=as<int> (p1a_val);
	size_t bq=1;
	size_t n=0,n1=0;
	size_t i=0,j=0;
  
	std::string yfile=as<std::string>(yfn);
	std::string cfile=as<std::string>(cfn);
	std::string mfile=as<std::string>(mfn);

  const std::vector <double> tbeta_val=as<std::vector<double> >(truebeta);
   const std::vector <double> ttheta_val=as<std::vector<double> >(truetheta);
  const std::vector <double> tgamma_val=as<std::vector<double> >(truegamma);
  const std::vector <double> reffect_val=as<std::vector<double> >(randeffect);

    if (tbeta_val.size() != p1)
    {
        Rprintf("Error in parameter p1 or truebeta.\n");
      return R_NilValue;	  
    }
   if (tgamma_val.size() != p2*2)
    {
        Rprintf("Error in parameter p2 (for two competing risks) or truegamma.\n");
      return R_NilValue;	  
    }
    /* allocate space for data */
    gsl_matrix *C = gsl_matrix_calloc(k,p2+2);
    gsl_matrix *FY= gsl_matrix_calloc(array_size, p1+p1a+bq+1);         
    gsl_vector *M1= gsl_vector_calloc(k);
    gsl_matrix_set_zero(FY);



    /* allocate space for help matrices */
    gsl_vector *RA= gsl_vector_calloc(p1a+1);
    gsl_vector *RI= gsl_vector_calloc(p1a+1);
    gsl_vector *S=gsl_vector_calloc(p1a+1);
    gsl_matrix *V=gsl_matrix_calloc(p1a+1,p1a+1);
    gsl_vector *W=gsl_vector_calloc(p1a+1);


    /* allocate space for true parameters */
    gsl_vector *tbeta = gsl_vector_calloc(p1),
               *ttheta = gsl_vector_calloc(K_num-1),
               *tvee = gsl_vector_calloc(g-1);
    gsl_matrix *tgamma = gsl_matrix_calloc(g,p2);
    gsl_matrix *tbeta2 = gsl_matrix_calloc(K_num-2,bq);

    gsl_matrix *VC= gsl_matrix_calloc(p1a+1,p1a+1);



  
	    /* assign the true value to parameters */
    
    double tsigmab0=reffect_val[0];
    double tsigmau=reffect_val[1];
    double	rho=0.9, sigmax=1.0, prob=0.5;
	  gsl_vector_set(tvee,0,reffect_val[2]);
    
    double tsigmab0u=sqrt(tsigmab0*tsigmau)*rho;

  for (size_t i=0;i<tbeta_val.size();i++) 
           gsl_vector_set(tbeta,i,tbeta_val[i]);
	   
  gsl_matrix_set(tbeta2,0,0,0);
  
  for (size_t i=0;i<ttheta_val.size();i++) 
           gsl_vector_set(ttheta,i,ttheta_val[i]);
  for (size_t i=0;i<g;i++)
        {   for(size_t j=0;j<p2;j++)
                {
       gsl_matrix_set(tgamma,i,j,tgamma_val[(i)*(p2)+j]);
				}
       }

     /* assign the true value to parameters */
	/*gsl_vector_set(tbeta,0,1);
    gsl_vector_set(tbeta,1,-1.5);
    gsl_vector_set(tbeta,2,-0.8);
    gsl_vector_set(ttheta,0,-0.5);
    gsl_vector_set(ttheta,1,1);

    gsl_matrix_set(tgamma,0,0,0.8);
    gsl_matrix_set(tgamma,0,1,-1);
    gsl_matrix_set(tgamma,1,0,0.5);
    gsl_matrix_set(tgamma,1,1,-1);

    gsl_vector_set(tvee,0,0.5);
*/

    gsl_matrix_set(VC,0,0,tsigmab0);
    gsl_matrix_set(VC,1,1,tsigmau);
    gsl_matrix_set(VC,1,0,tsigmab0u); 


    for(i=0;i<p1a+1;i++)
    {
        for(j=i+1;j<p1a+1;j++)    gsl_matrix_set(VC,i,j,gsl_matrix_get(VC,j,i));

    }


    gsl_linalg_SV_decomp(VC,V,S,W);

    gsl_matrix_set_zero(V);
    for(i=0;i<p1a+1;i++) gsl_matrix_set(V,i,i,sqrt(gsl_vector_get(S,i)));


    //double lamda01=0.15, lamda02=0.02;      
    double lamda01=0.15, lamda02=0.25;                       


    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);


    int point=0,t=0;
    double crate=0.0,rate1=0.0,rate2=0.0,max_censor=4.0;
    double temp=0.0,x1=0.0,x2=0.0,t1=0.0,t2=0.0,censor=0.0,u=0.0,b0=0.0,error=0.0,chisq_x=0.0;

   

    for(j=0;j<k;j++)
    {
        //startp=gsl_ran_flat(r,0,1);

        chisq_x=gsl_ran_chisq(r,degree_v)/degree_v;

        for(i=0;i<p1a+1;i++)    gsl_vector_set(RI,i,gsl_ran_gaussian(r,1));
        gsl_vector_scale(RI,1/sqrt(chisq_x));

        MulM(V,RI,RA);
        MulM(VC,RA,RI);
  
        b0=gsl_vector_get(RI,0);
        u=gsl_vector_get(RI,1);      

        x1=gsl_ran_gaussian(r,sqrt(sigmax))+2;
        x2=gsl_ran_bernoulli(r,prob);
        censor=gsl_ran_exponential(r,10);
        //censor=gsl_ran_exponential(r,40);
        gsl_matrix_set(C,j,2,x1);
        gsl_matrix_set(C,j,3,x2);
         
        temp=gsl_matrix_get(tgamma,0,0)*x1+gsl_matrix_get(tgamma,0,1)*x2;
        temp=exp(temp+u);
        temp=lamda01*temp;
        t1=gsl_ran_exponential(r, 1/temp);

        temp=gsl_matrix_get(tgamma,1,0)*x1+gsl_matrix_get(tgamma,1,1)*x2;
        temp=exp(temp+gsl_vector_get(tvee,0)*u);
        temp=lamda02*temp;
        t2=gsl_ran_exponential(r, 1/temp);


        if(t1<=t2 && t1<=censor && t1<=max_censor)
        {
            rate1+=1;
            gsl_matrix_set(C,j,0,t1);
            gsl_matrix_set(C,j,1,1);
            n=GetN(t1);               
        }

        if(t2<=t1 && t2<=censor && t2<=max_censor)
        {
            rate2+=1;
            gsl_matrix_set(C,j,0,t2);
            gsl_matrix_set(C,j,1,2);
            n=GetN(t2);               
        }

        if(Min(censor, max_censor)<=t1 && Min(censor, max_censor)<=t2)
        {
            crate+=1;
            gsl_matrix_set(C,j,0,Min(max_censor,censor));
            gsl_matrix_set(C,j,1,0);
            n=GetN(Min(max_censor,censor));
        }  


        for(i=point;i<point+n;i++)
        {
            gsl_matrix_set(FY,i,1,1);
            gsl_matrix_set(FY,i,2,x2-0.5);
            gsl_matrix_set(FY,i,3,(double)(i-point)*0.5);
            gsl_matrix_set(FY,i,4,x2-0.5);
            gsl_matrix_set(FY,i,5,(double)(i-point)*0.5*(x2-0.5));

            temp=gsl_ran_flat(r,0,1);
            error=log(temp/(1-temp));
            
            temp=gsl_vector_get(tbeta,0)*gsl_matrix_get(FY,i,3)+gsl_vector_get(tbeta,1)*(x2-0.5)
                 +gsl_vector_get(tbeta,2)*(x2-0.5)*gsl_matrix_get(FY,i,3);
            temp=temp+error+b0*gsl_matrix_get(FY,i,1);

            t=-1;
            do
            {
                t+=1;
                if(temp<gsl_vector_get(ttheta,t)) break;
            }while(t<K_num-2);

            gsl_matrix_set(FY,i,0,(double)(t+1));
            if(temp>gsl_vector_get(ttheta,K_num-2))  gsl_matrix_set(FY,i,0,(double)(K_num));
        }

        point+=n;
        gsl_vector_set(M1,j,(double)n);

    }

    n1 = point;

    gsl_matrix *Y = gsl_matrix_calloc(n1,p1+p1a+bq+1);      /* true matrix for longitudinal outcome Y **/
   
    for(i=0;i<n1;i++)
    {
        for(j=0;j<p1+p1a+bq+1;j++)  gsl_matrix_set(Y,i,j,gsl_matrix_get(FY,i,j));
    }




	 /***** ###### output data ############## *******/

  //write Y file
  
  FILE *output_F; 
  
  output_F=fopen(yfile.c_str(),"w");
  if (output_F == NULL)
  {Rprintf("Can't write Y file\n");
  }
    fprintf(output_F,"%s   ", "response time_RE intercept time x1 x2\n");

    for(i=0;i<Y->size1;i++)
    {
        for(j=0;j<Y->size2;j++)   fprintf(output_F,"%f   ", gsl_matrix_get(Y,i,j)); 
        if(i!=Y->size1-1)  fprintf(output_F,"\n");                      
    }
	
 //close writing files  
  fclose(output_F);
  
  //write C file
  
  output_F=fopen(cfile.c_str(),"w");
  if (output_F == NULL)
  {Rprintf("Can't write C file\n");
    //return R_NilValue;
  }
    fprintf(output_F,"%s   ", "surv  	comprisk	     x1         x2\n");

    for(i=0;i<C->size1;i++)
    {
        for(j=0;j<C->size2;j++)   fprintf(output_F,"%f   ", gsl_matrix_get(C,i,j)); 
        if(i!=C->size1-1)  fprintf(output_F,"\n");                      
    }
	
  //close writing files  
  fclose(output_F);
  
  //write M file
  
  output_F=fopen(mfile.c_str(),"w");
  if (output_F == NULL)
  {Rprintf("Can't write M file\n");
    // return R_NilValue;
  }
  
 for(i=0;i<M1->size;i++)
    {
        fprintf(output_F,"%f", gsl_vector_get(M1,i)); 
        if(i!=M1->size-1)  fprintf(output_F,"\n");                      
    }
	 //close writing files  
   fclose(output_F);
	

    gsl_matrix_free(FY);
    gsl_matrix_free(Y);
    gsl_matrix_free(C);
    gsl_vector_free(M1);

    gsl_matrix_free(tgamma);
    gsl_matrix_free(tbeta2);
    gsl_vector_free(tbeta);
    gsl_vector_free(ttheta);
    gsl_vector_free(tvee);


    gsl_matrix_free(VC);
    gsl_vector_free(RA);
    gsl_vector_free(RI);
    gsl_vector_free(S);
    gsl_matrix_free(V);
    gsl_vector_free(W);

   	//free random number memory
	gsl_rng_free(r);
	
  Rcpp::List ret;
  ret["censoring_rate"] = crate/(double)k;
  ret["rate1"] = rate1/(double)k;
  ret["rate2"] = rate2/(double)k;   
  ret["yfn"] = yfile.c_str(); 
  ret["cfn"] = cfile.c_str(); 
  ret["mfn"] = mfile.c_str(); 
  
  return ret;

}
