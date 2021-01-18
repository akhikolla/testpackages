#include "jmo.h"
namespace jmospace {
 int EM(
       gsl_vector *beta,
       gsl_matrix *beta2,
       gsl_vector *theta,
       gsl_matrix *gamma,
       gsl_vector *vee,
       gsl_matrix *H01,
       gsl_matrix *H02,
       gsl_matrix *sig,
       const gsl_matrix *Y,
       const gsl_matrix *C,
       const gsl_vector *M1,
       const int p1a,
       const int bq,
	   const int K_num,
	   const int j_max,
	   const int point,
	   const std::vector<double> xs,  
	   const std::vector<double> ws   
       )

{



    int p1=beta->size;
    int p2=gamma->size2;
    int g =gamma->size1;
    int a =H01->size2;
    int b =H02->size2;    
    
  
    int k = M1->size;

    int p,q,j,t,u,m;


    gsl_vector *FUNU=gsl_vector_calloc(k);
    gsl_vector *FUNUS=gsl_vector_calloc(k);
    gsl_matrix *FUNB=gsl_matrix_calloc(p1a,k);  
    gsl_matrix *FUNBS=gsl_matrix_calloc(p1a*(p1a+1)/2,k);  
    gsl_matrix *FUNBU=gsl_matrix_calloc(p1a,k);
   

    gsl_matrix *FUNE=gsl_matrix_calloc(g,k), 
               *FUNUSE=gsl_matrix_calloc(g-1,k),
               *FUNUE=gsl_matrix_calloc(g-1,k);


    gsl_matrix *FUNP=gsl_matrix_calloc(2*j_max,k), 
               *FUNPS=gsl_matrix_calloc(2*j_max,k),
               *FUNPR1=gsl_matrix_calloc(2*j_max,k),
               *FUNPR2=gsl_matrix_calloc(2*j_max,k);


    int status;


    status = GetE(FUNU,FUNUS,FUNB,FUNBS,FUNBU,FUNE,FUNUSE,FUNUE,FUNP,FUNPS,FUNPR1,FUNPR2,beta,beta2,theta,gamma,vee,H01,H02,sig,Y,C,M1,p1a,bq, K_num, j_max,point,xs,ws);

    if (status==100) return status;



    gsl_vector * SX = gsl_vector_calloc(p2);
    gsl_matrix * SXX = gsl_matrix_calloc(p2,p2);
    gsl_matrix * XX = gsl_matrix_calloc(p2,p2);
    gsl_vector * X = gsl_vector_calloc(p2);
    gsl_vector * gammai = gsl_vector_calloc(p2);
  
    double scalef;
  
    gsl_vector * Z = gsl_vector_calloc(p1);                       
    gsl_vector * SZ = gsl_vector_calloc(p1);
    gsl_matrix * SZZ = gsl_matrix_calloc(p1,p1);                  
    gsl_matrix * ZZ = gsl_matrix_calloc(p1,p1);      

    gsl_vector * D= gsl_vector_calloc(bq);                         
    gsl_vector * SD = gsl_vector_calloc(bq);
    gsl_matrix * SDD = gsl_matrix_calloc(bq,bq);                    
    gsl_matrix * DD = gsl_matrix_calloc(bq,bq);                            


   // gsl_vector * bi = gsl_vector_calloc(p1a);
   // gsl_matrix * bs = gsl_matrix_calloc(p1a,p1a);



    /* calculate beta, sig */

    gsl_matrix_set_zero(SZZ);
    gsl_vector_set_zero(SZ);

    p=0;
    gsl_matrix_set_zero(sig);

    for(j=0;j<k;j++)
    {
        u=(int)gsl_vector_get(M1,j);

        for(q=p;q<u+p;q++)
        {
            for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+bq+t));
            MulV(Z,ZZ);
           
            scalef=gsl_matrix_get(FUNP,2*(q-p),j)+gsl_matrix_get(FUNP,2*(q-p)+1,j)-1;
            gsl_vector_scale(Z,scalef);
            gsl_vector_add (SZ,Z);
 
            scalef=gsl_matrix_get(FUNP,2*(q-p),j)-gsl_matrix_get(FUNPS,2*(q-p),j)
                   +gsl_matrix_get(FUNP,2*(q-p)+1,j)-gsl_matrix_get(FUNPS,2*(q-p)+1,j);
            gsl_matrix_scale(ZZ,scalef);
            gsl_matrix_add(SZZ,ZZ);

        }
        
        p=p+u;

        gsl_matrix_set(sig,p1a,p1a,gsl_matrix_get(sig,p1a,p1a)+gsl_vector_get(FUNUS,j));
        for(t=0;t<p1a;t++)  gsl_matrix_set(sig,t,t,gsl_matrix_get(sig,t,t)+gsl_matrix_get(FUNBS,t,j));   
        for(t=0;t<p1a;t++)  gsl_matrix_set(sig,t,p1a,gsl_matrix_get(sig,t,p1a)+gsl_matrix_get(FUNBU,t,j));     
        
        for(q=1;q<p1a;q++)
        {
            for(t=0;t<p1a-q;t++)   gsl_matrix_set(sig,t,q+t,gsl_matrix_get(sig,t,q+t)+gsl_matrix_get(FUNBS,p1a+t+(q-1)*(p1a-1),j));
        }

    }

    gsl_matrix_scale(sig,1/(double)k);
    

    status=inv_matrix(SZZ);
    if(status==100)    return status;
    MulM(SZZ,SZ,Z);
  
    for(j=0;j<p1;j++)   gsl_vector_set(beta,j,gsl_vector_get(beta,j)+gsl_vector_get(Z,j));




    /* calculate beta2 */

    double scalefn,scalefd;


    for(t=2;t<K_num;t++)
    {
        gsl_vector_set_zero(SD);
        gsl_matrix_set_zero(SDD);

        p=0;
        for(j=0;j<k;j++)
        {
            u=(int)gsl_vector_get(M1,j);   
            for(q=p;q<u+p;q++)   
            {
                if(gsl_matrix_get(Y,q,0)==(double)t)
                {
                    for(m=0;m<bq;m++)   gsl_vector_set(D,m,gsl_matrix_get(Y,q,1+p1a+m));
                    MulV(D,DD);

                    scalefn=gsl_matrix_get(FUNPR1,2*(q-p),j);
                    scalefd=gsl_matrix_get(FUNPR2,2*(q-p),j);

                    gsl_vector_scale(D,0-scalefn);
                    gsl_matrix_scale(DD,scalefd);

                    gsl_vector_add(SD,D);
                    gsl_matrix_add(SDD,DD);
                }

                if(gsl_matrix_get(Y,q,0)==(double)t+1)  
                {

                    for(m=0;m<bq;m++)   gsl_vector_set(D,m,gsl_matrix_get(Y,q,1+p1a+m));
                    MulV(D,DD);

                    scalefn=0-gsl_matrix_get(FUNPR1,2*(q-p)+1,j);
                    scalefd=0-gsl_matrix_get(FUNPR2,2*(q-p)+1,j);

                    gsl_vector_scale(D,0-scalefn);
                    gsl_matrix_scale(DD,scalefd);

                    gsl_vector_add(SD,D);
                    gsl_matrix_add(SDD,DD);

                }
            }

            p=p+u;
 
        }

        if(gsl_matrix_get(SDD,0,0)<0)  
        {
            gsl_matrix_scale(SDD,-1);
            gsl_vector_scale(SD,-1);
        }

        status=inv_matrix(SDD);
        if(status==100)    return status;
        MulM(SDD,SD,D);
  
        for(j=0;j<bq;j++)   gsl_matrix_set(beta2,t-2,j,gsl_matrix_get(beta2,t-2,j)-gsl_vector_get(D,j));



    }



    /* calculate theta */

    gsl_vector *theta_n=gsl_vector_calloc(K_num-1);
    gsl_vector *theta_d=gsl_vector_calloc(K_num-1);
 
    gsl_vector_set_zero(theta_n);
    gsl_vector_set_zero(theta_d);

    p=0;

    for(j=0;j<k;j++)
    {
        u=(int)gsl_vector_get(M1,j);   
        
        for(q=p;q<u+p;q++)   
        {
            for(t=1;t<K_num;t++)
            {
                if(gsl_matrix_get(Y,q,0)==(double)t)  
                {
                    gsl_vector_set(theta_n,t-1,gsl_vector_get(theta_n,t-1)+gsl_matrix_get(FUNPR1,2*(q-p),j));
                    gsl_vector_set(theta_d,t-1,gsl_vector_get(theta_d,t-1)+gsl_matrix_get(FUNPR2,2*(q-p),j));
                }

                if(gsl_matrix_get(Y,q,0)==(double)t && t>1)  
                {
                    gsl_vector_set(theta_n,t-2,gsl_vector_get(theta_n,t-2)-gsl_matrix_get(FUNPR1,2*(q-p)+1,j));
                    gsl_vector_set(theta_d,t-2,gsl_vector_get(theta_d,t-2)-gsl_matrix_get(FUNPR2,2*(q-p)+1,j));
                }
            }

            if(gsl_matrix_get(Y,q,0)==K_num)  
            {
                gsl_vector_set(theta_n,K_num-2,gsl_vector_get(theta_n,K_num-2)-gsl_matrix_get(FUNPR1,2*(q-p)+1,j));
                gsl_vector_set(theta_d,K_num-2,gsl_vector_get(theta_d,K_num-2)-gsl_matrix_get(FUNPR2,2*(q-p)+1,j));
            }
 
        }

        p=p+u;
    }


    for(t=0;t<K_num-1;t++)   
    gsl_vector_set(theta,t,gsl_vector_get(theta,t)-gsl_vector_get(theta_n,t)/gsl_vector_get(theta_d,t));

    gsl_vector_free(theta_n);
    gsl_vector_free(theta_d);




    /* calculate H01 H02 */


    double dem, num;

    for(p=0;p<a;p++)
    {
        dem=0;
   
        for(j=0;j<k;j++)
        {
            if(gsl_matrix_get(C,j,0)>=gsl_matrix_get(H01,0,p))
            {   
                for(u=0;u<p2;u++)  
                {
                    gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                    gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
                }

                dem=dem+gsl_matrix_get(FUNE,0,j)*exp(MulVV(X,gammai));
                
            }
        }

        gsl_matrix_set(H01,2,p,gsl_matrix_get(H01,1,p)/dem);

    }

     
    for(p=0;p<b;p++)
    {
        dem=0;
   
        for(j=0;j<k;j++)
        {
            if(gsl_matrix_get(C,j,0)>=gsl_matrix_get(H02,0,p))
            {   
                for(u=0;u<p2;u++)  
                {
                    gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                    gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));
                }

                dem=dem+gsl_matrix_get(FUNE,1,j)*exp(MulVV(X,gammai));
                
            }
        }

        gsl_matrix_set(H02,2,p,gsl_matrix_get(H02,1,p)/dem);

    }


   
    /* calculate gamma */

    gsl_matrix_set_zero(SXX);
    gsl_vector_set_zero(SX);
  
  
    for(j=0;j<k;j++)
    {

        for(u=0;u<p2;u++)  
        {
            gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
            gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
        }
        
        MulV(X,XX);

        scalef=CH(H01,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,0,j);
        gsl_matrix_scale(XX,scalef); 

        if((int)gsl_matrix_get(C,j,1) == 1)
        {
            scalef=1-CH(H01,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,0,j);  
            gsl_vector_scale(X,scalef); 
        }

        if((int)gsl_matrix_get(C,j,1) != 1)
        {
            scalef=0-CH(H01,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,0,j);  
            gsl_vector_scale(X,scalef); 
        }

        gsl_matrix_add(SXX, XX);
        gsl_vector_add(SX, X);

    }
       

    status=inv_matrix(SXX);
    if(status==100)  return status;
    MulM(SXX,SX,X);
    
    for(j=0;j<p2;j++)   gsl_matrix_set(gamma,0,j,gsl_matrix_get(gamma,0,j)+gsl_vector_get(X,j));





    gsl_matrix_set_zero(SXX);
    gsl_vector_set_zero(SX);
  
    for(j=0;j<k;j++)
    {
        for(u=0;u<p2;u++)  
        {
            gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
            gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));
        }
        
        MulV(X,XX);

        scalef=CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,1,j);
        gsl_matrix_scale(XX,scalef); 

        if((int)gsl_matrix_get(C,j,1) == 2)
        {
            scalef=1-CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,1,j);  
            gsl_vector_scale(X,scalef); 
        }

        if((int)gsl_matrix_get(C,j,1) != 2)
        {
            scalef=0-CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNE,1,j);  
            gsl_vector_scale(X,scalef); 
        }

        gsl_matrix_add(SXX, XX);
        gsl_vector_add(SX, X);

    }
       

    status=inv_matrix(SXX);
    if(status==100)  return status;
    MulM(SXX,SX,X);
    
    for(j=0;j<p2;j++)   gsl_matrix_set(gamma,1,j,gsl_matrix_get(gamma,1,j)+gsl_vector_get(X,j));





    /* calculate vee */

    dem=0; num=0;
    for(j=0;j<k;j++)
    {
        for(u=0;u<p2;u++)  
        {
            gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
            gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));
        }

        dem=dem+CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNUSE,0,j);

        if((int)gsl_matrix_get(C,j,1)==2) 
        {
            num=num+gsl_vector_get(FUNU,j)-CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNUE,0,j);
        }

        if((int)gsl_matrix_get(C,j,1)!=2) 
        {
            num=num-CH(H02,gsl_matrix_get(C,j,0))*exp(MulVV(X,gammai))*gsl_matrix_get(FUNUE,0,j);
        }
 
    }

    gsl_vector_set(vee,0,gsl_vector_get(vee,0)+num/dem);



    gsl_vector_free(Z);
    gsl_vector_free(SZ);
    gsl_vector_free(X);
    gsl_vector_free(gammai);
    gsl_vector_free(SX);
    gsl_matrix_free(SZZ);
    gsl_matrix_free(ZZ);
    gsl_matrix_free(SXX);
    gsl_matrix_free(XX);

    gsl_vector_free(D);
    gsl_vector_free(SD);
    gsl_matrix_free(DD);
    gsl_matrix_free(SDD);


    gsl_vector_free(FUNU);
    gsl_vector_free(FUNUS);
    gsl_matrix_free(FUNB);
    gsl_matrix_free(FUNBS);
    gsl_matrix_free(FUNBU);
         

    gsl_matrix_free(FUNE); 
    gsl_matrix_free(FUNUSE); 
    gsl_matrix_free(FUNUE); 


    gsl_matrix_free(FUNP); 
    gsl_matrix_free(FUNPS); 
    gsl_matrix_free(FUNPR1);
    gsl_matrix_free(FUNPR2); 


    return 0;

}



 int GetCov(
           gsl_matrix *Cov,
           const gsl_vector *beta,
           const gsl_matrix *beta2,
           const gsl_vector *theta,
           const gsl_matrix *gamma,
           const gsl_vector *vee,
           const gsl_matrix *H01,
           const gsl_matrix *H02,
           const gsl_matrix *sig,
           const gsl_matrix *Y,
           const gsl_matrix *C,
           const gsl_vector *M1,
           const int p1a,
           const int bq,
		   const int K_num,
		   const int j_max,
		   const int point,
		   const std::vector<double> xs,  
	       const std::vector<double> ws
           )
{

    int p1=beta->size;
    int p2=gamma->size2;
    int g =gamma->size1;
    int a =H01->size2;
    int b =H02->size2;      


    int k = M1->size;


    int p,q,j,t,u,i,r,m;


    double temp,scalef;

    gsl_vector *S = gsl_vector_calloc(Cov->size1);
    gsl_matrix *SS= gsl_matrix_calloc(Cov->size1,Cov->size1);
    gsl_vector *TS = gsl_vector_calloc(Cov->size1);


    gsl_vector *FUNU=gsl_vector_calloc(k);
    gsl_vector *FUNUS=gsl_vector_calloc(k);
    gsl_matrix *FUNB=gsl_matrix_calloc(p1a,k);  
    gsl_matrix *FUNBS=gsl_matrix_calloc(p1a*(p1a+1)/2,k);  
    gsl_matrix *FUNBU=gsl_matrix_calloc(p1a,k);
 
   

    gsl_matrix *FUNE=gsl_matrix_calloc(g,k), 
               *FUNUSE=gsl_matrix_calloc(g-1,k),
               *FUNUE=gsl_matrix_calloc(g-1,k);


    gsl_matrix *FUNP=gsl_matrix_calloc(2*j_max,k), 
               *FUNPS=gsl_matrix_calloc(2*j_max,k),
               *FUNPR1=gsl_matrix_calloc(2*j_max,k),
               *FUNPR2=gsl_matrix_calloc(2*j_max,k);


    int status;
    status = GetE(FUNU,FUNUS,FUNB,FUNBS,FUNBU,FUNE,FUNUSE,FUNUE,FUNP,FUNPS,FUNPR1,FUNPR2,beta,beta2,theta,gamma,vee,H01,H02,
                  sig,Y,C,M1,p1a,bq,K_num, j_max,point,xs,ws);

    if (status==100) return status;


    gsl_vector * X = gsl_vector_calloc(p2);                         /* covariates for C */
    gsl_vector * RX= gsl_vector_calloc(p2);                         /* covariates for subjects in risk set */
    gsl_vector *SX = gsl_vector_calloc(p2);           
    gsl_vector *SRX= gsl_vector_calloc(p2);
    gsl_vector * gammai = gsl_vector_calloc(p2);

    gsl_vector * Z = gsl_vector_calloc(p1);                         /* covariates for Y */
    gsl_vector * SZ = gsl_vector_calloc(p1);
    gsl_matrix * SZZ = gsl_matrix_calloc(p1,p1);                    /* store sum of ZZ' */
    gsl_matrix * ZZ = gsl_matrix_calloc(p1,p1);                     /* store ZZ' */

    gsl_vector * D = gsl_vector_calloc(bq);                        
    gsl_vector * SD = gsl_vector_calloc(bq);

    gsl_vector * xtilde = gsl_vector_calloc(p1a);
    gsl_vector * xtilde1 = gsl_vector_calloc(p1a);
    gsl_vector * bi = gsl_vector_calloc(p1a);
    gsl_matrix * bs = gsl_matrix_calloc(p1a,p1a);
//    gsl_matrix * bu = gsl_matrix_calloc(p1a+1,p1a+1);

    gsl_matrix * VC = gsl_matrix_calloc(p1a+1,p1a+1);
    gsl_matrix * VI = gsl_matrix_calloc(p1a+1,p1a+1);
    gsl_matrix * HE = gsl_matrix_calloc(p1a+1,p1a+1);

    gsl_matrix_memcpy(VC,sig);

    for(i=0;i<p1a+1;i++)
    {
        for(j=0;j<i;j++)    gsl_matrix_set(VC,i,j,gsl_matrix_get(VC,j,i));
    }
  
    status=inv_matrix(VC);
    if(status==100) return 100;



    p=0;
    gsl_matrix_set_zero(Cov);
    gsl_vector_set_zero(TS);

    for(j=0;j<k;j++)
    {
    
        gsl_vector_set_zero(S);
       
        u=(int)gsl_vector_get(M1,j);

        for(t=0;t<p1a;t++)   gsl_vector_set(bi,t,gsl_matrix_get(FUNB,t,j));
        for(t=0;t<p1a;t++)   gsl_matrix_set(bs,t,t,gsl_matrix_get(FUNBS,t,j)); 
 
        for(q=1;q<p1a;q++)
        {
            for(t=0;t<p1a-q;t++)   gsl_matrix_set(bs,t,q+t,gsl_matrix_get(FUNBS,p1a+t+(q-1)*(p1a-1),j));
        }


        for(t=0;t<p1a;t++)
        {
            for(q=0;q<t;q++)   gsl_matrix_set(bs,t,q,gsl_matrix_get(bs,q,t));
        }



        /* calculate score for beta */

        gsl_vector_set_zero(SZ); 
   
        for(q=p;q<u+p;q++)
        {
            for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+bq+t));
           
            scalef=gsl_matrix_get(FUNP,2*(q-p),j)+gsl_matrix_get(FUNP,2*(q-p)+1,j)-1;
            gsl_vector_scale(Z,scalef);
            gsl_vector_add (SZ,Z);
        }

        for(q=0;q<p1;q++)  gsl_vector_set(S,q,gsl_vector_get(SZ,q));



        /* calculate score for beta2 */

        for(t=2;t<K_num;t++)
        {
            gsl_vector_set_zero(SD);
            for(q=p;q<u+p;q++)   
            {
                if(gsl_matrix_get(Y,q,0)==(double)t)
                {
                    for(m=0;m<bq;m++)   gsl_vector_set(D,m,gsl_matrix_get(Y,q,1+p1a+m));

                    scalef=gsl_matrix_get(FUNPR1,2*(q-p),j);
                    gsl_vector_scale(D,0-scalef);

                    gsl_vector_add(SD,D);
                }

                if(gsl_matrix_get(Y,q,0)==(double)t+1)  
                {
                    for(m=0;m<bq;m++)   gsl_vector_set(D,m,gsl_matrix_get(Y,q,1+p1a+m));

                    scalef=0-gsl_matrix_get(FUNPR1,2*(q-p)+1,j);
                    gsl_vector_scale(D,0-scalef);

                    gsl_vector_add(SD,D);

                }
            }
  
            for(m=0;m<bq;m++)   gsl_vector_set(S,p1+(t-2)*bq+m,gsl_vector_get(SD,m));

        }




        /* calculate score for theta */

        gsl_vector *theta_n=gsl_vector_calloc(K_num-1);
        gsl_vector_set_zero(theta_n);

        for(q=p;q<u+p;q++)   
        {
            for(t=1;t<K_num;t++)
            {
                if(gsl_matrix_get(Y,q,0)==(double)t)  
                {
                    gsl_vector_set(theta_n,t-1,gsl_vector_get(theta_n,t-1)+gsl_matrix_get(FUNPR1,2*(q-p),j));
                }

                if(gsl_matrix_get(Y,q,0)==(double)t && t>1)  
                {
                    gsl_vector_set(theta_n,t-2,gsl_vector_get(theta_n,t-2)-gsl_matrix_get(FUNPR1,2*(q-p)+1,j));
                }
            }

            if(gsl_matrix_get(Y,q,0)==K_num)  
            {
                gsl_vector_set(theta_n,K_num-2,gsl_vector_get(theta_n,K_num-2)-gsl_matrix_get(FUNPR1,2*(q-p)+1,j));
            }

        }

        for(q=0;q<K_num-1;q++)  gsl_vector_set(S,p1+(K_num-2)*bq+q,gsl_vector_get(theta_n,q));

        gsl_vector_free(theta_n);

        p=p+u;




        /* calculate score for gamma */
         

        /*  gamma11, gamma12 */

        for(u=0;u<p2;u++)   gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));

        gsl_vector_set_zero(SRX);
        for(r=0;r<a;r++)
        {
            if(gsl_matrix_get(H01,0,r)<=gsl_matrix_get(C,j,0))
            {                 
                gsl_vector_set_zero(SX);
                temp=0;                    

                for(q=0;q<k;q++)
                {
                    if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(H01,0,r))
                    {
                        for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                        temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q);
                        gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q));
                        gsl_vector_add(SX,RX);                 
                    }
                }

                gsl_vector_scale(SX, gsl_matrix_get(H01,1,r)/(temp*temp));
                gsl_vector_add(SRX,SX);
            }
        }

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        gsl_vector_scale(SX, CH(H01,gsl_matrix_get(C,j,0)));

        gsl_vector_sub(SRX,SX);

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        gsl_vector_scale(SRX, exp(MulVV(SX,gammai))*gsl_matrix_get(FUNE,0,j));


        if((int)gsl_matrix_get(C,j,1) != 1)
        {
            for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+(K_num-2)*bq+K_num-1,gsl_vector_get(SRX,q));
        }

        if((int)gsl_matrix_get(C,j,1) ==1)
        {
            
            gsl_vector_set_zero(SX);
            temp=0;

            for(q=0;q<k;q++)
            {
                if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(C,j,0))
                {
                    for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                    temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q);
                    gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q));
                    gsl_vector_add(SX,RX);
             
                }
            }
            gsl_vector_scale(SX, 1/temp);

            for(u=0;u<p2;u++)  gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));          
            gsl_vector_sub(X,SX); 
             
            gsl_vector_add(X, SRX);

            for(q=0;q<p2;q++)   gsl_vector_set(S,q+(K_num-2)*bq+p1+K_num-1,gsl_vector_get(X,q));

        }



        /*  gamma21, gamma22 */

        for(u=0;u<p2;u++)   gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));

        gsl_vector_set_zero(SRX);
        for(r=0;r<b;r++)
        {
            if(gsl_matrix_get(H02,0,r)<=gsl_matrix_get(C,j,0))
            {                 
                gsl_vector_set_zero(SX);
                temp=0;                    

                for(q=0;q<k;q++)
                {
                    if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(H02,0,r))
                    {
                        for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                        temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q);
                        gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q));
                        gsl_vector_add(SX,RX);                 
                    }
                }

                gsl_vector_scale(SX, gsl_matrix_get(H02,1,r)/(temp*temp));
                gsl_vector_add(SRX,SX);
            }
        }

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        gsl_vector_scale(SX, CH(H02,gsl_matrix_get(C,j,0)));

        gsl_vector_sub(SRX,SX);

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        gsl_vector_scale(SRX, exp(MulVV(SX,gammai))*gsl_matrix_get(FUNE,1,j));


        if((int)gsl_matrix_get(C,j,1) != 2)
        {
            for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+(K_num-2)*bq+K_num+p2-1,gsl_vector_get(SRX,q));
        }

        if((int)gsl_matrix_get(C,j,1) ==2)
        {
            
            gsl_vector_set_zero(SX);
            temp=0;

            for(q=0;q<k;q++)
            {
                if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(C,j,0))
                {
                    for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                    temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q);
                    gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q));
                    gsl_vector_add(SX,RX);
             
                }
            }
            gsl_vector_scale(SX, 1/temp);

            for(u=0;u<p2;u++)  gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));          
            gsl_vector_sub(X,SX); 
             
            gsl_vector_add(X, SRX);

            for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+(K_num-2)*bq+K_num+p2-1,gsl_vector_get(X,q));

        }




        /* calculate score for vee */

        double su,  sru; 



        /*  vee2 */

        for(u=0;u<p2;u++)   gsl_vector_set(gammai,u,gsl_matrix_get(gamma,1,u));

        sru=0;
        for(r=0;r<b;r++)
        {
            if(gsl_matrix_get(H02,0,r)<=gsl_matrix_get(C,j,0))
            {                 
                su=0;
                temp=0;                    

                for(q=0;q<k;q++)
                {
                    if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(H02,0,r))
                    {
                        for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                        temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q);
                        su+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNUE,0,q);               
                    }
                }

                sru=sru+su*gsl_matrix_get(H02,1,r)/(temp*temp);
            }
        }

        for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
        sru=sru*gsl_matrix_get(FUNE,1,j)*exp(MulVV(SX,gammai));

        su=CH(H02,gsl_matrix_get(C,j,0))*gsl_matrix_get(FUNUE,0,j)*exp(MulVV(SX,gammai));

        sru=sru-su;

        if((int)gsl_matrix_get(C,j,1) != 2)
        {
            gsl_vector_set(S,p1+(K_num-2)*bq+K_num+2*p2-1,sru);
        }

        if((int)gsl_matrix_get(C,j,1) ==2)
        {           
            su=0;
            temp=0;

            for(q=0;q<k;q++)
            {
                if(gsl_matrix_get(C,q,0)>=gsl_matrix_get(C,j,0))
                {
                    for(u=0;u<p2;u++)  gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));

                    temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,1,q);
                    su+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNUE,0,q);
                }
            }
            su=su/temp;

            su=gsl_vector_get(FUNU,j)-su;
            su=su+sru;

            gsl_vector_set(S,p1+(K_num-2)*bq+K_num+2*p2-1,su);

        }



        /*  Sigma matrix  */


        for(t=0;t<p1a;t++)
        {
            for(r=0;r<p1a;r++)   gsl_matrix_set(VI,t,r,gsl_matrix_get(bs,t,r));
        }

        gsl_matrix_set(VI,p1a,p1a,gsl_vector_get(FUNUS,j));
        for(t=0;t<p1a;t++)   gsl_matrix_set(VI,t,p1a,gsl_matrix_get(FUNBU,t,j));
        for(t=0;t<p1a;t++)   gsl_matrix_set(VI,p1a,t,gsl_matrix_get(FUNBU,t,j));        


        MulMM(VC,VI,HE);
        MulMM(HE,VC,VI);

        if(p1a==1)
        {
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num,(gsl_matrix_get(VI,0,0)-gsl_matrix_get(VC,0,0))/2);
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+1,(gsl_matrix_get(VI,1,1)-gsl_matrix_get(VC,1,1))/2);
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+2,gsl_matrix_get(VI,1,0)-gsl_matrix_get(VC,1,0));
        }

        if(p1a==2)
        {
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num,(gsl_matrix_get(VI,0,0)-gsl_matrix_get(VC,0,0))/2);
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+1,(gsl_matrix_get(VI,1,1)-gsl_matrix_get(VC,1,1))/2);
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+2,(gsl_matrix_get(VI,2,2)-gsl_matrix_get(VC,2,2))/2);
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+3,gsl_matrix_get(VI,0,1)-gsl_matrix_get(VC,0,1));
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+4,gsl_matrix_get(VI,1,2)-gsl_matrix_get(VC,1,2));
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+5,gsl_matrix_get(VI,0,2)-gsl_matrix_get(VC,0,2));
        }

        if(p1a==3)
        {
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num,(gsl_matrix_get(VI,0,0)-gsl_matrix_get(VC,0,0))/2);
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+1,(gsl_matrix_get(VI,1,1)-gsl_matrix_get(VC,1,1))/2);
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+2,(gsl_matrix_get(VI,2,2)-gsl_matrix_get(VC,2,2))/2);
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+3,(gsl_matrix_get(VI,3,3)-gsl_matrix_get(VC,3,3))/2);
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+4,gsl_matrix_get(VI,0,1)-gsl_matrix_get(VC,0,1));
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+5,gsl_matrix_get(VI,1,2)-gsl_matrix_get(VC,1,2));
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+6,gsl_matrix_get(VI,2,3)-gsl_matrix_get(VC,2,3));
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+7,gsl_matrix_get(VI,0,2)-gsl_matrix_get(VC,0,2));
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+8,gsl_matrix_get(VI,1,3)-gsl_matrix_get(VC,1,3));
            gsl_vector_set(S,p1+(K_num-2)*bq+2*p2+K_num+9,gsl_matrix_get(VI,0,3)-gsl_matrix_get(VC,0,3));
        }



        MulV(S,SS);
        gsl_matrix_add(Cov,SS);
        gsl_vector_add(TS,S);

    }

/*
 for(i=0;(unsigned)i<Cov->size1;i++) 
	 Rprintf("%f\n",gsl_vector_get(TS,i));*/


    gsl_vector_free(FUNU);
    gsl_vector_free(FUNUS);
    gsl_matrix_free(FUNB);
    gsl_matrix_free(FUNBS);
    gsl_matrix_free(FUNBU);
         

    gsl_matrix_free(FUNE); 
    gsl_matrix_free(FUNUSE); 
    gsl_matrix_free(FUNUE); 

    gsl_matrix_free(FUNP); 
    gsl_matrix_free(FUNPS); 
    gsl_matrix_free(FUNPR1);
    gsl_matrix_free(FUNPR2); 

    gsl_vector_free(D);                        
    gsl_vector_free(SD);
    gsl_vector_free(Z);                        
    gsl_vector_free(SZ);
    gsl_matrix_free(ZZ);
    gsl_matrix_free(SZZ);
    gsl_vector_free(X);
    gsl_vector_free(RX);                       
    gsl_vector_free(SX);           
    gsl_vector_free(SRX);
    gsl_vector_free(gammai);

    gsl_vector_free(bi);
    gsl_vector_free(xtilde);
    gsl_vector_free(xtilde1);
    gsl_matrix_free(bs);


    gsl_matrix_free(VC);
    gsl_matrix_free(VI);
    gsl_matrix_free(HE);
    gsl_vector_free(S);                        
    gsl_matrix_free(SS);

    gsl_vector_free(TS);


    return 0;
}




 int GetE(
          gsl_vector *FUNU,
          gsl_vector *FUNUS,
          gsl_matrix *FUNB,
          gsl_matrix *FUNBS,
          gsl_matrix *FUNBU,
          gsl_matrix *FUNE,
          gsl_matrix *FUNUSE,
          gsl_matrix *FUNUE,
          gsl_matrix *FUNP,
          gsl_matrix *FUNPS,
          gsl_matrix *FUNPR1,
          gsl_matrix *FUNPR2,
          const gsl_vector *beta,
          const gsl_matrix *beta2,
          const gsl_vector *theta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int bq,
		  const int K_num,
		  const int j_max,
		  const int point,
		  const std::vector<double> xs,  
	      const std::vector<double> ws
          )

{


    int p1=beta->size;
    int p2=gamma->size2;
   
    int k = M1->size;

    int i,j,q,t,m,c;
    double dem,temp,mu,mu1,mu2,et;
    double cuh01,cuh02,haz01,haz02,xgamma1,xgamma2;

    double pk,pk1,pk_1;


    gsl_vector *Z = gsl_vector_calloc(p1),
               *X = gsl_vector_calloc(p2),
               *D = gsl_vector_calloc(bq),
               *xtilde = gsl_vector_calloc(p1a),
               *gammai = gsl_vector_calloc(p2);


    gsl_vector_set_zero(FUNU);
    gsl_vector_set_zero(FUNUS);
    gsl_matrix_set_zero(FUNB);
    gsl_matrix_set_zero(FUNBS);
    gsl_matrix_set_zero(FUNBU);
    gsl_matrix_set_zero(FUNE);
    gsl_matrix_set_zero(FUNUSE);
    gsl_matrix_set_zero(FUNUE);
    gsl_matrix_set_zero(FUNP);
    gsl_matrix_set_zero(FUNPS);
    gsl_matrix_set_zero(FUNPR1);
    gsl_matrix_set_zero(FUNPR2);


 
    int db0,db1,db2,du;


    gsl_vector *xi = gsl_vector_calloc(point);
    gsl_vector *wi = gsl_vector_calloc(point);



    for(i=0;i<point/2;i++)   gsl_vector_set(xi,i,xs[i]);
    for(i=0;i<point/2;i++)   gsl_vector_set(wi,i,ws[i]);

    for(i=0;i<point/2;i++)   gsl_vector_set(xi,point-i-1, 0-gsl_vector_get(xi,i));
    for(i=0;i<point/2;i++)   gsl_vector_set(wi,point-i-1, gsl_vector_get(wi,i));



    gsl_matrix *VC=gsl_matrix_calloc(p1a+1,p1a+1);
    gsl_vector *S=gsl_vector_calloc(p1a+1);
    gsl_matrix *V=gsl_matrix_calloc(p1a+1,p1a+1);
    gsl_matrix *VV=gsl_matrix_calloc(p1a+1,p1a+1);
    gsl_vector *W=gsl_vector_calloc(p1a+1);


    gsl_matrix_memcpy(VC,sig);


    for(i=0;i<p1a+1;i++)
    {
        for(j=0;j<i;j++)    gsl_matrix_set(VC,i,j,gsl_matrix_get(VC,j,i));
    }



    gsl_linalg_SV_decomp (VC,V,S,W);
    gsl_matrix_set_zero(V);
    for(i=0;i<p1a+1;i++)  gsl_matrix_set(V,i,i,sqrt(gsl_vector_get(S,i)));
    MulMM(VC,V,VV);


    gsl_matrix_scale(VV,sqrt((double)2));


    gsl_vector *ri=gsl_vector_calloc(p1a+1);
    gsl_vector *ti=gsl_vector_calloc(p1a+1);


    m=0;    
    for(j=0;j<k;j++)
    {


        dem=0;

        q=(int)gsl_vector_get(M1,j);

        cuh01=CH(H01,gsl_matrix_get(C,j,0));
        cuh02=CH(H02,gsl_matrix_get(C,j,0));

        haz01=HAZ(H01,gsl_matrix_get(C,j,0));
        haz02=HAZ(H02,gsl_matrix_get(C,j,0));


        for(i=0;i<p2;i++) 
        {
            gsl_vector_set(X,i,gsl_matrix_get(C,j,2+i));
            gsl_vector_set(gammai,i,gsl_matrix_get(gamma,0,i)); 
        }
        xgamma1=MulVV(X,gammai);

        for(i=0;i<p2;i++) gsl_vector_set(gammai,i,gsl_matrix_get(gamma,1,i));
        xgamma2=MulVV(X,gammai);



        if(p1a==1)
        {

            for(db0=0;db0<point;db0++)
            {
                for(du=0;du<point;du++)
                {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);

                    temp=exp((double)10);

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+bq+1));
                        for(t=0;t<bq;t++)  gsl_vector_set(D,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        et=gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0);

                        c=(int)gsl_matrix_get(Y,m+i,0);
                        if(c==1)
                        {
                            pk=exp(gsl_vector_get(theta,c-1)-mu-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-et));
                            pk1=0;
                        }

                        else if(c==K_num)
                        {
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=1;
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }
                        else if(c==2)
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-et));
                        }

                        else
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }

                        temp*=(pk-pk1);
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=gsl_vector_get(wi,db0)*gsl_vector_get(wi,du);

                    dem+=temp;

                    gsl_vector_set(FUNU,j,gsl_vector_get(FUNU,j)+temp*gsl_vector_get(ti,p1a));
                    gsl_vector_set(FUNUS,j,gsl_vector_get(FUNUS,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a)));
                    
                    for(i=0;i<p1a;i++)  
                    {
                        gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)+temp*gsl_vector_get(ti,i));
                        gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i)));
                        gsl_matrix_set(FUNBU,i,j,gsl_matrix_get(FUNBU,i,j)+temp*gsl_vector_get(ti,i)*gsl_vector_get(ti,p1a));
                    }



                    for(i=1;i<p1a;i++)
                    {
                        for(t=0;t<p1a-i;t++)   gsl_matrix_set(FUNBS,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j)
                                               +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i));
                    }


                    gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)+temp*exp(gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNE,1,j,gsl_matrix_get(FUNE,1,j)+temp*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));

                    gsl_matrix_set(FUNUSE,0,j,gsl_matrix_get(FUNUSE,0,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a))*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNUE,0,j,gsl_matrix_get(FUNUE,0,j)+temp*gsl_vector_get(ti,p1a)*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));


                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+bq+1));
                        for(t=0;t<bq;t++)  gsl_vector_set(D,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        et=gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0);

                        c=(int)gsl_matrix_get(Y,m+i,0);
                        if(c==1)
                        {
                            pk=exp(gsl_vector_get(theta,c-1)-mu-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-et));
                            pk1=0;
                        }

                        else if(c==K_num)
                        {
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=1;
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }
                        else if(c==2)
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-et));
                        }

                        else
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }


                        pk_1=(pk-pk1);
                          /*if(pk_1==0) return 100; */

                        if(pk_1!=0)
                        {
                        gsl_matrix_set(FUNP,i*2,j,gsl_matrix_get(FUNP,i*2,j)+temp*pk);
                        gsl_matrix_set(FUNP,i*2+1,j,gsl_matrix_get(FUNP,i*2+1,j)+temp*pk1);

                        gsl_matrix_set(FUNPS,i*2,j,gsl_matrix_get(FUNPS,i*2,j)+temp*pk*pk);
                        gsl_matrix_set(FUNPS,i*2+1,j,gsl_matrix_get(FUNPS,i*2+1,j)+temp*pk1*pk1);

                        gsl_matrix_set(FUNPR1,i*2,j,gsl_matrix_get(FUNPR1,i*2,j)+temp*((pk-pk*pk)/pk_1));
                        gsl_matrix_set(FUNPR1,i*2+1,j,gsl_matrix_get(FUNPR1,i*2+1,j)+temp*((pk1-pk1*pk1)/pk_1));

                        gsl_matrix_set(FUNPR2,i*2,j,gsl_matrix_get(FUNPR2,i*2,j)+
                                       temp*((pk-pk*pk)*(1-2*pk)*(pk-pk1)/(pk_1*pk_1)-(pk-pk*pk)*(pk-pk*pk)/(pk_1*pk_1)));

                        gsl_matrix_set(FUNPR2,i*2+1,j,gsl_matrix_get(FUNPR2,i*2+1,j)+
                                       temp*((pk1-pk1*pk1)*(1-2*pk1)*(pk-pk1)/(pk_1*pk_1)+(pk1-pk1*pk1)*(pk1-pk1*pk1)/(pk_1*pk_1)));
                       }

                        if(pk_1==0)
                        {
                        gsl_matrix_set(FUNP,i*2,j,gsl_matrix_get(FUNP,i*2,j)+temp*pk);
                        gsl_matrix_set(FUNP,i*2+1,j,gsl_matrix_get(FUNP,i*2+1,j)+temp*pk1);

                        gsl_matrix_set(FUNPS,i*2,j,gsl_matrix_get(FUNPS,i*2,j)+temp*pk*pk);
                        gsl_matrix_set(FUNPS,i*2+1,j,gsl_matrix_get(FUNPS,i*2+1,j)+temp*pk1*pk1);

                        gsl_matrix_set(FUNPR1,i*2,j,gsl_matrix_get(FUNPR1,i*2,j)+temp);
                        gsl_matrix_set(FUNPR1,i*2+1,j,gsl_matrix_get(FUNPR1,i*2+1,j)+temp);

                        gsl_matrix_set(FUNPR2,i*2,j,gsl_matrix_get(FUNPR2,i*2,j)+
                                       temp*(pk*(1-pk)));

                        gsl_matrix_set(FUNPR2,i*2+1,j,gsl_matrix_get(FUNPR2,i*2+1,j)+
                                       temp*(pk1*(1-pk1)));
                       }

                    }

                }
            }
        }




        if(p1a==2)
        {

            for(db0=0;db0<point;db0++)
            {
                for(db1=0;db1<point;db1++)
                {
                    for(du=0;du<point;du++)
                    {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,db1));
                    gsl_vector_set(ri,2,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);


                    temp=exp((double)10);

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+bq+1));
                        for(t=0;t<bq;t++)  gsl_vector_set(D,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        et=gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0);

                        c=(int)gsl_matrix_get(Y,m+i,0);
                        if(c==1)
                        {
                            pk=exp(gsl_vector_get(theta,c-1)-mu-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-et));
                            pk1=0;
                        }

                        else if(c==K_num)
                        {
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=1;
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }
                        else if(c==2)
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-et));
                        }

                        else
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }

                        temp*=(pk-pk1);
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=gsl_vector_get(wi,db0)*gsl_vector_get(wi,du);

                    dem+=temp;

                    gsl_vector_set(FUNU,j,gsl_vector_get(FUNU,j)+temp*gsl_vector_get(ti,p1a));
                    gsl_vector_set(FUNUS,j,gsl_vector_get(FUNUS,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a)));
                    
                    for(i=0;i<p1a;i++)  
                    {
                        gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)+temp*gsl_vector_get(ti,i));
                        gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i)));
                        gsl_matrix_set(FUNBU,i,j,gsl_matrix_get(FUNBU,i,j)+temp*gsl_vector_get(ti,i)*gsl_vector_get(ti,p1a));
                    }



                    for(i=1;i<p1a;i++)
                    {
                        for(t=0;t<p1a-i;t++)   gsl_matrix_set(FUNBS,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j)
                                               +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i));
                    }


                    gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)+temp*exp(gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNE,1,j,gsl_matrix_get(FUNE,1,j)+temp*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));

                    gsl_matrix_set(FUNUSE,0,j,gsl_matrix_get(FUNUSE,0,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a))*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNUE,0,j,gsl_matrix_get(FUNUE,0,j)+temp*gsl_vector_get(ti,p1a)*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));


                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+bq+1));
                        for(t=0;t<bq;t++)  gsl_vector_set(D,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        et=gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0);

                        c=(int)gsl_matrix_get(Y,m+i,0);
                        if(c==1)
                        {
                            pk=exp(gsl_vector_get(theta,c-1)-mu-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-et));
                            pk1=0;
                        }

                        else if(c==K_num)
                        {
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=1;
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }
                        else if(c==2)
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-et));
                        }

                        else
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }

                        pk_1=(pk-pk1);
                        if(pk_1==0) return 100;

                        gsl_matrix_set(FUNP,i*2,j,gsl_matrix_get(FUNP,i*2,j)+temp*pk);
                        gsl_matrix_set(FUNP,i*2+1,j,gsl_matrix_get(FUNP,i*2+1,j)+temp*pk1);

                        gsl_matrix_set(FUNPS,i*2,j,gsl_matrix_get(FUNPS,i*2,j)+temp*pk*pk);
                        gsl_matrix_set(FUNPS,i*2+1,j,gsl_matrix_get(FUNPS,i*2+1,j)+temp*pk1*pk1);

                        gsl_matrix_set(FUNPR1,i*2,j,gsl_matrix_get(FUNPR1,i*2,j)+temp*((pk-pk*pk)/pk_1));
                        gsl_matrix_set(FUNPR1,i*2+1,j,gsl_matrix_get(FUNPR1,i*2+1,j)+temp*((pk1-pk1*pk1)/pk_1));

                        gsl_matrix_set(FUNPR2,i*2,j,gsl_matrix_get(FUNPR2,i*2,j)+
                                       temp*((pk-pk*pk)*(1-2*pk)/pk_1-(pk-pk*pk)*(pk-pk*pk)/(pk_1*pk_1)));

                        gsl_matrix_set(FUNPR2,i*2+1,j,gsl_matrix_get(FUNPR2,i*2+1,j)+
                                       temp*((pk1-pk1*pk1)*(1-2*pk1)/(pk_1)+(pk1-pk1*pk1)*(pk1-pk1*pk1)/(pk_1*pk_1)));
                    }

                }
            }
        }
        }


        if(p1a==3)
        {

        for(db0=0;db0<point;db0++)
        {
            for(db1=0;db1<point;db1++)
            {
                for(db2=0;db2<point;db2++)
                {
                    for(du=0;du<point;du++)
                    {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,db1));
                    gsl_vector_set(ri,2,gsl_vector_get(xi,db2));
                    gsl_vector_set(ri,3,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);


                    temp=exp((double)10);

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+bq+1));
                        for(t=0;t<bq;t++)  gsl_vector_set(D,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        et=gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0);

                        c=(int)gsl_matrix_get(Y,m+i,0);
                        if(c==1)
                        {
                            pk=exp(gsl_vector_get(theta,c-1)-mu-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-et));
                            pk1=0;
                        }

                        else if(c==K_num)
                        {
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=1;
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }
                        else if(c==2)
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-et));
                        }

                        else
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }

                        temp*=(pk-pk1);
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=gsl_vector_get(wi,db0)*gsl_vector_get(wi,du);

                    dem+=temp;

                    gsl_vector_set(FUNU,j,gsl_vector_get(FUNU,j)+temp*gsl_vector_get(ti,p1a));
                    gsl_vector_set(FUNUS,j,gsl_vector_get(FUNUS,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a)));
                    
                    for(i=0;i<p1a;i++)  
                    {
                        gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)+temp*gsl_vector_get(ti,i));
                        gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i)));
                        gsl_matrix_set(FUNBU,i,j,gsl_matrix_get(FUNBU,i,j)+temp*gsl_vector_get(ti,i)*gsl_vector_get(ti,p1a));
                    }



                    for(i=1;i<p1a;i++)
                    {
                        for(t=0;t<p1a-i;t++)   gsl_matrix_set(FUNBS,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j)
                                               +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i));
                    }


                    gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)+temp*exp(gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNE,1,j,gsl_matrix_get(FUNE,1,j)+temp*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));

                    gsl_matrix_set(FUNUSE,0,j,gsl_matrix_get(FUNUSE,0,j)+temp*gsl_pow_2(gsl_vector_get(ti,p1a))*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    gsl_matrix_set(FUNUE,0,j,gsl_matrix_get(FUNUE,0,j)+temp*gsl_vector_get(ti,p1a)*exp(gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));


                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+bq+1));
                        for(t=0;t<bq;t++)  gsl_vector_set(D,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        et=gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0);

                        c=(int)gsl_matrix_get(Y,m+i,0);
                        if(c==1)
                        {
                            pk=exp(gsl_vector_get(theta,c-1)-mu-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-et));
                            pk1=0;
                        }

                        else if(c==K_num)
                        {
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=1;
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }
                        else if(c==2)
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-et));
                        }

                        else
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }

                        pk_1=(pk-pk1);
                        if(pk_1==0) return 100;

                        gsl_matrix_set(FUNP,i*2,j,gsl_matrix_get(FUNP,i*2,j)+temp*pk);
                        gsl_matrix_set(FUNP,i*2+1,j,gsl_matrix_get(FUNP,i*2+1,j)+temp*pk1);

                        gsl_matrix_set(FUNPS,i*2,j,gsl_matrix_get(FUNPS,i*2,j)+temp*pk*pk);
                        gsl_matrix_set(FUNPS,i*2+1,j,gsl_matrix_get(FUNPS,i*2+1,j)+temp*pk1*pk1);

                        gsl_matrix_set(FUNPR1,i*2,j,gsl_matrix_get(FUNPR1,i*2,j)+temp*((pk-pk*pk)/pk_1));
                        gsl_matrix_set(FUNPR1,i*2+1,j,gsl_matrix_get(FUNPR1,i*2+1,j)+temp*((pk1-pk1*pk1)/pk_1));

                        gsl_matrix_set(FUNPR2,i*2,j,gsl_matrix_get(FUNPR2,i*2,j)+
                                       temp*((pk-pk*pk)*(1-2*pk)/pk_1-(pk-pk*pk)*(pk-pk*pk)/(pk_1*pk_1)));

                        gsl_matrix_set(FUNPR2,i*2+1,j,gsl_matrix_get(FUNPR2,i*2+1,j)+
                                       temp*((pk1-pk1*pk1)*(1-2*pk1)/(pk_1)+(pk1-pk1*pk1)*(pk1-pk1*pk1)/(pk_1*pk_1)));
                    }

                }
            }
        }
        }
        }

        if(dem==0) return 100;


        gsl_vector_set(FUNU,j,gsl_vector_get(FUNU,j)/dem);
        gsl_vector_set(FUNUS,j,gsl_vector_get(FUNUS,j)/dem);

        for(i=0;i<p1a;i++)
        {
            gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)/dem);
            gsl_matrix_set(FUNBU,i,j,gsl_matrix_get(FUNBU,i,j)/dem);
        }

        for(i=0;i<p1a*(p1a+1)/2;i++)
        {
            gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)/dem);
        }


        gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)/dem);
        gsl_matrix_set(FUNE,1,j,gsl_matrix_get(FUNE,1,j)/dem);

        gsl_matrix_set(FUNUSE,0,j,gsl_matrix_get(FUNUSE,0,j)/dem);
        gsl_matrix_set(FUNUE,0,j,gsl_matrix_get(FUNUE,0,j)/dem);


        for(i=0;i<j_max;i++)
        {
            gsl_matrix_set(FUNP,2*i,j,gsl_matrix_get(FUNP,2*i,j)/dem);
            gsl_matrix_set(FUNP,2*i+1,j,gsl_matrix_get(FUNP,2*i+1,j)/dem);

            gsl_matrix_set(FUNPS,2*i,j,gsl_matrix_get(FUNPS,2*i,j)/dem);
            gsl_matrix_set(FUNPS,2*i+1,j,gsl_matrix_get(FUNPS,2*i+1,j)/dem);


            gsl_matrix_set(FUNPR1,2*i,j,gsl_matrix_get(FUNPR1,2*i,j)/dem);
            gsl_matrix_set(FUNPR1,2*i+1,j,gsl_matrix_get(FUNPR1,2*i+1,j)/dem);
       
            gsl_matrix_set(FUNPR2,2*i,j,gsl_matrix_get(FUNPR2,2*i,j)/dem);
            gsl_matrix_set(FUNPR2,2*i+1,j,gsl_matrix_get(FUNPR2,2*i+1,j)/dem);

        }


        m+=q;


    }

    gsl_vector_free(Z);
    gsl_vector_free(D);
    gsl_vector_free(X);
    gsl_vector_free(xtilde);
    gsl_vector_free(gammai);

    gsl_vector_free(xi);
    gsl_vector_free(ti);
    gsl_vector_free(ri);
    gsl_vector_free(wi);
    gsl_matrix_free(VC);  
    gsl_vector_free(S);
    gsl_matrix_free(V);
    gsl_matrix_free(VV);
    gsl_vector_free(W);  

    return 0;
}




 double Getloglik(
          const gsl_vector *beta,
          const gsl_matrix *beta2,
          const gsl_vector *theta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int bq,
		  const int K_num,
		  const int point,
		  const std::vector<double> xs,  
	      const std::vector<double> ws
          )
{
    int p1=beta->size;
    int p2=gamma->size2;
    
    int k = M1->size;

    int i,j,q,t,m,c;
    double temp,temp1,mu,mu1,mu2,et;
    double cuh01,cuh02,haz01,haz02,xgamma1,xgamma2;

    double pk,pk1;


    gsl_vector *Z = gsl_vector_calloc(p1),
               *D = gsl_vector_calloc(bq),
               *X = gsl_vector_calloc(p2),
               *xtilde = gsl_vector_calloc(p1a),
               *gammai = gsl_vector_calloc(p2);




    int db0,db1,db2,du;


    gsl_vector *xi = gsl_vector_calloc(point);
    gsl_vector *wi = gsl_vector_calloc(point);



    for(i=0;i<point/2;i++)   gsl_vector_set(xi,i,xs[i]);
    for(i=0;i<point/2;i++)   gsl_vector_set(wi,i,ws[i]);

    for(i=0;i<point/2;i++)   gsl_vector_set(xi,point-i-1, 0-gsl_vector_get(xi,i));
    for(i=0;i<point/2;i++)   gsl_vector_set(wi,point-i-1, gsl_vector_get(wi,i));



    gsl_matrix *VC=gsl_matrix_calloc(p1a+1,p1a+1);
    gsl_vector *S=gsl_vector_calloc(p1a+1);
    gsl_matrix *V=gsl_matrix_calloc(p1a+1,p1a+1);
    gsl_matrix *VV=gsl_matrix_calloc(p1a+1,p1a+1);
    gsl_vector *W=gsl_vector_calloc(p1a+1);


    gsl_matrix_memcpy(VC,sig);


    for(i=0;i<p1a+1;i++)
    {
        for(j=0;j<i;j++)    gsl_matrix_set(VC,i,j,gsl_matrix_get(VC,j,i));
    }



    gsl_linalg_SV_decomp (VC,V,S,W);
    gsl_matrix_set_zero(V);
    for(i=0;i<p1a+1;i++)  gsl_matrix_set(V,i,i,sqrt(gsl_vector_get(S,i)));
    MulMM(VC,V,VV);


    gsl_matrix_scale(VV,sqrt((double)2));


    gsl_vector *ri=gsl_vector_calloc(p1a+1);
    gsl_vector *ti=gsl_vector_calloc(p1a+1);

    double loglik=0;

    m=0;    
    for(j=0;j<k;j++)
    {
        q=(int)gsl_vector_get(M1,j);

        cuh01=CH(H01,gsl_matrix_get(C,j,0));
        cuh02=CH(H02,gsl_matrix_get(C,j,0));

        haz01=HAZ(H01,gsl_matrix_get(C,j,0));
        haz02=HAZ(H02,gsl_matrix_get(C,j,0));


        for(i=0;i<p2;i++) 
        {
            gsl_vector_set(X,i,gsl_matrix_get(C,j,2+i));
            gsl_vector_set(gammai,i,gsl_matrix_get(gamma,0,i)); 
        }
        xgamma1=MulVV(X,gammai);

        for(i=0;i<p2;i++) gsl_vector_set(gammai,i,gsl_matrix_get(gamma,1,i));
        xgamma2=MulVV(X,gammai);

        temp1=0;

        if(p1a==1)
        {

            for(db0=0;db0<point;db0++)
            {
                for(du=0;du<point;du++)
                {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);

                    temp=1;

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+bq+1));
                        for(t=0;t<bq;t++)  gsl_vector_set(D,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        et=gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0);

                        c=(int)gsl_matrix_get(Y,m+i,0);
                        if(c==1)
                        {
                            pk=exp(gsl_vector_get(theta,c-1)-mu-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-et));
                            pk1=0;
                        }

                        else if(c==K_num)
                        {
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=1;
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }
                        else if(c==2)
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-et));
                        }

                        else
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }

                        temp*=(pk-pk1);
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=1/sqrt(M_PI)*gsl_vector_get(wi,db0)*gsl_vector_get(wi,du);

                    temp1+=temp;

                }
            }
        }




        if(p1a==2)
        {

            for(db0=0;db0<point;db0++)
            {
                for(db1=0;db1<point;db1++)
                {
                    for(du=0;du<point;du++)
                    {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,db1));
                    gsl_vector_set(ri,2,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);


                    temp=1;

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+bq+1));
                        for(t=0;t<bq;t++)  gsl_vector_set(D,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        et=gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0);

                        c=(int)gsl_matrix_get(Y,m+i,0);
                        if(c==1)
                        {
                            pk=exp(gsl_vector_get(theta,c-1)-mu-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-et));
                            pk1=0;
                        }

                        else if(c==K_num)
                        {
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=1;
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }
                        else if(c==2)
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-et));
                        }

                        else
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }

                        temp*=(pk-pk1);
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=1/sqrt(M_PI)*gsl_vector_get(wi,db0)*gsl_vector_get(wi,du);

                    temp1+=temp;

                }
            }
        }
        }


        if(p1a==3)
        {

        for(db0=0;db0<point;db0++)
        {
            for(db1=0;db1<point;db1++)
            {
                for(db2=0;db2<point;db2++)
                {
                    for(du=0;du<point;du++)
                    {

                    gsl_vector_set(ri,0,gsl_vector_get(xi,db0));
                    gsl_vector_set(ri,1,gsl_vector_get(xi,db1));
                    gsl_vector_set(ri,2,gsl_vector_get(xi,db2));
                    gsl_vector_set(ri,3,gsl_vector_get(xi,du));

                    MulM(VV,ri,ti);


                    temp=1;

                    for(i=0;i<q;i++)
                    {
                        for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+bq+1));
                        for(t=0;t<bq;t++)  gsl_vector_set(D,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                        for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                        mu=MulVV(Z,beta);
                        et=gsl_vector_get(ti,0)*gsl_vector_get(xtilde,0);

                        c=(int)gsl_matrix_get(Y,m+i,0);
                        if(c==1)
                        {
                            pk=exp(gsl_vector_get(theta,c-1)-mu-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-et));
                            pk1=0;
                        }

                        else if(c==K_num)
                        {
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=1;
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }
                        else if(c==2)
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-et));
                        }

                        else
                        {
                            mu1=0;
                            for(t=0;t<bq;t++)    mu1+=gsl_matrix_get(beta2,c-2,t)*gsl_vector_get(D,t);
                            mu2=0;
                            for(t=0;t<bq;t++)    mu2+=gsl_matrix_get(beta2,c-3,t)*gsl_vector_get(D,t);

                            pk=exp(gsl_vector_get(theta,c-1)-mu-mu1-et)/(1+exp(gsl_vector_get(theta,c-1)-mu-mu1-et));
                            pk1=exp(gsl_vector_get(theta,c-2)-mu-mu2-et)/(1+exp(gsl_vector_get(theta,c-2)-mu-mu2-et));
                        }

                        temp*=(pk-pk1);
                    }

                    if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+gsl_vector_get(ti,p1a));
                    if(gsl_matrix_get(C,j,1)==2)  temp*=haz02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a));

                    temp*=exp(0-cuh01*exp(xgamma1+gsl_vector_get(ti,p1a))-cuh02*exp(xgamma2+gsl_vector_get(vee,0)*gsl_vector_get(ti,p1a)));
                    temp*=1/sqrt(M_PI)*gsl_vector_get(wi,db0)*gsl_vector_get(wi,du);

                    temp1+=temp;

                }
            }
        }
        }
        }

        if(temp1==0) return 100;
 
        loglik+=log(temp1);

        m+=q;


    }

    gsl_vector_free(Z);
    gsl_vector_free(X);
    gsl_vector_free(D);
    gsl_vector_free(xtilde);
    gsl_vector_free(gammai);

    gsl_vector_free(xi);
    gsl_vector_free(ti);
    gsl_vector_free(ri);
    gsl_vector_free(wi);
    gsl_matrix_free(VC);  
    gsl_vector_free(S);
    gsl_matrix_free(V);
    gsl_matrix_free(VV);
    gsl_vector_free(W);  

    return loglik;
}







 double DiffV(const gsl_vector *veca, const gsl_vector *vecb)
{
    int k=veca->size;
    int i;
    double diff=0;

    for(i=0;i<k;i++)
    {
        if(Abs(gsl_vector_get(veca,i),gsl_vector_get(vecb,i))>diff) diff=Abs(gsl_vector_get(veca,i),gsl_vector_get(vecb,i));
    }

    return (diff);   
}


 double DiffM(const gsl_matrix *matrixa, const gsl_matrix *matrixb)
{
    int nrow=matrixa->size1, ncol=matrixa->size2;
    int i, j;
    double diff=0;

    for(i=0;i<nrow;i++)
    {
        for(j=0;j<ncol;j++)
        {
            if(Abs(gsl_matrix_get(matrixa,i,j),gsl_matrix_get(matrixb,i,j))>diff) 
               diff=Abs(gsl_matrix_get(matrixa,i,j),gsl_matrix_get(matrixb,i,j));
        }
    }

    return (diff);
}



 double Abs(const double a, const double b)
{
     
    if (a>=b) return a-b;
    else return b-a;
}



 int inv_matrix(gsl_matrix *x_square)
{
   int i,j;
   int k = x_square->size1;

   int status;

   gsl_vector *temp_vector=gsl_vector_calloc(k),
              *solution=gsl_vector_calloc(k);
   gsl_matrix *out = gsl_matrix_calloc(k,k);

   for(i=0;i<k;i++)
   {
       for(j=0;j<k;j++) gsl_matrix_set(out,i,j,gsl_matrix_get(x_square,i,j));
   }

   status=gsl_linalg_cholesky_decomp(out);
   if(status==100) return status;

   for (i = 0; i < k; i++)
   {
       gsl_vector_set_all(temp_vector,0);
       gsl_vector_set(temp_vector,i,1);

       status=gsl_linalg_cholesky_solve(out, temp_vector, solution);
       if(status==100) return status;

       gsl_matrix_set_col(x_square,i,solution);
   }
   
   gsl_vector_free(temp_vector);
   gsl_vector_free(solution);
   gsl_matrix_free(out);

   return 0;

}



 void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta)
{
    int p = XX->size1;
    int q = XX->size2;

    int i,j;
    double temp;

    for(i=0;i<p;i++)
    {
        temp=0;
        for(j=0;j<q;j++)  temp+=gsl_matrix_get(XX,i,j)*gsl_vector_get(X,j);
        gsl_vector_set(beta,i,temp);
    }

}
         

 void MulV(const gsl_vector *Z,gsl_matrix *ZZ)
{
    int p = Z->size;
    int i,j;

    for(i=0;i<p;i++)
    {
        for(j=0;j<p;j++) gsl_matrix_set(ZZ,i,j,gsl_vector_get(Z,i)*gsl_vector_get(Z,j));
    }
}
   

 double MulVV(const gsl_vector *Z,const gsl_vector *beta)
{
    int p=Z->size;
    int i;
    double temp=0;

    for(i=0;i<p;i++)  temp+=gsl_vector_get(Z,i)*gsl_vector_get(beta,i);

    return (temp);
}


 void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB)
{
    int p=A->size1;
    int q=A->size2;
    int k=B->size2;

    int i,j,t;
    double temp;

    for(i=0;i<p;i++)  
    {
        for(j=0;j<k;j++)
        {
            temp=0;
            for(t=0;t<q;t++)  temp+=gsl_matrix_get(A,i,t)*gsl_matrix_get(B,t,j);
            gsl_matrix_set(AB,i,j,temp);
        }
    }

}


 double CH(const gsl_matrix *H, double t)
{
    int a=H->size2;
    int i;

    double ch;

    if(t<gsl_matrix_get(H,0,0)) ch=0;
    else
    {
        ch=0;
        i=0;
        do
        {
            ch+=gsl_matrix_get(H,2,i);
            i+=1;
        }while(i<=a-1 && t>=gsl_matrix_get(H,0,i));
    }

    return (ch);
}


 double HAZ(const gsl_matrix *H, double t)
{
    int a=H->size2;
    int i;
    double temp=0;

    for(i=0;i<a;i++)
    {
        if(t==gsl_matrix_get(H,0,i)) temp=gsl_matrix_get(H,2,i);
    }

    return (temp);
}       




 double Min(const double t1, const double t2)
{
    if(t1<t2) return t1;
    else return t2;
}


 int Diff(
         const gsl_vector *prebeta,
         const gsl_vector *beta,
         const gsl_matrix *prebeta2,
         const gsl_matrix *beta2,
         const gsl_vector *pretheta,
         const gsl_vector *theta,
         const gsl_matrix *pregamma,
         const gsl_matrix *gamma,
         const gsl_vector *prevee,
         const gsl_vector *vee,
         const gsl_matrix *preH01,
         const gsl_matrix *H01,
         const gsl_matrix *preH02,
         const gsl_matrix *H02,
         const gsl_matrix *presig,
         const gsl_matrix *sig
         )
{    

     double epsilon=0.0001;  
        
     if(DiffV(prebeta,beta)>epsilon || DiffM(prebeta2,beta2)>epsilon || DiffV(pretheta,theta)>epsilon || DiffM(pregamma,gamma)>epsilon 
        || DiffV(prevee,vee)>epsilon || DiffM1(preH01,H01)==1 || DiffM1(preH02,H02)==1 || DiffM(presig,sig)>epsilon)
 
     return 1;

     else return 0;

}





 int DiffM1(const gsl_matrix *matrixa, const gsl_matrix *matrixb)
{
    int nrow=matrixa->size1, ncol=matrixa->size2;
    int i, j;
    int diff=0;

    i=nrow-1;

        for(j=0;j<ncol;j++)
        {
            if(gsl_matrix_get(matrixa,i,j)-gsl_matrix_get(matrixb,i,j)>0.001 || gsl_matrix_get(matrixa,i,j)-gsl_matrix_get(matrixb,i,j)<-0.001)

                diff=1;

        }

    return diff;
}


 int DiffM2(const gsl_matrix *preH1,const gsl_matrix *H1,const gsl_matrix *preH2,const gsl_matrix *H2)
{
     
     if(DiffM1(preH1,H1)==1 || DiffM1(preH2,H2)==1 ) return 1;
     else return 0;
}


 void TransM(const gsl_matrix *A, gsl_matrix *B)
{
    int rowa = A->size1;
    int cola = A->size2;

    int i, j;

    for(i=0;i<rowa;i++)
    {
        for(j=0;j<cola;j++)
        {
            gsl_matrix_set(B,j,i,gsl_matrix_get(A,i,j));
        }
    }

}   


 void STAT(gsl_matrix *store,int i,double *mean,double *sd)
{
    int n=store->size1;
    int j;

    *mean=0, *sd=0;
    
    for(j=0;j<n;j++)   *mean+=gsl_matrix_get(store,j,i);
    *mean=*mean/(double)n;

    for(j=0;j<n;j++)   *sd+=gsl_pow_2(gsl_matrix_get(store,j,i)-*mean);
    *sd = sqrt(*sd/(double)(n-1));

}
        
     

  int GetN(double t)
{

    return (int)(t/0.5+1);

}

Rcpp::List jmo_cmain(int k, int n1,int p1,int p2, int p1a, int bq,int K_num, int j_max, int point, std::vector<double>xs,  std::vector<double> ws,std::vector<double>beta_val, std::vector<double>theta_val, int maxiterations,std::string yfile, std::string cfile, std::string mfile, int trace)
{ 
    int i,j,t,g=2;  
	
/*
    k=586; n1=1993; p1=5; p2=2; p1a=1; bq=3; 

*/
    
    /* allocate space for data */

    gsl_matrix *C = gsl_matrix_calloc(k,p2+2);
    gsl_matrix *Y= gsl_matrix_calloc(n1, p1+p1a+bq+1);    
    gsl_vector *M1= gsl_vector_calloc(k);

	/* read Y matrix  */
  {  			
    FILE * f = fopen(yfile.c_str(), "r");  
    
	if (f == NULL)
    {
        Rprintf("File %s does not exist.\n", yfile.c_str());
      return R_NilValue;	  
    }
	
  
    int nrows=0;
      // Extract characters from file and store in character c
	for (char c = fgetc(f); c != EOF; c = fgetc(f))
			if (c == '\n')  nrows = nrows + 1;
	nrows=nrows+1;
	if (n1==nrows)
	{   rewind(f);
		gsl_matrix_fscanf(f, Y);  
		fclose(f);
	}
	else
	{ 
	    Rprintf("Input oberservations is %d, but the number of rows in %s is %d",n1,yfile.c_str(),nrows);
		fclose(f);
		return R_NilValue;
	}	 
 

    
  } 
  /* read C matrix  */
  {  
    FILE * f = fopen(cfile.c_str(), "r");
	
	if (f == NULL)
    {
        Rprintf("File %s does not exist.\n", cfile.c_str());
      return R_NilValue;	  
    }
	
  
    int nrows=0;
      // Extract characters from file and store in character c
	for (char c = fgetc(f); c != EOF; c = fgetc(f))
			if (c == '\n')  nrows = nrows + 1;
	nrows=nrows+1;
	if (k==nrows)
	{   rewind(f);
	    gsl_matrix_fscanf(f, C);
		fclose(f);
	}
	else
	{ 
	    Rprintf("Input subjects is %d, but the number of rows in %s is %d",k,cfile.c_str(),nrows);
		fclose(f);
		return R_NilValue;
	}	  
	    
    
  }
  /* read M1 vector  */
  {  
    FILE * f = fopen(mfile.c_str(), "r");
    
	if (f == NULL)
    {
        Rprintf("File %s does not exist.\n", mfile.c_str());
      return R_NilValue;	  
    }
	
  
    int nrows=0;
      // Extract characters from file and store in character c
	for (char c = fgetc(f); c != EOF; c = fgetc(f))
			if (c == '\n')  nrows = nrows + 1;
	nrows=nrows+1;
	if (k==nrows)
	{   rewind(f);
	    gsl_vector_fscanf(f, M1);
		fclose(f);
	}
	else
	{ 
	    Rprintf("Input subjects is %d, but the number of rows in %s is %d",k,mfile.c_str(),nrows);
		fclose(f);
		return R_NilValue;
	}	 
		
    
  }


    gsl_matrix * Cov=gsl_matrix_calloc(p1+(K_num-2)*bq+K_num+(p2+1)*g+(p1a+1)*(p1a+2)/2-2,p1+(K_num-2)*bq+K_num+(p2+1)*g+(p1a+1)*(p1a+2)/2-2);

/* switch the two competing risks (e.g., risk 1 becomes risk 2).
    for(i=0;i<k;i++)
    {
        if(gsl_matrix_get(C,i,1)>0)   gsl_matrix_set(C,i,1,3-gsl_matrix_get(C,i,1));
    }
*/
    /* allocate space for estimated parameters */    
    gsl_matrix * gamma=gsl_matrix_calloc(g, p2);                /* coefficient for covariates in competing risks*/
    gsl_vector * vee=gsl_vector_calloc(g-1);                    /* coefficient for random effect in competing risks */
    gsl_vector * beta=gsl_vector_calloc(p1);                    /* coefficient for longitudinal outcome  */ 
    gsl_matrix * beta2=gsl_matrix_calloc(K_num-2,bq); 
    gsl_vector * theta=gsl_vector_calloc(K_num-1);   
    gsl_matrix *sig = gsl_matrix_calloc(p1a+1,p1a+1);




    /* allocate space for pre parameters */            
    
    gsl_vector * prebeta=gsl_vector_calloc(p1);  
    gsl_vector * pretheta=gsl_vector_calloc(K_num-1);   
    gsl_matrix * prebeta2=gsl_matrix_calloc(K_num-2,bq);             
    gsl_matrix * pregamma=gsl_matrix_calloc(g, p2);                
    gsl_vector * prevee=gsl_vector_calloc(g-1);                                          

    gsl_matrix *presig = gsl_matrix_calloc(p1a+1,p1a+1);



    /* allocate space to store parameters */            
    
    gsl_vector * vbeta=gsl_vector_calloc(p1);  
    gsl_vector * vbeta2=gsl_vector_calloc((K_num-2)*bq);   
    gsl_vector * vtheta=gsl_vector_calloc(K_num-1);            
    gsl_matrix * vgamma=gsl_matrix_calloc(g, p2);                
    gsl_vector * vvee=gsl_vector_calloc(g-1); 

    gsl_vector *vsig = gsl_vector_calloc((p1a+1)*(p1a+2)/2);


/*     FILE *output_F; 

    output_F=fopen("output.dat","w");
              if (output_F == NULL)
              {Rprintf("Can't open file\n");
              } */

    int iter,status;

    


    gsl_matrix * FH01 = gsl_matrix_calloc(2,array_size);                 /* fake baseline hazard function for competing risk 1 */
    gsl_matrix * FH02 = gsl_matrix_calloc(2,array_size);                 /* fake baseline hazard function for competing risk 2 */
    
    gsl_matrix_set_zero(FH01);
    gsl_matrix_set_zero(FH02);


    /* find # events for risk 1 */
    int u, v, a=0, b=0;

        for(j=0;j<k;j++)
        {
            u=0;
            if(gsl_matrix_get(C,j,1)==1) 
            {
                if(a>=1)
                {
                    for(t=0;t<a;t++)
                    {
                        if(gsl_matrix_get(FH01,0,t)==gsl_matrix_get(C,j,0))
                        {
                            gsl_matrix_set(FH01,1,t,gsl_matrix_get(FH01,1,t)+1);
                            u=1;
                        }
                    }
                    if(u==0)
                    {
                        t=-1;

                        do
                        {
                            t=t+1;
                        }while(t<=a-1 && gsl_matrix_get(FH01,0,t) < gsl_matrix_get(C,j,0));

                        if(t==a)
                        {
                            gsl_matrix_set(FH01,0,t,gsl_matrix_get(C,j,0));
                            gsl_matrix_set(FH01,1,t,1);
                        }

                        if(t<a)
                        {
                            for(v=a-1;v>=t;v--)
                            {
                                gsl_matrix_set(FH01,0,v+1,gsl_matrix_get(FH01,0,v));
                                gsl_matrix_set(FH01,1,v+1,gsl_matrix_get(FH01,1,v));
                            }
                                gsl_matrix_set(FH01,0,t,gsl_matrix_get(C,j,0));
                                gsl_matrix_set(FH01,1,t,1);
                        }
                      
                        a=a+1;
                    }
                }
               
                if(a==0)
                {
                    gsl_matrix_set(FH01,0,0,gsl_matrix_get(C,j,0));
                    gsl_matrix_set(FH01,1,0,1);
                    a=a+1;
                }
            }
        }


    if(a==0) 
    {
        Rprintf("No failure time information for risk 1; Program exits\n");
        return R_NilValue;
    }
     

    gsl_matrix * H01 = gsl_matrix_calloc(3,a);                 /* baseline hazard function for competing risk 1 */

    for(i=0;i<3;i++)
    {
        if(i<=1)
        {
            for(j=0;j<a;j++)    gsl_matrix_set(H01,i,j, gsl_matrix_get(FH01,i,j));
        }
        if(i==2)   
        {
            for(j=0;j<a;j++)    gsl_matrix_set(H01,i,j,0.0001);
        }
    }
                                     

     
    /* find # events for risk 2 */

        for(j=0;j<k;j++)
        {
            u=0;
            if(gsl_matrix_get(C,j,1)==2) 
            {
                if(b>=1)
                {
                    for(t=0;t<b;t++)
                    {
                        if(gsl_matrix_get(FH02,0,t)==gsl_matrix_get(C,j,0))
                        {
                            gsl_matrix_set(FH02,1,t,gsl_matrix_get(FH02,1,t)+1);
                            u=1;
                        }
                    }
                    if(u==0)
                    {
                        t=-1;

                        do
                        {
                            t=t+1;
                        }while(t<=b-1 && gsl_matrix_get(FH02,0,t) < gsl_matrix_get(C,j,0));

                        if(t==b)
                        {
                            gsl_matrix_set(FH02,0,t,gsl_matrix_get(C,j,0));
                            gsl_matrix_set(FH02,1,t,1);
                        }

                        if(t<b)
                        {
                            for(v=b-1;v>=t;v--)
                            {
                                gsl_matrix_set(FH02,0,v+1,gsl_matrix_get(FH02,0,v));
                                gsl_matrix_set(FH02,1,v+1,gsl_matrix_get(FH02,1,v));
                            }
                                gsl_matrix_set(FH02,0,t,gsl_matrix_get(C,j,0));
                                gsl_matrix_set(FH02,1,t,1);
                        }
                        b=b+1;
                    }
                }
               
                if(b==0)
                {
                    gsl_matrix_set(FH02,0,0,gsl_matrix_get(C,j,0));
                    gsl_matrix_set(FH02,1,0,1);
                    b=b+1;
                }
            }
        }
 

    if(b==0) 
    {
        Rprintf("No failure time information for risk 2; Program exits\n");
        return R_NilValue;
    }


    gsl_matrix * H02 = gsl_matrix_calloc(3,b);                 /* baseline hazard function for competing risk 2 */

    for(i=0;i<3;i++)
    {
        if(i<=1)
        {
            for(j=0;j<b;j++)    gsl_matrix_set(H02,i,j, gsl_matrix_get(FH02,i,j));
        }
        if(i==2)   
        {
            for(j=0;j<b;j++)    gsl_matrix_set(H02,i,j,0.0001);
        }
    }


    gsl_matrix * preH01 = gsl_matrix_calloc(3,a);               
    gsl_matrix * preH02 = gsl_matrix_calloc(3,b);   




    /* initialize the parameters */

    gsl_matrix_set_zero(gamma);
    gsl_vector_set_zero(vee);
    gsl_matrix_set_zero(beta2);
    gsl_vector_set_zero(theta);

    gsl_matrix_set_identity(sig); 

    for (size_t i=0;i<beta_val.size();i++) 
           gsl_vector_set(beta,i,beta_val[i]);

    for (size_t i=0;i<theta_val.size();i++) 
           gsl_vector_set(theta,i,theta_val[i]);   
 /* gsl_vector_set(beta,0,-1.13);
    gsl_vector_set(beta,1,1.56);
    gsl_vector_set(beta,2,-1.28);
    gsl_vector_set(beta,3,-1.51);
    gsl_vector_set(beta,4,-1.80);

    gsl_vector_set(theta,0,-1);
    gsl_vector_set(theta,1,0);
    gsl_vector_set(theta,2,1);
*/

    

    iter=0;
    do
    {
		
		//check user interrupt just in case 
        Rcpp::checkUserInterrupt();
        iter=iter+1;

        /* store the pre-information */

        gsl_vector_memcpy(prebeta, beta);
        gsl_matrix_memcpy(prebeta2, beta2);
        gsl_vector_memcpy(pretheta, theta);
        gsl_vector_memcpy(prevee, vee);  
        gsl_matrix_memcpy(pregamma, gamma);
        gsl_matrix_memcpy(preH01, H01);
        gsl_matrix_memcpy(preH02, H02); 
        gsl_matrix_memcpy(presig, sig); 

 
        
        /* get new parameter estimates */


        status = EM(beta,beta2,theta,gamma,vee,H01,H02,sig,Y,C,M1,p1a,bq, K_num, j_max,point,xs,ws);
		if(trace==1){	
		
           Rprintf("iter=%d   status=%d\n",iter,status);


            Rprintf("Beta = \n");
            for (i=0;i<p1;i++)
            {
                Rprintf("%f     ", gsl_vector_get(beta,i)); 
                      
            }
            Rprintf("\n");



            Rprintf("Alpha = \n");
            for (i=0;i<K_num-2;i++)
            {
                for(j=0;j<bq;j++)
                {
                    Rprintf("%f    ", gsl_matrix_get(beta2,i,j));
                }
                      
            }
            Rprintf("\n");


            Rprintf("Theta = \n");
            for (i=0;i<K_num-1;i++)
            {
                Rprintf("%f     ", gsl_vector_get(theta,i)); 
                      
            }
            Rprintf("\n");


            Rprintf("Gamma = \n");
            for (i=0;i<g;i++)
            {
                for(j=0;j<p2;j++)
                {
                    Rprintf("%f    ", gsl_matrix_get(gamma,i,j));
                }
                Rprintf("\n");
            }
    
            Rprintf("Vee = \n");
            for (i=0;i<g-1;i++)
            {
                Rprintf("%f     ", gsl_vector_get(vee,i)); 
                      
            }
            Rprintf("\n"); 



            Rprintf("Sig = \n");
            for (i=0;i<p1a+1;i++)
            {
                for(j=0;j<p1a+1;j++)
                {
                    Rprintf("%f    ", gsl_matrix_get(sig,i,j));
                }
                Rprintf("\n");
            }
		}


    }while(Diff(prebeta,beta,prebeta2,beta2,pretheta,theta,pregamma,gamma,prevee,vee,preH01,H01,preH02,H02,presig,sig)==1
           && status != 100 && iter<100000);

    if(status==100) 
    {
        Rprintf("program stops because of error\n");
        return R_NilValue;
    }
    

    if(iter==100000) 
    {
        Rprintf("program stops because of nonconvergence\n");
        return R_NilValue;
    }

  NumericMatrix vcmatrix(Cov->size1,Cov->size1) ;
  NumericVector betas(p1);
  NumericMatrix beta2matrix(K_num-2,bq);
  NumericVector thetas(K_num-1);
  NumericMatrix gamma_matrix(g,p2);
  NumericVector v_estimate(g-1);
  NumericMatrix sigma_matrix(p1a+1,p1a+1);  
  NumericVector se_betas(p1);
  NumericVector se_betas2((K_num-2)*bq);
  NumericVector se_thetas(K_num-1);
  NumericMatrix sd_gamma_matrix(g,p2); 
  NumericVector sd_v_estimate(g-1);
  NumericVector sd_sigma((p1a+1)*(p1a+2)/2);
  double loglike=0.0;
  


    if(status != 100 && iter<100000)
    {

        /* if algorithm coverges, compute the variance-covariance matrix of parameters ***/

        status = GetCov(Cov,beta,beta2,theta,gamma,vee,H01,H02,sig,Y,C,M1,p1a,bq, K_num,j_max,point,xs,ws);

        if(status==100) 
        {
            Rprintf("program stops because of error\n");
            return R_NilValue;
        }


        if(status != 100)
        {
            status=inv_matrix(Cov);

            if(status==100) 
            {
                Rprintf("program stops because of error\n");
                return R_NilValue;
            }

            if(status!=100)
            {

           // fRprintf(output_F,"Variance-covariance matrix for all the parameters: \n");
            for (i=0;(unsigned)i<Cov->size1;i++)
            {
                for(j=0;(unsigned)j<Cov->size1;j++)
                {
                   // fRprintf(output_F,"%f    ", gsl_matrix_get(Cov,i,j));
				   vcmatrix(i,j)=gsl_matrix_get(Cov,i,j);
                }
               // fRprintf(output_F,"\n");
            }

            for (i=0;i<p1;i++)  gsl_vector_set(vbeta,i,gsl_matrix_get(Cov,i,i));
            for (i=0;i<(K_num-2)*bq;i++) gsl_vector_set(vbeta2,i,gsl_matrix_get(Cov,p1+i,p1+i));
            for (i=0;i<K_num-1;i++) gsl_vector_set(vtheta,i,gsl_matrix_get(Cov,p1+(K_num-2)*bq+i,p1+(K_num-2)*bq+i));

            for (i=0;i<p2;i++)  gsl_matrix_set(vgamma,0,i,gsl_matrix_get(Cov,p1+(K_num-2)*bq+K_num+i-1,p1+(K_num-2)*bq+K_num+i-1));
            for (i=0;i<p2;i++)  gsl_matrix_set(vgamma,1,i,gsl_matrix_get(Cov,p1+(K_num-2)*bq+K_num+p2+i-1,p1+(K_num-2)*bq+K_num+p2+i-1));

            for (i=0;i<g-1;i++)  gsl_vector_set(vvee,i,gsl_matrix_get(Cov,p1+(K_num-2)*bq+K_num+2*p2+i-1,p1+(K_num-2)*bq+K_num+2*p2+i-1));

            for (i=0;i<(p1a+1)*(p1a+2)/2;i++)  gsl_vector_set(vsig,i,gsl_matrix_get(Cov,p1+(K_num-2)*bq+K_num+2*p2+g+i-2,p1+(K_num-2)*bq+K_num+2*p2+g+i-2));


          //  fRprintf(output_F,"Beta = \n");
            for (i=0;i<p1;i++)
            {
               // fRprintf(output_F,"%f     ", gsl_vector_get(beta,i));
               betas(i)=gsl_vector_get(beta,i);
            }
          //  fRprintf(output_F,"\n");


           // fRprintf(output_F,"Beta2 = \n");
            for (i=0;i<K_num-2;i++)
            {
                for(j=0;j<bq;j++)
                {
                   // fRprintf(output_F,"%f    ", gsl_matrix_get(beta2,i,j));
				    beta2matrix(i,j)=gsl_matrix_get(beta2,i,j);
                }  
               // fRprintf(output_F,"\n");
            }

         //   fRprintf(output_F,"Theta = \n");
            for (i=0;i<K_num-1;i++)
            {
               // fRprintf(output_F,"%f     ", gsl_vector_get(theta,i));
			   thetas(i)=gsl_vector_get(theta,i);

            }
         //   fRprintf(output_F,"\n");



          //  fRprintf(output_F,"Gamma = \n");
            for (i=0;i<g;i++)
            {
                for(j=0;j<p2;j++)
                {
                   // fRprintf(output_F,"%f    ", gsl_matrix_get(gamma,i,j));
				gamma_matrix(i,j) =gsl_matrix_get(gamma,i,j);  
                }  
             //   fRprintf(output_F,"\n");
            }

          //  fRprintf(output_F,"V = \n");
            for (i=0;i<g-1;i++)
            {
               // fRprintf(output_F,"%f     ", gsl_vector_get(vee,i));
			   v_estimate(i)=gsl_vector_get(vee,i);
 
            }
            //fRprintf(output_F,"\n");


          //  fRprintf(output_F,"Sigma = \n");
            for (i=0;i<p1a+1;i++)  
            {
                for (j=0;j<p1a+1;j++)  
				{//fRprintf(output_F,"%f     ", gsl_matrix_get(sig,i,j));
				sigma_matrix(i,j)=gsl_matrix_get(sig,i,j);
				}
				
               // fRprintf(output_F,"\n");
            }

/*

            fRprintf(output_F,"H01 = \n");
            for (i=0;i<a;i++)
            {
                for(j=0;j<3;j++)
                {
                    fRprintf(output_F,"%f    ", gsl_matrix_get(H01,j,i));
                }
                fRprintf(output_F,"\n");
            }

            fRprintf(output_F,"H02 = \n");
            for (i=0;i<b;i++)
            {
                for(j=0;j<3;j++)
                {
                    fRprintf(output_F,"%f    ", gsl_matrix_get(H02,j,i));
                }
                fRprintf(output_F,"\n");
            }

*/


            /*** output variance estimates ***/

     //       fRprintf(output_F,"\n***************************************\n");

       //     fRprintf(output_F,"standard error of beta = \n");
            for (i=0;i<p1;i++)
            {
               // fRprintf(output_F,"%f\n", sqrt(gsl_vector_get(vbeta,i)));
			   se_betas(i)=sqrt(gsl_vector_get(vbeta,i));

            }


          //  fRprintf(output_F,"standard error of beta2 = \n");
            for (i=0;i<(K_num-2)*bq;i++)
            {
               // fRprintf(output_F,"%f\n", sqrt(gsl_vector_get(vbeta2,i)));
			   se_betas2(i)=sqrt(gsl_vector_get(vbeta2,i));

            }


           // fRprintf(output_F,"standard error of theta = \n");
            for (i=0;i<K_num-1;i++)
            {
              //  fRprintf(output_F,"%f\n", sqrt(gsl_vector_get(vtheta,i)));
                  se_thetas(i)=sqrt(gsl_vector_get(vtheta,i));
            }



          //  fRprintf(output_F,"standard error of gamma = \n");
            for (i=0;i<g;i++)
            {
                for(j=0;j<p2;j++)
                {
                   // fRprintf(output_F,"%f    ", sqrt(gsl_matrix_get(vgamma,i,j)));
				   sd_gamma_matrix(i,j)=sqrt(gsl_matrix_get(vgamma,i,j));
				   
                }  
              //  fRprintf(output_F,"\n");
            }

          //  fRprintf(output_F,"standard error of V = \n");
            for (i=0;i<g-1;i++)
            {
               // fRprintf(output_F,"%f     ", sqrt(gsl_vector_get(vvee,i)));
                sd_v_estimate(i)=sqrt(gsl_vector_get(vvee,i));
            }
          //  fRprintf(output_F,"\n");


           // fRprintf(output_F,"standard error of Sigma = \n");
            for (i=0;i<(p1a+1)*(p1a+2)/2;i++)  
			{
				//fRprintf(output_F,"%f\n", sqrt(gsl_vector_get(vsig,i)));
				sd_sigma(i)=sqrt(gsl_vector_get(vsig,i));
		    }



            loglike=Getloglik(beta,beta2,theta,gamma,vee,H01,H02,sig,Y,C,M1,p1a,bq, K_num,point,xs,ws);

           // fRprintf(output_F,"loglikelihood = %f\n",loglike);

            }
        }      
    }

    gsl_matrix_free(FH01);
    gsl_matrix_free(FH02);
    gsl_matrix_free(H01);
    gsl_matrix_free(H02);
    gsl_matrix_free(preH01);
    gsl_matrix_free(preH02);

    gsl_matrix_free(Y);
    gsl_matrix_free(C);
    gsl_vector_free(M1);

    gsl_matrix_free(gamma);
    gsl_vector_free(beta);
    gsl_matrix_free(beta2);
    gsl_vector_free(theta);
    gsl_vector_free(vee);

    gsl_matrix_free(pregamma);
    gsl_vector_free(prebeta);
    gsl_matrix_free(prebeta2);
    gsl_vector_free(prevee);


    gsl_matrix_free(vgamma);
    gsl_vector_free(vbeta);
    gsl_vector_free(vbeta2);
    gsl_vector_free(vtheta);
    gsl_vector_free(vvee);
    gsl_matrix_free(sig);
    gsl_matrix_free(presig);
    gsl_vector_free(vsig);
    gsl_matrix_free(Cov); 
	
  Rcpp::List ret;
   
  ret["vcmatrix"] = vcmatrix;
  ret["betas"] = betas; 
  ret["se_betas"] = se_betas; 
  ret["alphamatrix"] = beta2matrix;
  ret["se_alphas"] = se_betas2; 
  ret["thetas"] = thetas;
  ret["se_thetas"] = se_thetas;  
  ret["gamma_matrix"] = gamma_matrix; 
  ret["se_gamma_matrix"] = sd_gamma_matrix; 
  ret["v_estimate"] = v_estimate; 
  ret["se_v_estimate"] = sd_v_estimate;   
  ret["sigma_matrix"] = sigma_matrix;  
  ret["se_sigma"] = sd_sigma; 
  ret["loglike"] = loglike;
  ret["Nlevels"] = K_num;

   return ret;

}


}
