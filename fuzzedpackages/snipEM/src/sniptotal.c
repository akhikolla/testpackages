#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

#include "sniptotal.h"

double **mtxalloc(int n, int p);
int **Imtxalloc(int n, int p);
void mtxfree(double **a, int n, int p);
void Imtxfree(int **a, int n, int p);

void copymtx(double **b, double **a, int n, int p);
void copyImtx(int **b, int **a, int n, int p);
void copyvec(double *b, double *a, int n);
void copyIvec(int *b, int *a, int n);
void zeroIvec(int *a, int n);

double square(double x);
double loss(double **x, int n, int d, int k, double **mu, int *clust, int *s, int **V);
void estmea(double **x, int *s, int **V, int *clust, double **mu, int d, int n, int k);
void estmemb(double **x, int *clust, int k, int n, int d, double **mu, int **V);
void sniptotal(double *xin, int *kin, int *itermaxin, int *Vin, int *s, int *clust, double *Din, int *nin, int *din, double *muin);


/*
 * Compute Snipped k-means clustering as described in 
 * Farcomeni, A. (2014), "Snipping for robust k-means clustering 
 * under component-wise contamination".
 * This current C implemention is by Alessio Farcomeni and Andy Leung.
 *
 * INPUT:
 * xin 	   		- input data matrix, nin by din
 * kin	   		- number of clusters, kin >= 1
 * intermaxin	- number of iterations of the algorithm
 * Vin		 	- binary matrix for snipping, nin by din. starting solution. number of zeros
 * s			- binary vector of size nin for trimming. starting solution. Number
 * 				  of zeros will be preserved and correspond to trimmed rows. NULL performs no trimming
 * clust		- vector of size nin containing values from 1 to k. starting solution for class labels.
 * Din			- tuning parameter for the fitting algorithm. Corresponds approximately to the maximal 
 *				  change in loss by switching two non outlying entries. Comparing different choices 
 *				  is recommended.
 * nin			- row dimension of xin
 * din			- col dimension of xin
 * muin			- initial estimate of location, vector length din 
 *
 * OUTPUT:
 * Vin		 	- binary matrix for snipping, nin by din. optimal solution. 
 * s			- binary vector of size nin for trimming. optimal solution. 
 * clust		- vector of size nin containing values from 1 to k. starting solution for class labels.
 * muin			- final estimate of location, vector length din 
 * 
 */ 
void sniptotal(double *xin, int *kin, int *itermaxin, int *Vin, int *s, int *clust, double *Din, int *nin, int *din, double *muin)
{

  int n = nin[0];
  int k = kin[0];
  int itermax=itermaxin[0];
  int d = din[0];
  double D = Din[0];
  int i, j, o, iter, w, w1, wc,wc1; 
  double ll, llc, pr, llopt; 

  int szero, sone, vzero, vone; 

  double **mu     = mtxalloc(k, d);
  double **x      = mtxalloc(n, d);
  double **muc    = mtxalloc(k, d);
  double **optmu  = mtxalloc(k, d);
 
  double *sv = (double *) R_alloc(n, sizeof(double));
  
  int **V     = Imtxalloc(n, d);
  int **optV  = Imtxalloc(n, d);
  int **Vc    = Imtxalloc(n, d);
  
  int *sc       = (int *) R_alloc(n, sizeof(int));
  int *opts     = (int *) R_alloc(n, sizeof(int));
  int *clustc   = (int *) R_alloc(n, sizeof(int));
  int *optclust = (int *) R_alloc(n, sizeof(int));

  copyIvec(sc,s,n);
  copyIvec(opts,s,n);

  o=0; 
  for(j=0; j<d; j++) {
    for(i=0; i<n; i++) {
        x[i][j]=xin[o];
        o++;
    }
  }

  o=0; 
  for(i=0; i<n; i++) {
    for(j=0; j<d; j++) {
        V[i][j]=Vin[o];
        o++;
    }
  }

  copyImtx(optV,V,n,d);
  copyImtx(Vc,V,n,d);
  
  o=0; 
  for(i=0; i<k; i++) {
    for(j=0; j<d; j++) {
        mu[i][j]=muin[o];
        o++;
	}
  }

  copymtx(muc,mu,k,d);
  copymtx(optmu,mu,k,d);

  ll = loss(x, n, d, k, mu, clust, s, V);
  llopt = ll*2;

  szero=0; 
  sone=0; 
  vzero=0;
  vone=0;
  for(i=0; i<n; i++) {
    if(s[i]==0) {
        szero++;
	}
    if(s[i]==1) {
        sone++;
	}

    for(j=0; j<d; j++) {
        if(V[i][j]==0) {
            vzero++;
		}
        if(V[i][j]==1) {
            vone++;
		}
	}
  }

  for(iter = 0; iter<itermax; iter++) {
    w=(vzero + 1)*unif_rand() +1;
    w1=(vone + 1)*unif_rand() +1;
    wc=0; 
    wc1=0; 

    for(i=0; i<n; i++) {
        for(j=0; j<d; j++) {
            if(V[i][j]==0) {
                wc1++; 
                if(wc1==w) { 
	                Vc[i][j]=1;
				} else{
                    Vc[i][j]=0;
				}
			}

			if(V[i][j]==1) {
				wc++; 
				if(wc==w1) { 
					Vc[i][j]=0;
				} else{
					Vc[i][j]=1;
				}
			}
		}
	}

	zeroIvec(sc,n);
	for(i=0; i<n;i++) {
		for(j=0; j<d; j++) {
			if(V[i][j]==1) {
				sc[i]=1;
			}
		}
	}

	estmea(x, sc, Vc, clust, muc, d, n, k);
	estmemb(x,clustc, k, n, d, muc, Vc);
	llc = loss(x, n, d, k, muc, clustc, sc, Vc);
	pr = -log(iter+2)*(llc-ll);
	pr = exp(pr*d/D);
	if(pr>1){
		pr=1;
	}

	if((double)  unif_rand() < pr) {
		ll=llc;
		for(i=0; i<n; i++) {
			clust[i]=clustc[i];
			s[i]=sc[i];
			for(j=0; j<d; j++) {
			V[i][j]=Vc[i][j];

			if(V[i][j]==0) {sv++;}
			}
		}
		copymtx(mu,muc,k,d);
	}

	if(llopt>ll) {
		llopt=ll;
		for(i=0; i<n; i++) {
			opts[i]=s[i];
			optclust[i]=clust[i];

			for(j=0; j<d; j++) {
				optV[i][j]=V[i][j];
			}
		}
		copymtx(optmu,mu,k,d);
	}
  }

  copyIvec(clust,optclust,n);
  copyIvec(s,opts,n);

  o=0;
  for(i=0; i<n; i++){
	for(j=0; j<d; j++){
		Vin[o]=optV[i][j];
		o++;
	}
  }

  mtxfree(mu, k, d);
  mtxfree(x, n, d);
  mtxfree(muc, k, d);
  mtxfree(optmu, k, d);
  Imtxfree(V, n, d);
  Imtxfree(optV, n, d);
  Imtxfree(Vc, n, d);
 
}



// Allocate storage for an nxp matrix (double and int version)
double **mtxalloc(int n, int p)
{
    int i;
    double **a  = (double **) Calloc(n, double *);		
    for(i=0; i<n; i++) 
		a[i] = (double *) Calloc(p, double);
    return a;
}
int **Imtxalloc(int n, int p)
{
    int i;
    int **a  = (int **) Calloc(n, int *);		
    for(i=0; i<n; i++) 
		a[i] = (int *) Calloc(p, int);
    return a;
}

// Free the allocated  storage for an nxp matrix (double and int version)
void mtxfree(double **a, int n, int p)
{
    int i;
    for(i=0; i<n; i++) 
	Free(a[i]);
    Free(a);
}
void Imtxfree(int **a, int n, int p)
{
    int i;
    for(i=0; i<n; i++) 
	Free(a[i]);
    Free(a);
}

/* Copy a matrix, from a to b */
void copymtx(double **b, double **a, int n, int p)
{
    int i, j;
    for(i=0; i<n; i++)
	for(j=0; j<p; j++)
	    b[i][j] = a[i][j];
}
void copyImtx(int **b, int **a, int n, int p){
    int i, j;
    for(i=0; i<n; i++)
	for(j=0; j<p; j++)
	    b[i][j] = a[i][j];
}

/* Copy an integer vector, from a to b */
void copyvec(double *b, double *a, int n)
{
    int i;
    for(i=0; i<n; i++)
        b[i] = a[i];
}
void copyIvec(int *b, int *a, int n)
{
    int i;
    for(i=0; i<n; i++)
        b[i] = a[i];
}

/* Zero an integer vector */
void zeroIvec(int *a, int n) {
    int i;
    for(i=0; i<n; i++)
      a[i] = 0;
}



double square(double x)
{
	double res; 
	res=x*x; 
	return res;
}

double loss(double **x, int n, int d, int k, double **mu, int *clust, int *s, int **V)
{
  double res; 
  int i; int j; 
  res = 0;
  for(i=0; i<n; i++){
    if(s[i]==1){
        for(j=0; j<d; j++){
             res+=(x[i][j]-mu[clust[i]-1][j])*(x[i][j]-mu[clust[i]-1][j])*V[i][j];
        }
    }
  }
  return res;
}


void estmea(double **x, int *s, int **V, int *clust, double **mu, int d, int n, int k)
{
  int o;
  int j; 
  int i; 
  double sc; 
 
  for(o=0; o<k; o++){
    for(j=0; j<d; j++){
        sc=0;
        for(i=0; i<n; i++){
            if(clust[i]==o+1 && s[i]==1 && V[i][j]==1){
                sc++;
            }
        }
        mu[o][j]=0; 
        for(i=0; i<n; i++){
            if(clust[i]==o+1 && s[i]==1 && sc>0){ 
                mu[o][j]+=x[i][j]*V[i][j]/sc;
            }
        }
    }
  }
}


void estmemb(double **x, int *clust, int k, int n, int d, double **mu, int **V)
{
  int i; 
  int j; 
  int o; 
  double dm; 
  double jnk; 

  for(i=0; i<n; i++){
    dm=0;
    for(j=0; j<d; j++){
      dm+=square(x[i][j]-mu[0][j])*V[i][j];
	}
    jnk=0; 
    clust[i]=1; 
    for(o=1; o<k; o++){
        jnk=0; 
        for(j=0; j<d; j++) {
            jnk+=square(x[i][j]-mu[o][j])*V[i][j];
	    }
        if(jnk<dm) {
            dm=jnk;
            clust[i]=o+1;
        }
    }
  }
}
 


