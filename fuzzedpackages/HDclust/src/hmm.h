/*==========================================================================================*/
/*                                                                                          */
/* Copyright (C) [Dec 2015]-[April 2017] Jia Li, Department of Statistics,                  */
/* The Pennsylvania State University, USA - All Rights Reserved                             */
/*                                                                                          */
/* Unauthorized copying of this file, via any medium is strictly prohibited                 */
/*                                                                                          */
/* Proprietary and CONFIDENTIAL                                                             */
/*                                                                                          */
/* NOTICE: All information contained herein is, and remains the property of The             */
/* Pennsylvania State University. The intellectual and technical concepts                   */
/* contained herein are proprietary to The Pennsylvania State University and may            */
/* be covered by U.S. and Foreign Patents, patents in process, and are protected            */
/* by trade secret or copyright law. Dissemination of this information or                   */
/* reproduction of this material is strictly forbidden unless prior written                 */
/* permission is obtained from Jia Li at The Pennsylvania State University. If              */
/* you obtained this code from other sources, please write to Jia Li.                       */
/*                                                                                          */
/*                                                                                          */
/* The software is a part of the package for                                                */
/* Clustering with Hidden Markov Models on Variable Blocks                                  */
/*                                                                                          */
/* Written by Jia Li <jiali@stat.psu.edu>, April 7, 2017                                    */ 
/*                                                                                          */
/*==========================================================================================*/
#ifndef HMM_H
#define HMM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "cluster.h"
#define EPSILON 1.0e-5
#define LAMBDA 0.999 
// EPSILON controls the convergence of the loops in Baum-Welch algorithm
// old value of EPSILON is 1.0e-5, then it was set to 1.0e-4
// LAMBDA contraols shrinkage of the covariance matrix to the common covariance 
// in Baum-Welch algorithm. Old value of LAMBDA is 0.999
// then it was changed to 1.0
extern int DIAGCOV;

typedef struct 
{
  int id;
  int value;
} SORT_INT;
 
typedef struct 
{
  int id;
  double value;
} SORT_DOUBLE;
 
typedef struct gaussmodel_struct
{
  int dim;
  int exist;
  int cls;
  double *mean;
  double **sigma;
  double **sigma_inv;
  double sigma_det_log;
} GaussModel;

typedef struct hmmmodel_struct
{
  int dim;
  int numst;   /* numst is the number of states */
  int prenumst; //number of state of previous block
  int *var;    // The original variable in each dimension
  GaussModel **stpdf;
  double **a; // a[prenumst][numst] may not be a square matrix
  double *a00;  /* pmf of states at the boundary when there's no neighbor */
} HmmModel;

typedef struct gmmmodel_struct
{
  int dim;
  int numst;   /* numst is the number of states */
  GaussModel **stpdf;
  double *p;  // prior of states
} GmmModel;

typedef struct condchain_struct
{
  int dim;
  int nb; //number of blocks the variables are divided into
  int *bdim; // dimension of each block [nb]
  int *cbdim; //cumulated starting position of the dimension of each block
  int **var; // original variables in each block [nb]
  int *numst; //number of states per block [nb], duplicate info from mds
  int *cnumst; //cumulated starting position of numst of each block [nb]
  int maxnumst; //maximum of numst[nb]
  HmmModel **mds; //list of HMM models for each block [nb]
} CondChain;

typedef struct component_struct
{
  int *st; //states st[nb]
  double *mu; //mu[dim]
  double *mode;
  double logpdf; //log of the density at the mode
} CompMode;


/*---------- functions ---------------*/
/*--------- estimate.c ---------------*/
/*------------------------------------*/

extern void forward(double *u, double *thetalog, CondChain *md, double *loglikehd);
extern void backward(double *u, double *betalog, CondChain *md);
extern void CompLm(double *thetalog, double *betalog, double **Lm, CondChain *md);
extern void CompHml(double *u, double *thetalog, double *betalog, double ***Hml, 
		    CondChain *md);
extern void updatepar_adder(double *u, double *thetalog, double *betalog, 
				 CondChain *md, double **musum, 
				 double ***mom2sum, double ***Hml, double **Lm);
extern void viterbi(CondChain *md, double *u, int *optst, double *inita, 
		    double *lastmerit);
extern void viterbi_mulseq(CondChain *md, double **u, int nseq, int **st);
extern void initialize(double *u, int nseq, int dim, HmmModel *md, int ranflag);
extern void initialize2(double *u, int nseq, int dim, HmmModel *md, double *cdbk);
extern void initial_ccm(double **u,int nseq, CondChain *md);
extern void initial_ccm1(double **u,int nseq, CondChain *md, int sd);
extern double comploglike(CondChain *md, double **u, int nseq, double *wt, double *logl);
extern double classlikehd(CondChain *md, double **u, int nseq, double ***cprob, double *wt);
extern int baumwelch(double **u, int nseq, CondChain *md, double *loglikehd, 
		     double *lhsumpt, double epsilon, double *wt);
extern void ordervar(double **u,int nseq, int nb, int *bdim, int **var);
extern void hmmfit(double **u, int nseq, int nb, int *bdim, int **var, int *numst, 
		   CondChain *md, double *loglikehd, double *lhsumpt, 
		   double epsilon, double *wt);
extern void hmmfit_minit(double **u, int nseq, int nb, int *bdim, int **var, int *numst, CondChain **md, 
			 double *loglikehd, double *lhsumpt, double epsilon, double *wt, int ninit0,
			 int ninit1, int ninit2, int randomseed);
extern void hmmfit_minit2(double **u, int nseq, int nb, int *bdim, int **var, int *numst, CondChain **md, 
			  double *loglikehd, double *lhsumpt, double epsilon, double *wt, int ninit0,
			  int ninit1, int ninit2, int randomseed);
extern void printvbinfo(int nb, int *bdim, int **var, double lhsum);
extern int computenp(int nb, int *bdim, int *numst);
extern void hmmfit_vb(double **u, int nseq, int dim, int *snb, int **sbdim, int ***svar, int nperm, int nperm0, int **vlist0, CondChain **md, double *loglikehd, double *lhsumpt, double epsilon, double *wt, int ninit0, int ninit1, int ninit2, int randomseed, int *Numst0, int maxdim, int mindim, int relaxsearch);

/*-------------------------------------*/
/*-------------- prob.c ---------------*/
/*-------------------------------------*/

extern int newgauss(GaussModel *md, int dim, int exist);
extern int cpgauss(GaussModel *md1, GaussModel *md2);
extern void newhmm(HmmModel *md, int dim, int numst, int prenumst);
extern void freehmm(HmmModel **md_pt);
extern void cphmm(HmmModel *md1, HmmModel *md2);
extern void newccm(CondChain *md, int nb, int *bdim, int **var, int *numst);
extern void freeccm(CondChain **md_pt);
extern double gauss_pdf_log(double *ft, GaussModel *gm);
extern double gauss_pdf(double *ft, GaussModel *gm);
extern double mix_gauss_pdf_log(double *ft, GaussModel **gmlist, double *prior, int ncmp);

/*-------------------------------------*/
/*------------- modelio.c -------------*/
/*-------------------------------------*/

extern unsigned char write_hmm(HmmModel *md, FILE *outfile);
extern unsigned char read_hmm(HmmModel *md, FILE *infile);
extern unsigned char print_hmm(HmmModel *md, FILE *outfile);
extern unsigned char write_ccm(CondChain *md, FILE *outfile);
extern unsigned char read_ccm(CondChain *md, FILE *infile);
extern unsigned char print_ccm(CondChain *md, FILE *outfile);
extern unsigned char readascii_ccm(CondChain *md, FILE *infile);
extern unsigned char readascii_hmm(HmmModel *md, FILE *infile);
extern unsigned char readascii2_ccm(CondChain *md, FILE *infile);
extern unsigned char readascii2_hmm(HmmModel *md, FILE *infile);
extern void printclsinfo(FILE *outfile, int *clsid, int nseq, int numcls);
extern void readrefcls(FILE *reffile, int *rncls, int *rnseq, double ***refmode,
		       double **refsigma, int ***refpath, int **refcls, int dim, int nb, int *nseq, int **ct);
extern void printrefcls(FILE *reffile, int rncls, int rnseq, double **refmode,
			double *refsigma, int **refpath, int *refcls, int dim, int nb, int nseq, int *ct);


/*-------------------------------------*/
/*------------- modal.c ---------------*/
/*-------------------------------------*/

extern void SortInt(int *org, int *buf, int *invid, int sz);
extern void SortDouble(double *org, double *buf, int *invid, int sz);
extern void SortLexigraphicInt(int **orgmat, int **bufmat, int *invid, int dim, int sz);
extern int FindEntry(int **mat, int *entry, int dim, int sz);
extern int CountDif(int *buf, int len); 
extern double l2sq(double *ft1, double *ft2, int dim); 
extern double distmaxdim(double *ft1, double *ft2, int dim, double *sigma); 
extern void standarddev(double **u, int nseq, int dim, double *sigma);
extern void wtsum_matrix(double *wt, double ***mat, int len, int nr, int nc, 
			 double **smat);
extern void wtsum_vec(double *wt, double **mat, int len, int nr,  double *smat);
extern void matvec_multiply(double **mat, double *vec, int nr, int nc, 
			    double *res);
extern void sigmainv_array(CondChain *md, double *****sigma_inv_pt, 
			   double ****sigmainvmu_pt);
extern void sigmainv_array_gmm(GmmModel *md, double ****sigma_inv_pt, 
			       double ***sigmainvmu_pt);
extern void OverallSigma(CondChain *md, double *sigma);
extern void OverallSigma_Gmm(GmmModel *md, double *sigma);
extern void DataSigma(double **u, double *sigmadat, int dim, int nseq);
extern double bwmem(CondChain *md, double *x, double *res);
extern double posterior(GmmModel *md, double *x, double *p);
extern double mem(GmmModel *md, double *x, double *res);
extern void SetCompMode(CompMode *cpm, int *optst, CondChain *md);
extern void freeCompMode(CompMode *cpm);
extern void SetCompLogprior(double *logprior, int *mypath, CondChain *md);
extern int FuseGauss(GaussModel *gmd, int **mypath, int len, CondChain *md);
extern void FindDifSeq(int **optst, int nseq, int nb, int ***newst_pt,
		       int *newnseq, int *newid);
extern void groupmode(double **mode, int dim, int num, int *cls, int *numcls, 
		      double *sigma, double threshold, int meandist);
extern int FindCluster(double *mode, int dim, int rncls, double **refmode, 
		       double *sigma, double threshold, int meandist);
#endif
