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

#include "hmm.h"
#include <Rcpp.h>

//static int CompFcn(SORT_INT *a, SORT_INT *b)
static int CompFcn(const void *aa, const void *bb)
{
  SORT_INT* a = (SORT_INT *) aa;
  SORT_INT* b = (SORT_INT *) bb;
  
  if (a->value > b->value)
    return (1);
  if (a->value < b->value)
    return (-1);
  return (0);
}

//static int CompFcnDb(SORT_DOUBLE *a, SORT_DOUBLE *b)
static int CompFcnDb(const void *aa, const void *bb)
{
  SORT_DOUBLE *a = (SORT_DOUBLE *) aa;
  SORT_DOUBLE *b = (SORT_DOUBLE *) bb;
	
  if (a->value > b->value)
    return (1);
  if (a->value < b->value)
    return (-1);
  return (0);
}

void SortInt(int *org, int *buf, int *invid, int sz)
{
  int j;
  SORT_INT *score;

  score=(SORT_INT *)calloc(sz,sizeof(SORT_INT));
  
  try {
    if (score==NULL) {         	
      throw std::range_error("Unable to allocate space in SortInt");
    }
  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  //if (score==NULL) {
    
    //Rcpp::stop("Unable to allocate space in SortInt");
  //}

  for (j=0;j<sz;j++) {
    score[j].id=j;
    score[j].value=org[j];
  }
  qsort((SORT_INT *)score, sz, sizeof(SORT_INT), CompFcn);

  for (j=0;j<sz;j++) {
    buf[j]=org[score[j].id];
    invid[j]=score[j].id; //store the original id for the new order
  }

  free(score);
}

//Assume the integer arrays have been sorted lexicographically by
//SortLexigraphicInt()
int FindEntry(int **mat, int *entry, int dim, int sz)
{
  int i,k;
  int stpos, edpos, stpos2, edpos2;

  stpos=0;
  edpos=sz;
  for (k=0;k<dim;k++) {
    stpos2=-1;
    edpos2=0;
    for (i=stpos;i<edpos;i++) {
      if (mat[i][k]==entry[k]) {
	if (stpos2<0) stpos2=i;
	edpos2=i+1;
      }
      else {
	if (mat[i][k]>entry[k]) {
	  break;
	}
      }
    }
    stpos=stpos2;
    edpos=edpos2;
    if (stpos<0) break;
  }

  return(stpos);
}

int Difseq(int *s1, int *s2, int dim)
{
  int i,k;
  for (i=0,k=0;i<dim;i++)
    if (s1[i]==s2[i]) k++;

  if (k==dim) return(0);
  else return(1);
}

void SortLexigraphicInt(int **orgmat, int **bufmat, int *invid, int dim, int sz)
{
  int *org, *buf, **buf2, *invid0, *invid2;
  int i,k, stpos, edpos,mm;

  org=(int *)calloc(sz,sizeof(int));
  buf=(int *)calloc(sz,sizeof(int));
  invid0=(int *)calloc(sz,sizeof(int));
  invid2=(int *)calloc(sz,sizeof(int));
  buf2=(int **)calloc(sz,sizeof(int *));

  for (i=0;i<sz;i++) {
    bufmat[i]=orgmat[i];
    invid[i]=i;
  }

  for (k=0;k<dim;k++){
    for (i=0;i<sz;i++) org[i]=bufmat[i][k];
    stpos=0;
    if (k==0) edpos=sz;
    else {
      edpos=sz;
      for (i=1;i<sz;i++) {
	if (Difseq(bufmat[i],bufmat[stpos],k)) {
	  edpos=i;
	  break;
	}
      }
    }

    mm=0;
    while (stpos<edpos) {
      if (edpos-stpos>1) {
	mm++;
	
	SortInt(org+stpos,buf+stpos,invid0+stpos,edpos-stpos);
	
	for (i=stpos;i<edpos;i++) {
	  buf2[i]=bufmat[i];
	  invid2[i]=invid[i];
	}
	for (i=stpos;i<edpos;i++) {
	  bufmat[i]=buf2[invid0[i]+stpos];
	  invid[i]=invid2[invid0[i]+stpos];
	}
      }
      
      stpos=edpos;
      if (k>0) {
	if (stpos<sz) {
	  edpos=sz;
	  for (i=stpos+1;i<sz;i++) {
	    if (Difseq(bufmat[i],bufmat[stpos],k)){
	      edpos=i;
	      break;
	    }
	  }
	}
      }
    }

    if (mm==0) {
      //fprintf(stderr, "Already sorted at k=%d\n",k);
      // Already sorted at this k, no need to go further
      break;
    }
  }

  free(org);
  free(buf);
  free(buf2);
  free(invid0);
  free(invid2);
}

void SortDouble(double *org, double *buf, int *invid, int sz)
{
  int j;
  SORT_DOUBLE *score;

  score=(SORT_DOUBLE *)calloc(sz,sizeof(SORT_DOUBLE));
  
  try {
    if (score==NULL) {
      free(score);
      throw std::range_error("Unable to allocate space in SortDouble");
    }
  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  //if (score==NULL) {
    //free(score);
    //Rcpp::stop("Unable to allocate space in SortDouble.\n");
  //}

  for (j=0;j<sz;j++) {
    score[j].id=j;
    score[j].value=org[j];
  }
  qsort((SORT_DOUBLE *)score, sz, sizeof(SORT_DOUBLE), CompFcnDb);

  for (j=0;j<sz;j++) {
    buf[j]=org[score[j].id];
    invid[j]=score[j].id; //store the original id for the new order
  }

  free(score);
}

int CountDif(int *buf, int len) //assume sorted array 
{
  int i,k;

  for (i=1,k=1;i<len;i++) {
    if (buf[i]>buf[i-1]) {k++;}
  }
  return(k);
}

double l2sq(double *ft1, double *ft2, int dim) 
{
  int i;
  double db1;

  db1=0.0;
  for (i=0;i<dim;i++) {
    db1+=(ft1[i]-ft2[i])*(ft1[i]-ft2[i]);
  }

  return(db1);
}

double distmean(double *ft1, double *ft2, int dim, double *sigma) 
{
  int i;
  double db1,db2;

  db1=0.0;
  for (i=0;i<dim;i++) {
    db2=fabs(ft1[i]-ft2[i])/sigma[i];//normalize wrt sigma
    //db2=fabs(ft1[i]-ft2[i]);//No normalization
    //fprintf(stdout, "%d ft1=%f, ft2=%f, sigma=%f, db2=%f\n", i,ft1[i],ft2[i],sigma[i], db2);
    db1+=db2;
  }
  db1/=(float)dim;

  //fprintf(stdout, "db1=%f\n",db1);
  //exit(0);

  
  return(db1);
}

double distmaxdim(double *ft1, double *ft2, int dim, double *sigma) 
{
  int i;
  double db1,db2;

  db1=0.0;
  for (i=0;i<dim;i++) {
    db2=fabs(ft1[i]-ft2[i])/sigma[i];//normalize wrt sigma
    if (db1 < db2) db1=db2;
  }

  return(db1);
}

void standarddev(double **u, int nseq, int dim, double *sigma)
{
  int i,j;
  double *mu;
  
  mu=(double *)calloc(dim,sizeof(double));
  for (i=0;i<dim;i++) { mu[i]=0.0; sigma[i]=0.0;}

  for (i=0;i<nseq;i++) {
    for (j=0;j<dim;j++){
      mu[j]+=u[i][j];
      sigma[j]+= (u[i][j]*u[i][j]);
    }
  }

  for (i=0;i<dim;i++) {
    mu[i]/=(double)nseq;
    sigma[i]/=(double)nseq;
    sigma[i]-= (mu[i]*mu[i]);
    sigma[i]=sqrt(sigma[i]);
  }

  free(mu);
}

void wtsum_matrix(double *wt, double ***mat, int len, int nr, int nc, 
		  double **smat)
{
  int i,j,m;

  for (i=0;i<nr;i++)
    for (j=0;j<nc;j++) 
      smat[i][j]=0.0;

  for (m=0;m<len;m++) {
    for (i=0;i<nr;i++)
      for (j=0;j<nc;j++) 
	smat[i][j]+=wt[m]*mat[m][i][j];
  }
}

void wtsum_matrix_diag(double *wt, double ***mat, int len, int dim,
		       double **smat, int diagonal)
{
  int i,j,m;

  if (diagonal!=1) {
    wtsum_matrix(wt,mat,len,dim,dim,smat);
    return;
  }

  for (i=0;i<dim;i++)
    for (j=0;j<dim;j++) 
      smat[i][j]=0.0;

  for (m=0;m<len;m++) {
    for (i=0;i<dim;i++)
      smat[i][i]+=wt[m]*mat[m][i][i];
  }
}

void wtsum_vec(double *wt, double **mat, int len, int nr,  double *smat)
{
  int i,m;

  for (i=0;i<nr;i++)
      smat[i]=0.0;

  for (m=0;m<len;m++) {
    for (i=0;i<nr;i++)
      smat[i]+=wt[m]*mat[m][i];
  }
}

void matvec_multiply(double **mat, double *vec, int nr, int nc, double *res)
{
  int i,j;

  //mat[nr][nc], vec[nc], res[nr]
  for (i=0;i<nr;i++) {
    res[i]=0.0;
    for (j=0;j<nc;j++)
      res[i]+=mat[i][j]*vec[j];
  }
}

void squarematvec_multiply(double **mat, double *vec, int dim, double *res, int diagonal)
{
  int i;

  if (diagonal!=1)
    matvec_multiply(mat,vec,dim,dim,res);
  else {
    for (i=0;i<dim;i++) 
      res[i]=vec[i]*mat[i][i];
  }
}


void sigmainv_array(CondChain *md, double *****sigma_inv_pt, 
		    double ****sigmainvmu_pt)
//sigma_inv[nb][numst[]][bdim[]][bdim[]]
//sigmainvmu[nb][numst[]][bdim[]]
{
  int i,j;
  int nb, *numst, *bdim;
  double ****sigma_inv, ***sigmainvmu;

  nb=md->nb; numst=md->numst; bdim=md->bdim;

  sigma_inv=(double ****)calloc(nb,sizeof(double ***));
  sigmainvmu=(double ***)calloc(nb,sizeof(double **));

  for (i=0;i<nb;i++) {
    sigma_inv[i]=(double ***)calloc(numst[i],sizeof(double **));
    sigmainvmu[i]=(double **)calloc(numst[i],sizeof(double *));
  }

  for (i=0;i<nb;i++)
    for (j=0;j<numst[i];j++) {
      matrix_2d_double(sigma_inv[i]+j, bdim[i],bdim[i]);
      vector_double(sigmainvmu[i]+j,bdim[i]);

      matrix_2d_cpy_double(sigma_inv[i][j],md->mds[i]->stpdf[j]->sigma_inv,bdim[i],bdim[i]);
      squarematvec_multiply(sigma_inv[i][j], md->mds[i]->stpdf[j]->mean, bdim[i], sigmainvmu[i][j], DIAGCOV);

    }
  
  *sigma_inv_pt=sigma_inv;
  *sigmainvmu_pt=sigmainvmu;
}

//For GMM
void sigmainv_array_gmm(GmmModel *md, double ****sigma_inv_pt, 
			double ***sigmainvmu_pt)
//sigma_inv[numst][dim][dim]
//sigmainvmu[numst][dim][dim]
{
  int j, numst, dim;
  double ***sigma_inv, **sigmainvmu;

  numst=md->numst; dim=md->dim;

  sigma_inv=(double ***)calloc(numst,sizeof(double **));
  sigmainvmu=(double **)calloc(numst,sizeof(double *));

  for (j=0;j<numst;j++) {
    matrix_2d_double(sigma_inv+j, dim,dim);
    vector_double(sigmainvmu+j,dim);
    
    matrix_2d_cpy_double(sigma_inv[j],md->stpdf[j]->sigma_inv,dim,dim);
    squarematvec_multiply(sigma_inv[j], md->stpdf[j]->mean, dim,sigmainvmu[j], DIAGCOV);
    
  }
  
  *sigma_inv_pt=sigma_inv;
  *sigmainvmu_pt=sigmainvmu;
}

void OverallSigma(CondChain *md, double *sigma)
{
  int i,j,k,m;

  for (i=0,m=0;i<md->nb;i++)
    for (j=0;j<md->bdim[i];j++) {
      sigma[m]=0.0;
      for (k=0;k<md->numst[i];k++) {
	sigma[m]+=md->mds[i]->a00[k]*md->mds[i]->stpdf[k]->sigma[j][j];
      }
      sigma[m]=sqrt(sigma[m]);
      m++;
    }
}

void DataSigma(double **u, double *sigmadat, int dim, int nseq)
{
  double *mv;
  int i,j;

  if (nseq==0) return;
  
  mv=(double *)calloc(dim,sizeof(double));
  for (i=0;i<dim;i++) {
    mv[i]=0.0;
    sigmadat[i]=0.0;
  }
  for (i=0;i<nseq;i++){
    for (j=0;j<dim;j++)
      mv[j]+=u[i][j];
  }
  for (i=0;i<dim;i++) mv[i]/=(double)nseq;

  for (i=0;i<nseq;i++){
    for (j=0;j<dim;j++)
      sigmadat[j]+= ((u[i][j]-mv[j])*(u[i][j]-mv[j]));
  }

  for (i=0;i<dim;i++) sigmadat[i]=sqrt(sigmadat[i]/(double)nseq);
  free(mv);
}

void OverallSigma_Gmm(GmmModel *md, double *sigma)
{
  int j,k;

  for (j=0;j<md->dim;j++) {
    sigma[j]=0.0;
    for (k=0;k<md->numst;k++) {
      sigma[j]+=md->p[k]*md->stpdf[k]->sigma[j][j];
    }
    sigma[j]=sqrt(sigma[j]);
  }
}

//Baum-Welch Modal EM algorithm
double bwmem(CondChain *md, double *x, double *res)
{
  int i,j,m, ite, minite=30, maxite=1000;
  int nb, *numst, *bdim,maxdim,dim;
  double *thetalog, *betalog, **Lm, oldloglike, loglike;
  double dif, ratio, threshold1=1.0e-6, threshold2=1.0e-4, *sigma,db1;
  double **wtsigmainv, **wtsigma, *wtmu;
  double *mode, *oldmode;
  double ****sigma_inv_pool, ***sigmainvmu_pool;

  nb=md->nb; numst=md->numst; bdim=md->bdim;
  maxdim=bdim[0];
  for (i=1;i<nb;i++) {if (bdim[i]>maxdim) maxdim=bdim[i];}

  matrix_2d_double(&wtsigmainv,maxdim,maxdim);
  matrix_2d_double(&wtsigma,maxdim,maxdim);
  wtmu=(double *)calloc(maxdim,sizeof(double));
  for (i=0,m=0;i<nb;i++) m+=numst[i];
  thetalog=(double *)calloc(m, sizeof(double));
  betalog=(double *)calloc(m, sizeof(double));

  Lm=(double **)calloc(nb,sizeof(double *));
  for (i=0;i<nb;i++) Lm[i]=(double *)calloc(numst[i],sizeof(double));

  for (i=0,dim=0;i<nb;i++)  dim+=bdim[i];
  mode=(double *)calloc(dim,sizeof(double));
  oldmode=(double *)calloc(dim,sizeof(double));

  sigmainv_array(md, &sigma_inv_pool, &sigmainvmu_pool);

  for (i=0;i<dim;i++) oldmode[i]=mode[i]=x[i];

  ite=0;
  oldloglike=1.0;
  //compute the proper sigma=sqrt(sigmasq);
  sigma=(double *)calloc(dim,sizeof(double));
  OverallSigma(md,sigma); //sigma of each dimension

  while (ite<maxite) {
    forward(mode, thetalog, md, &loglike);
    backward(mode,betalog,md);
    CompLm(thetalog, betalog, Lm, md);

    ratio=fabs((loglike-oldloglike)/oldloglike);
    dif=distmaxdim(mode, oldmode, dim,sigma);
    if (ratio<threshold1 && dif < threshold2 && ite > minite)
      break;

    oldloglike=loglike;
    for (j=0;j<dim;j++) oldmode[j]=mode[j];

    //Update mode
    for (i=0,m=0;i<nb;i++) {
      wtsum_matrix_diag(Lm[i], sigma_inv_pool[i], numst[i], bdim[i], wtsigmainv, DIAGCOV);
      mat_det_inv_diag_double(wtsigmainv, wtsigma,&db1, bdim[i], DIAGCOV);
      wtsum_vec(Lm[i], sigmainvmu_pool[i], numst[i], bdim[i], wtmu);
      squarematvec_multiply(wtsigma, wtmu, bdim[i], mode+m, DIAGCOV);
      m+=bdim[i];
    }
    
    ite++;
  }

  for (i=0;i<dim;i++) res[i]=mode[i];
  forward(res, thetalog, md, &loglike);

  free(thetalog); free(betalog);
  free_matrix_2d_double(&wtsigmainv,maxdim);
  free_matrix_2d_double(&wtsigma,maxdim);
  free(wtmu);
  free_matrix_2d_double(&Lm,nb);
  free(mode); free(oldmode); free(sigma);

  for (i=0;i<nb;i++) {
    for (j=0;j<numst[i];j++) {
      free_matrix_2d_double(sigma_inv_pool[i]+j,bdim[i]);
      free(sigmainvmu_pool[i][j]);
    }
    free(sigma_inv_pool[i]);
    free(sigmainvmu_pool[i]);
  }
  free(sigma_inv_pool);
  free(sigmainvmu_pool);

  return(loglike);
}

double posterior(GmmModel *md, double *x, double *p)
{
  int i,numst;
  double v1,v2;

  numst=md->numst;

  for (i=0;i<numst;i++) {
    if (md->p[i]>0.0)
      p[i]=log(md->p[i])+gauss_pdf_log(x,md->stpdf[i]);
    else p[i]=-HUGE_VAL;
  }

  for (i=1,v1=p[0];i<numst;i++) {
    if (p[i]>v1) v1=p[i];
  }

  for (i=0,v2=0.0;i<numst;i++) {
    p[i]=exp(p[i]-v1);
    v2+=p[i];
  }
  
  for (i=0;i<numst;i++) p[i]/=v2;
  return(log(v2)+v1);
}

// For GMM
double mem(GmmModel *md, double *x, double *res)
{
  int i,j, ite, minite=30, maxite=1000;
  int numst,dim;
  double *Lm, oldloglike, loglike, db1;
  double dif, ratio, threshold1=1.0e-6, threshold2=1.0e-4, *sigma;
  double **wtsigmainv, **wtsigma, *wtmu;
  double *mode, *oldmode;
  double ***sigma_inv_pool, **sigmainvmu_pool;

  numst=md->numst; dim=md->dim;

  matrix_2d_double(&wtsigmainv,dim,dim);
  matrix_2d_double(&wtsigma,dim,dim);
  wtmu=(double *)calloc(dim,sizeof(double));
  Lm=(double *)calloc(numst,sizeof(double));//posterior prob

  mode=(double *)calloc(dim,sizeof(double));
  oldmode=(double *)calloc(dim,sizeof(double));

  sigmainv_array_gmm(md, &sigma_inv_pool, &sigmainvmu_pool);

  for (i=0;i<dim;i++) oldmode[i]=mode[i]=x[i];

  ite=0;
  oldloglike=1.0;
  //compute the proper sigma=sqrt(sigmasq);
  sigma=(double *)calloc(dim,sizeof(double));
  OverallSigma_Gmm(md,sigma);

  while (ite<maxite) {
    loglike=posterior(md,mode,Lm);//Needs definition this function

    ratio=fabs((loglike-oldloglike)/oldloglike);
    dif=distmaxdim(mode, oldmode, dim,sigma);
    if (ratio<threshold1 && dif < threshold2 && ite > minite)
      break;

    oldloglike=loglike;
    for (j=0;j<dim;j++) oldmode[j]=mode[j];

    //Update mode
    wtsum_matrix_diag(Lm, sigma_inv_pool, numst, dim, wtsigmainv, DIAGCOV);
    mat_det_inv_diag_double(wtsigmainv, wtsigma,&db1, dim, DIAGCOV);
    wtsum_vec(Lm, sigmainvmu_pool, numst, dim, wtmu);
    squarematvec_multiply(wtsigma, wtmu, dim, mode, DIAGCOV);
    
    ite++;
  }

  for (i=0;i<dim;i++) res[i]=mode[i];
  loglike=posterior(md,res,Lm);

  free_matrix_2d_double(&wtsigmainv,dim);
  free_matrix_2d_double(&wtsigma,dim);
  free(wtmu);
  free(Lm);
  free(mode); free(oldmode); free(sigma);

  for (j=0;j<numst;j++) {
    free_matrix_2d_double(sigma_inv_pool+j,dim);
    free(sigmainvmu_pool[j]);
  }
  free(sigma_inv_pool); free(sigmainvmu_pool);

  return(loglike);
}

void SetCompMode(CompMode *cpm, int *optst, CondChain *md)
{
  int i,j,k,dim;

  cpm->st=(int *)calloc(md->nb,sizeof(int));
  for (i=0,dim=0;i<md->nb;i++) {
    cpm->st[i]=optst[i];
    dim+=md->bdim[i];
  }
  cpm->mu=(double *)calloc(dim,sizeof(double));
  cpm->mode=(double *)calloc(dim,sizeof(double));

  for (i=0,k=0;i<md->nb;i++) 
    for (j=0;j<md->bdim[i];j++) {
      cpm->mu[k]=md->mds[i]->stpdf[optst[i]]->mean[j];
      k++;
    }
}
void freeCompMode(CompMode *cpm)
{
  free(cpm->st);
  free(cpm->mu);
  free(cpm->mode);
}

void SetCompLogprior(double *logprior, int *mypath, CondChain *md)
{
  int i;
  *logprior=log(md->mds[0]->a00[mypath[0]]);
  for (i=1;i<md->nb;i++) {
    (*logprior) +=log(md->mds[i]->a[mypath[i-1]][mypath[i]]);
  }
}

int FuseGauss(GaussModel *gmd, int **mypath, int len, CondChain *md)
{
  int i,j,k,mm,n,j1,j2;
  double *logprior;
  double *prior, db, db2;
  int dim;
  double **sigma2, *mu, mysigma_det;

  dim=gmd->dim;

  logprior=(double *)calloc(len,sizeof(double));
  prior=(double *)calloc(len,sizeof(double));
  for (i=0;i<len;i++)
    SetCompLogprior(logprior+i, mypath[i],md);

  //normalize logprior and put in prior
  for (i=1, db=logprior[0];i<len;i++) {
    if (logprior[i]>db) db=logprior[i];
  }
  for (i=0, db2=0.0;i<len;i++) {
    logprior[i]=exp(logprior[i]-db);
    db2+=logprior[i];
  }
  for (i=0;i<len;i++) prior[i]=logprior[i]/db2;

  mu=(double *)calloc(dim,sizeof(double));
  sigma2=(double **)calloc(dim,sizeof(double *));
  for (i=0;i<dim;i++) sigma2[i]=(double *)calloc(dim,sizeof(double));

  //Initialization
  for (i=0;i<dim;i++) gmd->mean[i]=0.0;
  for (i=0;i<dim;i++)
    for (j=0;j<dim;j++) {
      gmd->sigma[i][j]=0.0;
      sigma2[i][j]=0.0;
    }

  for (n=0;n<len;n++) {
    for (i=0,k=0;i<md->nb;i++) 
      for (j=0;j<md->bdim[i];j++) {
	gmd->mean[k]+=(prior[n]*md->mds[i]->stpdf[mypath[n][i]]->mean[j]);
	k++;
      }
  }

  //Compute the between-component covariance
  for (n=0;n<len;n++) {
    for (i=0,k=0;i<md->nb;i++) 
      for (j=0;j<md->bdim[i];j++) {
	mu[k]=md->mds[i]->stpdf[mypath[n][i]]->mean[j];
	k++;
      }
    for (i=0;i<dim;i++) { mu[i]-=gmd->mean[i]; }//centerize, remove mean
    for (i=0;i<dim;i++) 
      for (j=i;j<dim;j++) { 
	sigma2[i][j]+=(prior[n]*mu[i]*mu[j]);
	sigma2[j][i]=sigma2[i][j];
      }
  }

  //Store the within-component covariance in gmd->sigma first
  for (n=0;n<len;n++) {
    for (i=0,k=0;i<md->nb;i++) {
      for (j1=0;j1<md->bdim[i];j1++) {
	for (j2=0;j2<md->bdim[i];j2++) {
	  gmd->sigma[k+j1][k+j2]+=(prior[n]*md->mds[i]->stpdf[mypath[n][i]]->sigma[j1][j2]);
	}
      }
      k+=md->bdim[i];
    }
  }

  //Get the overall sigma by adding the between-component covariance

  for (i=0;i<dim;i++)
    for (j=0;j<dim;j++) {
      gmd->sigma[i][j]+=sigma2[i][j];
    }

  mm=mat_det_inv_diag_double(gmd->sigma, gmd->sigma_inv, &(mysigma_det),dim, DIAGCOV);
  if (DIAGCOV!=1) { gmd->sigma_det_log=log(mysigma_det);}
  else {
    mm=1;
    gmd->sigma_det_log=0.0;
    for (i=0;i<dim;i++) {
      gmd->sigma_det_log+=log(gmd->sigma[i][i]);
      if (gmd->sigma[i][i]<=0.0) mm=2;
    }
  }


  free(logprior);
  free(prior);
  free(mu);
  for (i=0;i<dim;i++) free(sigma2[i]);
  free(sigma2);

  if (mm==2 || (DIAGCOV!=1 && mysigma_det<=0)) { //non-positive definite covariance matrix    
    return(0);
  }
  return(1); //non-singular covariance                                                            
}

//id[nseq] stores for each original sequence the new id of in the array
//with all different sequences
//The array st[nseq][nb] is assumed to be sorted already by SortLexigraphicInc()
//id[nseq] has space allocated already
int CountDifArray(int **st, int nseq, int nb, int *id)
{
  int i,j,n;
  int start;

  start=0;
  id[0]=0;
  for (i=start+1;i<nseq;i++) {
    for (j=0,n=0;j<nb;j++) {
      if (st[i][j]==st[start][j]) n++;
    }
    if (n==nb) {//find same sequence
      id[i]=id[start];
    }
    else {
      id[i]=id[start]+1;
      start=i;
    }
  }

  return(id[nseq-1]+1);
}

void FindDifSeq(int **optst, int nseq, int nb, int ***newst_pt,
		int *newnseq, int *newid)
{
  int i,j;
  int **newst;
  int **bufoptst, *invid, *id2;

  bufoptst=(int **)calloc(nseq,sizeof(int *));
  invid=(int *)calloc(nseq,sizeof(int));
  id2=(int *)calloc(nseq,sizeof(int));

  SortLexigraphicInt(optst, bufoptst, invid, nb,nseq);
	
  *newnseq=CountDifArray(bufoptst,nseq,nb,id2);
  
  newst=(int **)calloc(*newnseq,sizeof(int *));
  for (i=0;i<*newnseq;i++) 
    newst[i]=(int *)calloc(nb,sizeof(int));

  for (i=0;i<nseq;i++) {
    newid[invid[i]]=id2[i];
  }

  for (i=0;i<nseq;i++) {
    //newst[][] may be given multiple times with the same result
    for (j=0;j<nb;j++) 
      newst[id2[i]][j]=bufoptst[i][j];
  }

  *newst_pt=newst;

  free(bufoptst);
  free(invid);
  free(id2);
}

void groupmode(double **mode, int dim, int num, int *cls, int *numcls, 
	       double *sigma, double threshold, int meandist)
{
  int i,m;
  int st, index, *done;
  double db1;
  //double db2, db3;
  
  done=(int *)calloc(num,sizeof(int));
  for (i=0;i<num;i++) {
    done[i]=0;
    cls[i]=0;
  }
    
  st=0;
  index=0;
  while (st<num) {
    cls[st]=index;
    for (i=st+1; i<num; i++) {
      if (done[i]) continue;

      if (meandist)
	db1=distmean(mode[st], mode[i], dim, sigma);
      else
	db1=distmaxdim(mode[st], mode[i], dim, sigma);
      //db2=distmean(mode[st], mode[i], dim, sigma);
      //db3=distmaxdim(mode[st], mode[i], dim, sigma);
      //fprintf(stderr, "st=%d, i=%d, mean dist=%f, max dist=%f\n",st,i,db2,db3);
      
      if (db1<threshold) {
	cls[i]=index;
	done[i]=1;
      }
    }

    m=st;
    for (i=m+1; i<num; i++) {
      if (!done[i]) {
	st=i;
	break;
      }
    }

    if (st==m) st=num; // all done

    index++;
  }

  *numcls=index;

  free(done);
}

int FindCluster(double *mode, int dim, int rncls, double **refmode,
		double *sigma, double threshold, int meandist)
{
  int i;
  int index=-1;
  double db1;
  
  for (i=0;i<rncls;i++) {
    if (meandist)
      db1=distmean(mode,refmode[i],dim,sigma);
    else
      db1=distmaxdim(mode,refmode[i],dim,sigma);

    if (db1<threshold) {
      index=i;
      break;
    }
  }

  return(index);
}
