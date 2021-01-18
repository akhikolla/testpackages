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

#include "cluster.h"
#include <Rcpp.h>



void split(double *cdwd, double *newcdwd, int dim, double *stddev)
{
  double mult_offset=0.1;
  int i;

  /** if cdwd is way off zero and the variance is small     **/
  /** the offset may be out of the data range and generates **/
  /** empty cells.                                          **/
  /*****
  for (i=0; i<dim; i++) {
    newcdwd[i] = cdwd[i]*(1+mult_offset*drand48());
  }
  *****/

  /* set random range to [0.25, 0.75] */
  for (i=0; i<dim; i++) {
    //newcdwd[i] = cdwd[i]+stddev[i]*mult_offset*(0.25+drand48()/2.0);
    newcdwd[i] = cdwd[i]+stddev[i]*mult_offset*(0.25+R::runif(0,1)/2.0);
  }
}

// index - array with cluster id
void centroid(double *cdbk, int dim, int numcdwd, double *vc, 
	 int *index, int numdata)
{
  int i,j,k;
  int *ct;

  ct=(int *)malloc(numcdwd*sizeof(int)); // number of points in clusters

  if (index==NULL)
    {
      for (k=0; k<dim; k++)
	cdbk[k] = 0.0;
      for (i=0; i<numdata; i++) 
	for (k=0; k<dim; k++) 
	  cdbk[k]+=vc[i*dim+k];
      for (k=0; k<dim; k++)
	cdbk[k] /= ((double)numdata);
    }
  else
    {
      for (j=0; j<numcdwd; j++) {
	for (k=0; k<dim; k++)
	  cdbk[j*dim+k] = 0.0;
	ct[j] = 0;
      }
      
      for (i=0; i<numdata; i++) {
	for (k=0; k<dim; k++) 
	  cdbk[index[i]*dim+k]+=vc[i*dim+k];
	ct[index[i]]++;
      }
      
      for (j=0; j<numcdwd; j++) 
	for (k=0; k<dim; k++)
	  cdbk[j*dim+k] /= ((double)ct[j]);
    }

  free(ct);

}  

void cellstdv(double *cdbk, double *stddev, int dim, int numcdwd, double *vc,
	      int *index,  int numdata)
{
  int i,j,k;
  int *ct;

  ct=(int *)malloc(numcdwd*sizeof(int));

  for (j=0; j<numcdwd; j++) {
    for (k=0; k<dim; k++)
      stddev[j*dim+k] = 0.0;
    ct[j] = 0;
  }
      
  for (i=0; i<numdata; i++) {
    for (k=0; k<dim; k++) 
      stddev[index[i]*dim+k]+=((vc[i*dim+k]-cdbk[index[i]*dim+k])*
	(vc[i*dim+k]-cdbk[index[i]*dim+k]));
    ct[index[i]]++;
  }
  
  for (j=0; j<numcdwd; j++) { 
    if (ct[j]>0) {
      for (k=0; k<dim; k++) {
	stddev[j*dim+k] /= ((double)ct[j]);
	stddev[j*dim+k]=sqrt(stddev[j*dim+k]);
      }
    }
    else {
      for (k=0; k<dim; k++) stddev[j*dim+k]=1.0;
    }
  }

  free(ct);

}  


double mse_dist(double *cdwd, double *vc, int dim)
{
  double mse= 0.0;
  int i;

  for (i=0; i<dim; i++)
    mse += (vc[i]-cdwd[i])*(vc[i]-cdwd[i]);

  return(mse);
}

// find cluster id as id of the nearest cluster to datapoint
void encode(double *cdbk, int dim, int numcdwd, double *vc, int *code,
	    int numdata)
{
  int i,j;
  double *dist,minv;

  dist=(double *)calloc(numcdwd,sizeof(double));

  for (i=0; i<numdata; i++) {
    for (j=0; j<numcdwd;j++)
      dist[j]=mse_dist(cdbk+j*dim, vc+i*dim, dim);    
    code[i]=0;
    minv=dist[0];
    for (j=1;j<numcdwd;j++)
      if (dist[j]<minv) {
	minv=dist[j];
	code[i]=j;
      }
  }

  free(dist);
}


double lloyd(double *cdbk, int dim, int numcdwd, double *vc, int numdata, 
	    double threshold)
     // cdbk - coordinates of centroids
     // vc - array with data in variable block
     // numcdwd - number of states in variable block
     // dim - dimensionality of data in variable block
     // index - array with cluster id
     /* return the value of the mean squared distance, i.e., average */
     /* squared distance between a sample and its centroid           */
     /* threshold is for controling when to stop the loop */
{
  int i,j,k,m;
  int ite;
  double olddist, minmse, mse;
  double dist = 0;
  int min_iteration=2;
  /*double threshold = 0.005;*/
  int *index; // array with state id of sample points
  double *tp;
  int cdbksz2; // number of times clusters split
  int temp_int, new_cdwds, cdbksz;
  double *stddev; // standard deviation for appropriate split

  //srand48(0);//5/26/2017, to ensure identical result when rerun
  
  cdbksz2 = 0;
  temp_int = 1;
  while (temp_int < numcdwd) {
    cdbksz2++;
    temp_int += temp_int;
  }  


  index = (int *)calloc(numdata, sizeof(int));
  stddev = (double *)calloc(numcdwd*dim,sizeof(double));

  centroid(cdbk, dim, 1, vc, NULL, numdata); // centroid coordinate for single cluster

  /* compute standard deviation for each cell */
  for (i=0;i<numdata;i++) index[i]=0;
  cellstdv(cdbk,stddev,dim,numcdwd,vc,index,numdata);

  if (numcdwd==1) {
    dist = 0.0;
    for (i=0, k=0; i<numdata; i++)
      for (j=0; j<dim; j++, k++)
	dist += (cdbk[j]-vc[k])*(cdbk[j]-vc[k]);
    dist /= ((double)numdata);
  }

  cdbksz = 1;

  for (ite=0; ite<cdbksz2; ite++) {
    new_cdwds = ((numcdwd - 2*cdbksz) >= 0) ? cdbksz : numcdwd - cdbksz;

    for (k=0; k<new_cdwds; k++)
      split(cdbk+k*dim, cdbk+(cdbksz+k)*dim, dim, stddev+k*dim);

    cdbksz += new_cdwds;

    dist=HUGE_VAL;
    m = 0;

    while (m < min_iteration || 
	   (fabs((double)(olddist-dist))/olddist > threshold
	    && dist < olddist))
      {
	m++;
	olddist = dist;
	tp = vc;
	dist = 0.0;
	// assign cluster id as id of the nearest cluster
	for (i=0; i<numdata; i++)
	  {
	    minmse = mse_dist(cdbk, tp, dim);
	    index[i]= 0;
	    
	    for (j=1; j<cdbksz; j++)
	      {
		mse = mse_dist(cdbk+j*dim, tp, dim);
		if (mse<minmse)
		  {
		    minmse=mse;
		    index[i]=j;
		  }
	      }
	    
	    dist += minmse;
	    tp += dim;
	  }
	
	dist /= ((double)numdata);

	centroid(cdbk, dim, cdbksz, vc, index, numdata);
      }
    cellstdv(cdbk,stddev,dim,cdbksz,vc,index,numdata);

  }

  free(index);
  free(stddev);

  return(dist);

}	    


/** The kmeans algorithm is close to lloyd, except that the number **/
/** codewords is not preselected.  Instead, it is determined by    **/
/** the minimum number of codewords that lead to a mean squared    **/
/** error below a given threshold.  The number of codewords is     **/
/** upper bounded by the given maxnumcdwd.                         **/

double kmeans(double *cdbk, int dim, int maxnumcdwd, int *fnumcdwd, 
	     double *vc, int numdata, double threshold, double distthred)
{
  int i,j,k;
  double dist;
 
  int numcdwd = 1;

  centroid(cdbk, dim, numcdwd, vc, NULL, numdata);

  dist = 0.0;
  for (i=0, k=0; i<numdata; i++)
    for (j=0; j<dim; j++, k++)
      dist += (cdbk[j]-vc[k])*(cdbk[j]-vc[k]);
  dist /= ((double)numdata);
  if (dist<distthred) {
    *fnumcdwd=1;
    return(dist);
  }

  numcdwd=2;
  do {
    dist=lloyd(cdbk, dim, numcdwd, vc, numdata, threshold);
    numcdwd++;
    //fprintf(stderr, "numcdwd=%d, dist=%f\n", numcdwd,dist);
  } while (numcdwd<=maxnumcdwd && dist > distthred);
  
  *fnumcdwd=numcdwd-1;
  
  return(dist);
}	    

