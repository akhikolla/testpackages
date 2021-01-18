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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include <Rcpp.h>


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*-------------------------- General vector ------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

unsigned char vector_uchar(unsigned char **mt, int dim)
{
  unsigned char *tp;

  tp = (unsigned char *)calloc(dim, sizeof(unsigned char));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in vector_uchar\n";
    return(0);
  }

  *mt = tp;
  return(1);
}

unsigned char vector_float(float **mt, int dim)
{
  float *tp;

  tp = (float *)calloc(dim, sizeof(float));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in vector_float\n";
    return(0);
  }

  *mt = tp;
  return(1);
}

unsigned char vector_double(double **mt, int dim)
{
  double *tp;

  tp = (double *)calloc(dim, sizeof(double));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in vector_double\n";
    return(0);
  }

  *mt = tp;
  return(1);
}

unsigned char vector_int(int **mt, int dim)
{
  int *tp;

  tp = (int *)calloc(dim, sizeof(int));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in vector_int\n";
    return(0);
  }

  *mt = tp;
  return(1);
}


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*---------------- General 2 dimension matrix ----------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

unsigned char matrix_2d_uchar(unsigned char ***mt, int rows, int cols)
{
  unsigned char **tp;
  int i;

  tp = (unsigned char **)calloc(rows, sizeof(unsigned char *));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in matrix_2d_uchar\n";
    return(0);
  }

  if (cols==0) {
    *mt = tp;
    return(1);
  }

  for (i=0; i<rows; i++) {
    tp[i] = (unsigned char *)calloc(cols, sizeof(unsigned char));
    if (tp[i]==NULL) {
      Rcpp::Rcout << "Can't allocate space in matrix_2d_uchar\n";
      return(0);
    }
  }

  *mt = tp;
  return(1);
}
  
unsigned char matrix_2d_float(float ***mt, int rows, int cols)
{
  float **tp;
  int i;

  tp = (float **)calloc(rows, sizeof(float *));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in matrix_2d_float\n";
    return(0);
  }

  if (cols==0) {
    *mt = tp;
    return(1);
  }

  for (i=0; i<rows; i++) {
    tp[i] = (float *)calloc(cols, sizeof(float));
    if (tp[i]==NULL) {
      Rcpp::Rcout << "Can't allocate space in matrix_2d_float\n";
      return(0);
    }
  }

  *mt = tp;
  return(1);
}
  
unsigned char matrix_2d_double(double ***mt, int rows, int cols)
{
  double **tp;
  int i;

  tp = (double **)calloc(rows, sizeof(double *));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in matrix_2d_double\n";
    return(0);
  }

  if (cols==0) {
    *mt = tp;
    return(1);
  }

  for (i=0; i<rows; i++) {
    tp[i] = (double *)calloc(cols, sizeof(double));
    if (tp[i]==NULL) {
      Rcpp::Rcout << "Can't allocate space in matrix_2d_double\n";
      return(0);
    }
  }

  *mt = tp;
  return(1);
}
  
unsigned char matrix_2d_int(int ***mt, int rows, int cols)
{
  int **tp;
  int i;

  tp = (int **)calloc(rows, sizeof(int *));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in matrix_2d_int\n";
    return(0);
  }

  if (cols==0) {
    *mt = tp;
    return(1);
  }

  for (i=0; i<rows; i++) {
    tp[i] = (int *)calloc(cols, sizeof(int));
    if (tp[i]==NULL) {
      Rcpp::Rcout << "Can't allocate space in matrix_2d_int\n";
      return(0);
    }
  }

  *mt = tp;
  return(1);
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*------------ Generate 3 dimension matrix ---------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

unsigned char matrix_3d_uchar(unsigned char ****mt, int rows, 
			      int cols, int depth)
{
  unsigned char ***tp;
  int i, j;

  tp = (unsigned char ***)calloc(rows, sizeof(unsigned char **));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in matrix_3d_uchar\n";
    return(0);
  }

  if (cols==0) {
    *mt = tp;
    return(1);
  }

  for (i=0; i<rows; i++) {
    tp[i] = (unsigned char **)calloc(cols, sizeof(unsigned char *));
    if (tp[i]==NULL) {
      Rcpp::Rcout << "Can't allocate space in matrix_3d_uchar\n";
      return(0);
    }
    if (depth==0) continue;
    for (j=0; j<cols; j++) {
      tp[i][j] = (unsigned char *)calloc(depth, sizeof(unsigned char));
      if (tp[i][j]==NULL) {
	Rcpp::Rcout << "Can't allocate space in matrix_3d_uchar\n";
	return(0);
      }
    }
  }

  *mt = tp;
  return(1);
}
  
unsigned char matrix_3d_float(float ****mt, int rows, 
			      int cols, int depth)
{
  float ***tp;
  int i, j;

  tp = (float ***)calloc(rows, sizeof(float **));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in matrix_3d_float\n";
    return(0);
  }

  if (cols==0) {
    *mt = tp;
    return(1);
  }

  for (i=0; i<rows; i++) {
    tp[i] = (float **)calloc(cols, sizeof(float *));
    if (tp[i]==NULL) {
      Rcpp::Rcout << "Can't allocate space in matrix_3d_float\n";
      return(0);
    }
    if (depth==0) continue;
    for (j=0; j<cols; j++) {
      tp[i][j] = (float *)calloc(depth, sizeof(float));
      if (tp[i][j]==NULL) {
	Rcpp::Rcout << "Can't allocate space in matrix_3d_float\n";
	return(0);
      }
    }
  }

  *mt = tp;
  return(1);
}
  

unsigned char matrix_3d_double(double ****mt, int rows, 
			      int cols, int depth)
{
  double ***tp;
  int i, j;

  tp = (double ***)calloc(rows, sizeof(double **));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in matrix_3d_double\n";
    return(0);
  }

  if (cols==0) {
    *mt = tp;
    return(1);
  }

  for (i=0; i<rows; i++) {
    tp[i] = (double **)calloc(cols, sizeof(double *));
    if (tp[i]==NULL) {
      Rcpp::Rcout << "Can't allocate space in matrix_3d_double\n";
      return(0);
    }
    if (depth==0) continue;
    for (j=0; j<cols; j++) {
      tp[i][j] = (double *)calloc(depth, sizeof(double));
      if (tp[i][j]==NULL) {
	Rcpp::Rcout << "Can't allocate space in matrix_3d_double\n";
	return(0);
      }
    }
  }

  *mt = tp;
  return(1);
}
  
unsigned char matrix_3d_int(int ****mt, int rows, 
			      int cols, int depth)
{
  int ***tp;
  int i, j;

  tp = (int ***)calloc(rows, sizeof(int **));
  if (tp==NULL) {
    Rcpp::Rcout << "Can't allocate space in matrix_3d_int\n";
    return(0);
  }

  if (cols==0) {
    *mt = tp;
    return(1);
  }

  for (i=0; i<rows; i++) {
    tp[i] = (int **)calloc(cols, sizeof(int *));
    if (tp[i]==NULL) {
      Rcpp::Rcout << "Can't allocate space in matrix_3d_int\n";
      return(0);
    }
    if (depth==0) continue;
    for (j=0; j<cols; j++) {
      tp[i][j] = (int *)calloc(depth, sizeof(int));
      if (tp[i][j]==NULL) {
	Rcpp::Rcout << "Can't allocate space in matrix_3d_int\n";
	return(0);
      }
    }
  }

  *mt = tp;
  return(1);
}
  

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*--------------- Free 2 dimension matrix ---------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void free_matrix_2d_uchar(unsigned char ***mt, int rows)
{
  unsigned char **tp;
  int i;

  tp = *mt;
  for (i=0; i<rows; i++)
    free(tp[i]);
  
  free(tp);
  
  *mt = NULL;
}

void free_matrix_2d_float(float ***mt, int rows)
{
  float **tp;
  int i;

  tp = *mt;
  for (i=0; i<rows; i++)
    free(tp[i]);
  
  free(tp);
  
  *mt = NULL;
}

void free_matrix_2d_double(double ***mt, int rows)
{
  double **tp;
  int i;

  tp = *mt;
  for (i=0; i<rows; i++)
    free(tp[i]);
  
  free(tp);
  
  *mt = NULL;
}

void free_matrix_2d_int(int ***mt, int rows)
{
  int **tp;
  int i;

  tp = *mt;
  for (i=0; i<rows; i++)
    free(tp[i]);
  
  free(tp);
  
  *mt = NULL;
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*--------------- Free 3 dimension matrix ---------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void free_matrix_3d_uchar(unsigned char ****mt, int rows, int cols)
{
  unsigned char ***tp;
  int i,j;

  tp = *mt;
  for (i=0; i<rows; i++)
    for (j=0; j<cols; j++)
      free(tp[i][j]);
  for (i=0; i<rows; i++)
    free(tp[i]);

  free(tp);
  
  *mt = NULL;
}

void free_matrix_3d_float(float ****mt, int rows, int cols)
{
  float ***tp;
  int i,j;

  tp = *mt;
  for (i=0; i<rows; i++)
    for (j=0; j<cols; j++)
      free(tp[i][j]);
  for (i=0; i<rows; i++)
    free(tp[i]);

  free(tp);
  
  *mt = NULL;
}

void free_matrix_3d_double(double ****mt, int rows, int cols)
{
  double ***tp;
  int i,j;

  tp = *mt;
  for (i=0; i<rows; i++)
    for (j=0; j<cols; j++)
      free(tp[i][j]);
  for (i=0; i<rows; i++)
    free(tp[i]);

  free(tp);
  
  *mt = NULL;
}

void free_matrix_3d_int(int ****mt, int rows, int cols)
{
  int ***tp;
  int i,j;

  tp = *mt;
  for (i=0; i<rows; i++)
    for (j=0; j<cols; j++)
      free(tp[i][j]);
  for (i=0; i<rows; i++)
    free(tp[i]);

  free(tp);
  
  *mt = NULL;
}



float mat_det_float(float **mt, int dim)
{
  int i,j,n;
  float res;
  float **submt;
  float *ptr1, *ptr2;

  if (dim==1)
    return(**mt);
  
  
  try {
    if (!matrix_2d_float(&submt, dim-1, dim-1))  {
      throw std::range_error("Couldn't allocate memory in matrix_2d_float!");
    }
  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  //if (!matrix_2d_float(&submt, dim-1, dim-1)) 
	//Rcpp::stop("Couldn't allocate memory in matrix_2d_float!\n");
    //exit(1);

  for (i=1; i<dim; i++) {
    ptr1 = submt[i-1];
    ptr2 = mt[i]+1;
    for (j=1; j<dim; j++) {
      *(ptr1++) = *ptr2;
      ptr2++;
    }
  }

  n=1;
  res = 0.0;
 
  for (i=0; i<dim; i++) {
    res += (n*mt[i][0]*mat_det_float(submt, dim-1));
    n = -n;
    if (i==dim-1)
      continue;
    ptr1 = submt[i];
    ptr2 = mt[i]+1;
    for (j=1; j<dim; j++) {
      *(ptr1++) = *ptr2;
      ptr2++;
    }
  }

  free_matrix_2d_float(&submt, dim-1);

  return(res);
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*----------------------- Set Memory for Vector ---------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void memcpy_1d_uchar(unsigned char *mt, int dim, unsigned char tp)
{
  int i;

  for (i=0; i<dim; i++)
    *(mt++) = tp;
}

void memcpy_1d_int(int *mt, int dim, int tp)
{
  int i;

  for (i=0; i<dim; i++)
    *(mt++) = tp;
}

void memcpy_1d_float(float *mt, int dim, float tp)
{
  int i;

  for (i=0; i<dim; i++)
    *(mt++) = tp;
}

void memcpy_1d_double(double *mt, int dim, double tp)
{
  int i;

  for (i=0; i<dim; i++)
    *(mt++) = tp;
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*----------------------- Set Memory for 2D matrix ------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void memcpy_2d_uchar(unsigned char **mt, int rows, int cols, unsigned char tp)
{
  int i;

  for (i=0; i<rows; i++)
    memcpy_1d_uchar(mt[i],cols,tp);
}

void memcpy_2d_int(int **mt, int rows, int cols, int tp)
{
  int i;

  for (i=0; i<rows; i++)
    memcpy_1d_int(mt[i],cols,tp);
}

void memcpy_2d_float(float **mt, int rows, int cols, float tp)
{
  int i;

  for (i=0; i<rows; i++)
    memcpy_1d_float(mt[i],cols,tp);
}

void memcpy_2d_double(double **mt, int rows, int cols, double tp)
{
  int i;

  for (i=0; i<rows; i++)
    memcpy_1d_double(mt[i],cols,tp);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*----------------------- Set Memory for 3D matrix ------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void memcpy_3d_uchar(unsigned char ***mt, int rows, int cols, int depth, 
		     unsigned char tp)
{
  int i;

  for (i=0; i<rows; i++)
    memcpy_2d_uchar(mt[i],cols,depth,tp);
}

void memcpy_3d_int(int ***mt, int rows, int cols, int depth, int tp)
{
  int i;

  for (i=0; i<rows; i++)
    memcpy_2d_int(mt[i],cols,depth,tp);
}

void memcpy_3d_float(float ***mt, int rows, int cols, int depth, float tp)
{
  int i;

  for (i=0; i<rows; i++)
    memcpy_2d_float(mt[i],cols,depth,tp);
}

void memcpy_3d_double(double ***mt, int rows, int cols, int depth, double tp)
{
  int i;

  for (i=0; i<rows; i++)
    memcpy_2d_double(mt[i],cols,depth,tp);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*------------------------ Vector copy ------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void vector_cpy_uchar(unsigned char *output, unsigned char *input, int dim)
{
  int i;

  for (i=0; i<dim; i++)
    {
      *(output++)=*input;
      input++;
    }
}

void vector_cpy_int(int *output, int *input, int dim)
{
  int i;

  for (i=0; i<dim; i++)
    {
      *(output++)=*input;
      input++;
    }
}

void vector_cpy_float(float *output, float *input, int dim)
{
  int i;

  for (i=0; i<dim; i++)
    {
      *(output++)=*input;
      input++;
    }
}

void vector_cpy_double(double *output, double *input, int dim)
{
  int i;

  for (i=0; i<dim; i++)
    {
      *(output++)=*input;
      input++;
    }
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*------------------------ Matrix copy ------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void matrix_2d_cpy_uchar(unsigned char **output, unsigned char **input, 
		      int rows, int cols)
{
  int i;

  for (i=0; i<rows; i++)
    vector_cpy_uchar(output[i], input[i],cols);
}

void matrix_2d_cpy_int(int **output, int **input, int rows, int cols)
{
  int i;

  for (i=0; i<rows; i++)
    vector_cpy_int(output[i], input[i],cols);
}

void matrix_2d_cpy_float(float **output, float **input, int rows, int cols)
{
  int i;

  for (i=0; i<rows; i++)
    vector_cpy_float(output[i], input[i],cols);
}

void matrix_2d_cpy_double(double **output, double **input, int rows, int cols)
{
  int i;

  for (i=0; i<rows; i++)
    vector_cpy_double(output[i], input[i],cols);
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*------------------------ Print out Matrix -------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void print_matrix_uchar(unsigned char **mt, int rows, int cols)
{
  int i,j;

  for (i=0; i<rows; i++) {
    for (j=0; j<cols; j++) {
      Rcpp::Rcout << mt[i][j] << " ";
      if ((j+1)%8==0)
	Rcpp::Rcout <<  "\n";
    }
    Rcpp::Rcout <<  "\n";
  }
}

void print_matrix_int(int **mt, int rows, int cols)
{
  int i,j;

  for (i=0; i<rows; i++) {
    for (j=0; j<cols; j++) {
      Rcpp::Rcout << mt[i][j] << " ";
      if ((j+1)%8==0)
	Rcpp::Rcout <<  "\n";
    }
    Rcpp::Rcout <<  "\n";
  }
}

void print_matrix_float(float **mt, int rows, int cols)
{
  int i,j;

  for (i=0; i<rows; i++) {
    for (j=0; j<cols; j++) {
      Rcpp::Rcout << mt[i][j] << " ";
      if ((j+1)%8==0)
	Rcpp::Rcout <<  "\n";
    }
    Rcpp::Rcout <<  "\n";
  }
}

void print_matrix_double(double **mt, int rows, int cols)
{
  int i,j;

  for (i=0; i<rows; i++) {
    for (j=0; j<cols; j++) {
     Rcpp::Rcout << mt[i][j] << " ";
      if ((j+1)%8==0)
	Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
  }
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*---------------      Numerical Programs             ---------------*/
/*---------------    LU decomposition programs        ---------------*/
/*---------------      Calculate matrix inversion     ---------------*/
/*---------------      Calculate matrix determinant   ---------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

double mat_det_double(double **mt, int dim)
{
  int i,j,n;
  double res;
  double **submt;
  double *ptr1, *ptr2;

  if (dim==1)
    return(**mt);
  
  n = dim-1;
  
  try {
    if (!matrix_2d_double(&submt, dim-1, dim-1))  {
      throw std::range_error("Couldn't allocate memory in matrix_2d_double!");
    }
  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  //if (!matrix_2d_double(&submt, dim-1, dim-1)) 
	//Rcpp::stop("Couldn't allocate memory in matrix_2d_double!\n");
    //exit(1);

  for (i=1; i<dim; i++) {
    ptr1 = submt[i-1];
    ptr2 = mt[i]+1;
    for (j=1; j<dim; j++) {
      *(ptr1++) = *ptr2;
      ptr2++;
    }
  }

  n=1;
  res = 0.0;
 
  for (i=0; i<dim; i++) {
    res += (n*mt[i][0]*mat_det_double(submt, dim-1));
    n = -n;
    if (i==dim-1)
      continue;
    ptr1 = submt[i];
    ptr2 = mt[i]+1;
    for (j=1; j<dim; j++) {
      *(ptr1++) = *ptr2;
      ptr2++;
    }
  }

  free_matrix_2d_double(&submt, dim-1);

  return(res);
}

unsigned char ludcmp_float(float **a, int n, int *indx, float *d)
{
  int i,imax,j,k;
  float big, dum,sum,temp;
  float *vv;
  float TINY=1e-20;

  if (!vector_float(&vv, n))
    return(0);

  *d = 1.0;
  for (i=0; i<n; i++) {
    big =0.0;
    for (j=0; j<n; j++)
      if ((temp=fabs(a[i][j]))>big) 
	big = temp;
    if (big==0.0) {
      Rcpp::Rcout << "Singular matrix in ludcmp_float\n";
      free(vv);
      return(2); /* 2 stands for singular matrix */
    }
    vv[i]=1.0/big;
  }

  for (j=0; j<n; j++) {
    for (i=0; i<j; i++) {
      sum = a[i][j];
      for (k=0; k<i; k++)
	sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i=j; i<n; i++) {
      sum = a[i][j];
      for (k=0; k<j; k++)
	sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum=vv[i]*fabs(sum))>=big) {
	big = dum;
	imax = i;
      }
    }
    if (j!=imax) {
      for (k=0; k<n; k++) {
	dum = a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j]==0.0) a[j][j] = TINY;
    if (j!=n-1) {
      dum = 1.0/(a[j][j]);
      for (i=j+1; i<n; i++) 
	a[i][j] *= dum;
    }
  }

  free(vv);
  return(1);
}


unsigned char ludcmp_double(double **a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big, dum,sum,temp;
  double *vv;
  double TINY=1e-20;

  if (!vector_double(&vv, n))
    return(0);

  *d = 1.0;
  for (i=0; i<n; i++) {
    big =0.0;
    for (j=0; j<n; j++)
      if ((temp=fabs(a[i][j]))>big) 
	big = temp;
    if (big==0.0) {
      Rcpp::Rcout << "Singular matrix in ludcmp_double" << std::endl;
      free(vv);
      return(2); /* 2 stands for singular matrix */
    }
    vv[i]=1.0/big;
  }

  for (j=0; j<n; j++) {
    for (i=0; i<j; i++) {
      sum = a[i][j];
      for (k=0; k<i; k++)
	sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i=j; i<n; i++) {
      sum = a[i][j];
      for (k=0; k<j; k++)
	sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum=vv[i]*fabs(sum))>=big) {
	big = dum;
	imax = i;
      }
    }
    if (j!=imax) {
      for (k=0; k<n; k++) {
	dum = a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j]==0.0) a[j][j] = TINY;
    if (j!=n-1) {
      dum = 1.0/(a[j][j]);
      for (i=j+1; i<n; i++) 
	a[i][j] *= dum;
    }
  }

  free(vv);
  return(1);
}

void lubksb_float(float **a, int n, int *indx, float *b)
{
  int i, ii=-1, ip,j;
  float sum;
  
  for (i=0; i<n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (j=ii; j<i; j++) sum -= a[i][j]*b[j];
    else if (sum != 0.0) ii=i;
    b[i]=sum;
  }
  for (i=n-1; i>=0; i--) {
    sum=b[i];
    for (j=i+1; j<n; j++) sum-= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


void lubksb_double(double **a, int n, int *indx, double *b)
{
  int i, ii=-1, ip,j;
  double sum;
  
  for (i=0; i<n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (j=ii; j<i; j++) sum -= a[i][j]*b[j];
    else if (sum != 0.0) ii=i;
    b[i]=sum;
  }
  for (i=n-1; i>=0; i--) {
    sum=b[i];
    for (j=i+1; j<n; j++) sum-= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


unsigned char mat_inv_float(float **mt, float **y, int dim)
{
  float d, *col;
  int i,j,*indx;
  float **a;
  float *ptr1, *ptr2;

  if (!matrix_2d_float(&a, dim,dim))
    return(0);

  for (i=0; i<dim; i++) {
    ptr1 = mt[i];
    ptr2 = a[i];
    for (j=0; j<dim; j++) {
      *(ptr2++) = *ptr1;
      ptr1++;
    }
  }

  if (!vector_float(&col, dim))
    return(0);

  if (!vector_int(&indx, dim))
    return(0);

  ludcmp_float(a,dim,indx,&d);
  for (j=0; j<dim; j++) {
    for (i=0; i<dim; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_float(a,dim,indx,col);
    for (i=0; i<dim;i++) 
      y[i][j]=col[i];
  }

  free(col);
  free(indx);
  free_matrix_2d_float(&a, dim);
  return(1);
}

unsigned char mat_inv_double(double **mt, double **y, int dim)
{
  double d, *col;
  int i,j,*indx;
  double **a;
  double *ptr1, *ptr2;

  if (!matrix_2d_double(&a, dim,dim))
    return(0);

  for (i=0; i<dim; i++) {
    ptr1 = mt[i];
    ptr2 = a[i];
    for (j=0; j<dim; j++) {
      *(ptr2++) = *ptr1;
      ptr1++;
    }
  }

  if (!vector_double(&col, dim))
    return(0);

  if (!vector_int(&indx, dim))
    return(0);

  ludcmp_double(a,dim,indx,&d);
  for (j=0; j<dim; j++) {
    for (i=0; i<dim; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_double(a,dim,indx,col);
    for (i=0; i<dim;i++) 
      y[i][j]=col[i];
  }

  free(col);
  free(indx);
  free_matrix_2d_double(&a, dim);
  return(1);
}

float mat_det_ludcmp_float(float **mt, int dim)
{
  float d;
  int i,j, *indx;
  float **a;
  float *ptr1, *ptr2;

  if (!matrix_2d_float(&a, dim,dim))
    return(0);

  for (i=0; i<dim; i++) {
    ptr1 = mt[i];
    ptr2 = a[i];
    for (j=0; j<dim; j++) {
      *(ptr2++) = *ptr1;
      ptr1++;
    }
  }

  try {
    if (!vector_int(&indx, dim))  {
      throw std::range_error("Couldn't allocate memory in vector_int!");
    }
  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  //if (!vector_int(&indx, dim))
	//Rcpp::stop("Couldn't allocate memory in vector_int!\n");
    
    //exit(1);

  ludcmp_float(a,dim,indx,&d);
  for (j=0; j<dim; j++) d *= a[j][j];

  free(indx);
  free_matrix_2d_float(&a, dim);
  return(d);
}

double mat_det_ludcmp_double(double **mt, int dim)
{
  double d;
  int i, j, *indx;
  double **a;
  double *ptr1, *ptr2;

  if (!matrix_2d_double(&a, dim,dim))
    return(0);

  for (i=0; i<dim; i++) {
    ptr1 = mt[i];
    ptr2 = a[i];
    for (j=0; j<dim; j++) {
      *(ptr2++) = *ptr1;
      ptr1++;
    }
  }

  try {
    if (!vector_int(&indx, dim))  {
      throw std::range_error("Couldn't allocate memory in vector_int!");
    }
  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  //if (!vector_int(&indx, dim))
	//Rcpp::stop("Couldn't allocate memory in vector_int!\n");
    
    //exit(1);

  ludcmp_double(a,dim,indx,&d);

  for (j=0; j<dim; j++) d *= a[j][j];

  free(indx);
  free_matrix_2d_double(&a, dim);
  return(d);
}


/** compute the determinant and at the same provide the inverse matrix **/
/** combining these two operations save computation since they share   **/
/** the call to ludcmp_double                                          **/

unsigned char mat_det_inv_double(double **mt, double **y, double *det, 
				 int dim)
{
  double d, *col;
  int i,j,m,*indx;
  double **a;
  double *ptr1, *ptr2;

  /** initialize matrix determinant **/
  *det=0.0;

  if (!matrix_2d_double(&a, dim,dim))
    return(0);

  for (i=0; i<dim; i++) {
    ptr1 = mt[i];
    ptr2 = a[i];
    for (j=0; j<dim; j++) {
      *(ptr2++) = *ptr1;
      ptr1++;
    }
  }

  if (!vector_double(&col, dim))
    return(0);

  if (!vector_int(&indx, dim))
    return(0);

  m=ludcmp_double(a,dim,indx,&d);
  if (m==2) {  /** singular matrix **/
    Rcpp::Rcout << "Singular matrix in mat_det_inv_double" << std::endl;
    *det=0.0;
    free(col);
    free(indx);
    free_matrix_2d_double(&a, dim);
    return(2);
  }

  for (j=0; j<dim; j++) d *= a[j][j];
  *det=d;
  
  for (j=0; j<dim; j++) {
    for (i=0; i<dim; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_double(a,dim,indx,col);
    for (i=0; i<dim;i++) 
      y[i][j]=col[i];
  }

  free(col);
  free(indx);
  free_matrix_2d_double(&a, dim);
  return(1);
}

unsigned char mat_det_inv_diag_double(double **mt, double **y, double *det, 
				      int dim, int diagonal)
{
  int i,j;
 
  if (diagonal!=1) {
    return(mat_det_inv_double(mt,y,det,dim));
  }

  /** initialize matrix determinant **/
  *det=1.0;
  for (i=0;i<dim;i++) (*det) *= mt[i][i];
  for (i=0;i<dim;i++)
    for (j=0;j<dim;j++) y[i][j]=0.0;
  for (i=0;i<dim;i++)
    y[i][i]=1.0/mt[i][i];

  if (*det==0.0) return(2);
  else return(1);
}

unsigned char mat_det_inv_float(float **mt, float **y, float *det, 
				 int dim)
{
  float d, *col;
  int i,j,m,*indx;
  float **a;
  float *ptr1, *ptr2;

  /** initialize matrix determinant **/
  *det=0.0;

  if (!matrix_2d_float(&a, dim,dim))
    return(0);

  for (i=0; i<dim; i++) {
    ptr1 = mt[i];
    ptr2 = a[i];
    for (j=0; j<dim; j++) {
      *(ptr2++) = *ptr1;
      ptr1++;
    }
  }

  if (!vector_float(&col, dim))
    return(0);

  if (!vector_int(&indx, dim))
    return(0);

  m=ludcmp_float(a,dim,indx,&d);
  if (m==2) {  /** singular matrix **/
    *det=0.0;
    free(col);
    free(indx);
    free_matrix_2d_float(&a, dim);
    return(2);
  }

  for (j=0; j<dim; j++) d *= a[j][j];
  *det=d;

  for (j=0; j<dim; j++) {
    for (i=0; i<dim; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_float(a,dim,indx,col);
    for (i=0; i<dim;i++) 
      y[i][j]=col[i];
  }

  free(col);
  free(indx);
  free_matrix_2d_float(&a, dim);
  return(1);
}

