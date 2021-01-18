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

/*-------------------------- mat_simple.c --------------------------*/

/*------------------------------------------------------------------*/
/*-------------------------- General vector ------------------------*/
/*------------------------------------------------------------------*/

extern unsigned char vector_uchar(unsigned char **, int );
extern unsigned char vector_float(float **, int );
extern unsigned char vector_double(double **, int );
extern unsigned char vector_int(int **, int);

/*------------------------------------------------------------------*/
/*---------------- General 2 dimension matrix ----------------------*/
/*------------------------------------------------------------------*/

extern unsigned char matrix_2d_uchar(unsigned char ***, int , int );
extern unsigned char matrix_2d_float(float ***, int , int );
extern unsigned char matrix_2d_double(double ***, int , int );
extern unsigned char matrix_2d_int(int ***, int , int );

/*--------------------------------------------------------------------*/
/*------------ Generate 3 dimension matrix ---------------------------*/
/*--------------------------------------------------------------------*/

extern unsigned char matrix_3d_uchar(unsigned char ****, int ,int , int );
extern unsigned char matrix_3d_float(float ****, int ,int , int );
extern unsigned char matrix_3d_double(double ****, int ,int , int );
extern unsigned char matrix_3d_int(int ****, int , int , int );

/*-------------------------------------------------------------------*/
/*--------------- Free 2 dimension matrix ---------------------------*/
/*-------------------------------------------------------------------*/
extern void free_matrix_2d_uchar(unsigned char ***, int );
extern void free_matrix_2d_float(float ***, int );
extern void free_matrix_2d_double(double ***, int );
extern void free_matrix_2d_int(int ***, int );

/*-------------------------------------------------------------------*/
/*--------------- Free 3 dimension matrix ---------------------------*/
/*-------------------------------------------------------------------*/
extern void free_matrix_3d_uchar(unsigned char ****, int , int );
extern void free_matrix_3d_float(float ****, int , int );
extern void free_matrix_3d_double(double ****, int , int );
extern void free_matrix_3d_int(int ****, int , int );

/*-------------------------------------------------------------------*/
/*----------------------- Set Memory for Vector ---------------------*/
/*-------------------------------------------------------------------*/
extern void memcpy_1d_uchar(unsigned char *, int , unsigned char );
extern void memcpy_1d_int(int *, int , int );
extern void memcpy_1d_float(float *, int , float );
extern void memcpy_1d_double(double *, int , double );

/*-------------------------------------------------------------------*/
/*----------------------- Set Memory for 2D matrix ------------------*/
/*-------------------------------------------------------------------*/
extern void memcpy_2d_uchar(unsigned char **, int , int , unsigned char );
extern void memcpy_2d_int(int **, int , int , int );
extern void memcpy_2d_float(float **, int , int , float );
extern void memcpy_2d_double(double **, int , int , double );

/*-------------------------------------------------------------------*/
/*----------------------- Set Memory for 3D matrix ------------------*/
/*-------------------------------------------------------------------*/
extern void memcpy_3d_uchar(unsigned char ***, int, int, int,unsigned char);
extern void memcpy_3d_int(int ***, int , int , int , int);
extern void memcpy_3d_float(float ***, int , int , int , float);
extern void memcpy_3d_double(double ***, int , int , int , double);

/*-------------------------------------------------------------------*/
/*------------------------ Vector copy ------------------------------*/
/*-------------------------------------------------------------------*/
extern void vector_cpy_uchar(unsigned char *, unsigned char *, int );
extern void vector_cpy_int(int *, int *, int );
extern void vector_cpy_float(float *, float *, int );
extern void vector_cpy_double(double *, double *, int );

/*-------------------------------------------------------------------*/
/*------------------------ Matrix copy ------------------------------*/
/*-------------------------------------------------------------------*/
extern void matrix_2d_cpy_uchar(unsigned char **, unsigned char **, 
		      int , int );
extern void matrix_2d_cpy_int(int **, int **, int , int );
extern void matrix_2d_cpy_float(float **, float **, int , int );
extern void matrix_2d_cpy_double(double **, double **, int , int );

/*-------------------------------------------------------------------*/
/*------------------------ Print out Matrix -------------------------*/
/*-------------------------------------------------------------------*/
extern void print_matrix_uchar(unsigned char **, int , int );
extern void print_matrix_int(int **, int , int );
extern void print_matrix_float(float **, int , int );
extern void print_matrix_double(double **, int , int );

/*-------------------------------------------------------------------*/
/*---------------      Numerical Programs             ---------------*/
/*---------------    LU decomposition programs        ---------------*/
/*---------------      Calculate matrix inversion     ---------------*/
/*---------------      Calculate matrix determinant   ---------------*/
/*-------------------------------------------------------------------*/
extern float mat_det_float(float **mt, int dim);
extern double mat_det_double(double **mt, int dim);
extern unsigned char ludcmp_float(float **, int , int *, float *);
extern unsigned char ludcmp_double(double **, int , int *, double *);
extern void lubksb_float(float **, int , int *, float *);
extern void lubksb_double(double **, int , int *, double *);
extern unsigned char mat_inv_float(float **, float **, int );
extern unsigned char mat_inv_double(double **, double **, int );
extern float mat_det_ludcmp_float(float **, int );
extern double mat_det_ludcmp_double(double **, int );
extern unsigned char mat_det_inv_double(double **, double **, double *, int);
extern unsigned char mat_det_inv_diag_double(double **, double **, double *, int,int);
extern unsigned char mat_det_inv_float(float **, float **, float *, int);
