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

extern void split(double *cdwd, double *newcdwd, int dim, double *stdv);
extern void centroid(double *cdbk, int dim, int numcdwd, double *vc, 
		int *index, int numdata);
extern void cellstdv(double *cdbk, double *stdv, int dim, int numcdwd, 
		     double *vc, int *index,  int numdata);
extern double mse_dist(double *cdwd, double *vc, int dim);
extern void encode(double *cdbk, int dim, int numcdwd, double *vc, int *code,
		   int numdata);
extern double lloyd(double *cdbk, int dim, int numcdwd, double *vc, int numdata, 
		   double threshold);
extern double kmeans(double *cdbk, int dim, int maxnumcdwd, int *fnumcdwd, 
		    double *vc, int numdata, double threshold, double distthred);
