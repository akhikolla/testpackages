/*******************************************************************************
 * Copyright (c) 2020, College of William & Mary
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the College of William & Mary nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * PRIMME: https://github.com/primme/primme
 * Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
 *******************************************************************************
 * File: primme_ls.h
 * 
 * Purpose - Main header with the PRIMME LS C interface functions.
 * 
 ******************************************************************************/

#ifndef PRIMME_LS_H
#define PRIMME_LS_H

#include <stdio.h>

#include "primme.h" // cyclic

typedef struct primme_ls_stats {
   PRIMME_INT numIterations;
   PRIMME_INT numGlobalSum;         /* times called globalSumReal */
   PRIMME_INT numBroadcast;         /* times called broadcastReal */
   PRIMME_INT volumeGlobalSum;      /* number of SCALARs reduced by globalSumReal */
   PRIMME_INT volumeBroadcast;      /* number of SCALARs broadcast by broadcastReal */
   double flopsDense;               /* FLOPS done by Num_update_VWXR_Sprimme */
   double elapsedTime; 
   double timeMatvec;               /* time expend by matrixMatvec */
   double timePrecond;              /* time expend by applyPreconditioner */
   double timeGlobalSum;            /* time expend by globalSumReal  */
   double timeBroadcast;            /* time expend by broadcastReal  */
   double timeDense;                /* time expend by Num_update_VWXR_Sprimme */
} primme_ls_stats;


/*--------------------------------------------------------------------------*/
typedef struct primme_ls_params {

   /* The user must input at least the following two arguments */
   PRIMME_INT n;
   void (*matrixMatvec)
      ( void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
        struct primme_params *primme, int *ierr);
   primme_op_datatype matrixMatvec_type; /* expected type of x and y */

   /* Preconditioner applied on block of vectors (if available) */
   void (*applyPreconditioner)
      ( void *x, PRIMME_INT *ldx,  void *y, PRIMME_INT *ldy, int *blockSize,
        struct primme_params *primme, int *ierr);
   primme_op_datatype applyPreconditioner_type; /* expected type of x and y */

   /* input for the following is only required for parallel programs */
   int numProcs;
   int procID;
   PRIMME_INT nLocal;
   void *commInfo;
   void (*globalSumReal)
      (void *sendBuf, void *recvBuf, int *count, struct primme_params *primme,
       int *ierr );
   primme_op_datatype globalSumReal_type; /* expected type of sendBuf and recvBuf */
   void (*broadcastReal)(
         void *buffer, int *count, struct primme_params *primme, int *ierr);
   primme_op_datatype broadcastReal_type; /* expected type of buffer */

   /*Though primme_initialize will assign defaults, most users will set these */
   int numRHS;         /* Number of right-hand-side to solve */
   int initSize;       /* Number of initial columns given */ 

   /* the following will be given default values depending on the method */
   PRIMME_INT maxIterations;
   PRIMME_INT iseed[4];
   double aNorm;                 /* Approximate 2-norm of the problem matrix */
   double tolerance;             /* Relative tolerance to the solutions relative to the matrix problem norm */
   primme_op_datatype internalPrecision; /* force primme to work in that precision */

   int printLevel;
   FILE *outputFile;

   void *matrix;                 /* User pointer for matrixMatvec */
   void *preconditioner;         /* User pointer for applyPreconditioner */
   PRIMME_INT ldevecs;           /* Leading dimension of the RHS */
   PRIMME_INT ldOPs;             /* Leading dimension of the input/output vectors given in matrixMatvec and applyPreconditioner */

   struct primme_ls_stats stats; /* Output stats */

   void (*convTestFun)(void *solution, double *rNorm, int *isconv, 
         struct primme_ls_params *primme, int *ierr);
   primme_op_datatype convTestFun_type; /* expected type of solution */
   void *convtest;
   void (*monitorFun)(int *numConverged,
         int *inner_its, void *LSRes, const char *msg, double *time,
         primme_event *event, struct primme_ls_params *primme, int *err);
   primme_op_datatype monitorFun_type; /* expected type of float-point arrays */
   void *monitor;
   void *queue;      /* magma device queue (magma_queue_t*) */
   const char *profile; /* regex expression with functions to monitor times */
} primme_params;
/*---------------------------------------------------------------------------*/



#endif /* PRIMME_LS_H */
