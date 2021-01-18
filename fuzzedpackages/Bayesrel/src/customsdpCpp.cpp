/*
 *  This is an easy to call version of the sdp routine.  It takes as
 *  input a problem (n,k,C,a,constraints,constant_offset), and an
 *  initial solution (X,y,Z), allocates working storage, and calls sdp()
 *  to solve the problem.  The solution is returned in X,y,Z,pobj,dobj, and
 *  the return code from sdp is returned as the return value from easy_sdp.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern "C" {
#include "declarations.h"
}
#include "customsdp.h"


int custom_sdpCpp(
     int n,
     int k,
     blockmatrix& C,
     double *a,
     struct constraintmatrix *constraints,
     double constant_offset,
     double *ppobj,
     double *pdobj,
     const arma::cube& car,
     arma::dvec& out,
     const int printlevel)
{
  int ret;
  struct constraintmatrix fill;
  struct paramstruc params;
  struct blockmatrix work1;
  struct blockmatrix work2;
  struct blockmatrix work3;
  struct blockmatrix bestx;
  struct blockmatrix bestz;
  struct blockmatrix Zi;
  struct blockmatrix dZ;
  struct blockmatrix dX;
  struct blockmatrix cholxinv;
  struct blockmatrix cholzinv;
  double *workvec1;
  double *workvec2;
  double *workvec3;
  double *workvec4;
  double *workvec5;
  double *workvec6;
  double *workvec7;
  double *workvec8;
  double *diagO;
  double *Fp;
  double *O;
  double *dy;
  double *dy1;
  double *rhs;
  double *besty;
  int ldam;
  struct sparseblock **byblocks;
  struct sparseblock *ptr;
  struct sparseblock *oldptr;
  int i;
  int j;
  int blk;
  struct sparseblock *p;
  struct sparseblock *q;
  struct sparseblock *prev=NULL;
  int nnz;
  struct blockmatrix X, Z;
  double *y;

   /*
    *  Initialize the parameters.
    */

   params.axtol=1.0e-8;
   params.atytol=1.0e-8;
   params.objtol=1.0e-8;
   params.pinftol=1.0e8;
   params.dinftol=1.0e8;
   params.maxiter=100;
   params.minstepfrac=0.90;
   params.maxstepfrac=0.97;
   params.minstepp=1.0e-8;
   params.minstepd=1.0e-8;
   params.usexzgap=1;
   params.tweakgap=0;
   params.affine=0;
   params.perturbobj=1;
   params.fastmode=0;
//   printlevel=0;

  /*
   *  Allocate working storage
   */

  // allocate storage for X, y , Z, that was previously done in initsoln
  alloc_mat(C,&X);
  alloc_mat(C,&Z);
  y=(double *)malloc(sizeof(double)*(k+1));

  alloc_mat(C,&work1);
  alloc_mat(C,&work2);
  alloc_mat(C,&work3);
  alloc_mat_packed(C,&bestx);
  alloc_mat_packed(C,&bestz);
  alloc_mat_packed(C,&cholxinv);
  alloc_mat_packed(C,&cholzinv);

  besty=(double *)malloc(sizeof(double)*(k+1));
   if (besty == NULL)
     {
       return(10);
     };

   if (n > k)
     {
       workvec1=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec1=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec1 == NULL)
     {
       return(10);
     };

   if (n > k)
     {
       workvec2=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec2=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec2 == NULL)
     {
       return(10);
     };

   if (n > k)
     {
       workvec3=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec3=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec3 == NULL)
     {
       return(10);
     };

   if (n > k)
     {
       workvec4=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec4=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec4 == NULL)
     {
       return(10);
     };

   if (n > k)
     {
       workvec5=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec5=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec5 == NULL)
     {
       return(10);
     };

   if (n > k)
     {
       workvec6=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec6=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec6 == NULL)
     {
       return(10);
     };

   if (n > k)
     {
       workvec7=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec7=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec7 == NULL)
     {
       return(10);
     };

   if (n > k)
     {
       workvec8=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec8=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec8 == NULL)
     {
       return(10);
     };


   if (n > k)
     {
       diagO=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       diagO=(double *)malloc(sizeof(double)*(k+1));
     };
   if (diagO == NULL)
     {
       return(10);
     };



   rhs=static_cast<double*>(malloc(sizeof(double)*(k+1)));
   if (rhs == NULL)
     {
       return(10);
     };

   dy=static_cast<double*>(malloc(sizeof(double)*(k+1)));
   if (dy == NULL)
     {
       return(10);
     };

   dy1=static_cast<double*>(malloc(sizeof(double)*(k+1)));
   if (dy1 == NULL)
     {
       return(10);
     };

   Fp=static_cast<double*>(malloc(sizeof(double)*(k+1)));
   if (Fp == NULL)
     {
       return(10);
     };

   /*
    *  Work out the leading dimension for the array.  Note that we may not
    *  want to use k itself, for cache issues.
    */
   if ((k % 2) == 0)
     ldam=k+1;
   else
     ldam=k;

   O=static_cast<double*>(malloc(sizeof(double)*ldam*ldam));
   if (O == NULL)
     {
       return(10);
     };

   alloc_mat(C,&Zi);
   alloc_mat(C,&dZ);
   alloc_mat(C,&dX);

   /*
    *  Fill in lots of details in the constraints data structure that haven't
    *  necessarily been done before now.
    */

   /*
    * Set up the cross links used by op_o
    * While we're at it, determine which blocks are sparse and dense.
    */

   /*
    * Next, setup issparse and NULL out all nextbyblock pointers.
    */

   for (i=1; i<=k; i++)
     {
       p=constraints[i].blocks;
       while (p != NULL)
     {
       /*
        * First, set issparse.
        */
       if (((p->numentries) > 0.25*(p->blocksize)) && ((p->numentries) > 15))
         {
           p->issparse=0;
         }
       else
         {
           p->issparse=1;
         };

       if (C.blocks[p->blocknum].blockcategory == DIAG)
         p->issparse=1;

       /*
        * Setup the cross links.
        */

       p->nextbyblock=NULL;
       p=p->next;
     };
     };

   /*
    * Now, cross link.
    */
   for (i=1; i<=k; i++)
     {
       p=constraints[i].blocks;
       while (p != NULL)
     {
       if (p->nextbyblock == NULL)
         {
           blk=p->blocknum;

           /*
        * link in the remaining blocks.
        */
           for (j=i+1; j<=k; j++)
         {
           q=constraints[j].blocks;

           while (q != NULL)
             {
               if (q->blocknum == p->blocknum)
             {
               if (p->nextbyblock == NULL)
                 {
                   p->nextbyblock=q;
                   q->nextbyblock=NULL;
                   prev=q;
                 }
               else
                 {
                   prev->nextbyblock=q;
                   q->nextbyblock=NULL;
                   prev=q;
                 };
               break;
             };
               q=q->next;
             };
         };
         };
       p=p->next;
     };
     };

   /*
    * If necessary, print out information on sparsity of blocks.
    */


   /*
    * Allocate space for byblocks pointers.
    */

   byblocks=(struct sparseblock **)malloc((C.nblocks+1)*sizeof(struct sparseblock *));
   if (byblocks == NULL)
     {
       return(10);
     };

   for (i=1; i<=C.nblocks; i++)
     byblocks[i]=NULL;

   /*
    * Fill in byblocks pointers.
    */
   for (i=1; i<=k; i++)
     {
       ptr=constraints[i].blocks;
       while (ptr != NULL)
     {
       if (byblocks[ptr->blocknum]==NULL)
         {
           byblocks[ptr->blocknum]=ptr;
         };
       ptr=ptr->next;
     };
     };

   /*
    *  Compute "fill".  This data structure tells us which elements in the
    *  block diagonal matrix have corresponding elements in one of the
    *  constraints, and which constraint this element first appears in.
    *
    */

    makefill(k,C,constraints,&fill,work1,printlevel);

    /*
     * Compute the nonzero structure of O.
     */

    nnz=structnnz(n,k,C,constraints);

    /*
     * Sort entries in diagonal blocks of constraints.
     */

    sort_entries(k,C,constraints);

	for(i=0; i<car.n_rows; i++) {

		for (j=0; j<k-1; j++) for (int j2 = j+1; j2<k; j2++)
		{
			C.blocks[1].data.mat[j2 + k*j ] = -car(i, j2, j);
			C.blocks[1].data.mat[j  + k*j2] = -car(i, j2, j);
		}

        for (j=0; j<k; j++)
		{
			C.blocks[2].data.vec[j+1] = -car(i, j, j);
		}




        initArma(n,k,C,a,constraints,&X,&y,&Z);

        ret=sdp(n,k,C,a,constant_offset,constraints,byblocks,fill,X,y,Z,
           cholxinv,cholzinv,ppobj,pdobj,work1,work2,work3,workvec1,
           workvec2,workvec3,workvec4,workvec5,workvec6,workvec7,workvec8,
           diagO,bestx,besty,bestz,Zi,O,rhs,dZ,dX,dy,dy1,Fp,
           printlevel,params);

		out(i) = 0;
		for (j = 1; j < k+1; j++)
			out(i) += y[j];

		// now compute the actual glb
		double scv = arma::accu(car.row(i));
		double svars = 0.0;
		for (j = 0; j < k; j++)
			svars += car(i, j, j);

		out(i) = (scv - svars + out(i)) / scv;

    }




   /*
    *  Now, free up all of the storage.
    */

   free_mat(work1);
   free_mat(work2);
   free_mat(work3);
   free_mat_packed(bestx);
   free_mat_packed(bestz);
   free_mat_packed(cholxinv);
   free_mat_packed(cholzinv);

   free_mat(Zi);
   free_mat(dZ);
   free_mat(dX);

   free(besty);
   free(workvec1);
   free(workvec2);
   free(workvec3);
   free(workvec4);
   free(workvec5);
   free(workvec6);
   free(workvec7);
   free(workvec8);
   free(rhs);
   free(dy);
   free(dy1);
   free(Fp);
   free(O);
   free(diagO);
   free(byblocks);
// free up the memory for X, y, Z, that was previosuly done in the parent function with free_prob
   free(y);
   free_mat(X);
   free_mat(Z);

   /*
    * Free up the fill data structure.
    */

   ptr=fill.blocks;
   while (ptr != NULL)
     {
       free(ptr->entries);
       free(ptr->iindices);
       free(ptr->jindices);
       oldptr=ptr;
       ptr=ptr->next;
       free(oldptr);
     };


  /*
   * Finally, free the constraints array.
   */

   return(ret);

}



