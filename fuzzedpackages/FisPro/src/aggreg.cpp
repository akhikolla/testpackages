//***********************************************************************
//
//
//                              AGGREG.CPP
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : MF class functions used by FISPRO, part of library fispro

//**********************************************************************
#include "fis.h"

void AGGREGSUM::Aggregate(RULE ** r, int nr, FISOUT *O, double deg)
  //***************************************************
{
#ifdef _OPENMP
#pragma omp single
#endif
  O->InitTabRes();
#ifdef _OPENMP
  int NbP = O->GetNbPossibles();
  vector<double> MuInfer(NbP);
  vector<int> RuleInfer(NbP);
#pragma omp for 
  for(int i = 0; i < NbP; i++)
    O->RuleInfer[i]=0;
#pragma omp for schedule(dynamic)
  for(int i = 0; i < nr; i++)
    {
      if ( ! r[i]->IsActive() )  continue;
      r[i]->Weight *= r[i]->GetExpertWeight();
#pragma omp critical 
      {MuInfer[O->ConcInfer[i]] += r[i]->Weight;}
      RuleInfer[O->ConcInfer[i]] = i; 
    }

  for(int i = 0; i < NbP; i++)
    {
      //#pragma omp critical 
      { O->MuInfer[i]+=MuInfer[i];}
      if (O->RuleInfer[i] < RuleInfer[i]) 
        O->RuleInfer[i]=RuleInfer[i];  //is this needed?
    }
#else
  int i;
  O->InitTabRes();
  for(i = 0; i < nr; i++)
    {
      if ( ! r[i]->IsActive() )  continue;

      r[i]->Weight *= r[i]->GetExpertWeight();
      O->MuInfer[O->ConcInfer[i]] += r[i]->Weight;
      O->RuleInfer[O->ConcInfer[i]] = i;
    }
#endif
}

void AGGREGMAX::Aggregate(RULE ** r, int nr, FISOUT *O, double deg)
  //***************************************************
{
  int i;
#ifdef _OPENMP
#pragma omp single
#endif
  O->InitTabRes();
#ifdef _OPENMP
  vector<double> MuInfer(O->GetNbPossibles());
  vector<int> RuleInfer(O->GetNbPossibles());
#pragma omp for schedule(dynamic)
  for(i = 0; i < nr; i++)
    {
      if ( ! r[i]->IsActive() )  continue;
      r[i]->Weight *= r[i]->GetExpertWeight();
      if(r[i]->Weight > O->MuInfer[O->ConcInfer[i]])
	{
	  O->MuInfer[O->ConcInfer[i]] = r[i]->Weight;
	  O->RuleInfer[O->ConcInfer[i]] = i;
	}
    }
  //#pragma omp critical
  for(int i = 0; i < O->GetNbPossibles(); i++)
    if (O->MuInfer[i] < MuInfer[i]) 
      {
          O->MuInfer[i]=MuInfer[i];
          O->RuleInfer[i]=RuleInfer[i];  
      }
#else
  for(i = 0; i < nr; i++)
    {
      if ( ! r[i]->IsActive() )  continue;
      r[i]->Weight *= r[i]->GetExpertWeight();
      if(r[i]->Weight > O->MuInfer[O->ConcInfer[i]])
	{
	  O->MuInfer[O->ConcInfer[i]] = r[i]->Weight;
	  O->RuleInfer[O->ConcInfer[i]] = i; 
	}
    }
#endif
}

void AGGREGIMP::Aggregate(RULE ** r, int nr, FISOUT *O, double deg)
  //*****************************************************************
{
  int i;
  bool empty;
  O->InitTabRes();

  MFDPOSS * mfg = NULL;
  MFDPOSS * mfgold = NULL;
  MFDPOSS * tmpdposs = NULL;
  MF * tmp=NULL;
  double *pPoss = O->GetPossibles();

  i=0;
  empty=false;//use to test empty intersection
  while (i<nr)
     {
       if(r[i]->IsActive())
	 {
	   if(r[i]->Weight > EPSILON)
	     {
	       if (O->MfConc[i]!= NULL) delete O->MfConc[i];
	       O->MfConc[i] = NULL;
	       // compute inferred output for each active rule
	       tmpdposs = imp->ComputeDposs(O->GetMF((int)(pPoss[O->ConcInfer[i]]-1)), r[i]->Weight);
	       O->MfConc[i]= new MFDPOSS(*tmpdposs);
	       delete tmpdposs;
	       tmpdposs=NULL;
	       // compute the intersection of all inferred outputs
	       // first time: mfg is MfConc[i] (rule conclusion with rule weight ponderation)
	       if (mfg == NULL && !empty) mfg = new MFDPOSS(*O->MfConc[i]);
	       // next times: mfg is the intersection of MfConc[i] and of previous result
	       // Inter returns NULL if intersection is empty
	       else if (mfg != NULL)
		 {
		   mfgold = new MFDPOSS(*mfg);
		   delete mfg;
		   mfg = O->MfConc[i]->Inter(mfgold);
		   delete mfgold;
		 }
	       if (mfg==NULL) empty=true;//case of empty intersection
	       //added by bch January 23rd - Tnorm on rule conclusions
	       if (O->MfConc[i]!= NULL)
		 if(deg < 1-EPSILON)
		   {
		     tmpdposs=imp->ComputeTnorme(O->MfConc[i],deg);
		     delete O->MfConc[i];
		     O->MfConc[i]=NULL;
		     O->MfConc[i]=new MFDPOSS(*tmpdposs);
		     delete tmpdposs;
		     tmpdposs=NULL;
		   }
	       //printf("\n in aggregate rule=%d, rule concl. result:\n",i);
	       //O->MfConc[i]->Print();
	     }
	   // case rule is active but low weight<epsilon
	   else
	     {
	       //assign universal conclusion to rule for visual satisfaction (interface)
	       if (O->MfConc[i]!= NULL) delete O->MfConc[i];
	       O->MfConc[i] = NULL;
	       tmp=new MFDOOR(O->min(),O->max());
	       O->MfConc[i]=new MFDPOSS(tmp,0);
	       delete tmp;
	     }
	 }
       i++;
     }
  // end of while loop on rules
  // if intersection was empty mfg is NULL, so MFGlob will be NULL also
  if (O->MfGlob!= NULL) delete O->MfGlob;
  O->MfGlob = NULL;

  //if (mfg == NULL) fprintf(stdout, "\nINFERED MFDPOSS IS EMPTY\n");
  //else  // truncate the output result
  //  {
      if(deg < 1-EPSILON)
	{
	  O->MfGlob = imp->ComputeTnorme(mfg, deg);
	  delete mfg;
	}
      else O->MfGlob = mfg;
  //  }
}

MFDPOSS * IMPLIGD::ComputeDposs(MF * A, double degre)
  //**************************************************
  //for impli = godel
{
  if(degre < EPSILON) return NULL;

  // first value of dilatation of the kernel
  double lval = -1.0;
  // second value for the dilatation of the kernel
  double rval = -1.0;
  // MF params
  double param[10];

  MFDPOSS * tmpdposs = NULL;
  MFTRAP *tmp = NULL;

  // the MFDPOSS is equivalent to the original MF
  if(degre > 1.0-EPSILON)  
    {
      tmpdposs = new MFDPOSS(A,0);
      return tmpdposs;
    }

  //other cases: triangle, trapezoid, semi trapezoid, door, universal
   A->GetParams(param);
  if(strcmp(A->GetType(),"trapezoidal") == 0)
    {
      lval = param[0]*(1-degre) + param[1]*degre;
      rval = param[3]*(1-degre) + param[2]*degre;
      tmp = new MFTRAP(param[0], lval, rval, param[3]);
      tmpdposs = new MFDPOSS(tmp,degre);
      delete tmp;
      return tmpdposs;
    }
  if(strcmp(A->GetType(),"triangular") == 0)
    {
      lval = param[0]*(1-degre) + param[1]*degre;
      rval = param[2]*(1-degre) + param[1]*degre;
      tmp = new MFTRAP(param[0], lval, rval, param[2]);
      tmpdposs = new MFDPOSS(tmp, degre);
      delete tmp;
      return tmpdposs;
    }
  if(strcmp(A->GetType(),"SemiTrapezoidalInf") == 0)
    {
      lval = param[0];
      rval = param[2]*(1-degre) + param[1]*degre;
      tmp = new MFTRAP(param[0], lval, rval, param[2]);
      tmpdposs = new MFDPOSS(tmp, degre);
      delete tmp;
      return tmpdposs;
    }
  if(strcmp(A->GetType(),"SemiTrapezoidalSup") == 0)
    {
      lval = param[0]*(1-degre) + param[1]*degre;
      rval = param[2];
      tmp = new MFTRAP(param[0], lval, rval, param[2]);
      tmpdposs = new MFDPOSS(tmp, degre);
      delete tmp;
      return tmpdposs;
    }
  if(strcmp(A->GetType(),"universal") == 0)
    {
      tmpdposs = new MFDPOSS(A,degre);
      return tmpdposs;
    }
  if(strcmp(A->GetType(),"door") == 0)
    {
      lval = param[0];
      rval = param[1];
      tmp = new MFTRAP(param[0], lval, rval, param[1]);
      tmpdposs = new MFDPOSS(tmp,degre);
      delete tmp;
      return tmpdposs;
    }
  // other cases: non allowed MFs
  sprintf( ErrorMsg, "~OnlyTriangularOrTrapezoidalShapesOrDoorsOrUniversalMFsAreManaged%s",
	   "~InOutputPartitionsWithImplicativeRules");

  throw std::runtime_error( ErrorMsg );
    
}

MFDPOSS* IMPLIRG::ComputeDposs(MF *A, double degre)
//*************************************************
{

  if(degre < EPSILON) return NULL;

  // first value of dilatation of the kernel
  double lval = -1.0;
  // second value for the dilatation of the kernel
  double rval = -1.0;
  // MF params
  double param[10];

  MFDPOSS * tmpdposs = NULL;
  MFDOOR *tmp = NULL;

  // correction bch 14 May 09 - get param(MF) in all cases
  A->GetParams(param);

  if(strcmp(A->GetType(),"trapezoidal") == 0)
    {
      if(degre > 1.0-EPSILON)
	{
	  //the MFDPOSS  is the kernel of the original MF
	  tmp = new MFDOOR(param[1], param[2]);
	  tmpdposs = new MFDPOSS(tmp,0);
	}
      else
	{
	  //A->GetParams(param);
	  lval = param[0]*(1-degre) + param[1]*degre;
	  rval = param[3]*(1-degre) + param[2]*degre;
	  tmp = new MFDOOR(lval,rval);
	  tmpdposs = new MFDPOSS(tmp,0);
	}
      delete tmp;
      return tmpdposs;
    }

  if(strcmp(A->GetType(),"triangular") == 0)
    {
      if(degre > 1.0-EPSILON)
	{
	  //the MFDPOSS is the kernel of the original MF
	  tmpdposs = new MFDPOSS(param[1]);
	  return tmpdposs;
				 
	}
      else
	{
	  //A->GetParams(param);
	  lval = param[0]*(1-degre) + param[1]*degre;
	  rval = param[2]*(1-degre) + param[1]*degre;
	  tmp = new MFDOOR(lval,rval);
	  tmpdposs = new MFDPOSS(tmp,0);
	  delete tmp;
	  return tmpdposs;
	}
     
    }

  if(strcmp(A->GetType(),"SemiTrapezoidalInf") == 0)
    {
      if(degre > 1.0-EPSILON)
	{
	  //the MFDPOSS is the kernel of the original MF
	  tmp = new MFDOOR(param[0],param[1]);
	  tmpdposs = new MFDPOSS(tmp,0);
	}
      else
	{
	  lval = param[0];
	  rval = param[2]*(1-degre) + param[1]*degre;
	  tmp = new MFDOOR(lval,rval);
	  tmpdposs = new MFDPOSS(tmp,0);
	}
      delete tmp;
      return tmpdposs;
    }

  if(strcmp(A->GetType(),"SemiTrapezoidalSup") == 0)
    {
      if(degre > 1.0-EPSILON)
	{
	  //the MFDPOSS is the kernel of the original MF
	  tmp = new MFDOOR(param[1],param[2]);
	  tmpdposs = new MFDPOSS(tmp,0);
	}
      else
	{
	  lval = param[0]*(1-degre) + param[1]*degre;
	  rval = param[2];
	  tmp = new MFDOOR(lval,rval);
	  tmpdposs = new MFDPOSS(tmp,0);
	}
      delete tmp;
      return tmpdposs;
    }

  if(strcmp(A->GetType(),"universal") == 0)
    {
      tmpdposs = new MFDPOSS(A,0);
      return tmpdposs;
    }

  if(strcmp(A->GetType(),"door") == 0)
    {
      tmpdposs = new MFDPOSS(A,0);
      return tmpdposs;
    }
  
  sprintf( ErrorMsg, "~OnlyTriangularOrTrapezoidalShapesOrDoorsOrUniversalMFsAreManaged%s",
	   "~InOutputPartitionsWithImplicativeRules");
	   throw std::runtime_error( ErrorMsg );
}

MFDPOSS* IMPLIGG::ComputeDposs(MF * A, double degre)
  //***************************************************
  //for impli = goguen
{
  if(degre < EPSILON) 	return NULL;

  // left value for dilatation of the kernel
  double lval = -1.0;
  // right value for dilatation of the kernel
  double rval = -1.0;
  //value of the MF
  double param[10];

  MFDPOSS *tmpdposs = NULL;
  MFTRAP *tmp = NULL;

  // if max degree the MFDPOSS is equivalent to the original MF
  if(degre > 1.0-EPSILON) 
    {
      tmpdposs = new MFDPOSS(A,0);
      return tmpdposs;
    }
  // allowed cases: triangle, trapezoid, semi trapezoid, door, universal
  A->GetParams(param);

  if(strcmp(A->GetType(),"trapezoidal") == 0)
    {
      lval = param[0]*(1-degre) + param[1]*degre;
      rval = param[3]*(1-degre) + param[2]*degre;
      tmp = new MFTRAP(param[0],lval,rval,param[3]);
      tmpdposs = new MFDPOSS(tmp,0);
      delete tmp;
      return tmpdposs;
    }
          
  if(strcmp(A->GetType(),"triangular") == 0)
    {
      A->GetParams(param);

      lval = param[0]*(1-degre) + param[1]*degre;
      rval = param[2]*(1-degre) + param[1]*degre;

      tmp = new MFTRAP(param[0], lval, rval, param[2]);
      tmpdposs = new MFDPOSS(tmp,0);
      delete tmp;
      return tmpdposs;
    }
if(strcmp(A->GetType(),"SemiTrapezoidalInf") == 0)
    {
      lval = param[0];
      rval = param[2]*(1-degre) + param[1]*degre;
      tmp = new MFTRAP(param[0], lval, rval, param[2]);
      tmpdposs = new MFDPOSS(tmp, 0);
      delete tmp;
      return tmpdposs;
    }
  if(strcmp(A->GetType(),"SemiTrapezoidalSup") == 0)
    {
      lval = param[0]*(1-degre) + param[1]*degre;
      rval = param[2];
      tmp = new MFTRAP(param[0], lval, rval, param[2]);
      tmpdposs = new MFDPOSS(tmp, 0);
      delete tmp;
      return tmpdposs;
    }
  if(strcmp(A->GetType(),"universal") == 0)
    {
      tmpdposs = new MFDPOSS(A,0);
      return tmpdposs;
    }
  if(strcmp(A->GetType(),"door") == 0)
    {
      tmpdposs = new MFDPOSS(A,0);
      delete tmp;
      return tmpdposs;
    }
  

  // other input MFs not allowed
  sprintf( ErrorMsg, "~OnlyTriangularOrTrapezoidalShapesOrDoorsOrUniversalMFsAreManaged%s",
	   "~InOutputPartitionsWithImplicativeRules");
  throw std::runtime_error( ErrorMsg );

}

MFDPOSS* IMPLI::ComputeTnorme(MFDPOSS * mf, double deg)
//*****************************************************
// default Tnorme is min
// used by resher-gaines and godel's impli.
{
  if (mf == NULL) return NULL;
  return (mf->minTnorme(deg));
}

MFDPOSS* IMPLIGG::ComputeTnorme(MFDPOSS * mf, double deg)
//*******************************************************
// prodTnorme is used by goguen impli
{
  if (mf == NULL) return NULL;
  return (mf->prodTnorme(deg));
}

