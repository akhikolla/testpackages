//***********************************************************************
//
//
//                              MFDPOSS.CPP
//
//
// Author(s) :Brigitte CHARNOMORDIC, Serge GUILLAUME
// Author(s) :Lydie DESPERBEN, Hazael JONES
// FISPRO Version :???
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:March 2008
// File : MFDPOSS class functions used by FISPRO, part of library fispro

//**********************************************************************
#include "fis.h"

using namespace fispro;

LIST* MFDPOSS::createList(double ls, double rs, double lk, double rk,
                          double deg, double maxp)
{
  LIST* lst = new LIST;
  lst->home();

  lst->add(ls, 0.0);
  if (!DBL_EQUAL(ls, lk) && (deg>EPSILON)&&(deg<maxp-EPSILON)) lst->add(lk, deg);

  lst->add(lk, maxp);
  if(!DBL_EQUAL(lk, rk)) lst->add(rk, maxp);

  if (!DBL_EQUAL(rs, rk) && (deg>EPSILON)&&(deg<maxp-EPSILON)) lst->add(rk, deg);

  lst->add(rs, 0.0);
  return(lst);
}

double MFDPOSS::Kernel(double & left, double & right) const
{
  return AlphaKernel(left, right, 1.0);
}


double MFDPOSS::Support(double & left, double & right) const
{
  pL->home();
  left = pL->Get()->x;

  pL->end();
  right = pL->Get()->x;

  return left + (right - left) / 2 ;
}

double MFDPOSS::AlphaKernel(double & left, double & right, double alpha) const
{
	if(DBL_INF_EQUAL(alpha, 0.)) return Support(left,right);
	if(DBL_SUP(alpha, maxposs))  return EMPTYVALUE;

	POINT* interP =NULL;

	// find the first point which is higher than input degree
	pL->home();
	while((!pL->IsEnd()) && (pL->Get()->y < alpha-EPSILON))
		pL->next();

	// find the point where the MFDPOSS crosses the line y=alpha
	if(DBL_EQUAL(pL->GetP()->x, pL->Get()->x))
		interP = new POINT(pL->Get()->x,alpha);
	else
	{
		POINT *p0 =new POINT(pL->GetP()->x,alpha);
		POINT *p1 =new POINT(pL->Get()->x,alpha);
		interP = InterSeg(pL->GetP(), pL->Get(),p0,p1);
		delete p0;
		delete p1;
	}

	if (interP == NULL) return EMPTYVALUE; // Inconsistency
	left=interP->x;
	delete interP;

	// find the last point which is bigger than input degree
	pL->end();
	while ((!pL->IsHome()) && (pL->Get()->y < alpha-EPSILON))
		pL->prev();

	// find the point where the MFDPOSS crosses the line y=alpha
	if(DBL_EQUAL(pL->GetN()->x, pL->Get()->x))
		interP = new POINT(pL->Get()->x,alpha);
	else
	{
		POINT *p0 =new POINT(pL->Get()->x,alpha);
		POINT *p1 =new POINT(pL->GetN()->x,alpha);
		interP = InterSeg(pL->Get(), pL->GetN(), p0,p1);
		delete p0;
		delete p1;
	}

	if (interP == NULL) return EMPTYVALUE; // Inconsistency
	right=interP->x;
	delete interP;

	return left + (right - left) / 2;
}

double MFDPOSS::GetDeg(double value) const
{
  MFDPOSS *f = new MFDPOSS(value);
  MFDPOSS *fi = Inter(f);
  delete f;
  if (fi==NULL) return 0.;

  double ret = fi->maxposs;
  delete fi;
  return ret;
}

int MFDPOSS::GetPoint(double &x, double &y, long index)
{
  // index is out of range
  if((index >=pL->GetSize()) || (index <0))  return -1;

  long save_index = pL->curI();
  pL->GotoI(index);
  x = pL->Get()->x;
  y= pL->Get()->y;
  pL->GotoI(save_index);
  return 0;
}


// check if the next point is an intersection point and
// return the found point if any
 POINT* MFDPOSS::CheckI(LIST *iL, LIST *sL, LIST *cL, int Ncro) const
{
  POINT* res_inter = NULL;
  POINT* interP = iL->Get();
  long  indexcro= cL->curI(); // save the pointer of the crossed curve
  bool firsteq = false;

  // loop on the segments of the curve to be crossed
  while((cL->curI() < Ncro-1) && (res_inter == NULL))
    {
      // First step:  find the intersection point
      // the final points of checked segments are the same
      if(*sL->GetN() == *cL->GetN()) res_inter = new POINT(*sL->GetN());

      // the checked segments are colinear
      else if(aligned(cL->GetN(), interP,sL->GetN()) &&
	      aligned(cL->Get(), interP,sL->GetN()))
	{
	  if(within(cL->GetN(), interP,sL->GetN())) res_inter = new POINT(*cL->GetN());
	  else if(within(sL->GetN(),interP,cL->GetN())) res_inter = new POINT(*sL->GetN());
	}
      else // the checked segments are not colinear
	  res_inter = InterSeg(interP, sL->GetN(), cL->Get(), cL->GetN());

      // Second step:  update the pointers on the segments
      if (res_inter == NULL)   cL->next();

      // intersection point is the last point of the list for the first time
      else if((*res_inter == *interP) && firsteq == false)
	{
	  firsteq = true;
	  delete res_inter;
	  res_inter = NULL;
	  cL->next();
	}
      // intersection point is not equal to the last point of the list
      // or intersection is equal to the last point for a secund time
      // ((*res_inter != *interP) || ((*res_inter == *interP) && firsteq))
      else
	{
	  // intersection point is the final point of a segment
	  if (*cL->GetN()== *res_inter) cL->next();

	  if (*sL->GetN()== *res_inter)
	    {
	      sL->next();
	      if(sL->IsEnd()) break;

	      // COMPLEX CASE
	      // check if crossed segment and next segment
	      // have infinite slopes and opposite slopes
	      if ((cL->Get()->x == sL->Get()->x) && (cL->Get()->x == cL->GetN()->x)&&
		  (cL->Get()->y > cL->GetN()->y) && (sL->Get()->y < sL->GetN()->y))
		{
		  if (within(cL->Get(), res_inter, sL->GetN()))
		    {
		      iL->add(res_inter->x, res_inter->y);
		      delete res_inter;
		      res_inter = NULL;
		      if(sL->Get()->y==cL->GetN()->y) iL->add(cL->Get()->x, cL->Get()->y);
		      else res_inter = new POINT(*cL->Get());
		      sL->next();
		    }
		  else if (within(sL->GetN(), res_inter,cL->Get()))
		    {
		      iL->add(res_inter->x, res_inter->y);
		      delete res_inter;
		      res_inter = NULL;
		      if(sL->Get()->y==cL->GetN()->y) iL->add(sL->GetN()->x, sL->Get()->y);
		      else res_inter =  new POINT(*sL->GetN());
		      sL->next();
		    }
		}
	    }
	  indexcro = cL->curI(); // update crossed curve index
	}
    } // end of while()
  cL->GotoI(indexcro);  // recover thecrossed curve  index
  return(res_inter);
}

//return NULL if the MF does not intersect
MFDPOSS* MFDPOSS::Inter(MFDPOSS *tdp) const
{
    if ((tdp==NULL) || (NbParams()<3) || (tdp->NbParams()<3)) return NULL;

  double l1,r1,l2,r2;
  Support(l1,r1);
  tdp->Support(l2,r2);

  // intersection is empty
  if (!(withinDbl(l2, l1, r1) || withinDbl(l1,l2, r2))) return NULL;


  const MFDPOSS *sel=NULL;      // pointer on the selected  curve
  const MFDPOSS *cro=NULL;      // pointer on the curve to be crossed
  LIST* lst =  new LIST;  // search list

  // case intersection includes a point or a segment with infinite slope
  if (DBL_EQUAL(r1, l2) || DBL_EQUAL(r2, l1))
    {
      if(DBL_EQUAL(r1, l2))
	{
	  sel = this;
	  cro = tdp;
	}
      else
	{
	  cro = this;
	  sel = tdp;
	}

      sel->pL->end();
      cro->pL->home();

      // intersection reduces to a point
      if (! DBL_EQUAL(sel->pL->GetP()->x, cro->pL->GetN()->x))
	{
	  delete lst;
	  return NULL;
	}

      // case intersection is a segment with infinite slope
      lst->home();
      lst->add(sel->pL->Get()->x,0);
      if (sel->pL->GetP()->y < cro->pL->GetN()->y)
	lst->add(sel->pL->GetP()->x, sel->pL->GetP()->y);
      else
	lst->add(sel->pL->GetP()->x, cro->pL->GetN()->y);
      lst->add(sel->pL->Get()->x,0);
    }

  // case intersection includes at least two segments
  else
    {
      POINT *res_inter =NULL;

      // init the list pointers
      tdp->pL->home();
      pL->home();

      // start with the MF at the right of the other
      // case start with first point of tdp
      if(DBL_SUP(l2, l1) || (DBL_EQUAL(l1, l2) && (pL->Get()->y>tdp->pL->Get()->y)))
	{
	  sel = tdp;
	  cro =this;
	}
      else if(DBL_SUP(l1, l2) || (DBL_EQUAL(l1, l2) && (pL->Get()->y<tdp->pL->Get()->y)))
	{
	  sel = this;
	  cro = tdp;
	}

      else  // the two curves start at the same point
	{
	  res_inter = new POINT(*pL->Get());
	  sel =this;
	  cro=tdp;
	}


      // initialize the searched list
      if(res_inter == NULL) // the two curves don't start at the same point
	{
	  lst->add(sel->pL->Get()->x, sel->pL->Get()->y);
	  res_inter = CheckI(lst, sel->pL, cro->pL, cro->NbParams());
	}
      //loop until the last point added is the last point of a curve
      while((pL->curI() < NbParams()-1) && (tdp->pL->curI() < tdp->NbParams()-1))
	{
	  // handle the intersection point
	  if (!(res_inter == NULL))
	    {
	      lst->add(res_inter->x, res_inter->y);
	      delete res_inter;
	      res_inter =NULL;

	      double cmp = (pL->GetN()->y-lst->Get()->y) * (tdp->pL->GetN()->x-lst->Get()->x);
	      cmp-= (pL->GetN()->x-lst->Get()->x) * (tdp->pL->GetN()->y-lst->Get()->y);

	      if (cmp < 0)
		{
		  sel = this;
		  cro = tdp;
		}
	      if(cmp > 0)
		{
		  sel = tdp;
		  cro = this;
		}
	      if (sel->pL->IsEnd()) break;
	    }
	  // handle the final point of the selected segment
	  else
	    {
	      lst->add(sel->pL->GetN()->x, sel->pL->GetN()->y);
	      sel->pL->next();
	      if (sel->pL->IsEnd()) break;
	    }
	  // update the pointer for the curve to be crossed
	  while((cro->pL->curI()<cro->NbParams()-1)
		&& (lst->Get()->x-cro->pL->GetN()->x>EPSILON))
	    cro->pL->next();

	  // check if there is an intersection in the next segment
	  res_inter =CheckI(lst, sel->pL, cro->pL, cro->NbParams());
	}
      // check the last point has been handled
      sel->pL->end();
      if(*lst->Get() != *sel->pL->Get()) lst->add(sel->pL->Get()->x, sel->pL->Get()->y);
      if(res_inter) delete res_inter;
      // Occurs when the last point is also an intersection point
      // returned by CheckI, but as the indices have incremenated
      // this point is not handled by the while loop.
    }

  MFDPOSS *mfi = new MFDPOSS(lst);
  mfi->Simplify();
  delete lst;
  return mfi;
}


//return NULL if the MF does not intersect
MFDPOSS* MFDPOSS::Union(MFDPOSS *tdp)
//***********************************
{
  // case lists of less than three points
  if ((NbParams() < 3) && (tdp->NbParams()<3)) return NULL;
  if (NbParams() < 3) return tdp->Clone();
  if (tdp->NbParams()<3) return Clone();

  double l1,r1,l2,r2;

  Support(l1,r1);
  tdp->Support(l2,r2);

  // intersection is empty
  if (!(withinDbl(l2, l1, r1) || withinDbl(l1,l2, r2))) return NULL;

  MFDPOSS *sel=NULL;      // pointer on the selected  curve
  MFDPOSS *cro=NULL;      // pointer on the curve to be crossed
  LIST* lst =  new LIST;  // search list

  // intersection includes a point or a segment with infinite slope
  if (DBL_EQUAL(r1, l2) || DBL_EQUAL(r2, l1))
    {
      if(DBL_EQUAL(r1, l2))
	{
	  sel = this;
	  cro = tdp;
	}
      else
	{
	  cro = this;
	  sel = tdp;
	}

      sel->pL->end();
      cro->pL->home();

      // intersection reduces to a point
      if (! DBL_EQUAL(sel->pL->GetP()->x, cro->pL->GetN()->x))
	{
	  delete lst;
	  return NULL;
	}

      // intersection is a segment with infinite slope
      sel->pL->home();
      lst->home();

      // include sel
      while( sel->pL->curI() < sel->NbParams()-1 )
	{
	  lst->add(sel->pL->Get()->x, sel->pL->Get()->y);
	  sel->pL->next();
	}

      cro->pL->next(); // the common point
      if (DBL_EQUAL(sel->pL->GetP()->y, cro->pL->GetN()->y)) cro->pL->next();


     while( cro->pL->curI() < cro->NbParams()-1 )
	{
	  lst->add(cro->pL->Get()->x, cro->pL->Get()->y);
	  cro->pL->next();
	}
     lst->add(cro->pL->Get()->x, cro->pL->Get()->y);
    }

// case intersection includes at least two segments
  else
    {
      POINT *res_inter = NULL;

      tdp->pL->home();
      pL->home();

      // start with the MF at the left of the other
      if (DBL_INF(l1, l2) || (DBL_EQUAL(l1, l2) && (pL->Get()->y > tdp->pL->Get()->y)))
	{
	  sel = this;
	  cro = tdp;
	}
      else if (DBL_SUP(l1, l2) || (DBL_EQUAL(l1, l2) && (pL->Get()->y < tdp->pL->Get()->y)))
	{
	  sel = tdp;
	  cro =this;
	}
      else  // the two curves start at the same point
	{
	  res_inter = new POINT(*pL->Get());
	  sel = this;
	  cro = tdp;
	}


      if(res_inter == NULL)  // the two curves don't start at the same point
	{
	  // check that the starting point has a null ordinate
	  if (sel->pL->Get()->y != 0) lst->add(sel->pL->Get()->x,0.);

	  lst->add(sel->pL->Get()->x, sel->pL->Get()->y);
	  res_inter = CheckI(lst, sel->pL, cro->pL, cro->NbParams());
	}

      //loop until the last point added is the last point of a curve
      while((pL->curI() < NbParams()-1) && (tdp->pL->curI() < tdp->NbParams()-1))
	{
	  if (res_inter != NULL)  // handle the intersection point
	    {
	      lst->add(res_inter->x, res_inter->y);
	      delete res_inter; res_inter = NULL;

	      double cmp;
	      cmp = (pL->GetN()->y-lst->Get()->y) * (tdp->pL->GetN()->x-lst->Get()->x);
	      cmp -= (pL->GetN()->x-lst->Get()->x) * (tdp->pL->GetN()->y-lst->Get()->y);

	      if (cmp < 0)
		{
		  sel = tdp;
		  cro = this;
		}
	      else if(cmp > 0)
		{
		  sel = this;
		  cro = tdp;
		}
	      if (sel->pL->IsEnd()) break;
	    }

	  else // the final point of the selected segment
	    {
	      lst->add(sel->pL->GetN()->x, sel->pL->GetN()->y);
	      sel->pL->next();
	      if (sel->pL->IsEnd()) break;
	    }
	  // update the pointer for the curve to be crossed
	  while((cro->pL->curI()<cro->NbParams()-1) && DBL_SUP(lst->Get()->x,cro->pL->GetN()->x))
	    cro->pL->next();

	  res_inter = CheckI(lst, sel->pL, cro->pL, cro->NbParams());
	}

      while(!sel->pL->IsEnd())  // check  the last points have been handled
	{
	  lst->add(sel->pL->Get()->x, sel->pL->Get()->y);
	  sel->pL->next();
    	}
      lst->add(sel->pL->Get()->x, sel->pL->Get()->y);

      if(res_inter) delete res_inter;
      // Occurs when the last point is also an intersection point
      // returned by CheckI, but as the indices have incremenated
      // this point is not handled by the while loop.
    }
  MFDPOSS *mfu = new MFDPOSS(lst);
  mfu->Simplify();
  delete lst;
  return mfu;
}

std::list<MFDPOSS> *MFDPOSS::Union(std::list<MFDPOSS> *unL)
//*********************************************************
{

  std::list<MFDPOSS> *ouL = new std::list<MFDPOSS>;

  if ((unL== NULL) || (unL->size()== 0))
    {
      ouL->push_back(*this);
      return ouL;
    }

  MFDPOSS *mfres = this;
  MFDPOSS *mfunion = NULL;
  for (std::list<MFDPOSS>::iterator it =unL->begin(); it !=unL->end(); it++)
    {
      mfunion = (*it).Union(mfres);
      if (mfunion == NULL) ouL->push_back(*it);
      else
	{
	  mfres = mfunion->Clone();
	  delete mfunion;
	}
    }
  ouL->push_back(*mfres);
  if(mfres != this) delete mfres;

  return ouL;
}


MFDPOSS* MFDPOSS::Join(MFDPOSS *tdp)
// *********************************
// the MFDPOSS has at least 3 points
// the MFDPOSS is increasing from 0 and then is decreasing to 0
// this method computes the convex enveloppe of this MF with the input MF
{
  if ((NbParams()<3) || (tdp->NbParams()<3))return NULL;

	// case lists have not the same maximum
  if(fabs(maxposs - tdp->maxposs) > EPSILON) return NULL;

  double l1,r1,l2,r2,l,r;
  Support(l1,r1);
  tdp->Support(l2,r2);

  LIST* lst =  new LIST;  // search list

  MFDPOSS *mfunion = Union(tdp);
  MFDPOSS *left = NULL;
  MFDPOSS *right = NULL;

  // case no union is possible
  if (mfunion == NULL)
	{
	  if(l1<l2)
	    {
	      left = this;
	      right = tdp;
	    }
	  else
	    {
	      right = this;
	      left = tdp;
	    }
	}
  else
    {
      left = mfunion;
      right = mfunion;
    }
  // include the  points on the left of the left MFDPOSS
  left->AlphaKernel(l,r,left->maxposs);
  POINT *lP = new POINT(l, left->maxposs);
  lst->home();
  left->pL->home();
  while ((!left->pL->IsEnd()) && (*left->pL->Get() != *lP))
    {
      lst->add(left->pL->Get()->x, left->pL->Get()->y);
      left->pL->next();
    }
  lst->add(left->pL->Get()->x, left->pL->Get()->y);
  // move the pointer of the right MFDPOSS on the kernel
  right->AlphaKernel(l,r,right->maxposs);
  POINT *rP = new POINT(r, right->maxposs);
  right->pL->end();
  while ((!right->pL->IsHome()) && (*right->pL->Get() != *rP)) right->pL->prev();

  if (*rP == *lst->Get()) right->pL->next();
  // include the  points on the right of the right MFDPOSS
  while (!right->pL->IsEnd())
    {
      lst->add(right->pL->Get()->x, right->pL->Get()->y);
      right->pL->next();
    }
  lst->add(right->pL->Get()->x, right->pL->Get()->y);

  MFDPOSS *mf_out = new MFDPOSS(lst);
  mf_out->Simplify();
  if (mfunion != NULL) delete mfunion;
  delete lst;
  delete rP;
  delete lP;
  return mf_out;
}

MFDPOSS* MFDPOSS::minTnorme(double deg)
// ************************************
{
  if(DBL_SUP_EQUAL(deg, maxposs))    return Clone();
  if(DBL_INF_EQUAL(deg, 0.))    return NULL;

  double l,r;  // abcissa of the intersection points
  if (AlphaKernel(l,r,deg) == EMPTYVALUE) return NULL; // Inconsistency

  LIST* lst =  new LIST;
  lst->home();
  pL->home();
  // include the  points on the left of the alpha kernel
  while ((!pL->IsEnd()) && (pL->Get()->y < deg-EPSILON))
    {
      lst->add(pL->Get()->x, pL->Get()->y);
      pL->next();
    }
  // include the intersection points
  lst->add(l,deg);
  if (! DBL_EQUAL(l, r))   lst->add(r,deg);

  // skip the  points in the alpha kernel
  while ((!pL->IsEnd()) && (pL->Get()->y >= deg-EPSILON)) pL->next();
  // include the  points on the right of the alpha kernel
  while (!pL->IsEnd())
    {
      lst->add(pL->Get()->x, pL->Get()->y);
      pL->next();
    }
  lst->add(pL->Get()->x, pL->Get()->y);

  MFDPOSS *res = new MFDPOSS(lst);
  delete lst;
  return res;
}

MFDPOSS* MFDPOSS::prodTnorme(double deg)
// *************************************
{
  if(DBL_SUP_EQUAL(deg, 1.)) return Clone();

  LIST* lst =  new LIST;

  pL->home();
  lst->home();
  while (!pL->IsEnd())
    {
      lst->add(pL->Get()->x, pL->Get()->y * deg);
      pL->next();
    }
  lst->add(pL->Get()->x, pL->Get()->y * deg);

  MFDPOSS *res = new MFDPOSS(lst);
  delete lst;
  return res;

}

MFDPOSS* MFDPOSS::translate(double val, double vmin, double vmax)
// **************************************************************
{

  LIST* lst =  new LIST;
  lst->home();
  pL->home();
  while (!pL->IsEnd())
    {
      lst->add(pL->Get()->x+val, pL->Get()->y);
      pL->next();
    }
  lst->add(pL->Get()->x+val, pL->Get()->y);

  MFDPOSS *mftr = new MFDPOSS(lst);
  delete lst;

  // reduce the MFDPOSS to the range[vmin,vmax]
  ACUT *mask = new ACUT(vmin, vmax, maxposs);
  MFDPOSS *mfmask = new MFDPOSS(mask);
  delete mask;

  MFDPOSS *mfinter = mftr->Inter(mfmask);
  delete mftr;

  if (mfinter == NULL)  return mfmask;
  delete mfmask;
  return mfinter;
}


void MFDPOSS::Simplify(void)
//**************************
{
  if (pL->GetSize() < 3) return;

  // if 2 points are equal, remove the first one
  pL->home();
  while(! pL->IsEnd())
    {
      if(*pL->Get() == *pL->GetN()) pL->RemD();
      pL->next();
    }

  // if 3 points are aligned, erase the middle one
  pL->home();
  pL->next();
  while(! pL->IsEnd())
    {
      if(InSegment(pL->Get(), pL->GetP(), pL->GetN())) pL->RemD();
      pL->next();
    }
}

double MFDPOSS::computeArea()
//***************************
{
  double deltax, sumy;
  double area2 =0.;

  pL->home();
  while (!pL->IsEnd())
    {
      deltax =fabs(pL->Get()->x-pL->GetN()->x);
      if(deltax > EPSILON)
    	{
	  sumy = pL->Get()->y + pL->GetN()->y;
	  area2 += deltax*sumy;
    	}
      pL->next();
    }
  return (0.5*area2);
}

void MFDPOSS::DecompAcut(int nb)
//******************************
{
  if (nb < 1) return;
  double ld, rd, alpha,maxp;
  int i;

  maxp=min(maxposs,1.0);
  acuts = new ACUT[nb];
  for(i = 1; i <= nb; i++)
    {
      alpha = maxp*double(i)/double(nb);
      AlphaKernel(ld, rd, alpha);
      acuts[i-1].l = ld;
      acuts[i-1].r = rd;
      acuts[i-1].alpha = alpha;
    }
}


