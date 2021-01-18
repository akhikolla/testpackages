#include "common.h"
#include "pt.h"


//**********  Fonctions utilisant la classe POINT **********

// Test if points p, p1 and p2 are located on the same line
int aligned(const fispro::POINT *p, const fispro::POINT *p1,const fispro::POINT * p2)
{
  double align = (p->x-p1->x) * (p2->y - p1->y); 
  align -= (p->y - p1->y) * (p2->x - p1->x); 
  
  if(DBL_EQUAL(align, 0.)) return 1;
  return 0; 	
}

// Test if val belongs to  intervall [b1,b2] or  [b2,b1] 
int withinDbl(double val, double b1, double b2)
{
  double bmin, bmax;
  if (b1 < b2)
    {
      bmin =b1;
      bmax =b2;
    }
  else
    {
      bmax =b1;
      bmin =b2;
    }
  if ((val>=bmin-EPSILON) && (val<=bmax+EPSILON)) return 1;
  return 0;
}  

// Test if point p belongs to  segment [p1,p2] or  [p2,p1] 
// knowing that p,p1, p2 are located on the same line 
int within(const fispro::POINT* p, const fispro::POINT* p1, const fispro::POINT* p2)
{
  // case with infinite slope
  if(DBL_EQUAL(p1->x, p2->x)) return withinDbl(p->y, p1->y, p2->y);
  else return withinDbl(p->x, p1->x, p2->x);
}

// Test if point p belongs to  segment [p1,p2] or  [p2,p1] 
int InSegment(const fispro::POINT* p, const fispro::POINT* p1, const fispro::POINT* p2)
{
  if(aligned( p, p1,p2)) return(within(p, p1, p2));
  return 0;
}
   
// return the intersection point if the segments are not parallel and 
//        do intersect with the parallel flag set to 0
// return NULL if segments are parallel with 
//        the parallel flag set to 1
// return NULL if segments are not parallel and 
//        do not intersect with the parallel flag set to 0 
 
fispro::POINT* InterSeg(const fispro::POINT* pt1, const fispro::POINT* pt2, 
                           const fispro::POINT* pt3, const fispro::POINT* pt4)
{
  fispro::POINT *tmp = NULL;
  double va,vb,vc,vd;
  
  // case lines have finite slopes
  if(!(DBL_EQUAL(pt1->x, pt2->x)) && !(DBL_EQUAL(pt3->x, pt4->x)))
    {		
      // compute taps for the lines going through pt1,pt2, pt3,pt4
      va = (pt2->y - pt1->y) / (pt2->x - pt1->x);
      vb = (pt1->y * pt2->x - pt2->y * pt1->x) / (pt2->x - pt1->x);
      vc = (pt4->y - pt3->y) / (pt4->x - pt3->x);
      vd = (pt3->y * pt4->x - pt4->y * pt3->x) / (pt4->x - pt3->x);

      if(DBL_EQUAL(va, vc)) return NULL;   // parallel lines
      else 
	tmp = new fispro::POINT((vd-vb)/(va-vc),(vd*va- vb*vc)/(va-vc));
    }
  // case first line has an infinite slope		
  else if(DBL_EQUAL(pt1->x, pt2->x) && !(DBL_EQUAL(pt3->x, pt4->x)))
    {
      vc = (pt4->y - pt3->y) / (pt4->x - pt3->x);
      vd = (pt3->y * pt4->x - pt4->y * pt3->x) / (pt4->x - pt3->x);		
      tmp = new fispro::POINT(pt1->x, vc* pt1->x + vd);
    }
  // case second line has an infinite slope	
  else if(!(DBL_EQUAL(pt1->x, pt2->x)) && DBL_EQUAL(pt3->x, pt4->x))
    {
      va = (pt2->y - pt1->y) / (pt2->x - pt1->x);
      vb = (pt1->y * pt2->x - pt2->y * pt1->x) / (pt2->x - pt1->x);
      tmp = new fispro::POINT(pt3->x, va* pt3->x + vb);
    }
  // case both lines have an infinite slope
  else return NULL;
	 	
  // check that the found intersection belongs to [p1,p2] and [p3,p4]
  if (within(tmp,pt1,pt2) && within(tmp,pt3,pt4)) return (tmp);				
  else
    {
      delete tmp;
      return NULL;	
    }
}    

