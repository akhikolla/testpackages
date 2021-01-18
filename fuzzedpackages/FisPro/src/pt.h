#ifndef __POINT_H
#define __POINT_H

namespace fispro {

//*******************   POINT and LIST  **************************
class POINT
{
  public :
    double x,y;

  POINT(){ x=0.; y=0.;};
  POINT(double xp, double yp){x=xp; y=yp;};
  //  POINT(const POINT *p) { x = p->x; y = p->y; };

  double Norme(){ return sqrt(x*x+y*y);};

  void   Print()
  {
    cout << "POINT " << x << " " << y << " " << endl;
    return;
  };

  bool operator == (const POINT &pt) const
  {
    return (DBL_EQUAL(x, pt.x) && DBL_EQUAL(y, pt.y));
  }

  bool operator!= (POINT  &pt)
  {
    return (!(*this==pt)) ;
  }


  POINT* operator= (POINT* pt)
    {
      x = pt->x;
      y = pt->y;
      return this;
    }
};

} // namespace fispro


//**********  Fonctions based on class POINT **********

// Test if points p, p1 and p2 are located on the same line
int aligned(const fispro::POINT *p, const fispro::POINT *p1,const fispro::POINT * p2);

// Test if val belongs to  intervall [b1,b2] or  [b2,b1] 
int withinDbl(double val, double b1, double b2);

// Test if point p belongs to  segment [p1,p2] or  [p2,p1] 
// knowing that p,p1, p2 are located on the same line 
int within(const fispro::POINT* p, const fispro::POINT* p1, const fispro::POINT* p2);

// Test if point p belongs to  segment [p1,p2] or  [p2,p1] 
int InSegment(const fispro::POINT* p, const fispro::POINT* p1, const fispro::POINT* p2);
 
// return the intersection point if the segments are not parallel and 
//        do intersect with the parallel flag set to 0
// return NULL if segments are parallel with 
//        the parallel flag set to 1
// return NULL if segments are not parallel and 
//        do not intersect with the parallel flag set to 0 
 
fispro::POINT* InterSeg(const fispro::POINT* pt1, const fispro::POINT* pt2, 
		const fispro::POINT* pt3, const fispro::POINT* pt4);



//*********************************************************************

class LIST
{
 private:
  struct data
  {
    fispro::POINT* pt;
    data*  next;
    data*  prev;
  };
  
  data*	head;		//Always points to first element in list
  data*	tail;		//Always points to last element in list
  data*	cur;		//Always points to cur element in list
  int	size;		//Holds the true(based on 1) size of list
  long	index;		//This is the current index
  
 public:

  LIST()
    {
      head = tail = cur = 0;
      size = 0;
      index = -1;
    }
  
  ~LIST()
    {
 	home();
  	while(!IsEmpty()) RemD();
    }
	
  void Print(FILE * f)
  {
    fprintf(f, "\n Head: %p, Tail:  %p, Cur:  %p, Size:  %d", 
	   (void *) head,  (void *) tail,  (void *) cur, size);
  }

  // Move functions	
  int next()
    {
      if(!IsEmpty())
	{
	  if(cur->next)
	    {
	      cur = cur->next;
	      index++;
	      return 1;
	    }
	}
      return 0;	
    }
  
  int prev()
    {
      if(!IsEmpty())
	if(cur->prev)
	  {
	    cur = cur->prev;
	    index--;
	    return 1;		
	  }
      return 0;
    }
 	
  int home()
    {
      if(IsEmpty()) return 0;

      cur = head;
      index = 0;
      return 1;
    }
	
  int end()
    {
      if(IsEmpty()) return 0;

      cur = tail;
      index = size-1;
      return 1;
    }
	
  void GotoI(long where)
    {
      if(where == index) return;
      
      if( where > index ) while( where>index && next() );
      else  while( where<index && prev() );
    }
	
	
  // Data access functions

  void add(double x, double y)
    {
      data* temp = new data;
      
      temp->next = 0;
      temp->prev = 0;
      temp->pt = new fispro::POINT(x, y);

      if(end())
	{
	  cur->next = temp;
	  temp->prev = cur;
	}
      else  head = temp;

      size++;
      index = size-1;
      cur = temp;
      tail = cur;

      temp = NULL;
    }
  

  fispro::POINT* GetP()
    {
      if( (!IsEmpty()) && (cur) && (cur->prev!=NULL) )
		return(cur->prev->pt);
      else return NULL;
    }
	
  fispro::POINT* GetN()
    {
      if( (!IsEmpty()) && (cur) && (cur->next!=NULL) )
		return(cur->next->pt);
      else  return NULL;
    }

  fispro::POINT* Get()
    {
      if((!IsEmpty()) && (cur)) return(cur->pt);
      else return NULL;
    }
	
  void RemD()
    {
      data* temp;
      if(!IsEmpty())
	{
	  if(IsHome())
	    {
	      head = head->next;
	      if(head) head->prev = 0;
	      delete cur->pt;
	      delete cur;
	      home();
	    }
	  else
	    {
	      temp = cur->prev;
	      temp->next = cur->next;
	      if(!IsEnd()) cur->next->prev = temp;
	      else tail = temp;
	      delete cur->pt;
	      delete cur;
	      cur = temp;	
	      index--;
	    }
	  size--;
	}
    }
	

  // Query functions
  int IsEmpty()
    {
      if(head != 0) return 0;
      return 1;
    }
  
  int IsHome()
    {
      if(!IsEmpty() && (cur == head)) return 1;
      return 0;
    }
  	
  int IsEnd()
    {
      if(!IsEmpty() && (cur == tail)) return 1;
      return 0;
    }
  	
  long curI()
    {
    return index;
    }
	
  int GetSize()
    {
      return size;
    }
};

#endif
//**************************  POINT.H  **********************************
