//***********************************************************************
//
//
//                              FIS.H
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  December 1, 2007
// File : FIS class functions used by FISPRO, part of library fispro


//**********************************************************************


//******************************************************************
//
//                  DEFINITION OF THE  NECESSARY CLASSES
//
//                  FOR THE FUZZY INFERENCE SYSTEM
//
//
//******************************************************************


#ifndef __FIS_H
#define __FIS_H


#include "common.h"


#ifdef _OPENMP
#include <omp.h>
#endif

//! Maximum number of MF for a given input
#define MAX_MF 999

//! Maximum number of classes a given output
#define MAX_CLASS 15

#define OUT_DEFAUT     -1.0

//! Maximum number of breakpoints
#define MAX_BP         10000

//! Performance result file header labels
#define OBSERVED  "OBS"
#define INFERRED  "INF"
#define ALARM     "Al"
#define CLASINF   "ClINF"
#define CLASALARM "CLAl"
#define FIS_ERROR     "Err"
#define BLANK     "Bl"
#define CUMUL     "CErr2"
#define INFERREDSYMBMATCH "MATCH"
#define MINKERNEL "MINK"
#define MAXKERNEL "MAXK"
#define MINSUPPORT "MINS"
#define MAXSUPPORT "MAXS"
#define MUTHIMPLI 0.667




//! Merge two rules base coming from different FIS into one.
//! Fis1 and Fis2 are the two Fis file name.
//! Merge is the merge file name.
//! Occur is the file name containing the rule occurences.
//! tabconc is an array containing the rule conclusions,
//! if distinct conclusions are not taken into account in the rule comparison
//! conc is a flag, true if distinct conclusions yield distinct rules, false otherwise.
//! In this case, the mean and the standard deviation of conclusions for a given
//! premise are computed.int MergeRules(const char * Fis1, const char *Fis2, const char * Merge, const char *Occur, double **&tabconc, int conc);

//! Contains the MergeRules function, make the comparison between multiple fis
//! firstpart is the firstpart of the FIS file name
//! n is the number of fis to compare
//! lastfispart is the eventual last part of the fis file name
//! Resulting Fis file name: "firstfispart"."number" (if lastfispart is NULL)
//! or "firstfispart"."number"."lastfispart"
//! nR, meanO and stdO are final results (respectively the total number of
//! distinct rules, the mean number of rule occurrence in the rule base,
//! and the standard deviation of occurrences.
//! conc is a flag, true if distinct conclusions yield distinct rules, false otherwise.
//! In this case, the mean and the standard deviation of conclusions for a given
//! premise are computed.
int StableRules(char * firstfispart, int n, char * lastfispart, char * res, int & nR, double & meanO, double & stdO, int conc);

// int ReadResult(char * in, int nrow, char sep, char *out = NULL);


//! For Kernel and Support, left and right
struct Trapeze
{
  double lk, rk, ls, rs;
};


// For segment left and right
struct CSeg
{
  double l, r;
};

//! Rule base information structure
//! Includes Print and WriteHeader functions
//! maxR   : maximal number of rules according to the active inputs
//!          and the FIS structure
//! nR     : number of rules
//! maxVr  : maximum number of variables in a rule
//! meanVr : mean number of variables per rule
//! nVar   : number of distinct variables used in rules
//! meanMF : mean number of MF per used variable
//! nClass : number of classes or MF in the fuzzy output
//!          (0 if crisp output and not classif)
//! nRc    : number of rules per class
//! nIn    : number of inputs
//! nOut   : number of outputs
//! out    : the output number
//! The number of MF in negative when si input is inactive
struct InfoRB
{
private:
  InfoRB& operator=(const InfoRB&);
  InfoRB(const InfoRB&);

public:
  int maxR, nR, maxVr, nVar, nClass, *nRc, *nMf;
  double meanVr, meanMF, *labels;
  int nIn, nOut, out;

  InfoRB()
  {
    maxR = nR = maxVr = nVar =  -1;
    meanVr = meanMF = -1.;
    nClass = 0;
    nRc = NULL;
    nMf = NULL,
    labels = NULL;
    nIn = 0;
    nOut = 0;
    out = 0;
  }

  ~InfoRB()
  {
    if(nRc) delete [] nRc;
    if(nMf) delete [] nMf;
    if(labels) delete [] labels;
  }

  void Print(FILE * f) const
  {
    int i;

    for(i = 0; i < nIn; i++)
      fprintf(f, "%d & ", nMf[i]);
    for(i = 0; i < nOut; i++)
      fprintf(f, "%d & ", nMf[i+nIn]);

    fprintf(f, "%d & %d & %d & %d & %f & %d & %f ", out+1, maxR, nR, maxVr, meanVr, nVar, meanMF);
    if(nClass && nRc && labels)
      for(int i = 0; i < nClass; i++)
	fprintf(f, "& (%f) & %d ", labels[i],  nRc[i]);
  }

  void WriteHeader(FILE * f) const
  {
    int i;

    for(i = 0; i < nIn; i++)
      fprintf(f, "In %d & ", i+1);
    for(i = 0; i < nOut; i++)
      fprintf(f, "Out %d & ",i+1);
    fprintf(f, " Out  &   maxR  &   nR  &   maxVr &   meanVr &  nVar &  meanMF ");
    if(nClass && nRc && labels)
      for(int i = 0; i < nClass; i++)
	fprintf(f, "& (class/MF)  &  nRc  ");
  }
};


struct ACUT
{
  double l, r, alpha;

  ACUT()
  {
    l=r=alpha=0.0;
  }
  ACUT(double a,double b,double degre)
  {
    l=a;
    r=b;
    alpha=degre;
  }
};



//*******************************************************************
//
//
//                         CLASS for membership function
//
//
//*******************************************************************


//! An abstract class to handle various types of membership functions:
//! 'triangular', 'trapezoidal', 'SemiTrapezoidalInf', 'SemiTrapezoidalSup',
//! 'gbell', 'gaussian', 'discrete', 'universal', 'door'.

class MFPWLinear;

class MF
//*******
{
 MF& operator=(const MF&);

 public:
  char *Name;
  ACUT *acuts;

  MF()
    {
      Name = new char [1];
      Name[0] = 0;
      acuts = NULL;
    }


  MF( const MF & sef ) { Name = NULL; SetName( sef.Name ); acuts = NULL; }


  virtual ~MF() { delete [] Name; delete [] acuts;}


  int operator != (const MF &);
  virtual int NbParams() const { return 0; }
  virtual void GetParams( double *params ) const {}
  virtual void Update( double *params) =0;//!!!pb
  //! Returns the type of the MF
  virtual const char *GetType() const = 0;
  virtual MF* Clone() const = 0;
  //! Normalize a MF, param's are entry bound
  virtual void Normalize(double,double){};
  //! Unnormalize a MF, param's are entry bound before normalisation.
  virtual void UnNormalize(double,double){};
  void SetName(const char *name );


  //! The main function of the class
  //! Compute the membership degree of a given value to the MF
  //! Has to be describe for each type of MF
  virtual double GetDeg(double value) const = 0;

  virtual double GetLeftx(double degre) {return EMPTYVALUE; };

  virtual double GetRightx(double degre) {return EMPTYVALUE; };

  //! Computes the degree of similarity of the MF T with the current one
  double MFMatchDeg(MF * T) const;

  //! Return:  middle of the kernel
  double Kernel(CSeg &c) const { return Kernel(c.l, c.r); }
  virtual double Kernel(double & left, double & right) const = 0;

  //! Used by FIS::Crisp2Fuz to turn a crisp output into a fuzzy one.
  //! It is defined as Kernel() except for SemiTrapezoidal shapes.
  virtual double Corner() const { double left, right; return Kernel(left, right); }

  //! Return:  middle of the support
  void Support(CSeg &c) { Support(c.l, c.r); }
  virtual double Support(double & left, double & right) const = 0;

  //! Compute the alpha cut kernel (nivel alpha)
  //! Return:  middle of the segment
  double AlphaKernel(CSeg &c, double al) const { return AlphaKernel(c.l, c.r, al); }
  virtual double AlphaKernel(double & left, double & right, double alpha) const = 0;


  virtual void Print(FILE *f) const
    { fprintf(f, "\nMF : %s\tType : %s\t", Name, GetType()); }


  //! Print membership function  configuration.
  //! Function used to store a FIS configuration
  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
    {
      fprintf(f, "MF%d=%c%s%c%c%c%s%c%c", num+1, STRING_SEP, Name, STRING_SEP,
	      SEPARE, STRING_SEP, GetType(), STRING_SEP, SEPARE);
    }

  void PrintAlphaCuts(int nbalf, FILE *f)
  {
    for(int i=0; i<nbalf;i++)
      fprintf(f, "AlphaCut %d %8.3f, %8.3f, %8.3f\n",
	     i, acuts[i].l,acuts[i].r, acuts[i].alpha);
  }


  //! Calculates the centroid in x, centre,  and the area, mass,
  //! of the specified alpha cut.
  //! The alpha-cut is modelled as a trapeze defined by
  //! the MF support and the alpha-cut kernel.
  //! The correponding coordinates are stored in 'coord'.
  void Centroid(double a_coupe, double & centre, double & masse, Trapeze * coord) const;
   
  //Decompose the input in n alpha-cuts
  //Fill the acuts array
  void DecompAcut(int n)
  {
    double ld, rd;
    acuts = new ACUT[n];
    for(int i = 1; i <= n;i++)
      {
	AlphaKernel(ld,rd,double(i-1)/double(n));
	acuts[i-1].l = ld;
	acuts[i-1].r = rd;
	acuts[i-1].alpha = double(i)/double(n);
      }
  }

  virtual MFPWLinear getMFPWLinear() const;
}; // end class MF



#include "mf.h"  // Inherited classes

  //! Return a fuzzy number, either MFTRI or MFTRAP
  //! v is the center 
  //! kw is the kernel width
  //! sw is the support width
MF * FuzNumber(double v, double kw, double sw);

struct Vec: public std::vector<double>
{
  double& operator[](size_t i) {
    assert(i<size());
    return std::vector<double>::operator[](i);
  }
  double operator[](size_t i) const {
    assert(i<size());
    return std::vector<double>::operator[](i);
  }
};

namespace {
  bool maxSupLess(const MF* x, const MF* y)
  {
    double xMax, yMax, dummy;
    x->Support(dummy,xMax);
    y->Support(dummy,yMax);
    return xMax<yMax;
  }
}


//******************************************************************
//
//
//                             CLASS FISIN
//
//
//******************************************************************


//! Base class to handle input and output variables.
class FISIN
//**********
{
  //void ldLinMFs();//moved into public for OPENMP


   protected :


    //! Range bounds
    double ValInf, ValSup;

  //! The number of MF of the input
  int Nmf;
  //! The MF's description
  MF ** Fp;
  //! The input is taken into account (active is true) or not
  int active;

  //! Read MF parameters from a buffer
  void ReadMf(char * ligne, int NumSef); //throw(runtime_error);

  //! Read input parameters from a file
  void Init( ifstream &f, int bufsize, int num ); //throw(runtime_error);

  //! Set all values to default and pointers to NULL
  void Init();
  /* For GetDegs, we store a list of PW linear MFs, and then a list of remainders */
  std::vector<std::pair<int,MFPWLinear> > linMF;
  std::vector<std::pair<int, MF*> > nonLinMF;

#ifdef _OPENMP
  // for OpenMP use, we need a separate Mfdeg vector for each thread
  std::vector<Vec> _Mfdeg;
#else
  //  std::vector<double> _Mfdeg;
  Vec _Mfdeg;
#endif

  void destroy()
    { 
      delete [] Name;
      if(Nmf>0)
	{
	  if(Fp !=NULL)
	    {
	      for(int i = 0 ; i < Nmf ; i++) delete Fp[i];
	      delete [] Fp; 
	      Fp=NULL;
	    }
	  //delete [] Mfdeg;
	  //Mfdeg=NULL;
	}
      if(dPart) delete [] dPart;
      dPart=NULL;
    }

 public:

  void ldLinMFs();// for OPENMP

   FISIN& operator=(const FISIN& x) {
    ValInf=x.ValInf;
    ValSup=x.ValSup;
    Nmf=x.Nmf;
    active=x.active;
    MinObserve=x.MinObserve;
    MaxObserve=x.MaxObserve;
    OLower=x.OLower;
    OUpper=x.OUpper;
    // now handle heap objects
    destroy();
    
    Fp=new MF*[Nmf];
    for (int i=0; i<Nmf; i++)
      Fp[i]=x.Fp[i]->Clone();

    Name=new char[strlen(x.Name)+1];
    strcpy(Name,x.Name);
    return *this;
  }
  //array of MFDOOR for partioning decomposition
  MFDOOR *dPart;
  int nPart;

  char *Name;
  //! Observed Range Bounds in the learning data set
  double MinObserve, MaxObserve;
  //! This array contain the last computed values
  //double * Mfdeg;
  //  std::vector<double>& Mfdeg() {
  Vec& Mfdeg() {
#ifdef _OPENMP
    if (privMfdeg)
      return _Mfdeg[omp_get_thread_num()];
    else
      return _Mfdeg[0];
#else
    return _Mfdeg;
#endif
  }
  const Vec& Mfdeg() const {return const_cast<FISIN*>(this)->Mfdeg();}

  //! Indicate whether Mfdeg should be thread private, or shared
  bool privMfdeg;


  //! Save Range bounds for normalisation unormalisation
  double OLower, OUpper;// OLower > OUpper to be sure that it has been initialized

  //for interface of fuzzy input
  double Kw; //kernelWeight
  double Sw; //supportWeight

  FISIN( const FISIN & );


  FISIN(ifstream &f, int bufsize, int num) : privMfdeg(false) //throw(runtime_error)
    {
      Init();
      Init( f, bufsize, num );
    }


  FISIN(int n = 0, int flagactive = true) : privMfdeg(false)
    {
      Init();
      Nmf = n;
      active = flagactive;
      if( Nmf != 0 )
        {
	  Fp = new MF *[Nmf];
	  //Mfdeg = new double [Nmf];
        }
    }


  //! Generates a standardized fuzzy partition:
  //! regular grid with n  triangle MFs between  min and max.
  //! If 'tri' is false the extreme MF's are trapezoidal shaped,
  //! otherwise they are  triangles with one infinite breakpoint
  FISIN(int n, double min, double max, int tri = false); //throw(runtime_error);


  //! Generates a standardized fuzzy partition:
  //! regular grid with n  triangle MFs between  min and max.
  //! MF centers are in the 't' array.
  //! If sort is true, the center array is first sorted in a ascending way.
  FISIN(double *t, int n, double min, double max, int sort = true); //throw(runtime_error);

  //! Generates a standardized fuzzy partition with trapezoidal MF.
  //! Range limits: min and max. MF kernel limits are in the 't' array.
  //! n is the 't' array size. NbMF=(n/2)+1. 
  FISIN(int n, double *t, double min, double max); //throw(runtime_error);
  
  // QQ
  FISIN(double *, int *, int, double , double, double, double,int );


  /* virtual ~FISIN()
    {
      int i ;
      delete [] Name;
      if(Nmf>0)
	{
	  if(Fp !=NULL)
	    {
	      for(i = 0 ; i < Nmf ; i++) delete Fp[i];
	      delete [] Fp; 
	      Fp=NULL;
	    }
	   delete [] Mfdeg;
	   Mfdeg=NULL;
	}
      if(dPart) delete [] dPart;
      dPart=NULL;
    }
  */
  virtual ~FISIN() {destroy();}


  int operator != (const FISIN &);
  void SetName(const char *name );
  void SetStdMfNames(void);

  void SetRange( double range_inf, double range_sup ); //throw(runtime_error);
  void SetRangeOnly( double range_inf, double range_sup ); //throw(runtime_error);
  MF *GetMF( int sef_number )
  {
	  return Fp[sef_number];
  }
  virtual const char *GetType() const { return "Input"; };
  //void AddMF( MF *sef );
  void AddMF( MF *sef) {AddMF(sef,Nmf);}
  void AddMF( MF *sef, int insert_before ); ///< add MF before \a insert_before    
  void RemoveMF( int sef_number );
  void ReplaceMF( int sef_number, MF *new_sef );
  void MoveMF( int sef_number, int move );
  void SortMFs() {std::sort(Fp,Fp+Nmf,maxSupLess);}

  // initialize normalisation flags
  void initNormalize();
  void Normalize();
  void UnNormalize();


    //! Returns  1 if the whole vector is equal to zero otherwise 0
  //! and fills the array Mfdeg
  int GetDegs(double v);
  void GetDegsV(double v);  ///<same as GetDegs, but don't calc return flag
  
  //! Sets all degrees as equal with the sum being 0.5
  int SetEqDegs(double v);

  //! Gets a random value (FisRand) and then computes its membership degrees
  int GetRandDegs(double v);

  //! Returns the membership degree of 'v' to the MF #'s'
  double GetADeg(int s, double v);

  //! Computes the degree of similarity of the MF T with the MF's of the partition
  //! and fills the array Mfdeg
  double MFMatchDegs(MF * T);

  //! Computes the degree of similarity of the MF T with the MF #'s'
  double MFMatchADeg(int s, MF * T);

  //! Calculates the centroid, Kernel and  Support  of a MF
  //! Calculates the centroid in x, c,  and the area, m ,
  //! of the MF '#s'  alpha cut, 'level'.
  //! The alpha-cut is modelled as a trapeze defined by the MF support and the alpha-cut kernel.
  //! The correponding coordinates are stored in 't'.
  void Centroid(int s, double level, double & c, double & m, Trapeze * t)
    {
      if( (s >= 0) && (s < Nmf) )
	Fp[s]->Centroid(level, c, m, t);
      else
	m = 0.;
    }

  //! Return:  middle of the kernel of MF #'s'
  double Kernel(int s, double &deb, double &fi)
    {
      if( (s >= 0) && (s < Nmf) )
	return Fp[s]->Kernel(deb, fi);
      else return FisMknan();
    }


  //! Return:  middle of the support of MF #'s'
  double Support(int s, double &deb, double &fi)
    {
      if( (s >= 0) && (s < Nmf) )
	return Fp[s]->Support(deb, fi);
      else return FisMknan();
    }

  //! Fills the 'p' array with the Nmf centers.
  //! The pointer 'p' has to be previoulsly allocated.
  void GetMfCenters(double *p);
  //! Fills the 'p' array with the PFF characteristic points (trapezoidal PFF)
  void GetSFPparams(double *&p, int* &typemf, int & size, FILE *display);
 
  //! This function filled the array 'Bp' with the coordinates of the break points
  //! of the partition. 'np' is the number of breakpoints.
  //! The allocation of space memory for storing 'Bp' is done within the function.
  //! A break point delimites the region of influence for a rule (or a combination of rules).
  //! For a n-MF partition, at most n+n-1 break points are defined:  the middle of the kernel
  //! of each MF and the intersection, if exists, for two adjacent MFs, of the lines which
  //! join the support to the kernel limits. Considered neighbours for ith MF are MF i-1 and
  //! MF i+1:  there is no MF sorting done by this function.
  void GetBreakPoints(double * & Bp, int & np);

  //! Return current value of active
  int IsActive(void) { return active; }
  //! Activate the input :   active = true
  void Activate(void) { active = true; }
  //! Deactivate the input :  active = false
  void Deactivate(void) { active = false; }


  //! Return the number of MF of the partition
  int GetNbMf(void) { return Nmf; }


  //! Return lower bound of range
  double min(void) { return ValInf; }
  //! Return upper bound of range
  double max(void) { return ValSup; }


  // Different printing functions
  void PrintDeg(FILE * f) const
    {
      int i;
      fprintf(f, "MF degrees for input : %s\n", Name);
      for(i = 0; i < Nmf; i++)  fprintf(f, "\t%8.3f", Mfdeg()[i]);
      fprintf(f, "\n");
    }

  virtual void Print(FILE * f) const
    {
      int i;

      fprintf(f, "\n%s : %s   Active (oui = 1) : %d", GetType(), Name, active);
      fprintf(f, "\nRange : %8.3f%c%8.3f", ValInf, SEPARE, ValSup);
      fprintf(f, "\nNmf : %d", Nmf);
      for(i = 0; i < Nmf; i++)    Fp[i]->Print(f);
      if(!strcmp(GetType(), "Input")) fprintf(f, "\n");
    }



  //! For output type handling the tag must be isolated
  //! printing the tag
  //! Function used to store a FIS configuration
  void PrintCfgTag(int num, FILE * f) const
    {
      fprintf(f, "\n[%s%d]\n", GetType(), num);    // Tag
    }


  //! Print common part of input and output configuration.
  //! Function used to store a FIS configuration
  virtual void PrintCfgCont(FILE * f, const char *fd = FORMAT_DOUBLE) const
    {
      int i;
      char StrActive[4];

      if(active) sprintf(StrActive, "yes");
      else sprintf(StrActive, "no");

      fprintf(f, "Active=%c%s%c\n", STRING_SEP, StrActive, STRING_SEP);
      fprintf(f, "Name=%c%s%c\n", STRING_SEP, Name, STRING_SEP);
      fprintf(f, "Range=%c", START_NB);
      fprintf(f, fd, ValInf);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, ValSup);
      fprintf(f, "%c\n", END_NB);
      fprintf(f, "NMFs=%d\n", Nmf);
      for(i = 0; i < Nmf; i++)    Fp[i]->PrintCfg(i, f, fd);
    }


  //! Print input configuration.
  //! Function used to store a FIS configuration
  virtual void PrintCfg(int num, FILE * f, const char *fd = FORMAT_DOUBLE) const
    {
      PrintCfgTag(num, f);
      PrintCfgCont(f, fd);
    }

  void PrintDPart(FILE * f)
  {
    for(int i = 0; i < nPart; i++)
      {
	fprintf(f, "\ndoor %d  ",i);
	dPart[i].Print(f);
      }
    }

    //Decompose input partition
  void DecomposePart(FILE *display);
  void DecomposePart(std::list<double> *dL);

  // Find the intervals of the partition which intersect
  // the alpha-cut door, fill the arrays of bounds
  // and return their number.
  int getIntersect(ACUT &door, double *inf, double *sup);

   //! Check if the input is a Strong Fuzzy Partition.
  //! If s is non NULL, the MF have re-ordered according to s.
  bool IsSfp(int *& s);

 //for GUI of fuzzy input
  void SetTemplate(double kw, double sw);
  double GetKernelWeightTemplate();
  double GetSupportWeightTemplate();

  // This function only works for normalized data between 0 and 1
  // and for a standard fuzzy partition
  // Use CheckFuzDist before and UnNormalize after.
  // norm=0: unnormalized distance
  // norm=1: normalized distance/(nmf-1)
  // d=1 debug mode
  // Return the distance between x and y
  double Distance(double x, double y, int norm=0, int d = 0);

  void CheckFuzDist(void);

  // Compute the partition coefficient and the partition entropy 
  // dat are are the data of size n
  // pc and pe are set by the function  
  void PcPe(double * dat, int n, double &pc, double &pe);
  void Tri2Trap(void);
}; // Fin classe FISIN

class FISOUT;   // For cross references
class RULE;
//******************************************************************
//
//
//                             CLASS IMPLI
//
//
//******************************************************************

//! Abstract class to define implication method.
//! Three are proposed 'rg', 'gg', 'gd'.
//! Each AGGREGIMP include an IMPLI object.
class IMPLI
//**********
{
  protected :

  public :
  //! compute the possibility function
  virtual MFDPOSS* ComputeDposs(MF *A, double degree)=0;
  virtual MFDPOSS* ComputeTnorme(MFDPOSS *A, double degree);

  IMPLI()   {}
  virtual ~IMPLI() {}
};


class IMPLIRG : public IMPLI
//*****************************
{
  protected :

    public :
    virtual MFDPOSS* ComputeDposs(MF *A, double degree);

  IMPLIRG() : IMPLI() {}

  ~IMPLIRG() { }
};

class IMPLIGD : public IMPLI
//*****************************
{
  protected :

    public :
    virtual MFDPOSS* ComputeDposs(MF *A, double degree);


  IMPLIGD() : IMPLI() {}

  ~IMPLIGD() { }
};

class IMPLIGG : public IMPLI
//*****************************
{
  protected :

    public :
    virtual MFDPOSS* ComputeDposs(MF *A, double degree);
    virtual MFDPOSS* ComputeTnorme(MFDPOSS *A, double degree);

  IMPLIGG() : IMPLI() {}

  ~IMPLIGG() { }
};
//******************************************************************
//
//
//                             CLASS AGGREG
//
//
//******************************************************************

//! Abstract class to define the rule aggregation method.
//! Three are proposed 'sum', 'max', 'imp'.
//! Each FISOUT include an AGGREG object.
class AGGREG
//**********
{
  protected :

    public :
    //! Fills the 'MuInfer' and 'RuleInfer' arrays (members of FISOUT)
    virtual void Aggregate(RULE ** r, int nr, FISOUT * O, double deg=1.0) = 0;

  AGGREG()   {}

  virtual ~AGGREG() {}
};


class AGGREGSUM : public AGGREG
//*****************************
{
  protected :


    public :
    virtual void Aggregate(RULE ** r, int nr, FISOUT * O, double deg =1.0);

  AGGREGSUM() : AGGREG() {}

  ~AGGREGSUM() { }
};

class AGGREGMAX : public AGGREG
//*****************************
{
  public :

    virtual void Aggregate(RULE ** r, int nr, FISOUT * O, double deg =1.0);

  AGGREGMAX() : AGGREG() {}

  ~AGGREGMAX() { }
};

class AGGREGIMP : public AGGREG
//*****************************
{
 protected :

 public :
  IMPLI * imp;

  virtual void Aggregate(RULE ** r, int nr, FISOUT * O, double deg =1.0);

  AGGREGIMP()  { imp= new IMPLIGD();}

  AGGREGIMP(IMPLI* impl) { imp = impl; }

  ~AGGREGIMP() { if(imp!=NULL) delete imp;}
};

//******************************************************************
//
//
//                             CLASS DEFUZ
//
//
//******************************************************************


//! Abstract class to define the defuzzification method.
//! Each FISOUT include a DEFUZ object.
//! The available types for a 'fuzzy' output are: 'area', 'sugeno', and 'MeanMax'.
//! Two can be used for a 'crisp' output:  'sugeno' and 'MaxCrisp'.

class DEFUZ
//*********
{
  protected :

    int NbClas;       // Number of classes

  //! This value, whose semantic depends on the operator type (see Alarm to get
  //! a detailled description), is set by default in
  //! order to not modify the existing constructors. It can be managed through the
  //! GetThres and SetThres functions.
  double Thres;


  public :

    DEFUZ()
    {
      NbClas = 0;
      Thres=0.0;
      Alarm = 0;
    }

  virtual ~DEFUZ() {}

  //! Various types of Alarm can be set when inferring depending on the
  //! defuzzification operator. They are managed by the EvalOut function and
  //! the current value can be print in the result file.
  //! The available values are:
  //! NOTHING             0 - All types. - No comment.
  //! NO_ACTIVE_RULE      1 - All types. - No comment.
  //! AMBIGUITY           2 - SugenoClassif and MaxCrisp - The difference
  //! between the two main classes is less than a threshold (default value
  //! AMBIGU_THRES        0.2).
  //! NON_CONNEX_AREA     3 - WeArea - Set when the area defined by the fired
  //! fuzzy sets (threshold  set to MIN_THRES = 0.1 by default) is non connex.
  //! NON_CONNEX_SEGMENT  4 - MeanMax - Set when the max corresponds to two
  //! fuzzy sets (with a tolerance threshold set by default to EQUALITY_THRES = 0.1)
  //! and the resulting segment is non connex.
  int Alarm;

  //! Compute the inferred output using the RuleInfer and MuInfer arrays filled by AGGREG::Aggregate().
  //! 'fa' is the archive file for Performance.
  //! It is handled independently from the display flag 'aff'.
  virtual double EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa = NULL, FILE *display = NULL) = 0;

  //! This function is called by FIS::WriteHeader, it fills the comment line
  //! with the labels of the columns written by each defuz operator.
  virtual void WriteHeader(FILE *p, FISOUT *O) const = 0;

  //! Filled the values max and i_max with the level of activity and the number of
  //! output MF. The initial values are all set to -1. When the difference between
  //! max1 and max2 is less than Thres, both are filled, otherwise only the first one.
  //! This function is used to detect ambiguity for MaxCrisp
  //! as well as to detect a non connex segment for MeanMax.
  void GetMax(FISOUT * O, double & max1, double & max2, int & i_max1, int & i_max2) const;


  int NbClasses(void) const { return NbClas; }
  double GetThres(void) const { return Thres; }
  void SetThres(double t) { Thres = t; }

};

#include "defuz.h"  // Inherited class definitions

//*****************************************************************
//
//
//                             CLASS FISOUT
//
//
//******************************************************************

//! This abstract class is inherited from FUSIN, as an output description includes an input one.
//! Two types of FISOUT are to be defined 'OUT_CRISP' and 'OUT_FUZZY'.
class FISOUT : public FISIN
//**************************
{
  protected :
    //!   Defuzzification operator
    char * Defuz;
  //!   Aggregation operator
  char * Disj;   //!   Aggregation operator
  //! Output Default value when no rule is activated
  double Default;

  //! True if the output is a class label
  int Classif;

  //! Number of distinct conclusions
  int NbPossibles;
  //! Values of distinct conclusions (size 'NbPossibles')
  double * Possibles;

  //! Set all values to default and pointers to NULL
  //! except the explicit parameters.
  void Init( const char *OpDefuz,  const char *OpDisj, int ccl, double VDef)
    {
      //warning! inits to NULL must be done before functions calls that could fail
      Defuz = NULL;
      Disj = NULL;
      Default = VDef;
      Classif = ccl;
      Possibles = MuInfer = NULL;
      RuleInfer = ConcInfer = NULL;
      NbPossibles = 0;
      Def = NULL;
      Ag = NULL;
      SetOpDefuz( OpDefuz );
      SetOpDisj( OpDisj );
      MfConc= NULL;     // Only used by implicative rules
      MfGlob = NULL;    // Only used by implicative rules
    }

  //! Read output parameters from a file
  void Init( ifstream &f, int bufsize, int num, const char *OpDefuz, const char *OpDisj, int ccl, double VDef );

  AGGREG * Ag;

 public:

  DEFUZ *Def;

  //! For each distinct output value, the degree of match for the current example
  //! cumulated (Sum) or max (Max).
  //! This array is fitted by AGGREG::Aggregate called by FISOUT::InferOut (size 'NbPossibles').
  //! Attention:  there is a shift of 1 between the informations
  //! displayed (the first element is 1) and the machine implementation (the first element is 0).

  double * MuInfer;

  //! For each distinct output value, the rule number corresponding to
  //! the max (if disjunction max) or the last (if disjunction sum)
  //! degree of match for the current example.
  //! This array is fitted by AGGREG::Aggregate called by FISOUT::InferOut (size 'NbPossibles').
  //! Attention:  there is a shift of 1 between the informations
  //! displayed (the first element is 1) and the machine implementation (the first element is 0).
  int * RuleInfer;

  //! For each rule, each conclusion is labelled by its index in the array Possibles.
  //! This table is suitable to save time in the infer function (size NbRules).
  int * ConcInfer;

  //! the inferred MF for the current example
  //! Only used by implicative rules
  MFDPOSS * MfGlob;

  //! For each rule the MF for the current example
  //! This array is fitted by FISOUT::InferOut (size NbRules)
  //! Only used by implicative rules
  MFDPOSS ** MfConc;

  void init() {
	  Def = NULL;
	  MuInfer = NULL;
	  RuleInfer = NULL;
	  ConcInfer = NULL;
	  MfGlob = NULL;
	  MfConc = NULL;
  }

  FISOUT(int n=0) : FISIN(n)
    {
	  init();
    };

  FISOUT(int n, double min, double max, int tri=false)
    : FISIN(n, min, max, tri)
    {
	  init();
    };


  FISOUT(double *t, int n, double min, double max, int sort = true)
    : FISIN(t, n, min, max, sort)
    {
	  init();
    };

  FISOUT(double *AD,int *CT, int i1, double d1, double d2,double b1,double b2,int i2)
    : FISIN(AD,CT,i1,d1,d2,b1,b2,i2 )
    {
	  init();
    };

  //! The copy constructor does not copy the 'MuInfer', 'RuleInfer' and 'ConcInfer' arrays.
  //! The size of 'ConcInfer' array, the number of rules, is unknown here.
  FISOUT( FISOUT &sortie ) : FISIN( sortie ) {
	  init();
  }

  virtual ~FISOUT()
    {
      if (Defuz!=NULL) delete [] Defuz;
      if (Disj != NULL) delete [] Disj;
      if (Def != NULL) delete Def;
      if (Ag != NULL) delete Ag;
      DeletePossibles();
    }


  int operator != (const FISOUT &);
  //! This function will be overloaded in inherited class to control the authorized
  //! defuzzification operator types according to the output nature (crisp or fuzzy).
  virtual void SetOpDefuz( const char *op_defuz );
  //! This function will be overloaded in inherited class to control the authorized
  //! aggragation operator types according to the output nature (crisp or fuzzy).
  virtual void SetOpDisj( const char *op_disj );
  const char *GetDisj() { return Disj; };
  //! This value is inferred when none of the rules are activated.
  void DefaultValue( double default_value ) { Default = default_value; }
  void Classification( int classif );
  const char *GetType() const { return "Output"; };

  char * Defuzzify(void) { return Defuz; }
  char * Disjunct(void) { return Disj; }
  double DefaultValue(void) { return Default; }
  int Classification(void) { return Classif; }

  MFDPOSS *GetMFGlob() { return MfGlob; }
 MFDPOSS *GetMFConc( int rule_number ) { return MfConc[rule_number]; }

  //! Returns either 'crisp' or 'fuzzy'.
  virtual const char *GetOutputType() const = 0;
  virtual FISOUT *Clone() = 0;
  //! check MF types for implicative output
  void  CheckImpliMFs();
  void  CheckImpliMF(MF* mfcandidate);
  //! Update the output configuration defuzzification and aggregation operators,
  //! size of 'RuleInfer' and 'MuInfer' arrays, and fills the 'ConcInfer' correspondance table.
  //! This function should be called before the inference after a change in the rule base.
  //! During the init process it cannot be within FISOUT constructor because the
  //! OUTPUT(s) are set up previously to RULEs.
  //! 'nr' is the number of rules,
  //! 'nb' is the output number.
  void InitPossibles(RULE ** r, int nr, int nb);
  //! Update Possible array following call to Rule[changedrulenum] ->SetConc
  //! to be called by learning modules (fisopt ...)
  //! or from the Java interface
  void UpdatePossibles(RULE ** r, int nr, int changedrulenum, int nb);
  //! Free the pointers allocated by InitPossibles. Used by InitPossibles and ~FISOUT.
  void DeletePossibles();//called by FISOUT destructor
  void DeletePossibles(int nbrul);//called by InitPossibles
  void DeleteMFConc (int nbrul) const;
  void DeleteMFConcArray();
  //! Returns the number of distinct possible conclusions
  int GetNbPossibles(void) { return NbPossibles; }

  //! Returns a pointer on the values of the distinct possible conclusions
  double * GetPossibles(void) { return Possibles; }

  //! Filled 'MuInfer' and 'RuleInfer' arrays with default values.
  //! Used by InitPossibles() and AGGREG::Aggregate().
  void InitTabRes(void)
    {
      for(int j = 0; j < NbPossibles; j++)
	{ MuInfer[j] = 0.; RuleInfer[j] = -1;}
    }


  //! Shift the rule numbers to set the the first rule to 1 instead of 0.
  //! Used by InferOut.
  void ShiftRuleNbs(void)
    {
      for(int j = 0; j < NbPossibles; j++)
	RuleInfer[j] ++;
    }

  virtual void Print(FILE * f) const
    {
      char ChaineClassif[4];
      int i;
      if(Classif) sprintf(ChaineClassif, "yes");
      else sprintf(ChaineClassif, "no");

      FISIN::Print(f);
      fprintf(f, "\nOutput %s   Defuzzification : %s Classification : %s\n",
	      GetOutputType(), Defuz, ChaineClassif);
      fprintf(f, "\nDefault value : %11.3f", Default);
      fprintf(f, "\nNb of possibles conclusions: %d (", NbPossibles);
      for(i = 0; i < NbPossibles; i ++) fprintf(f, "%11.3f", Possibles[i]);
      fprintf(f, ")\n");
      fprintf(f, "\nMuInfer et RuleInfer : ");
      for(i = 0; i < NbPossibles; i ++)
	fprintf(f, "%11.3f %d", MuInfer[i], RuleInfer[i]);
    }

  //! Encapsulation of aggregation and defuzzification processes.
  //! Returns the inferred value. Called by FIS::Infer().
  double InferOut(RULE ** r, int nr, FILE * fic, FILE *display, double deg =1.0)
    {
      double val;

      Ag->Aggregate(r, nr, this, deg);
      val = Def->EvalOut(r, nr, this, fic, display);
      ShiftRuleNbs();

      return val;
    }

  //! Print output configuration.
  //! Function used to store a FIS configuration
  virtual void PrintCfg(int num, FILE * f, const char *fd = FORMAT_DOUBLE) const
    {
      char ChaineClassif[4];
      if(Classif) sprintf(ChaineClassif, "yes");
      else sprintf(ChaineClassif, "no");

      FISIN::PrintCfgTag(num, f);
      fprintf(f, "Nature=%c%s%c\n", STRING_SEP, GetOutputType(), STRING_SEP);
      fprintf(f, "Defuzzification=%c%s%c\n", STRING_SEP, Defuz, STRING_SEP);
      fprintf(f, "Disjunction=%c%s%c\n", STRING_SEP, Disj, STRING_SEP);
      fprintf(f, "DefaultValue=");
      fprintf(f, fd, Default);
      fprintf(f, "\n");
      fprintf(f, "Classif=%c%s%c \n", STRING_SEP, ChaineClassif, STRING_SEP);
      FISIN::PrintCfgCont(f, fd);
    }
  // SPECIFIC TESTS REQUIRED FOR MF TYPES (IMPLICATIVE CASE)
  void ReplaceMF(int numMF, MF* mfnew);
  void AddMF(MF* mfnew);

}; // Fin classe FISOUT


//! Inherited from the abstract class FISOUT to implement a crisp output.
class OUT_CRISP : public FISOUT
//*****************************
{

  public :

  OUT_CRISP(ifstream &f, int bufsize, int num, const char *OpDefuz, const char *OpDisj, int ccl, double VDef = OUT_DEFAUT )
    : FISOUT()
    {
      FISOUT::Init( f, bufsize, num, OpDefuz, OpDisj, ccl, VDef );
      if( Nmf != 0 )
	{
	  sprintf( ErrorMsg, "~Output~%d~:~NoMfAllowedForCrispOutput~\n", num);
	  throw std::runtime_error( ErrorMsg );
	}
    }

  OUT_CRISP(double inf, double sup, const char *OpDefuz, const char *OpDisj, int ccl, double VDef = OUT_DEFAUT )
    : FISOUT()
    {
      FISOUT::Init(OpDefuz, OpDisj, ccl, VDef );
      SetRange(inf, sup);
    }

  OUT_CRISP() : FISOUT()
    // Nmf = 0 par defaut, pas de classif
    {
      FISOUT::Init( (char *)SugenoDefuz(), (char *)DisjSum(), 0, OUT_DEFAUT );
    }

  OUT_CRISP( OUT_CRISP &sortie ) : FISOUT(sortie)
    {
      FISOUT::Init( sortie.Defuzzify(), sortie.Disjunct(), sortie.Classification(), sortie.DefaultValue() );
    }

  OUT_CRISP( FISOUT &sortie ) : FISOUT(sortie)
    {
      if( strcmp( sortie.GetOutputType(), OUT_CRISP::OutputType() ) == 0 )
	FISOUT::Init( sortie.Defuzzify(), sortie.Disjunct(), sortie.Classification(),
		      sortie.DefaultValue() );
      else
	FISOUT::Init( OUT_CRISP::SugenoDefuz(), OUT_CRISP::DisjSum(),
		      sortie.Classification(), sortie.DefaultValue() );
      // if fuzzy output delete mfs
      while( Nmf != 0 ) RemoveMF( 0 );
    }


  virtual ~OUT_CRISP() {}


  FISOUT *Clone() { return new OUT_CRISP(*this); }
  const char *GetOutputType() const { return OutputType(); }
  static const char *OutputType() { return "crisp"; }
  static const char *SugenoDefuz() { return "sugeno"; }
  static const char *MaxCrispDefuz() { return "MaxCrisp"; }
  //! 'sugeno' and 'MaxCrisp' are authorized.
  void SetOpDefuz( const char *op_defuz );

  //! Both 'max' and 'sum' are allowed.
  void SetOpDisj( const char *op_disj );
  static const char *DisjSum() { return "sum"; }
  static const char *DisjMax() { return "max"; }
};



class OUT_FUZZY : public FISOUT
//*****************************

{
  public :
    OUT_FUZZY(ifstream &f, int bufsize, int num, const char *OpDefuz, const char *OpDisj, int ccl, double VDef, int Cover = false)
    : FISOUT()
    {
      FISOUT::Init( f, bufsize, num, OpDefuz, OpDisj, ccl, VDef );
      if( Nmf == 0 )
        {
        sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~NumberOfMFInOutput~ = 0", num );
        throw std::runtime_error( ErrorMsg );
        }
      if(Cover && Nmf > 1 && !(strcmp(Fp[0]->GetType(), MFTRAPINF::Type())) &&
	 !(strcmp(Fp[Nmf-1]->GetType(), MFTRAPSUP::Type())) )
	OutCoverage();
    }

  //! creates a fuzzy output with discrete MFs for values in t array
  OUT_FUZZY(int discrete, double *t, int n, double min, double max, const char *OpDefuz,
	    const char *OpDisj, int classif, double VDef = OUT_DEFAUT) :  FISOUT()
    {
      FISOUT::Init(OpDefuz, OpDisj, classif, VDef );
      InitDiscrete(t,n, min, max);
    }

  //! creates a fuzzy output with a standardized fuzzy partition
  OUT_FUZZY(double *t, int n, double min, double max, int sort, const char *OpDefuz,
	    const char *OpDisj, int ccl, double VDef = OUT_DEFAUT) : FISOUT(t, n, min, max, sort)
    {
      FISOUT::Init( OpDefuz, OpDisj, ccl, VDef );
      if(n > 1) OutCoverage();
    }

  OUT_FUZZY(double *t, int * tt, int n, double min, double max, double omin,  double omax,
	    int sort, const char *OpDefuz, const char *OpDisj, int ccl, double VDef = OUT_DEFAUT)
    : FISOUT(t, tt, n, min, max, omin, omax, sort)
    {
      FISOUT::Init( OpDefuz, OpDisj, ccl, VDef );
      if(n > 1) OutCoverage();
    }

  OUT_FUZZY(int n, double min, double max, const char *OpDefuz, const char *OpDisj, int ccl,
	    double VDef = OUT_DEFAUT, int tri=false) : FISOUT(n, min, max, tri)
    {
      FISOUT::Init( OpDefuz, OpDisj, ccl, VDef );
      if(n > 1) OutCoverage();
    }

  OUT_FUZZY( OUT_FUZZY &sortie ) : FISOUT(sortie)
    {
      FISOUT::Init( sortie.Defuzzify(), sortie.Disjunct(), sortie.Classification(), sortie.DefaultValue() );
    }

  OUT_FUZZY( FISOUT &sortie ) : FISOUT(sortie)
    {
      if( strcmp( sortie.GetOutputType(), OUT_FUZZY::OutputType() ) == 0 )
	FISOUT::Init( sortie.Defuzzify(), sortie.Disjunct(), sortie.Classification(), sortie.DefaultValue() );
      else
	FISOUT::Init( OUT_FUZZY::AreaDefuz(), OUT_FUZZY::DisjSum(), sortie.Classification(), sortie.DefaultValue() );
    }


  virtual ~OUT_FUZZY() {}

  void InitDiscrete(double *t, int n, double min, double max);
  FISOUT *Clone() { return new OUT_FUZZY(*this); }
  const char *GetOutputType() const { return OutputType(); }
  static const char *OutputType() { return "fuzzy"; }
  static const char *AreaDefuz() { return "area"; }
  static const char *MeanMaxDefuz() { return "MeanMax"; }
  static const char *SugenoDefuz() { return "sugeno"; }
  static const char *ImpFuzzyDefuz() { return "impli"; }
  //! Authorized oprators:  'area', 'sugeno', 'MeanMax'.
  void SetOpDefuz( const char *op_defuz );

  //! 'max'  'sum' 'imp' are allowed.
  void SetOpDisj( const char *op_disj ); //throw( runtime_error );
  static const char *DisjSum() { return "sum"; }
  static const char *DisjMax() { return "max"; }
  static const char *DisjIrg() { return "irg"; }
  static const char *DisjIgg() { return "igg"; }
  static const char *DisjIgd() { return "igd"; }
  //! Only use with partitions ended by semi-trapezoidal shaped fuzzy
  //! sets, when the range limits belong to the scorresponding kernels.
  //! Modify the extreme bouds of each semi-trapezoidal MF to ensure
  //! the range limits be inferred.
  void OutCoverage(void);
  // symbolic match between inferred and observed outputs
  double SymbMatch(double * muObs, double * muInf,int nmf, double mutObs, double mutInf, FILE *display=NULL);

  //! Convert a Strong Fuzzy Partition to a Quasi Strong Partition
  //! used by implicative rules.
  //! If s is non NULL, the MF have re-ordered according to s.
  //! Partition limits are set to input range bounds.
  //! Sfp2Qsp returns -1 if Nmf=0, -2 if output partition is not SFP, 0 otherwise
  int  Sfp2Qsp(int *& s);
  //! Convert a Quasi Strong Partition used by implicative rules
  //! to a Strong Fuzzy Partition.
  //! If s is non NULL, the MF have re-ordered according to s.
  //! Call the OutCoverage function to ensure the range limits be inferred.
  //! if argument onlytestQSP is true, no transformation is done, and the function
  //! returns true if transformation into SFP  was possible, i.e. initial partion was QFP
  //! false otherwise
  //! if argument onlytestQSP is false, the transformation is done if possible, and the function
  //! returns true, or false otherwise
  bool Qsp2Sfp(int *& s, bool onlytestQSP=false);
  //! returns true if Qsp partition, false otherwise
  //! calls Qsp2Sfp(s,true)
  bool IsQsp();
};


//******************************************************************
//
//
//                             CLASS RULE
//
//
//******************************************************************


#include "rule.h" //! PREMISE and CONCLUSION classes


//! The class RULE is made up of two distinct objects.  The PREMISE deals with the
//! input part of the rule and compute the the degree of match to the rule of a given input
//! by the function MatchDeg. The PREMISE class is an abstract one:  three opretaors are
//! implemented:  'min', 'prod', 'luka'.
//! The MatchDeg fonction uses the MfDeg array of each input, which have to be previously filled
//! by calling the fonction FISIN::GetDegs.
//! The CONCLUSION class is used only to described the conclusion part of the rule. It could
//! be used to implement different implication operators.
//! This class contains different comparison methods and private data handling functions.
class RULE
//*********
{
  protected :
    PREMISE * Prem ;
  CONCLUSION * Conclu;
  //! A rule may be activated or deactivated using this flag.
  int Active;
  //! Expert rules can be manually weighted
  double ExpertWeight;

 public:

  //! The matching degree, which also depends on the conjunction operator
  double Weight;
  //! To store the cumulated weight of the rule for a training set for instance.
  double CumWeight;
  //! Used by some learning methods, to store the number of examples which activate the rule.
  int NbOccur;


  //! Val corresponds to a line in the configuration file.
  //! It is made up of a number of integers egal to the number of inputs and
  //! followed by a number of doubles egal to the number of outputs.
  //! The field separator is SEPARE.
  RULE(int nI, FISIN ** E, int nO, FISOUT ** S, char * cConj, char * Val); //throw(runtime_error);


  RULE() { Init(); } //! Empty constructor  necessary for derived classes


  RULE(RULE &, FISIN ** E, FISOUT ** S);

  RULE(RULE &, FISIN ** E);


  virtual ~RULE()
    {
      delete Prem;
      delete Conclu;
    }


  int operator != (const RULE &);


  void Init()
    {
      Prem = NULL;
      Conclu = NULL;
      ExpertWeight = 1.0;
    }


  //! types de premisse
  static const char *PremiseProd() { return "prod"; }
  static const char *PremiseMin() { return "min"; }
  static const char *PremiseLuka() { return "luka"; }


  void SetPremise( int nI, FISIN ** E, char * cConj );
  void SetConclusion( int nO, FISOUT ** S );


  //! The implicit argument of Prem->MatchDeg is the membership vector
  //! of each input (MfDeg). It must be correctly initialized.
  double MatchDeg(void)
    {
      return Prem->MatchDeg();
    }

  //! Set the rule weight to the matching degree.
  //! The implicit argument of Prem->MatchDeg is the membership vector
  //! of each input (MfDeg). It must be correctly initialized.
  void ExecRule(void)
    {
      if(Prem) Weight = Prem->MatchDeg();
    }


  // Miscellaneous printing functions
  virtual void Print(FILE * f)
    {
      Prem->Print(f);
      Conclu->Print(f);
      if(! Active) fprintf(f, "  Inactive ");
      else  fprintf(f, "          ");
      fprintf(f, "\n");
    }


  //! Print membership function  configuration.
  //! Function used to store a FIS configuration

  virtual void PrintCfg(FILE * f, const char *fd = FORMAT_DOUBLE, bool ExpW = false) const
    {
      Prem->Print(f);
      Conclu->Print(f, fd);
      if(ExpW) fprintf(f, fd, ExpertWeight);
      fprintf(f, "\n");
    }

  // Handling private data
  virtual void PrintProps(FILE * f) const
    { Prem->PrintProps(f); }

  // Print premises with separation (printprops pritn premises without separation)
  virtual void PrintPrems(FILE * f) const
    { Prem->Print(f);}


  int IsActive(void) { return Active; }
  void Activate(void) { Active = true; }
  void Deactivate(void) { Active = false; }


  int GetNbProp(void) { return Prem->GetNbProp(); }
  int GetNbConc(void) { return Conclu->GetNbConc(); }


  //! Returns the number of different labels in the  premise.
  //! 'diff' contains the number of the first different label
  //!  and -1 if the rules are identical.
  //! 'indiff' contains the number of differences between a given factor and
  //! an  indifferent factor (the partial distances are equal to zero).
  int Distance(RULE *R, int & diff, int & indiff)

    { return Prem->Distance(R->Prem, diff, indiff); }


  //! Returns 0 if the rules are the same, otherwise 1.
  int ComparePremises(RULE *R)
    { return Prem->Compare(R->Prem); }


  //!Returns 0 if the rules are the same, otherwise 1.
  int Compare(RULE *R)
    {
      if(! Prem->Compare(R->Prem)) return Conclu->Compare(R->Conclu);
      return 1;
    }


  int SetAProps(RULE *R, int alloc = true)
    { return Prem->SetAProps(R->Prem, alloc); }


  void SetAProps(int *Tab)
    { Prem->SetAProps(Tab); }


  int SetAProp(int f, int num)
    { return Prem->SetAProp(f, num); }


  void GetProps(int *Tab)
    { Prem->GetProps(Tab); }


  int GetAProp(int &f, int num)
    { return Prem->GetAProp(f, num); }


  void SetConcs(double * Tab)
    { Conclu->SetConcs(Tab); }

  void Normalize ( int nOut ,double lowerb , double upperb )
    {
      Conclu->Normalize ( nOut ,lowerb , upperb );
    }

  void UnNormalize ( int nOut ,double lowerb , double upperb )
    {
      Conclu->UnNormalize ( nOut ,lowerb , upperb );
    }
  void SetAConc(int num, double Val)
    { Conclu->SetAConc(num, Val);  }

  double GetAConc(int k)
    { return Conclu->GetAConc(k); }

  void SetExpertWeight(double ew) 
  { 
    if(ew > EPSILON) ExpertWeight = ew; 
    else
      {
	sprintf( ErrorMsg, "~ExpertWeight~MustBePositive~: %f\n", ew);
	throw std::runtime_error( ErrorMsg );
      }
  }

  double GetExpertWeight() 
    { return ExpertWeight; }

}; // End of  RULE class




//******************************************************************
//
//
//                             CLASS FIS
//
//
//******************************************************************

int CmpCumDec(const void * a, const void *b);
int CmpCumInc(const void * a, const void *b);

// This class handle a Fuzzy Inference System.

class FIS
//*******
{
 private:
  void Init(); //throw();


 protected:


  char * cConjunction, * strMissingValues, *strErrorIndex;
  int NbIn, NbOut, NbRules, NbExceptions, NbActRules;


  //! Fonction called by constructors to read a FIS configuration file.
  virtual void InitSystem(const char * fichier, int Cover = false);
  //! Fonction called by InitSystem to read the header of the configuration file.
  virtual void ReadHdr(ifstream &f, int bufsize); //throw(runtime_error);
  //! Fonction called by InitSystem to read an input description from the configuration file.
  virtual void ReadIn(ifstream &f, int bufsize, int num); //throw(runtime_error);
  //! Fonction called by InitSystem to read an output description from the configuration file.
  virtual void ReadOut(ifstream &f, int bufsize, int num, int Cover = false);
  //! Fonction called by InitSystem to read rule description from the configuration file.
  virtual void ReadRules(ifstream &f, int bufsize); //throw(runtime_error);
  //! Fonction called by InitSystem to read exce4ption description from the configuration file.
  virtual void ReadExcep(ifstream &f, int bufsize);



 public:

  //! Pointer on the output array. Attention:  there is a shift of 1 between the informations
  //! displayed (the first element is 1) and the machine implementation (the first element is 0).
  FISOUT ** Out;

  //! Pointer on the input array. Attention:  there is a shift of 1 between the informations
  //! displayed (the first element is 1) and the machine implementation (the first element is 0).
  FISIN ** In;

  //! Pointer on the rule array. Attention:  there is a shift of 1 between the informations
  //! displayed (the first element is 1) and the machine implementation (the first element is 0).
  RULE ** Rule;

  char * Name;
  //! Infered values computed by Infer.
  double * OutValue,
    //! Difference between observed and infered values.
    * OutErr;

  // Error indices for regression problems computed by the Perf function
  double PIn, RMSE, MAE;

  FIS(const char *fichier, int Cover = false) //throw(runtime_error)
    {
      Init();
      InitSystem(fichier, Cover);
    }


  FIS();


  FIS(const FIS &);

  FIS& operator=(const FIS&);

  int operator != (const FIS &) const;


  virtual ~FIS()
    {
      int i,irule;
      // check that all ptrs are non NULL before calling delete
      if (In != NULL)
	{
	  for( i=0 ; i<NbIn ; i++)
	    if (In[i]!= NULL) delete In[i];
	  delete [] In;
	}

      if (Out != NULL)
	{
	  for( i=0 ; i<NbOut ; i++)
	    {
	      if (Out[i]!= NULL)
		{
		  if( Out[i]->MfConc != NULL )
		    {
		      for (irule=0;irule<NbRules;irule++)
			delete Out[i]->MfConc[irule];
		      delete [] Out[i]->MfConc;
		      Out[i]->MfConc=NULL;
		    }
		  if (Out[i]->MfGlob!=NULL)
		    {
		      delete Out[i]->MfGlob;
		      Out[i]->MfGlob=NULL;
		    }
		  delete Out[i];
		  Out[i]=NULL;
		}
	    }
	  delete [] Out;
	  Out=NULL;
	}
      if (Rule != NULL)
	{
	  for( i=0 ; i<NbRules ; i++)
	    if (Rule[i]!= NULL)  delete Rule[i];
	  delete [] Rule;
	  Rule=NULL;
	}
      if (OutValue != NULL) delete [] OutValue;
      OutValue=NULL;
      if (OutErr != NULL) delete [] OutErr;
      OutErr=NULL;
      if (Name != NULL) delete [] Name;
      Name=NULL;
      if (cConjunction != NULL) delete [] cConjunction;
      cConjunction=NULL;
      if (strMissingValues != NULL) delete [] strMissingValues;
      strMissingValues=NULL;
      if (strErrorIndex != NULL) delete [] strErrorIndex;
      strErrorIndex=NULL;
    }

  void SetName( const char *name );
  void SetConjunction( const char *conjunction );
  void SetMissingValues( const char *missing_value );
  void SetErrorIndex( const char *index );
  void AddInput( FISIN *entree );
  void AddOutput( FISOUT *sortie );
  void RemoveInput( int input_number );
  void ReplaceInput( int input_number, FISIN *new_input );
  void RemoveOutput( int output_number );
  void ReplaceOutput( int output_number, FISOUT *new_output );
  void RemoveMFInInput( int input_number, int mf_number );
  void RemoveMFInOutput( int output_number, int mf_number );
  // change output to classif type - calls initdefuz to change defuz object
  void AddRule( RULE *regle );
  void RemoveRule( int rule_number );
  void RemoveAllRules(void);
  void DeleteMFConc(int output_number) const;
  void DeleteMFConcArray(int output_number);
    //! Returns the array index  corresponding to the first identical rule
  //! and -1 if the rule is not in the rule base.
  //! Research can be restarted on the base by modifying the starting index value.
  //! The equality test is made on the premise part only, if 'conc' is true,
  //! the test includes consequent part comparison.
  virtual int RulePos(RULE *R, int depart = 0, int conc = false);

  //! Compute the matching degree of input data 'values' for all
  //! the active rules and fill the 'weights' array.
  //! The 'weights' array has to be previously allocated.
  void RuleWeights (double * values, double * weights);

  // fis normalisation and data normalisation if u give one.
  void Normalize(double **SampleData=NULL,int nbrow=0);
  // fis unnormalisation and data unnormalisation if u give one.
  void UnNormalize(double **SampleData=NULL,int nbrow=0);

  static const char *RandomMissingValues() { return "random"; }
  static const char *MeanMissingValues() { return "mean"; }

  static const char *PiErrorIndex() { return "PI"; }
  static const char *RmseErrorIndex() { return "RMSE"; }
  static const char *MaeErrorIndex() { return "MAE"; }

  //! To turn a fuzzy output to a crisp one, 'o' being the output number
  void Fuz2Crisp(int o);

  //! To turn a crisp output to a fuzzy one, 'o' being the output number.
  //! 'c' is the array with the fuzzy set centers and 'nc' is its dimension.
  //! Default (c = NULL) : The number of MF in the fuzzy output corresponds
  //! to the number of distinct conclusions in the rule base
  void Crisp2Fuz(int o, const char * DefuzType, double * c = NULL, int nc = 0);
  // FIS2Qsp urns the strong fuzzy partition (SFP) into a
  // quasi strong fuzzy partition (QSP)
  // If needed re-order the meembership functions according to
  // their kernel in an ascending order
  // updates the rules
  // Return value:
  // The MF may have been re-ordered - see sorted argument in Sfp2Qsp
//In that case if needed re-order the membership functions according to
//             their kernel in an ascending order
//             Return value:
//             2: output partition is already Qsp
//             -1: Invalid output number
//             -3: the output is not fuzzy
//             -4:  Nmf=0
//             -5 : Problem in function Sfp2Qsp-initial output partition is not Sfp
//              0 if Qsp2Sfp has been done OK
//              1 if Qsp2Sfp has been done OK and if MFs have been re-sorted
//              2 if output partition is already a Qsp partition

  int FIS2Qsp(int o, const char * DisjType);
  // Convert a conjunctive fuzzy output into an implicative one
  // no change in the ouptut partition if transfPart=false
  // calls Conj2ImpQFP and changes the ouptut partition if transfPart=true
  // if transfPart=true calls Sfp2Qsp and
//                   turns the  strong fuzzy partition (SFP) into a quasi strong fuzzy partition (QSP)
//In that case if needed re-order the membership functions according to
//             their kernel in an ascending order
//             Return value:
//             -1: Invalid output number
//             -2: Already an implicative output
//             -3: the output is not fuzzy
//             -4:  Nmf=0
//             -5 : Problem in function Sfp2Qsp - initial output partition is not Sfp
//              0 if Qsp2Sfp has been done OK
//              1 if Qsp2Sfp has been done OK and if MFs have been re-sorted
//              2 if output partition is already a Qsp partition
// if transfPart=false returns 1 if partition is already Sfp, 0 otherwise
  int Conj2Imp(int o, const char * DisjType, bool transfPart);
  // FIS2Sfp turns a Qsp partition into a Sfp one
  // updates the rules
  // it returns
  // -1 if invalid outputnumber
  // -3 the output is not fuzzy
  // -4 if  Nmf=0
  // -5 if transf. into Sfp does not work
  // 0 if  transf. into Sfp has been done OK
  // 1 if transf. into SfpQsp has been done OK and if MFs have been re-sorted
  // 2 if output partition is already an Sfp partition
  int FIS2Sfp(int o, const char * DefuzType, const char * DisjType);

  // Convert an implicative output into a conjunctive fuzzy one
// if transfPart=true calls Qsp2Sfp and
//                   turns the quasi strong fuzzy partition (QSP) into a strong fuzzy partition (SFP)
//In that case if needed re-order the membership functions according to
//             their kernel in an ascending order
//             Return value:
//             2: output partition is already Sfp
//             -1: Invalid output number
//             -2: Already a conjunctive output
//             -3: the output is not fuzzy
//             -4:  Nmf=0
//             -5 : Problem in function Qsp2Sfp-output partition is not Qsp
//              0 if Qsp2Sfp has been done OK
//              1 if Qsp2Sfp has been done OK and if MFs have been re-sorted
//              2 if output partition is already an Sfp partition
// if transfPart=false returns 1 if partition is already Sfp, 0 otherwise
  int Imp2Conj(int o, const char * DefuzType, const char * DisjType, bool transfPart);

  //! Infer the output values according to the input one and the FIS configuration.
  //! The inference mechanism is decribed in the user documentation.
  //! 'v' is the item description (input values and optionally output ones) using crisp numbers,
  //! 'out_number' is the output to infer, if -1 all outputs are inferred,
  //! 'fic' is the file in which the result is printed,
  //! 'display' is the file in which intermediate result are displayed.
  //! Returns the maximum degree of match of the example over the rule base.
  //! The inferred value is stored in the 'OutValue' array.
  //! When the observed output value is part of 'v', the difference between
  //! observed and inferred is stored in the 'OutErr' array.

  virtual double Infer(double * v, int out_number = -1, FILE * fic = NULL, FILE *display = NULL, double deg = 1.) const;
  //! Infer the output values according to the input one and the FIS configuration.
  //! The inference mechanism is decribed in the user documentation.
  //! 'v' is the item description using fuzzy sets,
  //! the input values could be fuzzy sets or crisp values
  //! (use the MFDISCRETE constructor to turn a crisp value into a fuzzy set),
  //! 'out_number' is the output to infer, if -1 all outputs are inferred,
  //! 'fic' is the file in which the result is printed,
  //! 'display' is the file in which intermediate result are displayed.
  //! Returns the maximum degree of match of the example over the rule base.
  //! The inferred value is stored in the 'OutValue' array.
  //! When the observed output value is part of 'v', the difference between
  //! observed and inferred is stored in the 'OutErr' array.
  double Infer(MF ** v, int out_number, FILE * fic=NULL, FILE *display = NULL) const;

  //! Measures the performance of the FIS for the 'NumS' output, for the data set 'fdata'.
  //! 'Coverage' is the coverage level.
  //! 'MaxErr' is the maximum of the error over the data set.
  //! 'MuThresh' is a threshold, the examples whose cumulated weight over the whole rule base
  //! is less than 'MuSeuil' are not taken into account in the performance index.
  //! 'fres' is the name of the output file.
  //! The output formatted for a spreadsheet, the number of columns depends on the configuration
  //! (defuzzification operator, flag classif, observed output in the data set)
  //! 'display' is the file in which intermediate result are displayed.
  //! Returns -2 en case of inconsistency (not enough input variables or not active output),
  //! -1  when the observed output is not part of the sample file, the performance index otherwise:
  //! the number of misclassified for a crisp output with classif, or the root mean square of errors.
  //! The outfile file, 'fres', structure is as follow, a line corresponds to an example:
  //! Obs (if available) - Defuz (see below) - Err (if available) - Blank - Perf.
  //! 'Obs' is the observed output value, if it is available 'Err' indicates the difference between
  //! observed and inferred output values, 'Blank' is an integer whose value is 1 if the cumulated
  //! weight is less than 'MuSeuil', else 'Blank' is set to 0, the last colum is the cumulated error
  //!  or performance index.
  //! The fields filled by the defuzzification operator are:
  //! For a crisp output :
  //! 'sugeno' and classif= 'no' : the inferred value - Alarm
  //! 'sugeno' and classif= 'yes': the inferred value - Alarm - the assigned class - Alarm
  //! 'MaxCrisp' : the inferred value - Alarm
  //! For a fuzzy output :
  //! 'sugeno' and classif= 'no' : the inferred value - Alarm
  //! 'sugeno' and classif= 'yes': the inferred value - Alarm -
  //! the membership degree of the inferred value to each of the output fuzzy sets.
  //! 'area' and classif= 'no' : the inferred value - Alarm
  //! 'area' and classif= 'yes': the inferred value  - Alarm -
  //! the membership degree of the inferred value to each of the output fuzzy sets.
  //! 'MeanMax' and classif= 'no' : the inferred value - Alarm
  //! 'MeanMax' and classif= 'yes': the inferred value - Alarm  -
  //! the membership degree of the inferred value to each of the output fuzzy sets, if and only if the example activate at least one rule.
  virtual double Performance(int NumS, const char * fdata, double & Coverage, double & MaxErr, double MuThresh,  const char * fres = NULL, FILE *display = NULL);


  //! Measures the performance of the FIS for the 'NumS' output, for the data array Data
  //! number of lines of this array 'NbEx'
  //! other  arguments required for classification only, and not modified:
  //! MisClassified and Lab
  //! must be preallocated (see ResClassifAlloc)
  // The function modifies the following arguments:
  //! 'Coverage' is the coverage level.
  //! 'MaxErr' is the maximum of the error over the data set.
  //! 'MuThresh' is a threshold, the examples whose cumulated weight over the whole rule base
  //! is less than 'MuSeuil' are not taken into account in the performance index.
  //! optional: 'fres' is the name of the output file.
  //! The output formatted for a spreadsheet, the number of columns depends on the configuration
  //! (defuzzification operator, flag classif, observed output in the data set)
  //! 'display' is the file in which intermediate result are displayed.
  //! Returns -2 en case of inconsistency (not enough input variables or not active output),
  //! -1  when the observed output is not part of the sample file, the performance index otherwise:
  //! the number of misclassified for a crisp output with classif, or the root mean square of errors.
  //! The outfile file, 'fres', structure is as follow, a line corresponds to an example:
  //! Obs (if available) - Defuz (see below) - Err (if available) - Blank - Perf.
  //! 'Obs' is the observed output value, if it is available 'Err' indicates the difference between
  //! observed and inferred output values, 'Blank' is an integer whose value is 1 if the cumulated
  //! weight is less than 'MuSeuil', else 'Blank' is set to 0, the last colum is the cumulated error
  //!  or performance index.
  //! The fields filled by the defuzzification operator are:
  //! For a crisp output :
  //! 'sugeno' and classif= 'no' : the inferred value - Alarm
  //! 'sugeno' and classif= 'yes': the inferred value - the assigned class - Alarm
  //! 'MaxCrisp' : the inferred value - Alarm
  //! For a fuzzy output :
  //! 'sugeno' and classif= 'no' : the inferred value - Alarm
  //! 'sugeno' and classif= 'yes': the inferred value - Alarm -
  //! the membership degree of the inferred value to each of the output fuzzy sets.
  //! 'area' and classif= 'no' : the inferred value - Alarm
  //! 'area' and classif= 'yes': the inferred value  - Alarm -
  //!  the membership degree of the inferred value to each of the output fuzzy sets.
  //! 'MeanMax' and classif= 'no' : the inferred value - Alarm
  //! 'MeanMax' and classif= 'yes': the inferred value - Alarm  -
  //! the membership degree of the inferred value to each of the output fuzzy sets,
  //! if and only if the example activate at least one rule.

  virtual double Perf(int NumS, double ** Data, int NbEx, double & Coverage, double & MaxErr,
		      double MuThresh, int * MisClassified, double * Lab, int ref = true,
		      FILE * fres = NULL, FILE *display = NULL);

  //! Write a comment line in the performance function result file '*p'.
  //! 'o' is the output number.
  //! 'obs' is true when the observed output is part of the data.
  void WriteHeader(int o, FILE *p, int obs) const;

  //! Check the current system consistency and update the outputs arrays.
  //! Return values:  0 when the system is consistent, -100 when the number of inputs
  //! is different of the number of propositions in the premise part of the first rule
  //! -101 + 'number of the input' (the first one being 0) when the premise part
  //! of a rule uses an unknown input label,  -200 when the number of outputs is different of
  //! the number of propositions in the conclusion part of the first rule,
  //! -201 - 'number of the ouput' (the first one being 0) when the conclusion part
  //! of a rule uses an unknown output label.
  //! updates NbAct
  int CheckConsistency(void);

  //! Check system consistency and calls InitClassLabels() before calling Infer().
  double InferCheck(double * v, double ** Val = NULL, int nb = 0, int out_number=-1, FILE * fic=NULL, FILE *display=NULL);
  double InferCheck(MF ** v, double ** Val = NULL, int nb = 0, int out_number=-1, FILE * fic=NULL, FILE *display=NULL);
  int InferFatiCheck(MFDPOSS ** v, int out_number, int nalphacut=5, double ** Val = NULL, int nb=0, FILE * fic = NULL, FILE *display=NULL);
  //***
  //! Reallocate Classes array from data or from rule conclusions
  //! for all the outputs.
  //! 'Val' is the data array, of size 'nb',  used to initialize
  //! the class labels in case of classification. If (Val == NULL) the
  //! labels are initilaized from rule conclusions.
  void InitClassLabels(double **Val, int nb);

  //! Initialize the class labels of output 'num', if classification,
  //! 'Val' is an array of size 'nb' containig the labels.
  void SetClassLabels(int num, double *Val, int nb);

  //! rebuild classes array before calling Performace on data array
  //! Data array must be allocated before call to ClassCheck
  //! calls ClassifCheck and ResClassifAlloc
  //! which allocates MisClassified and Lab arrays
  int ClassCheck(int *& MisClassified, double  *& Lab, double ** Data, int NbEx, int numS);

  //! rebuild classes array before calling Performance on data array
  //! Data array must be allocated before call to ClassCheck
  //! does not allocate MisClassified array nor Lab array
  int ClassCheckNoAllocResClassif(double ** Data, int NbEx, int NumS);

  //! deletes and recreates Classes array
  //! only in classification crisp case
  //! does nothing otherwise
  //! does not allocate MisClassified array nor Lab array
  int ClassifCheck(double ** Data, int NbEx, int NumS) const;

  //! allocates resclassif and Lab arrays (called by ClassifCheck)
  int ResClassifAlloc(int *& MisClassified, double  *& Lab, int numS) const;

  //! Sort the rules of the current rule base according to their influence in a data file.
  //! 'dat' is the data array of size 'n' (number of samples).
  //! The order is according to the rule cumulated weight,
  //! decreasing if 'order' is any positive value, increasing if 'order' is any negative value.
  //! None sorting is done when order is zero.
  void SortRules(double **dat, int n, int order);

  //! Print all the possible combinations of the input breakpoints in a the file whose name is 'archive'.
  //! Call the function FISIN::GetBreakPoints.
  //! A break point delimites the region of influence for a rule (or a combination of rules).
  //! For a n-MF partition, at most n+n-1 break points are defined:  the middle of the kernel
  //! of each MF and the intersection, if exists, for two adjacent MFs, of the lines which
  //! join the support to the kernel limits. Considered neighbours for ith MF are MF i-1 and
  //! MF i+1:  there is no MF sorting done by this function.
  //! The return value is 0 when the file has been generated or the number of breakpoint
  //! combinations (the number of lines in the file) when it is greater than 'NbMax'
  //! (default=10000)
  int GetBreakPoints(char * archive, int NbMax=MAX_BP);

  //! Recursive fonction used by GetBreakPoints
  void GenereCombi(int i, FILE *f, int *NbBpG, int *cBpG, double ** BpG);

  //! Fonction used by GetBreakPoints
  void PrintBreakPoints(FILE *f, int *cBpG, double ** BpG) const;


  virtual void Print(FILE *f) const;

  //! Produces a configuration file for the current FIS
  virtual void PrintCfg(FILE *f, const char * fd = FORMAT_DOUBLE) const;

  //! pi: performance index, ci: coverage index, errmax: max of error, nout: output number
  //! Calls AnalyzeRB and print all the results in file 'f'.
  //! Return value : n out when nout greater than the number of outputs, 0 otherwise.
  int PerfRB(double pi, double ci, double errmax, int nout, FILE * f);

   //! Return value : n out when nout greater than the number of outputs, 0 otherwise.
 int WriteHeaderPerfRB(int nout, FILE * f);


  //! Analyzes the rule base and fills the structure InfoRB
  //! n is the output number (default:  0).
  //! Return value : n when n greater than the number of outputs, 0 otherwise.
  //! Calls the InitClassLabels function.
  int AnalyzeRB(InfoRB & i, int n = 0, double ** Varray = NULL, int nb = 0);

  // Access to private data
  int GetNbIn(void)     { return NbIn;  }  // Activates or not
  int GetNbOut(void)     { return NbOut; }
  int GetNbRule(void)      { return NbRules; }
  int GetNbExcept(void)      { return NbExceptions; }
  char *TypeConj(void)   const { return cConjunction; }
  char *MissingValues(void) const { return strMissingValues; }
  char *ErrorIndex(void) const { return strErrorIndex; }
  int GetNbActRule(void)     { return NbActRules; }

  //! recalculates the number of active rules and updates fis.NbActRules
  int ComputeNbActRule(void);

  virtual FIS* Clone() { return new FIS(*this); }

  //! Reduce the vocabulary used in the conclusions.
  //! NumS is the number of the Output, Data the array containing the data.
  //! MuSeuil is the activity Threshold for an example not to be considered blank.
  //! The reduction can be made in three different ways:
  //! - Classif Output : Conclusions are set to the nearest class label.
  //! - NConc != 0 : The monodimensional k-means is made one time, with NConc centers,
  //! Conclusions are set to their nearest center.
  //! - NConc = 0, regression type : Loop incrementing the number of centers resulting
  //! from the monodimensional k-means.
  //! The loop end when the number of centers equal the number of samples, or when the percentage
  //! of Performance Loss between the original system and the reduced system is less than PerfLoss.
  //! return the number of new item in the vocabulary.
  //! NB : except for the classif case, resulting conclusions are adjusted so that upper
  //! and lower bounds of output variation domain are part of them.
  //! ExtVoc=true: centres are built from rule conclusion values else from learning output data.
  //! return the new number of items in the vocabulary.
  int VocReduc(int NumS, double **Data, int NbEx, double MuSeuil, double PerfLoss, int NConc, int ExtVoc);

  //! replace Old Conclusions of the Conc array by the nearest one of the centres array
  //! Conc is the array of the old conclusions, centres the array of the new ones
  //! Conc is the number of new conclusions, or the size of the centres array
  void NewConc(double *&Conc, double *Centres, int nconc);

  //! Include the performance and the ComputeWeightedPerf functions, return the weighted performance.
  //! NPart is the number of parts in the local performance (used only in a crisp regression case,
  //! automatically computed in a classif and fuzzy output case).
  //! DomBreakPoints is a character chain containing the new frontier points
  //! (used only in a crisp regression output case, chain format : [ , , , ] with NPart-1 values).
  //! PartWeight is a character chain containing the weight of each part
  //!(equal to NPart, number of classes or output fuzzy subsets, chain format : [ , , , ]) .
  //! WeightedCov is the WeightedCoverage
  //! MaxError is the maximum error encounted
  //! MuSeuil is the activation threshold
  //! fres is the name of the resulting file name
  double WeightedPerf(int NumS, char * fdata, int NPart, char *DomBreakpoints, char *PartWeight,  double &WeightedCov, double &MaxError, double MuSeuil, int ErrorType, char *fres, FILE *display);

  //! Given an Array containing the results and the coverage index each part,
  //! compute and return the wieghted performance.
  //! PartWeight is the character chain containing the Wieghts of each part.
  //! Resultab and Couvert are the array containing respectively the performance
  //! index and the coverage index of each part (last item : global index).
  //! WieghtedCov is the WeightedCoverage.
  double ComputeWeightedPerf(char *PartWeight, int NPart, double *&ResultTab, double *&Couvert, double &WeightedCov);

  //! Compute the performance and coverage index of each specified part
  //! fData is the file name containing the... data
  //! NPart and DomBreakPoints are only used in a crisp regression case
  //! ResultTab contains the performance indexes of each part
  //! Couvert contains the coverage indexes of each part
  //! MaxErr contains the maximum error encounted in each part
  //! NItem contains the number of samples for each part
  //! Museuil is the activation threshold for a sample to be active
  //! ErrorType is the computed error type
  //! (1 : Fispro error, 2 : MSE, 3 RMSE, 4 MAE (mean error), 5 RMAE (relative mean error)).
  //! fres and affiche are the same as for the global perf
  //! for Each Tab, the last Item is the global index
  int Performance(int NumS, char * fdata, int NPart, char *DomBreakPoints, double *&ResultTab, double *&Couvert, double *&MaxErr, double *&NItem, double MuSeuil, int ErrorType, char * fres, FILE *display);

  //! Init the weights in the Weights array and normalize them
  //! NPart is the number of part
  //! PartWeights is the character chain containing the weights
  //! If PartWeights is NULL, weights are equireparted and normalized.
  void InitWeights(int NPart, char *PartWeights, double *&Weights);

  //! Init the frontier points in the BreakPoints array
  //! NPart is the number of part
  //! DomBreakPoints is the character chain containing the weights
  //! If DomBreakPoints is NULL, each part is of eaqul size
  //! This fucntion is used in a crisp regression case only
  void InitBreakPoints(int NumS, int Npart, char *DomBreakPoints, double *&BreakPoints);

  //! Same as Performance, but Lab and Misclassified arrays must be initialized,
  //! and Data is the data array, not a file name.
  //! NbEx is the number of samples
  int Perf(int NumS, double ** Data, int NbEx, int NPart, double *&ResultTab,  double *&Couvert, double *&MaxErr, double *&NItem, double MuSeuil, int ErrorType, double *BreakPoints, int * MisClassified, double * Lab, int ref, FILE * f, FILE *display);

  // FATI Inference for implicative rules and fuzzy inputs
  void InferFatiPrep(int nout);
  void KinkPoints(std::list<double> ** dl, int nout);
  void UpdatePartList(int iout, std::list<double> **dL, double mposs, int m1, int m2);

  void BuildFuzIn(double * val, MFDPOSS ** tpl, MFDPOSS ** fuzval);
  void BuildFuzIn(double * crispIn, double * KW, double * SW, MFDPOSS ** &fuzval, double maxposs=1.0);
  MFDPOSS * InferAcut(double *binf, double *bsup, int iout, FILE *fg, double alpha, FILE *display);
  MFDPOSS * InferFatiAlpha(MFDPOSS ** v, int a, int nout, FILE * f, FILE *display);
  MFDPOSS * InferFati(MFDPOSS ** v, int nalf, int nout, FILE * fic, FILE *display);

  double** dist(char * datafile, int &nbex, char*dissoutfile, double *numc,bool norm,double mink,bool display,bool wordless);
  //!computes distances of data elements in datafile, based on Fuzzy Partitions in FIS
  //! writes distance matrix in dissoutfile
  //! numc, for multidimensional distance, integer array, for each position
  //! warning: numc must have a length equal to NbIn (active or not)
  //! 0 = Euclidean distance, 1 = fuzzy partition-based distance in the corresponding dimension
   //! norm=true, normalized distance (/f-1), false, unnormalized distance
  //! mink=Minkowski power
  //! returns the distance matrix
  //! 
  double ** distWithNormedData(double **normdata,int nbex,int nbcol,char*dissoutfile, double* numc,bool norm,double mink,bool display,bool wordless);
    //!like dist function, but works with an array or normalized data
  double ** NormCheckDataDist(char * datafile,int &nbex,int &nbcol,bool display,bool wordless);

}; // End of FIS class


#endif



//*********************      End of FIS.H      ***********************




