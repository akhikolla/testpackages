//***********************************************************************

//
//
//                              MF.H
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : MF class functions used by FISPRO, part of library fispro

//**********************************************************************
//***********************************************************************
//*** CLASS: MFDPOSS                                                  ***
//*** AUTHOR: Hazaël JONES (jones@ensam.inra.fr)                      ***
//*** Release Date: August 22, 2005                                   ***
//***********************************************************************


#ifndef __MF_H
#define __MF_H

#if defined(__MINGW32__) && defined(__STRICT_ANSI__) || defined(_MSC_VER)
#define M_PI		3.14159265358979323846
#endif

#include <list>
#include "fis.h"
#include <assert.h>
#include "pt.h"


class MFUNIV:public MF
//********************
{
  protected :

    double a, b;

 public:

  MFUNIV() : MF(), a(0), b(0) {}

  MFUNIV(const MFUNIV &s) : MF(s)
    {
      a = s.a;
      b = s.b;
    }

  MFUNIV(const MFUNIV *s) : MF(*s)
    {
      a = s->a;
      b = s->b;
    }

  MFUNIV(double left, double right) : MF()
    {
      a = left; b = right;

      if(DBL_SUP_EQUAL(a, b))  // a >= b
	throw std::runtime_error( "~S2~MustBeHigherThan~S1~" );
    }

  MF* Clone() const { return new MFUNIV(*this); }

  void operator=(MFUNIV  &T)
  {
    if(this == &T) return;
    a = T.a; b = T.b;
    SetName(T.Name);
  }

  const char *GetType() const { return Type(); }
  static const char *Type() { return "universal"; }

  void ValParam(double &aa, double &bb) const
    { aa = a;  bb = b; }

  void Update(double *Values)
    {
      a = Values[0];
      b = Values[1];
    }

  int NbParams() const { return 2; };
  void GetParams( double *params ) const { ValParam( params[0], params[1] ); }

  double GetDeg(double value) const
    {
      return 1.;
    }

  double Kernel(double & left, double & right) const
    { left = a; right = b; return left + (right - left) / 2.;}
  double Support(double & left, double & right) const
    { left = a; right = b;  return left + (right - left) / 2 ;}
  double AlphaKernel(double & left, double & right, double alpha) const
    {
      return Kernel(left, right);
    }

  void Print(FILE * f) const
    {
      MF::Print(f);
      fprintf(f, "%8.3f%c%8.3f", a, SEPARE, b);
    }

  //! Print input configuration.
  //! Function used to store a FIS configuration
  virtual void PrintCfg(int num, FILE *f, const char * fd = FORMAT_DOUBLE) const
    {
      MF::PrintCfg(num, f, fd);
      //    fprintf(f, "%c%8.3f%c%8.3f%c\n", START_NB, a, SEPARE, b, END_NB);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, a);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, b);
      fprintf(f, "%c\n", END_NB);
    }

  void SetInf( double inf ) { a = inf; }
  void SetSup( double sup ) { b = sup; }
  virtual MFPWLinear getMFPWLinear() const;

};


class MFTRI:public MF
//*******************
{
  protected :

    double a, b, c;

 public:

  MFTRI() : MF(), a(0), b(0), c(0) {}

  MFTRI(const MFTRI &s) : MF(s)
    {
      a = s.a;
      b = s.b;
      c = s.c;
    }

  MFTRI(const MFTRI *s) : MF(*s)
    {
      a = s->a;
      b = s->b;
      c = s->c;
    }

  MFTRI(double centre, double delta) : MF()  // Triangle symetrique
    {
      if( delta < EPSILON )
	throw std::runtime_error( "~ValueMustBePositive~" );

      b = centre;
      a = centre - delta;
      c = centre + delta;
    }

  MFTRI(double aa, double bb, double cc) : MF()
    {
      a = aa; b = bb; c = cc;

      if(DBL_SUP(a, b)) // a >= b   DBL_SUP _EQUAL
	throw std::runtime_error( "~S2~MustBeHigherThan~S1~" );
      if(DBL_SUP_EQUAL(a, c))
	throw std::runtime_error( "~S3~MustBeHigherThan~S1~" );
      if(DBL_SUP(b, c))   // DBL_SUP_EQUAL
	throw std::runtime_error( "~S3~MustBeHigherThan~S2~" );
    }

  MF* Clone() const { return new MFTRI(*this); }

  //! Function used by the optimization module
  void Normalize ( double lowerb , double upperb )
    {
      double ratio;
      ratio=upperb - lowerb;
      if (fabs(ratio) >EPSILON)
	{
	  a = (a - lowerb)/ratio;
	  b = (b - lowerb)/ratio;
	  c = (c - lowerb)/ratio;
	}

    }

  //! Function used by the optimization module
  void UnNormalize( double lowerb , double upperb )
    {
      double ratio;
      ratio=upperb - lowerb;
      a = a*ratio + lowerb;
      b = b*ratio + lowerb;
      c = c*ratio + lowerb;
    }

  void operator=(MFTRI  &T)
  {
    if(this == &T) return;
    a = T.a; b = T.b; c = T.c;
    SetName(T.Name);
  }

  const char *GetType() const { return Type(); }
  static const char *Type() { return "triangular"; }

  void ValParam(double &aa, double &bb, double &cc) const
    { aa = a;  bb = b; cc = c; }

  void Update(double *Values)
    {
      a = Values[0];
      b = Values[1];
      c = Values[2];
    }

  int NbParams() const { return 3; };
  void GetParams( double *params ) const { ValParam( params[0], params[1], params[2] ); }

  double GetDeg(double value) const
    {
      if ( (value < a) || (value > c) ) return  0.0;
      // modif May 15 2009 - case of perfect right-angled triangles: a=b or b=c
      if (value==b) return 1.0;
      if (value <b) return ( (value - a) / (b - a) );
      return ( (c - value) / (c - b) );
    }

  double Kernel(double & left, double & right) const
    { left = right = b ; return left;}
  double Support(double & left, double & right) const
    { left = a; right = c;  return left + (right - left) / 2 ;}
  double AlphaKernel(double & left, double & right, double alpha) const
    {
      left = alpha * b + (1.0 - alpha) * a ;
      right = alpha * b + (1.0 - alpha) * c ;
      if(right == left) return right;
      return left + (right - left) / 2 ;
    }

  void Print(FILE * f) const
    {
      MF::Print(f);
      fprintf(f, "%8.3f%c%8.3f%c%8.3f", a, SEPARE, b, SEPARE, c);
    }

  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
    {
      MF::PrintCfg(num, f, fd);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, a);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, b);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, c);
      fprintf(f, "%c\n", END_NB);
    }

    double GetLeftx(double degre)
    {
    	return a*(1-degre)+b*degre;
    }

    double GetRightx(double degre)
    {
    	return c*(1-degre) + b * degre;
    }

    virtual MFPWLinear getMFPWLinear() const;

};



class MFTRAP:public MF
//********************
{
  protected :

    double a, b, c, d;

 public:

  MFTRAP() : MF(), a(0), b(0), c(0), d(0) {}
  MFTRAP(double aa, double bb, double cc, double dd) : MF()
    {
      a = aa; b = bb; c = cc; d = dd;

       if(DBL_SUP(a, b)) // a >= b   DBL_SUP _EQUAL
	throw std::runtime_error( "~S2~MustBeHigherThan~S1~" );
      if(DBL_SUP(b, c)) // b > c
	throw std::runtime_error( "~S3~MustBeHigherThan~S2~" );
      if( a - d > EPSILON )
	throw std::runtime_error( "~S4~MustBeHigherThan~S1~" );
      if(DBL_SUP_EQUAL(b, d))
	throw std::runtime_error( "~S4~MustBeHigherThan~S2~" );
      if(DBL_SUP(c, d))   // DBL_SUP_EQUAL
	throw std::runtime_error( "~S4~MustBeHigherThan~S3~" );
    }

  MFTRAP(const MFTRAP &s) : MF(s)
    {
      a = s.a;
      b = s.b;
      c = s.c;
      d = s.d;
    }

  MF* Clone() const { return new MFTRAP(*this); }

  //! Function used by the optimization module
  void Normalize ( double lowerb , double upperb )
    {
      double ratio;
      ratio=upperb - lowerb;
      if( ratio >EPSILON)
	{
	  a = (a - lowerb)/ratio;
	  b = (b - lowerb)/ratio;
	  c = (c - lowerb)/ratio;
	  d = (d - lowerb)/ratio;
	}
    }

  //! Function used by the optimization module
  void UnNormalize( double lowerb , double upperb )
    {
      double ratio;
      ratio=upperb - lowerb;
      if( ratio >EPSILON)
	{
	  a = a*ratio + lowerb;
	  b = b*ratio + lowerb;
	  c = c*ratio + lowerb;
	  d = d*ratio + lowerb;
	}
    }

  const char *GetType() const { return Type(); }
  static const char *Type() { return "trapezoidal"; }

  void ValParam(double &aa, double &bb, double &cc, double &dd) const
    { aa = a;  bb = b; cc = c; dd = d; }

  void Update(double *Values)
    {
      a = Values[0];
      b = Values[1];
      c = Values[2];
      d = Values[3];
    }

  int NbParams() const { return 4; };
  void GetParams( double *params ) const { ValParam( params[0], params[1], params[2] , params[3] ); }

  double GetDeg(double value) const
    {
      if ( (value < a) || (value > d) ) return  0.0;
      // modif May 15 2009 - case of perfect right-angled trapezoids : a=b or c=d
      if (value==b || value==c) return 1.0;
      if (value <b) return ( (value - a) / (b - a) );
      if (value <c) return 1.0;
      return ( (d - value) / (d - c) );
    }

  double Kernel(double & left, double & right) const
    {
      left = b; right = c ;
      if(right == left) return right;
      return left + (right - left) / 2 ;
    }
  double Support(double & left, double & right) const
    { left = a; right = d;  return left + (right - left) / 2 ; }

  double AlphaKernel(double & left, double & right, double alpha) const
    {
      left = alpha * b + (1.0 - alpha) * a ;
      right = alpha * c + (1.0 - alpha) * d ;
      if(right == left) return right;
      return left + (right - left) / 2 ;
    }

  void Print(FILE * f) const
    {
      MF::Print(f);
      fprintf(f, "%8.3f%c%8.3f%c%8.3f%c%8.3f", a, SEPARE, b, SEPARE, c, SEPARE, d);
    }

  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
    {
      MF::PrintCfg(num, f, fd);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, a);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, b);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, c);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, d);
      fprintf(f, "%c\n", END_NB);
    }

    double GetLeftx(double degre)
    {
    	return a*(1-degre)+b*degre;
    }

    double GetRightx(double degre)
    {
    	return d*(1-degre) + c * degre;
    }

    virtual MFPWLinear getMFPWLinear() const;

};


class MFTRAPINF:public MF
//***********************
{
  double a, b, c;

 public:

  MFTRAPINF(double aa, double bb, double cc) : MF()
    {
      a = aa; b = bb; c = cc;

      if(DBL_SUP(a, b)) // a >= b
	throw std::runtime_error( "~S2~MustBeHigherThan~S1~" );
      if(DBL_SUP_EQUAL(b, c)) //
	throw std::runtime_error( "~S3~MustBeHigherThan~S2~" );
    }

  MFTRAPINF(const MFTRAPINF &s) : MF(s)
    {
      a = s.a;
      b = s.b;
      c = s.c;
    }

  MF* Clone() const { return new MFTRAPINF(*this); }

  //! Function used by the optimization module
  void Normalize ( double lowerb , double upperb )
    {
      double ratio;
      ratio=upperb - lowerb;

      if (fabs(ratio) >EPSILON)
	{
	  a = (a - lowerb)/ratio;
	  b = (b - lowerb)/ratio;
	  c = (c - lowerb)/ratio;
	}
    }

  //! Function used by the optimization module
  void UnNormalize( double lowerb , double upperb )
    {
      double ratio;
      ratio=upperb - lowerb;
      if (fabs(ratio) >EPSILON)
	{
	  a = a*ratio + lowerb;
	  b = b*ratio + lowerb;
	  c = c*ratio + lowerb;
	}
    }

  const char *GetType() const { return Type(); }
  static const char *Type() { return "SemiTrapezoidalInf"; }

  void ValParam(double &aa, double &bb, double &cc) const
    { aa = a;  bb = b; cc = c; }

  void Update(double *Values)
    {
      a = Values[0];
      b = Values[1];
      c = Values[2];
    }

  int NbParams() const { return 3; };
  void GetParams( double *params ) const { ValParam( params[0], params[1], params[2] ); }
  void SetInf( double inf ) { a = inf; }

  double GetDeg(double value) const
    {
      if (value <= b)  return  1.0;
      else if (value >=c) return 0.0;
      return ( (c - value) / (c - b) );
    }

  double Kernel(double & left, double & right) const
    {
      left = a;
      right = b ;
      if(right == left) return right;
      return  left + (right - left) / 2;
    }

  double Corner() const { return b; }

  double Support(double & left, double & right) const
    { left = a; right = c;  return left + (right - left) / 2 ; }

  double AlphaKernel(double & left, double & right, double alpha) const
    {
      left = a;
      right = alpha * b + (1.0 - alpha) * c ;
      if(right == left) return right;
      return left + (right - left) / 2 ;
    }

  void Print(FILE * f) const
    {
      MF::Print(f);
      fprintf(f, "%8.3f%c%8.3f%c%8.3f", a, SEPARE, b, SEPARE, c);
    }

  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
    {
      MF::PrintCfg(num, f, fd);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, a);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, b);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, c);
      fprintf(f, "%c\n", END_NB);
    }
  virtual MFPWLinear getMFPWLinear() const;

};


class MFTRAPSUP:public MF
//***********************
{
  double a, b, c;

 public:

  MFTRAPSUP(double aa, double bb, double cc) : MF()
    {
      a = aa; b = bb; c = cc;

      if(DBL_INF_EQUAL(b, a)) // b <= a
	throw std::runtime_error( "~S2~MustBeHigherThan~S1~" );
      if(DBL_INF(c, b)) // c<b
	if( b -c > EPSILON2 )
	  throw std::runtime_error( "~S3~MustBeHigherThan~S2~" );
    }

  MFTRAPSUP(const MFTRAPSUP &s) : MF(s)
    {
      a = s.a;
      b = s.b;
      c = s.c;
    }

  MF* Clone() const { return new MFTRAPSUP(*this); }

  //! Function used by the optimization module
  void Normalize ( double lowerb , double upperb )
    {
      double ratio;
      ratio= upperb - lowerb;
      if (fabs(ratio) >EPSILON)
	{
	  a = (a - lowerb)/ratio;
	  b = (b - lowerb)/ratio;
	  c = (c - lowerb)/ratio;
	}
    }

  //! Function used by the optimization module
  void UnNormalize( double lowerb , double upperb )
    {
      double ratio;
      ratio= upperb - lowerb;
      if (fabs(ratio) >EPSILON)
	{
	  a = a*ratio + lowerb;
	  b = b*ratio + lowerb;
	  c = c*ratio + lowerb;
	}
    }

  const char *GetType() const { return Type(); }
  static const char *Type() { return "SemiTrapezoidalSup"; }

  void ValParam(double &aa, double &bb, double &cc) const
    { aa = a;  bb = b; cc = c;}

  void Update(double *Values)
    {
      a = Values[0];
      b = Values[1];
      c = Values[2];
    }

  int NbParams() const { return 3; };
  void GetParams( double *params ) const { ValParam( params[0], params[1], params[2] ); }
  void SetSup( double sup ) { c = sup; }

  double GetDeg(double value) const
    {
      if (value <= a)  return  0.0;
      else if (value >= b) return 1.0;
      return ( (value -  a) / (b - a) );
    }

  double Kernel(double & left, double & right) const
    {
      left = b;
      right = c;
      if(right == left) return right;
      return left + (right - left) / 2 ;
    }

  double Corner() const { return b; }

  double Support(double & left, double & right) const
    { left = a; right = c; return left + (right - left) / 2 ; }

  double AlphaKernel(double & left, double & right, double alpha) const
    {
      left = alpha * b + (1.0 - alpha) * a ;
      right = c;
      if(right == left) return right;
      return left + (right - left) / 2 ;
    }

  void Print(FILE * f) const
    {
      MF::Print(f);
      fprintf(f, "%8.3f%c%8.3f%c%8.3f", a, SEPARE, b, SEPARE, c);
    }

  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
    {
      MF::PrintCfg(num, f, fd);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, a);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, b);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, c);
      fprintf(f, "%c\n", END_NB);
    }

  virtual MFPWLinear getMFPWLinear() const;

};

class MFGBELL:public MF
//*********************
{
  double a, b, c;

 public:

  MFGBELL(double aa, double bb, double cc) : MF()
    {
      a  = aa;  b = bb; c = cc;
    }

  MFGBELL(const MFGBELL &s) : MF(s)
    {
      a = s.a;
      b = s.b;
      c = s.c;
    }

  MF* Clone() const{ return new MFGBELL(*this); }
  const char *GetType() const { return Type(); }
  static const char *Type() { return "gbell"; }

  void ValParam(double &aa, double &bb, double &cc) const
    { aa = a;  bb = b; cc = c; }

  void Update(double *Values)
    {
      a = Values[0];
      b = Values[1];
      c = Values[2];
    }

  int NbParams() const { return 3; };
  void GetParams( double *params ) const { ValParam( params[0], params[1], params[2] ); }

  double GetDeg(double value) const
    {
      return (1 / pow( (1 + fabs( (value - c)/a) ), b*2 ) );
    }

  double Kernel(double & left, double & right) const
    { left = c; right = c; return left; }
  double Support(double & left, double & right) const
    { left = c - 3*a; right = c + 3*a; return left + (right - left) / 2 ; }

  double AlphaKernel(double & left, double & right, double alpha) const
    {
      double ec = a * (exp(log(alpha) / (-2.0 * b)) - 1.0);
      left = c - ec ;
      right =  c + ec ;
      if(right == left) return right;
      return left + (right - left) / 2 ;
    }

  void Print(FILE * f) const
    {
      MF::Print(f);
      fprintf(f, "%8.3f%c%8.3f%c%8.3f", a, SEPARE, b, SEPARE, c);
    }

  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
    {
      MF::PrintCfg(num, f, fd);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, a);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, b);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, c);
      fprintf(f, "%c\n", END_NB);
    }

};



class MFGAUSS:public MF
//*********************
{
  double mean, std;

 public:

  MFGAUSS(double a, double b) : MF()
    {
      if(a <= 0)
	throw std::runtime_error( "~StandardDeviation~MustBePositive~" );
      std = a;  mean = b; // for compatibility with MATLAB
    }

  MFGAUSS(const MFGAUSS &s) : MF(s)
    {
      std = s.std;
      mean = s.mean;
    }

  MF* Clone() const { return new MFGAUSS(*this); }
  const char *GetType() const { return Type(); }
  static const char *Type() { return "gaussian"; }

  void ValParam(double &a, double &b) const
    { a = mean;  b = std; }

  void Update(double *Values)
    {
      mean = Values[0];
      std = Values[1];
    }

  int NbParams() const { return 2; };
  void GetParams( double *params ) const { ValParam( params[0], params[1] ); }

  double GetDeg(double value) const
    {
      return (exp(-1.0 * (value - mean)* (value - mean) / (2.0 * std * std) ));
    }

  double Kernel(double & left, double & right) const
    { left = right = mean; return left; }
  double Support(double & left, double & right) const
    {
      left = mean - 3*std;
      right = mean + 3*std;
      return left + (right - left) / 2 ;
    }

  double AlphaKernel(double & left, double & right, double alpha) const
    {
      double ec = sqrt(log(alpha) * -2.0 * std * std);
      left = mean - ec ;
      right =  mean + ec ;
      if(right == left) return right;
      return left + (right - left) / 2 ;
    }

  void Print(FILE * f) const
    {
      MF::Print(f);
      fprintf(f, "%8.3f%c%8.3f", std, SEPARE, mean);
    }

  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
    {
      MF::PrintCfg(num, f, fd);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, std);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, mean);
      fprintf(f, "%c\n", END_NB);
    }

};


class MFDISCRETE:public MF
//************************
{
  double *MfValues;
  int Nb;
  double min, max;

 public:

  MFDISCRETE(double v) : MF()
    {
      Nb = 1;
      MfValues = new double [1];
      MfValues[0] = v;
      min = max = v;
    }

  MFDISCRETE(int n, double * tab) : MF()
    {
      int i;

      Nb = n;
      MfValues = new double[Nb];

      min = max = tab[0];
      for(i = 0; i < Nb; i ++)
	{
	  MfValues[i] = tab[i];
	  if(MfValues[i] < min) min = MfValues[i];
	  if(MfValues[i] > max) max = MfValues[i];
	}
    }

  MFDISCRETE(const MFDISCRETE &s) : MF(s)
    {
      Nb = s.Nb;
      min = s.min;
      max = s.max;
      MfValues = new double[Nb];
      for( int i=0 ; i<Nb ; i++ )
	MfValues[i] = s.MfValues[i];
    }

  virtual ~MFDISCRETE() { delete [] MfValues; }

  MF* Clone() const { return new MFDISCRETE(*this); }
  const char *GetType() const { return Type(); }
  static const char *Type() { return "discrete"; }

  int NbParams() const { return Nb; };
  void GetParams( double *params ) const { ValParam( params ); }

  void ValParam(double *t) const
    {
      for( int i=0; i<Nb; i++ )      t[i] = MfValues[i];
    }

  void Update(double *Values)
    {
      for( int i=0; i<Nb; i++ )  MfValues[i]=Values[i];
    }

  double GetDeg(double value) const
    {
      for( int i=0; i<Nb; i++ )
	{
	  if(FisIsnan(value)) return value;   //  NaN
	  if( value == MfValues[i])   return 1.0;
	}
      return 0.0;
    }

  double Kernel(double & left, double & right) const
    {
      left = min;
      right = max;
      if(right == left) return right;
      return left + (right - left) / 2 ;
    }

  double Support(double & left, double & right) const
    {
      return Kernel(left, right);
    }

  double AlphaKernel(double & left, double & right, double alpha) const
    {
      left = alpha;
      return Kernel(left, right);
    }
  //alpha-cut is meaningless for this MF type
  // the first instruction is for eliminating a compilation warning

  void Print(FILE * f) const
    {
      int i;
      MF::Print(f);
      for(i = 0; i < Nb; i ++)    fprintf(f, "%8.3f%c", MfValues[i], SEPARE) ;
    }

  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
    {
      int i;

      MF::PrintCfg(num, f, fd);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, MfValues[0]);
      for(i = 1; i < Nb; i ++)
	{
	  fprintf(f, "%c", SEPARE);
	  fprintf(f, fd, MfValues[i]);
	}
      fprintf(f, "%c\n", END_NB);
    }

};

class MFDOOR:public MF
//********************
{
  protected :

    double a, b;

 public:

  MFDOOR() : MF(), a(0), b(0) {}

  MFDOOR(const MFDOOR &s) : MF(s)
    {
      a = s.a;
      b = s.b;
    }

  MFDOOR(const MFDOOR *s) : MF(*s)
    {
      a = s->a;
      b = s->b;
    }

  MFDOOR(double left, double right) : MF()
    {
      a = left; b = right;

      if(DBL_SUP(a, b))  // a > b
	throw std::runtime_error( "~S2~MustBeHigherThan~S1~" );
    }

  MF* Clone() const { return new MFDOOR(*this); }

  void operator=(MFDOOR  &T)
  {
    if(this == &T) return;
    a = T.a; b = T.b;
    SetName(T.Name);
  }

  const char *GetType() const { return Type(); }
  static const char *Type() { return "door"; }

  void ValParam(double &aa, double &bb) const
    { aa = a;  bb = b; }

  void Update(double *Values)
    {
      a = Values[0];
      b = Values[1];
    }

  int NbParams() const { return 2; };
  void GetParams( double *params ) const { ValParam( params[0], params[1] ); }

  double GetDeg(double value) const
    {
      if((value < a) || (value > b)) return 0.0;
      return 1.;
    }

  double Kernel(double & left, double & right) const
    { left = a; right = b; return left + (right - left) / 2.;}
  double Support(double & left, double & right) const
    { left = a; right = b;  return left + (right - left) / 2 ;}
  double AlphaKernel(double & left, double & right, double alpha) const
    {
      return Kernel(left, right);
    }

  void Print(FILE * f) const
    {
      MF::Print(f);
      fprintf(f, "%8.3f%c%8.3f", a, SEPARE, b);
    }

  //! Print input configuration.
  //! Function used to store a FIS configuration
  virtual void PrintCfg(int num, FILE *f, const char * fd = FORMAT_DOUBLE) const
    {
      MF::PrintCfg(num, f, fd);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, a);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, b);
      fprintf(f, "%c\n", END_NB);
    }

  void SetInf( double inf ) { a = inf; }
  void SetSup( double sup ) { b = sup; }

    double GetLeftx(double degre)
    {
    	return a;
    }

    double GetRightx(double degre)
    {
    	return b;
    }

    virtual MFPWLinear getMFPWLinear() const;

};



/*
  A class representing the piecewise linear membership functions with inlinable
  Degree function
*/

#include <limits>

class MFPWLinear
{
  double a,b,c,d;
  double onba, ondc;

public:
  static double inf()
  {
//    const long long v=0x7ff0000000000000LL;
//    return *(double*)&v;
//IEEE inf appears to foul up g++'s optimiser!!!
    return std::numeric_limits<double>::max();
  }

  std::vector<double> getValues() const {
    std::vector<double> r(4); r[0]=a; r[1]=b; r[2]=c; r[3]=d; return r;
  }
  double GetDeg(double value) const {
    if (value <=a || value>=d)
      return 0;
    else if (value>=b && value <=c)
      return 1;
    else if (value<b)
      return (value-a)*onba;
    else
      return (d-value)*ondc;
  }
  MFPWLinear(const MFUNIV& x) {a=b=-inf(); c=d=inf();}
  MFPWLinear(const MFTRI& x) {
    x.ValParam(a,b,d); c=b; onba=1/(b-a); ondc=1/(d-c);
    assert(validParms());
  }
  MFPWLinear(const MFTRAP& x) {
    x.ValParam(a,b,c,d); onba=1/(b-a); ondc=1/(d-c);
    assert(validParms());
  }
  MFPWLinear(const MFTRAPINF& x) {
    x.ValParam(b,c,d); a=b=-inf(); ondc=1/(d-c);
    assert(validParms());
  }
  MFPWLinear(const MFTRAPSUP& x) {
    x.ValParam(a,b,c); c=d=inf(); onba=1/(b-a);
    assert(validParms());
  }
  MFPWLinear(const MFDOOR& x) : onba(0), ondc(0) {
    x.ValParam(b,c); a=b; d=c;
    assert(validParms());
  }
  MFPWLinear(const MF& x) : onba(0), ondc(0) {a=b=c=d=0;}  /* a==d => invalid PWlinear MF */
  bool valid() const {return a!=d;}
  bool validParms() const {
    if (a>b)
      printf("in validParms a=%g, b=%g\n",a,b);
    //    if (b>c)
    if(DBL_SUP(b, c))
      printf("in validParms b=%g, c=%g\n",b,c);
    if(c>d)
      printf("in validParms c=%g, d=%g\n",c,d);
    return a<=b && DBL_INF_EQUAL(b, c) && c<=d;}
};

inline MFPWLinear MF::getMFPWLinear() const {return MFPWLinear(*this);}
inline MFPWLinear MFUNIV::getMFPWLinear() const {return MFPWLinear(*this);}
inline MFPWLinear MFTRI::getMFPWLinear() const {return MFPWLinear(*this);}
inline MFPWLinear MFTRAP::getMFPWLinear() const {return MFPWLinear(*this);}
inline MFPWLinear MFTRAPINF::getMFPWLinear() const {return MFPWLinear(*this);}
inline MFPWLinear MFTRAPSUP::getMFPWLinear() const {return MFPWLinear(*this);}
inline MFPWLinear MFDOOR::getMFPWLinear() const {return MFPWLinear(*this);}


class MFDPOSS:public MF
//********************
{
 protected :
  LIST* pL;
  double maxposs;

 public:

 MFDPOSS() : MF() {pL = new LIST;maxposs=0.0;}

 MFDPOSS(const MFDPOSS &T) : MF(T)
   {
     if(this == &T) return;

     pL = new LIST;
     // check if T not empty
     if(T.NbParams() >0)
       {
	 T.pL->home();
	 pL->home();
	 maxposs = T.maxposs;
	 pL->add(T.pL->Get()->x, T.pL->Get()->y);

	 while( !T.pL->IsEnd() )
	   {
	     pL->next();
	     T.pL->next();
	     pL->add(T.pL->Get()->x, T.pL->Get()->y);
	   }

	 SetName(T.Name);
       }
   }

 MFDPOSS(double v) : MF()
    {
      pL = createList(v, v, v, v, 0., 1.);
      maxposs=1.0;
    }

 MFDPOSS(ACUT *d) : MF()
    {
      maxposs =d->alpha;
      pL = createList(d->l, d->r, d->l, d->r, 0., maxposs);
    }

 MFDPOSS(ACUT &d) : MF()
    {
      maxposs = d.alpha;
      pL = createList(d.l, d.r, d.l, d.r, 0., maxposs);
    }

 MFDPOSS(MF &mfi, double degree) : MF()
   {
	 double lk=0.0, rk=0.0;
	 double ls=0.0, rs=0.0;
	 mfi.Support(ls,rs);
	 mfi.Kernel(lk,rk);
	 maxposs = 1.0;
	 pL = createList(ls, rs, lk, rk, degree, maxposs);
   }

 MFDPOSS(MF *mfi, double degree) : MF()
   {
     if (mfi != NULL)
       {
	 double lk=0.0, rk=0.0;
	 double ls=0.0, rs=0.0;
	 mfi->Support(ls,rs);
	 mfi->Kernel(lk,rk);
	 maxposs = 1.0;
	 pL = createList(ls, rs, lk, rk, degree, maxposs);

       }
     else
       pL = new LIST;
   }

 MFDPOSS(LIST *tmplist) : MF()
   {
     pL = new LIST;
     if(tmplist->GetSize() >0)
       {
	 pL->home();
	 tmplist->home();

	 pL->add(tmplist->Get()->x, tmplist->Get()->y);
	 maxposs = tmplist->Get()->y;
	 while( !tmplist->IsEnd() )
	   {
	     pL->next();
	     tmplist->next();
	     pL->add(tmplist->Get()->x, tmplist->Get()->y);
	     if(tmplist->Get()->y>maxposs) maxposs = tmplist->Get()->y;
	   }
	 Simplify();
       }
   }

     MFDPOSS(const char * fic): MF()
    {
      FILE *f;
      char string[30];
      double x,y, oldx, oldy;

      f = fopen(fic, "rt");
      if(f == NULL)
	{
	  sprintf( ErrorMsg, "~CannotOpenFISFile~: %.100s~", fic);
	  throw std::runtime_error( ErrorMsg );
	}
      if(fgets(string, 30, f) == NULL) // check the first line
	{
	  sprintf( ErrorMsg, "~FirstLineEmptyInFile~: %.100s~", fic);
	  throw std::runtime_error( ErrorMsg );
	}

      pL = new LIST;
      pL->home();

      sscanf( string, "%lg  %lg",&x,&y);
      oldx = x;
      oldy = y;
      maxposs = y;

      if (y > EPSILON) // check that first y = 0
	{
	  printf(" WARNING added start value, %8.3f %8.3f \n", x,0.0);
	  pL->add(x,0);
	}
      pL->add(x,y);

      while(fgets(string, 30, f)!=NULL)
	{
	  sscanf( string, "%lg  %lg",&x,&y);
	  // check that the inputs are in increasing x order
	  if(fabs(x - oldx)> EPSILON)
	    {
	      pL->add(x,y);
	      if(y > maxposs) maxposs = y;
	      oldx = x;
	      oldy = y;
	    }
	  else
	    {
	      if(fabs(oldy - y) > EPSILON)
		{
		  pL->add(x,y);
		  if(y > maxposs) maxposs = y;
		  oldx =x+2*EPSILON;
		  oldy =y;
		}
	      else printf(" WARNING skipped value(x too small) %8.3f %8.3f \n", x,y);
	    }
	}

      if (y >EPSILON)  // check that last y = 0
	{
	  printf(" WARNING added end value, %8.3f %8.3f \n", x,0.0);
	  pL->add(x,0);
	}
      fclose(f);
    }

  virtual ~MFDPOSS()
    {
      delete pL;
    }

  MFDPOSS* Clone() const
  {
    return new MFDPOSS(this->pL);
  }

  void Update(double *Values)
    { //???????
    }

  void operator=(MFDPOSS  &T)
    {
      if(this == &T) return;

      if(pL) delete pL;
      pL = new LIST;
      if(T.NbParams() >0)
	{
	  T.pL->home();
	  pL->home();
	  maxposs = T.maxposs;
	  pL->add(T.pL->Get()->x, T.pL->Get()->y);

	  while( !T.pL->IsEnd() )
	    {
	      pL->next();
	      T.pL->next();
	      pL->add(T.pL->Get()->x, T.pL->Get()->y);
	    }
	  SetName(T.Name);
	}
    }

  const char *GetType() const { return Type(); }
  static const char *Type() { return "possibility_distribution"; }

  int NbParams() const { return pL->GetSize(); };
  double getMaxposs() const { return maxposs; };

  void Print(FILE * f)
  {
    long save_index = pL->curI();

    // pL->Print();  // Debug info
    pL->home();
    fprintf(f, "%8.3f%c%8.3f\n", pL->Get()->x, ESPACE,
	    pL->Get()->y);
    while( !pL->IsEnd() )
      {
	pL->next();
	fprintf(f, "%8.3f%c%8.3f\n", pL->Get()->x, ESPACE,
		pL->Get()->y);
      }
    pL->GotoI(save_index);
  }

  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
  {

    long save_index = pL->curI();
    MF::PrintCfg(num, f, fd);
    fprintf(f, "%c", START_NB);
    pL->home();
    fprintf(f, "\n%8.3f%c%8.3f\n", pL->Get()->x, ESPACE,  pL->Get()->y);

    while( !pL->IsEnd() )
      {
	pL->next();
	fprintf(f, fd, pL->Get()->x);
	fprintf(f, "%c", ESPACE);
	fprintf(f, fd, pL->Get()->y);
	fprintf(f,"\n");
      }
    pL->GotoI(save_index);
  }

  LIST* createList(double ls, double rs, double lk, double rk,
                 double deg, double maxp);

  double Kernel(double & left, double & right) const;
  double Support(double & left, double & right)const ;
  double AlphaKernel(double & left, double & right, double alpha) const;
  double GetDeg(double value) const;

  int  GetPoint(double &x, double &y, long index);

  fispro::POINT* CheckI(LIST *iL, LIST *sL, LIST *cL, int Ncro) const;
  MFDPOSS* Inter(MFDPOSS *tdp) const;
  MFDPOSS* Union(MFDPOSS *tdp);
  MFDPOSS* Join(MFDPOSS *tdp);

  std::list<MFDPOSS> *Union(std::list<MFDPOSS> *unL);

  MFDPOSS* minTnorme(double deg);
  MFDPOSS* prodTnorme(double deg);

  MFDPOSS* translate(double val, double vmin, double vmax);
  void Simplify(void);
  double computeArea();
  void DecompAcut(int nbalf);
};

class MFSINUS:public MF
//***********************
{
protected:
  double a, b;

public:

  MFSINUS(double aa, double bb) : MF()
    {
    if( (fabs(bb-aa)) < EPSILON )  // a > b
    	throw std::runtime_error( "~S2~MustBeDifferentfrom~S1~" );
    if( (bb-aa) < EPSILON )  // a > b
    	throw std::runtime_error( "~S2~MustBeHigherThan~S1~" );

    a = aa; b = bb;
    }

  MFSINUS(const MFSINUS &s) : MF(s)
    {
      a = s.a;
      b = s.b;
    }

  MF* Clone() const { return new MFSINUS(*this); }

  //! Function used by the optimization module
  void Normalize ( double lowerb , double upperb )
    {
      double ratio;
      ratio = upperb - lowerb;
      if (fabs(ratio) >EPSILON)
	{
	  a = (a - lowerb)/ratio;
	  b = (b - lowerb)/ratio;
	}

    }

  //! Function used by the optimization module
  void UnNormalize( double lowerb , double upperb )
    {
      double ratio;
     ratio=upperb - lowerb;
      if( fabs(ratio) >EPSILON)
	{
	  a = a*ratio + lowerb;
	  b = b*ratio + lowerb;
	}

    }

  const char *GetType() const { return Type(); }
  static const char *Type() { return "sinus"; }

  void ValParam(double &aa, double &bb) const
    { aa = a;  bb = b; }

  void Update(double *Values)
    {
      a = Values[0];
      b = Values[1];
    }

  int NbParams() const { return 2; };
  void GetParams( double *params ) const { ValParam( params[0], params[1]); }
  void SetInf( double inf ) { a = inf; }

  double GetDeg(double value) const {
      double val=0.0;
	  if (value < a) return(0.0);
	  if (value > b) return(0.0);
	  value=M_PI*((value-a)/(b-a));
	  val=sin(value);
      if (val < 0) return(0.0);
      return(val);
    }

  double Kernel(double & left, double & right) const{
     left=(a+b)*0.5;
     right=left;
     return left;
  }

  double Corner() const { return b; }

  double Support(double & left, double & right) const
    {
      left = a; right = b;
      return left + (right - left) / 2 ;
    }

  double AlphaKernel(double & left, double & right, double alpha) const {
      double aker=0.0;
      left=right=0.0;
      aker=asin(alpha);
      aker=a+((b-a)*aker/M_PI);
      left = aker;
      right = a+b-aker;
      return aker;
  }

  void Print(FILE * f) const
    {
      MF::Print(f);
      fprintf(f, "%8.3f%c%8.3f", a, SEPARE, b);
    }

  virtual void PrintCfg(int num, FILE *f, const char *fd = FORMAT_DOUBLE) const
    {
      MF::PrintCfg(num, f, fd);
      fprintf(f, "%c", START_NB);
      fprintf(f, fd, a);
      fprintf(f, "%c", SEPARE);
      fprintf(f, fd, b);
      fprintf(f, "%c\n", END_NB);
    }
};

class MFSINUSINF : public MFSINUS {

public:
	MFSINUSINF(double a, double b) : MFSINUS(a, b) {}

	MF* Clone() const { return new MFSINUSINF(*this); }

	const char *GetType() const { return Type(); }

	static const char *Type() { return "SinusInf"; }

	double GetDeg(double value) const {
		double val=0.0;
		if (value < a) return(1.0);
		if (value > b) return(0.0);
		value=M_PI*0.5*((value-a)/(b-a));
		val=cos(value);
		if (val < 0) return(0.0);
		return(val);
    }

	double Kernel(double & left, double & right) const {
		left=a;
	    right=left;
	    return left;
	}

	double AlphaKernel(double & left, double & right, double alpha) const {
		double aker=0.0;
		left=right=0.0;
		aker=acos(alpha);
		aker=a+((b-a)*aker*2.0/M_PI);
		left = a;
		right = aker;
		return aker;
	}
};

class MFSINUSSUP : public MFSINUS {

public:
	MFSINUSSUP(double a, double b) : MFSINUS(a, b) {}

	MF* Clone() const { return new MFSINUSSUP(*this); }

	const char *GetType() const { return Type(); }

	static const char *Type() { return "SinusSup"; }

	double GetDeg(double value) const {
		double val=0.0;
		if (value < a) return(0.0);
		if (value > b) return(1.0);
		value=M_PI*0.5*((value-a)/(b-a));
		val=sin(value);
		if (val < 0) return(0.0);
		return(val);
	}

	double Kernel(double & left, double & right) const {
		left=b;
	    right=left;
	    return left;
	}

    double AlphaKernel(double & left, double & right, double alpha) const {
    	double aker=0.0;
    	left=right=0.0;
    	aker=asin(alpha);
    	aker=a+((b-a)*aker*2.0/M_PI);
    	left = aker;
    	right = b;
    	return aker;
    }
};

#endif

//**************************  MF.H  **********************************

