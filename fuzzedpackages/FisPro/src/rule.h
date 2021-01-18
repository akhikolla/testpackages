//***********************************************************************
//
//
//                              RULE.H
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : MF class functions used by FISPRO, part of library fispro

//**********************************************************************
#include "fis.h"

#ifndef __RULE_H
#define __RULE_H


//! To handle indifferent factors, corresponds to an inactive variable for the rule.
//! 0 is a convenient value as, from the interface side, the first MF is the number 1.
#define ANY      0



//******************************************************************
//
//
//                         CLASSE PREMISE
//
//
//
//******************************************************************

class PREMISE
//************
{
  protected :
    //! The number of proposisitions in the opremise, one per input
    int NbProps;
  //! The values of the proposisitions, one per input
  int * Props;
  //! Pointer on the input array
  FISIN ** Input;

  public :

    PREMISE(int ni, FISIN ** E, int Alloc = true)
    {
      NbProps = ni;
      Input = E;
      Props = NULL;
      if(NbProps && Alloc)
	Props = new int [NbProps];
      RstProps();
    }

  PREMISE()
    {
      NbProps = 0;
      Props = NULL;
      Input = NULL;
    }

  virtual PREMISE* Clone(FISIN ** E) = 0;

  PREMISE(PREMISE &prem, FISIN ** E)
    {
      NbProps = prem.NbProps;
      Input = E;
      Props = new int [NbProps];
      for( int i=0 ; i<NbProps ; i++ )
	Props[i] = prem.Props[i];
    }

  virtual ~PREMISE()
    {  delete [] Props;  }

  int operator != (const PREMISE &prem)
  {
    if( NbProps != prem.NbProps )
      return 1;
    for( int i=0 ; i<NbProps ; i++ )
      if( Props[i] != prem.Props[i] )
	return 1;
    return 0;
  }

  //! Virtual function to be defined within all the inherited classes.
  //! Returns the degree of match of the current the current example to the rule.
  //! This degree is computed using the public array of each of the input variables, Mfdeg.
  virtual double MatchDeg(void) = 0;


  //! Returns the number of different labels in the premise
  //! and -1 if the rules have a different number of factors.
  //! 'd' contains the number of the first not ANY label different
  //! from its  homologous in 'P' et and -1 when identical rules.
  //! 'ind' contains the number of differences between a given value factor
  //! and an indifferent factor (the corresponding partial distances are zero).
  int Distance(PREMISE * P, int & d, int & ind)
    {
      d = -1;
      ind = 0;
      if(NbProps != P->NbProps) return -1;

      int i, r;
      r = 0;
      for(i = 0; i < NbProps; i ++)
	{
	  // Both are ANY
	  if( (Props[i] == ANY) && (P->Props[i] == ANY) ) continue;

	  // Only one is ANY
	  if( (Props[i] == ANY) || (P->Props[i] == ANY) )
	    {
	      ind ++;
	      // if((P->Props[i] == ANY) && (d == -1)) d = i;
	      if(d == -1)  d = i;
	      continue;
	    }

	  if(Props[i] != P->Props[i])
	    {
	      if(d == -1 || r == 0) d = i;
	      // Priority is given to this configuration
	      r++;
	    }
	}
      return r;
    }

  int Compare(PREMISE * P)
    {
      if(NbProps != P->NbProps) return 1;

      int i;
      for(i = 0; i < NbProps; i ++)
	{
	  if( (Props[i] == ANY) || (P->Props[i] == ANY) ) continue;
	  if(Props[i] != P->Props[i])  return 1;
	}
      return 0;
    }

  virtual void Print(FILE *f)
    {
      int i;
      for(i = 0; i < NbProps; i++)
	fprintf(f, "%d%c ", Props[i], SEPARE);
    }

  void PrintProps(FILE * f)
    {
      int i;
      for(i = 0; i < NbProps; i++)
	fprintf(f, "%d", Props[i]);
    }

  // Factor handling functions from rule class
  int GetNbProp(void) { return NbProps; }

  int SetAProp(int value, int numero)
    {
      if( value > Input[numero]->GetNbMf() )
        ThrowFactorError( value, numero );
      if(numero >= 0 && numero < NbProps)
	{
	  Props[numero] = value;
	  return 0;
	}
      return 1;
    }

  void SetAProps(int * Valeurs)
    {
      int i;
      for(i = 0; i < NbProps; i ++)
        {
	  if( Valeurs[i] > Input[i]->GetNbMf() )
	    ThrowFactorError( Valeurs[i], i );
	  Props[i] = Valeurs[i];
        }
    }

  void RstProps()
    {
      for( int i=0; i<NbProps; i++ )
        Props[i] = 0;
    }
  void ThrowFactorError( int factor, int input_number )
    {
      char error_msg[100];
      sprintf( error_msg, "~RuleFactor~: %d >~NumberOfMFInInput~%d", factor, input_number+1 );
      throw std::runtime_error( error_msg );
    }

  void SetAProps(double * Valeurs)
    {
      int i;
      for(i = 0; i < NbProps; i ++)
        {
	  if( (int)Valeurs[i] > Input[i]->GetNbMf() )
	    ThrowFactorError( (int)Valeurs[i], i );
	  Props[i] = (int)Valeurs[i];
        }
    }


  // for rule generation, allocation is done little by little when needed
  int SetAProps(PREMISE *P, int alloc = true)
    {
      int i;
      if(NbProps != P->NbProps)  return 1;
      if(alloc)
	{
	  Props = new int [NbProps];
	}
      for(i = 0; i < NbProps; i ++)
	{
	  if( P->Props[i] > Input[i]->GetNbMf() )
	    ThrowFactorError( (int)P->Props[i], i );
	  Props[i] = (int)P->Props[i];
        }
      return 0;
    }

  void GetProps(int * v)
    {
      for(int i=0 ; i<NbProps ; i++ )
        v[i] = Props[i];
    }

  int GetAProp(int &v, int numero)
  {
	  if(numero >= 0 && numero < NbProps)
	  {
		  v = Props[numero];
		  return v;
	  }
	  return -1;
  }
};

//! Inherited class to implement Zadeh T-norm:  V(a,b)=min(a,b)
class PREMISE_MIN : public PREMISE
//**********************************
{
  public :

    PREMISE_MIN(int ni, FISIN ** E, int mem = true) : PREMISE(ni, E, mem) {}

  PREMISE_MIN(PREMISE_MIN &prem, FISIN ** E) : PREMISE(prem, E) {}

  virtual ~PREMISE_MIN() {}

  virtual PREMISE* Clone(FISIN ** E) { return new PREMISE_MIN(*this, E); }

  double MatchDeg(void)
    {
      int i;
      double degre;
      int valid;

      degre = 1.0;
      valid = false;

      for(i = 0; i < NbProps; i++)
	{
	  if(Input[i]->IsActive() == false)  continue;

	  valid = true;
	  if(Props[i] != ANY && Input[i]->Mfdeg()[Props[i] - 1] < degre)
	    degre = Input[i]->Mfdeg()[Props[i] - 1];
	}
      if(valid) return degre;
      return 0.0;
    }
};


//! Inherited class to implement the algebraic or probabilistic T-norm:  V(a,b)=a*b
class PREMISE_PROD : public PREMISE
//***********************************
{

  public :

    PREMISE_PROD(int ni, FISIN ** E, int mem = true) : PREMISE(ni, E, mem) {}

  PREMISE_PROD(PREMISE_PROD &prem, FISIN ** E) : PREMISE(prem, E) {}

  virtual ~PREMISE_PROD() {}

  virtual PREMISE* Clone(FISIN ** E) { return new PREMISE_PROD(*this, E); }

  double MatchDeg(void)
    {
      int i;
      double val;
      int valid;

      valid = false;
      val = 1.0;

      for(i = 0; i < NbProps; i++)
	{
	  if(Input[i]->IsActive() == false) continue;

	  valid = true;
	  if(Props[i] != ANY) val *= Input[i]->Mfdeg()[Props[i] - 1];
	}
      if(valid) return val;
      return 0.0;
    }
};

//! Inherited class to implement the Lukasiewicz T-norm:  V(a,b)=max(0, a + b - 1)
class PREMISE_LUKA : public PREMISE
//***********************************
{
  public :

    PREMISE_LUKA(int ni, FISIN ** E, int mem = true) : PREMISE(ni, E, mem){}

  PREMISE_LUKA(PREMISE_LUKA &prem, FISIN ** E) : PREMISE(prem, E) {}

  virtual ~PREMISE_LUKA() {}

  virtual PREMISE* Clone(FISIN ** E) { return new PREMISE_LUKA(*this, E); }

  double MatchDeg(void)
    {
      int i;
      double val;
      int valid;

      valid = false;

      val = 1.0 - NbProps;

      for(i = 0; i < NbProps; i++)
	{
	  if(Input[i]->IsActive() == false) val += 1.0;
	  else
	    {
	      if(Props[i] == ANY) val += 1.0;
	      else val += Input[i]->Mfdeg()[Props[i] - 1];
	      valid = true;
	    }
	}
      if(valid && val > 0) return val;
      return 0.;
    }
};


//******************************************************************
//
//
//                          CONCLUSION CLASS
//
//
//
//******************************************************************

class CONCLUSION
//**************
{
  protected :
    //! The number of rule conclusions, one per output
    int NbConcs;
  //! The values of the conclusions, one per output
  double * Concs;
  //! Pointer on the output array
  FISOUT **  Output;

  public :

    CONCLUSION(int no, FISOUT ** S)
    {
      NbConcs = no;
      Output = S;
      Concs = NULL;
      if(NbConcs)
	Concs = new double [NbConcs];
      RstConcs();
    }

  CONCLUSION()
    {
      NbConcs = 0;
      Concs = NULL;
      Output = NULL;
    }

  CONCLUSION(CONCLUSION &c, FISOUT ** S)
    {
      NbConcs = c.NbConcs;
      Output = S;
      Concs = new double[NbConcs];
      for( int i=0 ; i<NbConcs ; i++ )
	Concs[i] = c.Concs[i];
    }

  virtual ~CONCLUSION()
    { delete [] Concs; }

  int operator != (const CONCLUSION &c)
  {
    if( NbConcs != c.NbConcs )
      return 1;
    for( int i=0 ; i<NbConcs ; i++ )
      if( Concs[i] != c.Concs[i] )
	return 1;
    return 0;
  }

  void RstConcs()
    {
      for( int i=0 ; i<NbConcs ; i++ )
	Concs[i] = 0;
    }

  //! Function used by the optimization module
  void Normalize ( int nOut ,double lowerb , double upperb )
    {
      double ratio;
      ratio = upperb - lowerb;
      SetAConc( nOut , ( GetAConc ( nOut ) - lowerb)/ratio );
    }

  //! Function used by the optimization module
  void UnNormalize( int nOut, double lowerb , double upperb )
    {
      double ratio,result;
      ratio = ( double ) (upperb - lowerb);
      result =  ( double ) GetAConc ( nOut );
      result = result *ratio+lowerb;
      SetAConc( nOut , result );
    }

  void ThrowConcError( int action, int output_number )
    {
      char error_msg[100];
      sprintf( error_msg, "~RuleConc~: %d >~NumberOfMFInOutput~%d", action, output_number+1 );
      throw std::runtime_error( error_msg );
    }

  int Compare(CONCLUSION * P)
    {
      if(NbConcs != P->NbConcs) return 1;

      int i;
      for(i = 0; i < NbConcs; i ++)
	if(Concs[i] != P->Concs[i])  return 1;

      return 0;
    }

  virtual void Print(FILE *f, const char *fd = FORMAT_DOUBLE)
    {
      int i;
      for(i = 0; i < NbConcs; i++)
	{
	  fprintf(f, fd, Concs[i]);
	  fprintf(f, "%c", SEPARE);
	}
    }

  // Functions to handle actions from the rule class
  void SetConcs(double * Valeurs)
    {
      int i;
      for(i = 0; i < NbConcs; i ++)
        {
	  if( strcmp( Output[i]->GetOutputType(), OUT_FUZZY::OutputType() ) ==0 )
	    if( (int)Valeurs[i] > Output[i]->GetNbMf()  || (int)Valeurs[i] < 1 )
	      ThrowConcError( (int)Valeurs[i], i );
	  Concs[i] = Valeurs[i];
        }
    }

  void SetAConc(int i, double Valeur)
    {
      // change a rule conclusion
      // depending on situation call classcheck or ModPossibles
      // to update Possibles and Classes arrays
      if( strcmp( Output[i]->GetOutputType(), OUT_FUZZY::OutputType() ) ==0 )
        if( (int)Valeur > Output[i]->GetNbMf() || (int)Valeur < 1 )
          ThrowConcError( (int)Valeur, i );
      if(i >= 0  && i < NbConcs)
        Concs[i] = Valeur;
    }

  double GetAConc(int i)
    {
      if(i >= 0 && i < NbConcs)   	return Concs[i];
      return FisMknan();
    }

  int GetNbConc(void) { return NbConcs; }

};



#endif

//************************** RULE.H *********************************


