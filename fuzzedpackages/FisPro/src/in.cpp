//***********************************************************************

//
//
//                              IN.CPP
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : FISIN class functions used by FISPRO, part of library fispro

//**********************************************************************

#include "fis.h"

//*********     Membership degree functions ***************************

double FISIN::GetADeg(int s, double v)
  //**********************************
  // Returns the membership degree of 'v' to fuzzy set 's'
{
  if((s >= 0) && (s < Nmf))
    return (Fp[s]->GetDeg(v));
  return FisMknan();
}

void FISIN::ldLinMFs()
{
#ifdef _OPENMP
#pragma omp critical(ldLinMFs)
  // must retest the entry condition to prevent a race condition
  if (linMF.size()+nonLinMF.size()!=(unsigned) Nmf)  
#endif
    {
      linMF.clear();
      nonLinMF.clear();
      for(int i = 0; i < Nmf; i++)
        {
          MFPWLinear m=Fp[i]->getMFPWLinear();
          if (m.valid())
            linMF.push_back(make_pair(i,m));
          else
            nonLinMF.push_back(make_pair(i,Fp[i]));
        }
#ifdef _OPENMP
      _Mfdeg.resize(omp_get_max_threads());
      for (unsigned int i=0; i<_Mfdeg.size(); i++)
        _Mfdeg[i].resize(Nmf);
#else
      _Mfdeg.resize(Nmf);
#endif
    }
}

/*
  int FISIN::GetDegs(double v)
{
  int i;
  int retour = 1;
  for(i = 0; i < Nmf; i++)
    if((Mfdeg[i] = Fp[i]->GetDeg(v)) != 0) retour = 0;
  return retour;
}
*/

int FISIN::GetDegs(double v)  
  //************************
  // Returns:  1 if all degrees are zero else 0
{
  unsigned int i;
  int retour = 1;

//  for(i = 0; i < Nmf; i++)
//    if((Mfdeg[i] = Fp[i]->GetDeg(v)) != 0) retour = 0;

/* 
   instead of the previous commented out code, we copy all piece-wise linear membership 
   functions into a new non-polymorphic data structure that has an inlinable GetDeg() function.
   This allows better vectorisation of the GetDeg() loop, which is impeded by the virtual 
   function call. Nonlinear MFs need to be handled by the usual virtual function
*/


  if (linMF.size()+nonLinMF.size()!=(unsigned) Nmf) ldLinMFs();
          
  unsigned int sz=linMF.size();
#ifdef __ICC
#pragma ivdep
#pragma vector always
#endif
  for (i=0; i<sz; i++)
    {
      retour &= (Mfdeg()[linMF[i].first]=linMF[i].second.GetDeg(v)) == 0;
      //assert(fabs(Mfdeg()[linMF[i].first]-Fp[linMF[i].first]->GetDeg(v))<1e-5);
    }
  for (i=0; i<nonLinMF.size(); i++)
    retour &= (Mfdeg()[nonLinMF[i].first]=nonLinMF[i].second->GetDeg(v)) == 0;

  return retour;
}

// version of GetDegs with return value optimised away
void FISIN::GetDegsV(double v)  
{
  unsigned int i;

  if (linMF.size()+nonLinMF.size()!= (unsigned) Nmf) ldLinMFs(); 
          
  unsigned int sz=linMF.size();
#ifdef __ICC
#pragma ivdep
#pragma vector always
#endif
  for (i=0; i<sz; i++)
    {
      Mfdeg()[linMF[i].first]=linMF[i].second.GetDeg(v);
      //assert(fabs(Mfdeg()[linMF[i].first]-Fp[linMF[i].first]->GetDeg(v))<1e-5);
    }
  for (i=0; i<nonLinMF.size(); i++)
    Mfdeg()[nonLinMF[i].first]=nonLinMF[i].second->GetDeg(v);
}



int FISIN::SetEqDegs(double v)
  //**************************
{
  int i;
  v = 0.5 / (double) Nmf;
  Mfdeg().resize(Nmf);
  for(i = 0; i < Nmf; i++)
    Mfdeg()[i] = v;
  return 0;
}

int FISIN::GetRandDegs(double v)
  //****************************
  // Returns:  1 if all degrees are zero else 0
{
  v = ValInf + FisRand() * (ValSup - ValInf);
  return GetDegs(v);
}

double FISIN::MFMatchADeg(int s, MF * T)
  //************************************
  // Returns the membership degree of 'v' to fuzzy set 's'
{
  if((s >= 0) && (s < Nmf))
    return (Fp[s]->MFMatchDeg(T));
  return FisMknan();
}

double FISIN::MFMatchDegs(MF * T)
  //*****************************
  // Fill the Mfdeg() vector
{
  int i;
  int retour = 1;
   Mfdeg().resize(Nmf);
   for(i = 0; i < Nmf; i++)
    if((Mfdeg()[i] = Fp[i]->MFMatchDeg(T)) != 0) retour = 0;
  return retour;
}

void FISIN::GetMfCenters(double *p)
  //*******************************
{
  int i;
  char * type;
  double t[20];  // The number of params for MFDISCRETE is unknown

  for(i = 0; i < Nmf; i++)
    {
      type = (char *) Fp[i]->GetType();
      Fp[i]->GetParams(t);
      if(!(strcmp(type, MFTRI::Type())))
	p[i] = t[1];
      else if(!(strcmp(type, MFTRAPINF::Type())))
	p[i] = t[1];
      else if(!(strcmp(type, MFTRAPSUP::Type())))
	p[i] = t[1];
      else if(!(strcmp(type, MFTRAP::Type())))
	p[i] = (t[1]+t[2])/2;
      else if(!(strcmp(type, MFGAUSS::Type())))
	p[i] = t[0];
      else if(!(strcmp(type, MFGBELL::Type())))
	p[i] = t[1];
      else if(!(strcmp(type, MFDISCRETE::Type())))
	p[i] = t[0];
      else if(!(strcmp(type, MFDOOR::Type())))
	p[i] = (t[0]+t[1])/2;
      else if(!(strcmp(type, MFUNIV::Type())))
	p[i] = (t[0]+t[1])/2;
      else if(!(strcmp(type, MFSINUS::Type())))
	{
	  if (fabs(t[2]) <EPSILON)
	    p[i] = (t[0]+t[1])/2;
	  if  (fabs(t[2]-90.0) <EPSILON)
	    p[i] = t[0];
	  if (fabs(t[2]+90.0) <EPSILON)
	    p[i] = t[1];
	}
    }
}
void FISIN::GetSFPparams(double *&p, int* &typemf, int & size, FILE *display)
  //***********************************************************
{
  int i,k;
  char * type;
  double t[4];
  bool sfp=false;
  int * sorted=NULL;
  // test SFP-exit if not
  sfp=IsSfp(sorted);
  delete [] sorted;
  if (!sfp)
    throw std::runtime_error("Input partition is not sfp");

  // case of 2 MFs
  k=0;
  
  if (Nmf<2)
    {
      sprintf( ErrorMsg, "~Nmf~must~be~>=2~");
      throw std::runtime_error( ErrorMsg );
    }
  typemf=new int[Nmf];

  if (Nmf==2)
    {
      // 2 MF SFP
      size=2;
      p=new double[size];
      typemf[0]=0;
      typemf[1]=0;
      Fp[0]->GetParams(t);
      p[0]=t[0];
      p[1]=t[1];
    }
  else 
    //2 passages-first one to determine tri and trap
    {
      size=2;//for first and last parameters
      for(i = 1; i < Nmf-1; i++)
	{
	  type = (char *) Fp[i]->GetType();
	  if(!(strcmp(type, MFTRAP::Type())))
	    {
	      typemf[i]=1;
	      size=size+2;
	    }
	  else if (!(strcmp(type, MFTRI::Type())))
	    {
	      typemf[i]=2;
	      size=size+1;
	    }
	  else
	    {
	      sprintf( ErrorMsg, "~only~tri~or~trap~MFs~allowed~");
	      throw std::runtime_error( ErrorMsg );
	    }
	}
      //allocate arrays with the right size
      typemf[0]=0;
      p=new double[size];
      Fp[0]->GetParams(t);
      p[0]=t[1];
      k=1;
      //second passage to store parameters according to type
      for(i = 1; i < Nmf-1; i++)
	{
	  type = (char *) Fp[i]->GetType();
	  Fp[i]->GetParams(t);
	  if(typemf[i]==1)//trap
	    {
	      p[k] = t[1];
	      p[k+1]=t[2];
	      k=k+2;
	    }
	  else//tri
	    {
	      p[k] = t[1];
	      k=k+1;
	    }
	}
      //last param=semi trap sup mid parameter
      typemf[Nmf-1]=0;
      Fp[Nmf-1]->GetParams(t); 
      p[k]=t[1];
      if (display)
	{
	  fprintf(display, "in GetSFPParams k=%d,size=%d, parameters:",k,size);
	  for (int ii=0;ii<size;ii++)
	    fprintf(display, "%g ",p[ii]);
	  fprintf(display, "\n");
	}
    }
}

void FISIN::GetBreakPoints(double * & Bp, int & np)
  //***********************************************
{
  double *tmp;
  double sl, sr, stl, str, kl, kr, ktl, ktr, v;
  int i, n;

  Bp = NULL;
  np = 0;
  if(! Nmf) return;

  tmp = new double [Nmf * 2 - 1];
  i = n = 0;

  tmp[n++] = Fp[i]->Kernel(kl, kr);
  Fp[i]->Support(sl, sr);

  for(i = 1; i < Nmf; i++)
    {
      v = Fp[i]->Kernel(ktl, ktr);
      Fp[i]->Support(stl, str);
      if(stl < sr) // Intersection
	tmp[n++] = (stl*(sr-kr)+sr*(ktl-stl)) / ((sr-kr)+(ktl-stl));

      tmp[n++] = v;
      kl = ktl; kr = ktr;
      sl = stl; sr = str;
    }

  np = n;
  if(np == Nmf * 2 - 1)
    {
      Bp = tmp;
      return;
    }

  Bp = new double [np];
  for(i = 0; i < np; i++) Bp[i] = tmp[i];
  delete [] tmp;
  return;
}

//*********     Constructors and init functions *********************

void FISIN::Init()
  //**************
{
  // Pointers
  Fp = NULL;
  Name = NULL;
  //Mfdeg = NULL;
#ifdef _OPENMP
  _Mfdeg.resize(omp_get_max_threads());
#endif
  dPart =NULL;
  nPart =0;
  // Data
  Nmf = 0;
  ValInf = 0.;
  ValSup = 1.;
  SetName( "" );
  initNormalize();
  Kw = 0.0;
  Sw = 0.0;
}

FISIN::FISIN( const FISIN & entree ): privMfdeg(entree.privMfdeg)
  //********************************
{
  int i;
  Init();

  SetName(entree.Name);
  SetRange( entree.ValInf, entree.ValSup);
  OLower = entree.OLower;
  OUpper = entree.OUpper;
  active = entree.active;
  Nmf = entree.Nmf;
  if( Nmf != 0 )
    {
      Fp = new MF * [Nmf];
      //Mfdeg = new double [Nmf];
      for( i=0 ; i<Nmf ; i++ )
	Fp[i] = NULL;
      for( i=0 ; i<Nmf ; i++ )
	{
	  Fp[i] = entree.Fp[i]->Clone();
	  Fp[i]->SetName(entree.Fp[i]->Name);//added by bch March 8,2010
	}
    }
}

FISIN::FISIN(int n, double min, double max, int tri): privMfdeg(false) //throw(std::runtime_error)
  //************************************************
{
  double dyn, gauche, centre, droite;
  int i;

  Init();
  SetRange( min, max );
  Nmf = n;
  active = true;

  if(Nmf == 0)
    return;
  Fp = new MF *[Nmf];
  for( int i=0 ; i<Nmf ; i++ )
    Fp[i] = NULL;
  //Mfdeg = new double [Nmf];

  // Generate the grid
  dyn = ValSup - ValInf;

  if(Nmf == 1)
    {
      Fp[0] = new MFTRI(-INFINI, dyn/2., INFINI);
      return;
    }

  dyn /= (Nmf - 1);  // no center incrementing

  for(i = 0; i < Nmf; i++)
    {
      if(i == 0) gauche = -INFINI;
      else gauche = ValInf + (dyn * (i-1));

      centre = ValInf + (dyn * i);

      if(i == Nmf - 1) droite = INFINI;
      else droite = ValInf + (dyn * (i+1));

      if((i == 0) &&  (! tri))
	Fp[i] = new MFTRAPINF(ValInf, centre, droite);
      else if((i == Nmf - 1)  &&  (! tri))
	Fp[i] = new MFTRAPSUP(gauche, centre, ValSup);
      else Fp[i] = new MFTRI(gauche, centre, droite);
    }
} // End of the Grille constructor

FISIN::FISIN(double *t, int n, double min, double max, int sort): privMfdeg(false) //throw(runtime_error);
  //************************************************************
{
  int i;
  double left, right, centre;

  Init();
  SetRange(min, max);
  Nmf = n;
  active = true;

  if(Nmf == 0)
    return;
  Fp = new MF *[Nmf];
  for(i = 0 ; i < Nmf ; i++ )
    Fp[i] = NULL;
  //Mfdeg = new double [Nmf];

  if(sort)
    qsort(t, n, sizeof(double), CmpDblAsc);

  for(i = 0; i < Nmf; i++)
    {
      if(i == 0) left = -INFINI;
      else left = t[i-1];

      centre = t[i];

      if(i == Nmf - 1) right = INFINI;
      else right = t[i+1];

      if(i == 0)
	Fp[i] = new MFTRAPINF(ValInf, centre, right);
      else if(i == Nmf - 1)
	Fp[i] = new MFTRAPSUP(left, centre, ValSup);
      else Fp[i] = new MFTRI(left, centre, right);
    }
}


FISIN::FISIN(int n, double *t, double min, double max) //throw(runtime_error);
//****************************************************
// Build a strong fuzzy partition with trapezoidal MF
{
  if(!n || n%2) 
    {
      sprintf(ErrorMsg, "~EvenNumberOfPointsNeededFor~TrapezoidalSFP~(n=%d)" , n);
      throw std::runtime_error(ErrorMsg);
    }
  int i;
  Init();
  SetRange(min, max);
  active = true;
  Nmf = (n/2)+1;
  Fp = new MF *[Nmf];
  for(i = 0 ; i < Nmf ; i++ ) Fp[i] = NULL;

  Fp[0] = new MFTRAPINF(ValInf, t[0], t[1]);
  Fp[Nmf-1] = new MFTRAPSUP(t[n-2], t[n-1], ValSup);

  for(i = 1; i < Nmf-1; i++)
    Fp[i] = new MFTRAP(t[i*2-2], t[i*2-1], t[i*2], t[i*2+1]);
}


FISIN::FISIN(double *ArrayCenter, int *CenterType, int n, double lowerb, double upperb,
	     double olower,double oupper,int indPFF=1) : privMfdeg(false)//throw(std::runtime_error)
  //***********************************************************************************
  // Generates a standardized fuzzy partition if indPFF = 1
  // Generates a regular input if indPFF=0
{
  if ( n == 0  ||  n <  2  ) return;
  if ( lowerb > upperb  ) return;
  Init();
  active = true;
  int cpt2=0;
  Nmf = n;
  Fp = new MF *[Nmf];
  for( int i=0 ; i<Nmf ; i++ )
    Fp[i] = NULL;

  //Mfdeg = new double [Nmf];
  // Generates the grid
  if(Nmf == 1) {
    Fp[0] = new MFTRI(lowerb, ArrayCenter[0] , upperb);
    return;
  };

  for( int cpt = 0; cpt < Nmf; cpt++)
    {
      if ( indPFF == 1 )
	{
	  if (cpt == 0) {
	    Fp[cpt] = new MFTRAPINF( ArrayCenter[cpt2], ArrayCenter[cpt+1] , ArrayCenter[cpt+2] );
	    cpt2+=2;
	  }
	  else {
	    if (cpt == Nmf - 1) {
	      Fp[cpt] = new MFTRAPSUP(  ArrayCenter[cpt2-1] , ArrayCenter[cpt2] , ArrayCenter[cpt2+1] ) ;
	      cpt2+=2;
            }
	    else {
	      if ( CenterType[cpt] == 1)
		{
		  Fp[cpt] = new MFTRI( ArrayCenter[cpt2-1] , ArrayCenter[cpt2] , ArrayCenter[cpt2+1] );
		  cpt2++;
		};
	      if ( CenterType[cpt] == 2)
		{
		  Fp[cpt] = new MFTRAP( ArrayCenter[cpt2-1] ,
					ArrayCenter[cpt2] , ArrayCenter[cpt2+1],
					ArrayCenter[cpt2+2]
					);
		  cpt2+=2;
		};

	    };
	  };
	}
      else
	// case of not standardized partition
	{
	  switch (  CenterType[cpt] )
	    {
	    case 1 :
	      Fp[cpt] = new MFTRI( ArrayCenter[cpt2] , ArrayCenter[cpt2+1] , ArrayCenter[cpt2+2] );
	      cpt2+=3;
	      break;
	    case 2 :
	      Fp[cpt] = new MFTRAP( ArrayCenter[cpt2] , ArrayCenter[cpt2+1] , ArrayCenter[cpt2+2], ArrayCenter[cpt2+3] );
	      cpt2+=4;
	      break;
	    case 3 :
	      Fp[cpt] = new MFTRAPINF( ArrayCenter[cpt2] , ArrayCenter[cpt2+1], ArrayCenter[cpt2+2] );
	      cpt2+=3;
	      break;
	    case 4 :
	      Fp[cpt] = new MFTRAPSUP( ArrayCenter[cpt2] , ArrayCenter[cpt2+1] , ArrayCenter[cpt2+2]) ;
	      cpt2+=3;
	      break;
	    case 5 :
	      Fp[cpt] = new MFGBELL( ArrayCenter[cpt2] , ArrayCenter[cpt2 + 1] , ArrayCenter[cpt2 + 2] ) ;
	      cpt2+=3;
	      break;
	    case 6 :
	      Fp[cpt] = new MFGAUSS( ArrayCenter[cpt2] , ArrayCenter[cpt2+1] ) ;
	      cpt2+=2;
	      break;
	    case 7 :
	      // loop must be done to reconstruct the double array for MFDISCRETE ...
	      //          Fp[cpt2] = new MFDISCRETE( ArrayCenter[n-2] , ArrayCenter[n-1] , upperb) ;
	      cpt2+=cpt2;
	      break;
	    case 8 :
	      Fp[cpt] = new MFUNIV( ArrayCenter[cpt2] , ArrayCenter[cpt2+1] ) ;
	      cpt2+=2;
	      break;
	    default :
	      break;
	    }
	};
    };

  SetRangeOnly( lowerb, upperb );
  OLower=olower;
  OUpper=oupper;

  return;

} // End of the Grille constructor



void FISIN::initNormalize()
  //***********************
{
  OUpper=0;
  OLower=1;
  return;
}

void FISIN::Normalize()
  //*******************
{ 
// save values for unnormalizing later
 OUpper=ValSup;
 OLower=ValInf;
 for ( int cptMF =0 ; cptMF < Nmf ; cptMF++ ) 
   GetMF(cptMF)->Normalize(OLower,OUpper);
 
 SetRangeOnly(0,1);
 return;
}

void FISIN::UnNormalize() //throw(std::runtime_error)
  //*********************
{
  if(OLower > OUpper)
    {
      sprintf( ErrorMsg, "~NotPossibleTheFISWasNotNormalized~");
      throw std::runtime_error( ErrorMsg );
    }
  for ( int cptMF =0 ; cptMF < Nmf ; cptMF++ ) GetMF(cptMF)->UnNormalize(OLower,OUpper);
  SetRangeOnly(OLower,OUpper);
  return;
}

void FISIN::Init(ifstream &f, int bufsize, int num) //throw(std::runtime_error)
  //***********************************************
{
  char *tmp=NULL, *buf=NULL;
  double *Tab=NULL;
  int i;

  try
    {
      tmp = new char[bufsize];
      buf = new char[bufsize];

      do{                           // Skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      // Active = Activate
      sprintf( tmp, "Active=" );
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %-3d\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", GetType(), num, tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if(SearchStr(buf, tmp))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %-3d\n~StringSeparatorNotFoundInString~: %.50s~", GetType(), num, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if( ! strncmp(tmp, "no", 4) )  active = false;   //  inactive input
      else  if( ! strncmp(tmp, "yes", 4) ) active = true;
      else
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %-3d\n~ExpectedString~: Activate=yes or no\n~ReadString~: %.50s~", GetType(), num, tmp );
	  throw std::runtime_error( ErrorMsg );
	}

      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %
      sprintf( tmp, "Name=" );
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %-3d\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", GetType(), num, tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if(SearchStr(buf, tmp))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %-3d\n~StringSeparatorNotFoundInString~: %.50s~", GetType(), num, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      SetName(tmp);

      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %


      //  Range
      sprintf( tmp, "Range=" );
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %.50s\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", GetType(), Name, tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      Tab = new double [2];
      if((strlen(buf + strlen(tmp)) == 0) || (*(buf + strlen(tmp))==0x0D))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~EmptyString~");
	  throw std::runtime_error( ErrorMsg );
	}

      try { SearchNb(buf, Tab, 2, SEPARE, START_NB, END_NB);}
      catch( std::exception &e )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %.50s (~Range~)\n%.80s", GetType(), Name, e.what() );
	  throw std::runtime_error( ErrorMsg );
	}

      try { SetRange( Tab[0], Tab[1] ); }
      catch( std::exception &e )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %.50s\n%.80s", GetType(), Name, e.what() );
	  throw std::runtime_error( ErrorMsg );
	}

      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf( tmp, "NMFs=" );
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %.50s\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", GetType(), Name, tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      i =  atoi(buf + strlen(tmp));
      Nmf = i;

      if(Nmf != 0)
	{
	  Fp = new MF *[Nmf];
	  for(i = 0; i < Nmf; i++) Fp[i] = NULL;

	  //Mfdeg = new double [Nmf];

	  for(i = 0; i < Nmf; i++)
	    {
	      do{                           // skip empty lines
		f.getline(buf, bufsize);
	      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
	      // ignores the ^M of MS-DOS files read in Linux-Unix
	      // ignores comment lines beginning with # or %
	      ReadMf(buf, i + 1);  // tags are in  base 1
	    }
	}
      delete [] Tab;
      delete [] tmp;
      delete [] buf;
    }
  catch( std::exception &e )
    {
      if(Tab) delete [] Tab;
      delete [] tmp;
      delete [] buf;
      throw;
    }
} // end of FISIN(file)

void FISIN::ReadMf(char * ligne, int NumSef) //throw(std::runtime_error)
  //****************************************
{
  char * tmp=NULL, *nom=NULL, *type=NULL;
  int offset;
  double * bornes=NULL;
  int len, retour;

  try
    {
      len = strlen(ligne);
      tmp = new char[len];
      nom  = new char[len];
      type = new char[len];

      sprintf(tmp,"MF%d=",NumSef);
      if( strncmp(tmp, ligne, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %.50s\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", GetType(), Name, tmp, ligne );
	  throw std::runtime_error( ErrorMsg );
	}

      offset = strlen(tmp);
      if(SearchStr(ligne + offset, nom))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %.50s\n~MF~: %-3d\n~StringSeparatorNotFoundInString~: %.50s~", GetType(), Name, NumSef, ligne );
	  throw std::runtime_error( ErrorMsg );
	}

      char * index = strchr(ligne, SEPARE);
      if(SearchStr(index, type))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %.50s\n~MF~: %-3d\n~StringSeparatorNotFoundInString~: %.50s~", GetType(), Name, NumSef, ligne );
	  throw std::runtime_error( ErrorMsg );
	}

      // read bounds
      offset = index - ligne + 2;
      index = strchr(ligne + offset, SEPARE);
      if((strlen(index) == 0) ||  (index[0]==0x0D))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~~EmptyString~");
	  throw std::runtime_error( ErrorMsg );
	}

      int nb_bornes = CntNbs(index, SEPARE, START_NB, END_NB);
      bornes = new double [nb_bornes];
      try { retour = SearchNb(index,  bornes, nb_bornes, SEPARE, START_NB, END_NB);}
      catch(std::exception &e)
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %.50s\n~MF~: %-3d\n~FormatErrorInMFParameters~: %.50s\n%.50s~", GetType(), Name, NumSef, index, e.what() );
	  throw std::runtime_error( ErrorMsg );
	}

      try
	{
	  if( (!strcmp(type, MFTRI::Type())) && (retour == 3))
	    Fp[NumSef - 1] = new MFTRI(bornes[0], bornes[1], bornes[2]);
	  else if( (!strcmp(type, MFTRI::Type())) && (retour == 2))
	    Fp[NumSef - 1] = new MFTRI(bornes[0], bornes[1]);
	  else if( (!strcmp(type, MFTRAP::Type())) && (retour == 4))
	    Fp[NumSef - 1] = new MFTRAP(bornes[0], bornes[1], bornes[2], bornes[3]);
	  else if( (!strcmp(type, MFTRAPINF::Type())) && (retour == 3))
	    Fp[NumSef - 1] = new MFTRAPINF(bornes[0], bornes[1], bornes[2]);
	  else if( (!strcmp(type, MFTRAPSUP::Type())) && (retour == 3))
	    Fp[NumSef - 1] = new MFTRAPSUP(bornes[0], bornes[1], bornes[2]);
	  else if( (!strcmp(type, MFGBELL::Type())) && (retour == 3))
	    Fp[NumSef - 1] = new MFGBELL(bornes[0], bornes[1], bornes[2]);
	  else if( (!strcmp(type, MFGAUSS::Type())) && (retour == 2))
	    Fp[NumSef - 1] = new MFGAUSS(bornes[0], bornes[1]);
	  else if( (!strcmp(type, MFUNIV::Type())) && (retour == 2))
	    Fp[NumSef - 1] = new MFUNIV(bornes[0], bornes[1]);
	  else if(!strcmp(type, MFDISCRETE::Type()))
	    Fp[NumSef - 1] = new MFDISCRETE(retour, bornes);
	  else if( (!strcmp(type, MFDOOR::Type())) && (retour == 2))
	    Fp[NumSef - 1] = new MFDOOR(bornes[0], bornes[1]);
	  else if( (!strcmp(type, MFSINUS::Type())) && (retour == 2))
	    Fp[NumSef - 1] = new MFSINUS(bornes[0], bornes[1]);
	  else if( (!strcmp(type, MFSINUSINF::Type())) && (retour == 2))
	    Fp[NumSef - 1] = new MFSINUSINF(bornes[0], bornes[1]);
	  else if( (!strcmp(type, MFSINUSSUP::Type())) && (retour == 2))
	    Fp[NumSef - 1] = new MFSINUSSUP(bornes[0], bornes[1]);
	  else
	    {
	      sprintf( ErrorMsg, "UnknownMFType~: %.50s~\n~Or~\n~IncorrectNumberOfBounds~: %-3d", type, retour );
	      throw std::runtime_error( ErrorMsg );
	    }
	  Fp[NumSef - 1]->SetName(nom);  // assigns the  MF name
	}
      catch( std::exception &e )  // intercepts the exceptions in the MF constructors
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~%.50s~: %.50s\n~MF~: %-3d\n~%.50s~", GetType(), Name, NumSef, e.what() );
	  throw std::runtime_error( ErrorMsg );
	}
      delete [] nom;
      delete [] type;
      delete [] tmp;
      delete [] bornes;
    }
  catch( std::exception &e )
    {
      delete [] nom;
      delete [] type;
      delete [] tmp;
      if(bornes) delete [] bornes;
      throw;
    }
} // End of  ReadMf()

int FISIN::operator != (const FISIN &entree)
  //****************************************
  // Return value: 0 when equal else 1
{
  if( (strcmp(Name, entree.Name) != 0) ||
      (ValInf != entree.ValInf) ||
      (ValSup != entree.ValSup) ||
      (active != entree.active) ||
      (Nmf != entree.Nmf) )
    return 1;

  for( int i=0 ; i<Nmf ; i++ )
    if( *Fp[i] != *entree.Fp[i] )
      return 1;

  return 0;
}

void FISIN::SetName(const char *name)
  //*********************************
{
  delete [] Name;
  Name = new char[strlen(name)+1];
  sprintf( Name, "%s", name );
}

void FISIN::SetStdMfNames(void)
  //***************************
{
  char MfName[15];
  int i;

  for(i = 0; i < Nmf; i ++)
    {
      sprintf(MfName, "MF%d", i+1);
      Fp[i]->SetName(MfName);
    }
}
//*********     Membership function management  *********************


void FISIN::SetRangeOnly( double range_inf, double range_sup )  //throw(std::runtime_error)
  //**********************************************************
{
  if( range_inf >= range_sup )
    throw std::runtime_error( "~Range~Upper~MustBeHigherThan~Range~Lower~" );
  ValInf = range_inf;
  ValSup = range_sup;
}

void FISIN::SetRange( double range_inf, double range_sup )  //throw(std::runtime_error)
  //******************************************************
{
  if( range_inf >= range_sup )
    {
      sprintf(ErrorMsg, "~Range~Upper~(%8.3f)~MustBeHigherThan~Range~Lower~(%8.3f)" , range_inf, range_sup);
   throw std::runtime_error(ErrorMsg);
    }
  ValInf = range_inf;
  ValSup = range_sup;
  // modifying MS' bounds for semitrapezoidal inf and sup
  for( int i=0 ; i<Nmf ; i++ )
    {
      if( strcmp(Fp[i]->GetType(), MFTRAPINF::Type()) == 0 )
	((MFTRAPINF *)Fp[i])->SetInf( ValInf );
      if( strcmp(Fp[i]->GetType(), MFTRAPSUP::Type()) == 0 )
	((MFTRAPSUP *)Fp[i])->SetSup( ValSup );
      if( strcmp(Fp[i]->GetType(), MFUNIV::Type()) == 0 )
	{
	  ((MFUNIV *)Fp[i])->SetInf( ValInf );
	  ((MFUNIV *)Fp[i])->SetSup( ValSup );
	}
    }
}


void FISIN::ReplaceMF( int sef_number, MF *new_sef )
  //************************************************
{
  if (sef_number>=0 && sef_number<Nmf)
    {
      delete Fp[sef_number];
      Fp[sef_number] = new_sef;
    }
}

void FISIN::MoveMF( int sef_number, int move )
  //*******************************************
{
  int new_position = sef_number + move;
  if( (new_position >= 0) && (new_position < Nmf) )
    {
      MF *temp = Fp[new_position];
      Fp[new_position] = Fp[sef_number];
      Fp[sef_number] = temp;
    }
}

void FISIN::AddMF( MF *sef , int insert_before = -1 )
  //************************
{
  // warning: sef argument must not be deleted
  // new MF (sef) is added after all existing MFs
  if (insert_before<0) insert_before=Nmf;
  MF **temp = new MF * [Nmf+1];
  //double *MfdegTemp = new double [Nmf+1];
  for( int i=0 ; i< insert_before; i++ )
    temp[i] = Fp[i]->Clone();
  temp[insert_before] = sef;
  for (int i=insert_before;i<Nmf;i++)
     temp[i+1] = Fp[i];
  delete [] Fp;
  //delete [] Mfdeg;

  Nmf++;
  Fp = temp;
  //Mfdeg = MfdegTemp;
  Mfdeg().resize(Nmf);
}

void FISIN::RemoveMF( int sef_number )
  //**********************************
{

  // PROTECTION
  if ((sef_number>=0) && (sef_number<Nmf))
    {
      MF **temp = new MF * [Nmf-1];
      //double *MfdegTemp = new double [Nmf-1];


      for( int i=0, j=0 ; j<Nmf ; j++ )
	if( j != sef_number )
	  {
	    temp[i] = Fp[j]->Clone();
	    i++;
	  }
      for (int i=0;i<Nmf;i++)
	delete Fp[i];
      delete [] Fp;
      //delete [] Mfdeg;

      Nmf--;
      Fp = temp;
      //Mfdeg = MfdegTemp;
      Mfdeg().resize(Nmf);
    }
  // else no change
}

void FISIN::DecomposePart(std::list<double> *dL)
//**********************************************
{
  double l,r;
  int j;
  for (j=0; j< Nmf; j++)
    {
      GetMF(j)->Kernel(l,r);
      dL->push_back(l);
      dL->push_back(r);
      GetMF(j)->Support(l,r);
      dL->push_back(l);
      dL->push_back(r);
    }
  dL->sort();
  dL->unique();

  // create the doors of the partition
  dPart = new MFDOOR[dL->size()-1];
  double inf = 0;
  double sup = 0;
  nPart = 0;
  for(std::list<double>::iterator it = dL->begin(); it!=dL->end(); it++)
    {
      if (it == dL->begin()) inf = *it;
      else
	{
	  sup =*it;
	  if (fabs(sup-inf)>EPSILON)
	    {
	      dPart[nPart].SetInf(inf);
	      dPart[nPart].SetSup(sup);
	      nPart++;
	    }
	  inf = sup;
	}
    }
}

// Find the intervals of the partition which intersect
// the alpha-cut door, fill the arrays of bounds
// and return their number.
int FISIN::getIntersect(ACUT &door, double *inf, double *sup)
//***********************************************************
{
  MFDPOSS *mfd, *mfi, *mfdp;
  int i, n=0;

  mfd = new MFDPOSS(door);
  for( i = 0; i < nPart; i++)
    {
      mfdp =new MFDPOSS(dPart[i],0);
      mfi = mfd->Inter(mfdp);
      if (mfi!=NULL)
	{
	  mfi->Support(inf[n],sup[n]);
	  delete mfi;
	  n++;
	}
      delete mfdp;
    }
  delete mfd;
  return(n);
}

void FISIN::DecomposePart(FILE *display)  // Not used
//*****************************
{
	double lk,rk;
	dPart = new MFDOOR[(2*Nmf-1)];

	Fp[0]->Kernel(lk,rk);
	if(display){
		fprintf(display, "Nmf %d\n",Nmf);
		fprintf(display, "i 0, lk %8.3f, rk %8.3f\n",lk, rk);
	}
	//kernel door
	dPart[0].SetInf(lk);
	dPart[0].SetSup(rk);
	int j=1;
	for( int i=1; i<Nmf ; i++ )
	{
		//intersection door
		dPart[j].SetInf(rk);

		Fp[i]->Kernel(lk,rk);
		if(display)
			fprintf(display, "i %d, lk %8.3f, rk %8.3f\n",i,lk, rk);
		dPart[j++].SetSup(lk);

		//kernel door
		dPart[j].SetInf(lk);
		dPart[j++].SetSup(rk);
	}
	nPart = j;
}
//
double *kG=NULL;
int CmpKAsc(const void * a, const void *b)
  //**************************************
{
  if(kG[*(unsigned int *)a] > kG[*(unsigned int *)b]) return 1;
  if(kG[*(unsigned int *)a] < kG[*(unsigned int *)b]) return  -1;
  return 0;
}
//


bool FISIN::IsSfp(int *& s)
//*************************
{
  int i, j;
  char *t;
  double *par, *ppar;
  bool ret = true;
  bool order = false;
  double inf, sup;
  t = NULL;
  inf = sup = 0.;

  if(Nmf == 0) return false;
  if(Nmf == 1) return true;

  for(i = 0; i < Nmf; i++)
    {
      t = (char *) Fp[i]->GetType();
      if(strcmp(t, MFTRI::Type()) && strcmp(t, MFTRAP::Type())
	 && strcmp(t, MFTRAPINF::Type()) && strcmp(t, MFTRAPSUP::Type()))
	ret = false;
      if(i && Fp[i]->Kernel(inf, sup) <  Fp[i-1]->Kernel(inf, sup)) order = true;
    }

  if(ret == false) return ret;

  if(order)
    {
      s = new int [Nmf];
      kG = new double [Nmf];

      for(i = 0; i < Nmf; i++)
	{
	  kG[i]=Fp[i]->Kernel(inf, sup);
	  s[i] = i;
	}
      qsort(s, Nmf, sizeof(int),CmpKAsc);
      delete [] kG;

      MF **tmp = new MF * [Nmf];
      for(i = 0; i < Nmf; i++)
	{
	  tmp[i] = Fp[s[i]]->Clone();
	  tmp[i]->SetName("");
	}
      for(i = 0; i < Nmf; i++)
	{
	  delete Fp[i];
	  Fp[i] = NULL;
	}
      delete [] Fp;
      Fp = tmp;
     }

  if(strcmp(Fp[0]->GetType(), MFTRAPINF::Type())) ret = false;
  if(strcmp(Fp[Nmf-1]->GetType(), MFTRAPSUP::Type())) ret = false;

  par = new double[4];
  ppar = new double [4];

  Fp[0]->GetParams(ppar);
  for(i = 1; i < Nmf; i++)
    {
      Fp[i]->GetParams(par);
      if(i > 1 && !strcmp(Fp[i-1]->GetType(), MFTRAP::Type()))
	{ if((par[0] != ppar[2]) || (par[1] != ppar[3])) ret = false;}
      else if((par[0] != ppar[1]) || (par[1] != ppar[2])) ret = false;

     for(j = 0; j < 4; j++) ppar[j] = par[j];
    }

  delete [] par; delete [] ppar;
  return ret;
}

void FISIN::SetTemplate(double kernelWeight, double supportWeight)
//****************************************************************
//for GUI of fuzzy input
{
	Kw = kernelWeight;
	Sw = supportWeight;
}
double FISIN::GetKernelWeightTemplate()
//**************************************
{return Kw;}

double FISIN::GetSupportWeightTemplate()
//**************************************
{return Sw;}
//**************************************


// This function only works for normalized data between 0 and 1
// and for a standard fuzzy partition. 
// Use CheckFuzDist before and UnNormalize after
// d=1 debug mode
// Return the distance between x and y

double FISIN::Distance(double x, double y, int norm, int d)
//*********************************************************
{
  double mux,muy,dist;
  int k,mfx,mfy;
	
  if (fabs(x-y) < EPSILON) 
    {
 //     if(d) fprintf(stdout, "\nNull distance\n");
      return 0.;
    }
       
  mux = muy = 0.0;
  mfx = mfy = -1;

  GetDegs(x);
  for(k = 0; k < Nmf ; k++) 
    if(Mfdeg()[k]>0.)
      {
	mux=Mfdeg()[k];
	mfx=k;
	break;
      }

  GetDegs(y);
  for(k = 0; k < Nmf ; k++) 
    if(Mfdeg()[k]>0.)
      {
	muy=Mfdeg()[k];
	mfy=k;
	break;
      }
	
  dist=fabs((mux-muy+mfy-mfx));
  if (norm)
    dist=dist/(Nmf-1);

  //if (d)
  //  fprintf(stdout, "\nx:%f y: %f  mfx:%d mfy:%d  mux: %f muy:%f  dist:%f ",x, y, mfx, mfy, mux, muy, dist);

  return dist;
}

void FISIN::CheckFuzDist(void)
//****************************
{
  bool sfp=false;
  int * sorted=NULL;

  sfp=IsSfp(sorted);
  delete [] sorted;
  if (!sfp)
    throw std::runtime_error("Input partition is not sfp");
  
  Normalize();

}

void FISIN::PcPe(double * dat, int n, double &pc, double &pe)
//***********************************************************
{ 
  double deg, sdeg = 0.;
  int i,j;

  pc = pe = 0.;
  for(i = 0; i < n; i++)
    {
      GetDegsV(dat[i]); 
      for(j = 0; j < Nmf; j++)
	{
	  deg = Mfdeg()[j];
	  sdeg += deg;
	  pc += deg * deg;
	  if(deg > EPSILON && deg < (1.0 - EPSILON)) 
	    pe += deg * log(deg);
	}
    }
  pc /= sdeg;
  pe /= sdeg;
  pe *= -1.;
}
void FISIN::Tri2Trap(void)
//************************
{ 
  int i;
  double *par;
  par = new double [3];
  char * oldName=NULL;
  int len;

  for(i = 0; i < Nmf; i++)
    {
     if(! strcmp(Fp[i]->GetType(), MFTRI::Type())) 
       {
	 Fp[i]->GetParams(par); 
	 len=strlen(Fp[i]->Name);
	 oldName=new char[len+1];
	 strcpy(oldName,Fp[i]->Name);
	 delete Fp[i];
	 Fp[i] = new MFTRAP(par[0], par[1], par[1], par[2]);
	 Fp[i]->SetName(oldName);
	 delete [] oldName;
       } 
    } 
  delete [] par;
}
//**************************    IN.CPP   *******************************
