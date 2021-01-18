//***********************************************************************

//
//
//                              MF.CPP
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC 
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : MF class functions used by FISPRO, part of library fispro

//**********************************************************************

#include "fis.h"

void MF::Centroid(double a_coupe, double & centre, double & masse, Trapeze * coord) const
  //*******************************************************************************
{
  double lkb, lke, lsb, lse, a_cut;
  double m1, m2, m3;
  double c1, c2, c3;

  if (a_coupe < EPSILON) { masse = 0.0; centre = 0.0; return; }
  if (a_coupe > 1.)  a_cut = 1.;
  else a_cut = a_coupe;

  lkb = lke = lsb = lse = -1.;

  AlphaKernel(lkb, lke, a_cut);
  Support(lsb, lse);
  // These four points define one trapeze whose area is computed from 
  // one rectangle and two triangles.
  coord->lk = lkb;  coord->rk = lke; coord->ls = lsb;  coord->rs = lse; 

  // Rectangle
  m1 = (lke - lkb) * a_cut;
  c1 = lkb + (lke - lkb) / 2;

  // Triangles
  m2 = (lkb - lsb) * a_cut / 2;
  c2 = lsb + (lkb - lsb) * 2 / 3;

  m3 = (lse - lke) * a_cut / 2;
  c3 = lke + (lse - lke) / 3;

  masse = m1 + m2 + m3;
  if( masse )
    centre = (m1 * c1 + m2 * c2 + m3 * c3) / masse;
  else
    centre = c1;
}

double MF::MFMatchDeg(MF * T) const
  //*************************
  // Fill the Mfdeg vector
{

  double sl, sr, stl, str, kl, kr, ktl, ktr, intersect;

  sl =  sr = stl = str = kl = kr = ktl = ktr = intersect = 0;

  // Case 1
  // return 0. if the supports do not intersect
  Support(sl, sr);
  T->Support(stl, str);
  if(sr < stl || str < sl) return 0.;

  // Case 2
  // return 1. if the kernels intersect
  Kernel(kl, kr);
  T->Kernel(ktl, ktr);
  if(ktl <= kr && kl < ktr) return 1.;

  // Case 3
  // T kernel is before this kernel
  if(ktr < kl) 
    intersect = (sl*(str-ktr)+str*(kl-sl)) / ((str-ktr)+(kl-sl));

  // T kernel is after this kernel
  else
    intersect = (stl*(sr-kr)+sr*(ktl-stl)) / ((sr-kr)+(ktl-stl));

  return GetDeg(intersect);
}

int MF::operator != (const MF &sef)
  //*******************************
{
  if( (strcmp( Name, sef.Name) != 0) || (strcmp( GetType(), sef.GetType()) != 0) || (NbParams() != sef.NbParams()) )
    return 1;

  double *params = new double[NbParams()];
  GetParams(params);
  double *sef_params = new double[NbParams()];
  sef.GetParams(sef_params);
  int retour = 0;
  for( int i=0 ; i<NbParams() ; i++ )
    if( params[i] != sef_params[i] )
      {
	retour = 1;   
	break;
      }
  delete [] params;
  delete [] sef_params;
  return retour;
}
  
void MF::SetName(const char *name)
  //******************************
{
  delete [] Name;
  Name = new char[strlen(name)+1];
  sprintf( Name, "%s", name );
}

MF * FuzNumber(double v, double kw, double sw)
//********************************************
{
  if(sw < 0. || kw < 0.)
    throw std::runtime_error("~SupportWidth~and~KernelWidth~MustBeNonNegative~");
  if(DBL_INF_EQUAL(sw, kw)) 
    throw std::runtime_error("~SupportWidth~MustBeHigherThan~KernelWidth~");
  
  if(DBL_EQUAL(kw, 0)) return new MFTRI(v, sw/2);
  return new MFTRAP(v-sw/2, v-kw/2, v+kw/2, v+sw/2);
}
