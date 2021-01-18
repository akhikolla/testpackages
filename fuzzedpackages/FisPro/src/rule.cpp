//***********************************************************************
//
//
//                              RULE.CPP
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC 
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : MF class functions used by FISPRO, part of library fispro

//**********************************************************************

#include "fis.h"


RULE::RULE(int nI, FISIN ** E, int nO, FISOUT ** S, char * cConj, char * Val) //throw(std::runtime_error)
  //*************************************************************************
{
  // Val corresponds to one line in the file
  // as many integers as the number of inputs
  // as many doubles as the number of outputs separated by SEPARE.
  // The alloc flag can be false. In this case the memory for the factor array
  // is not allocated. This option is used for generating 
  // the rules.
  int ret;

  Init();

  if((strlen(Val) == 0) ||  (Val[0] == 0x0D))
    {
      sprintf( ErrorMsg, "~EmptyStringInRuleConstructor~\n");
      throw std::runtime_error( ErrorMsg );
    }

  try
    { 
      SetPremise( nI, E, cConj ); 
    }
  catch( std::exception &e )
    {
      sprintf( ErrorMsg, "~ErrorInFISFile~\n%.100s", e.what() );
      throw std::runtime_error( ErrorMsg );
    }

  SetConclusion( nO, S );

  double *Tab = new double [nI + nO + 1];
  try { ret = SearchNb(Val, Tab, nI + nO + 1);} 
  catch( std::exception &e )
    {
      delete [] Tab;
      sprintf( ErrorMsg, "~ErrorInFISFile~\n~ErrorInRuleValues~: %.50s\n%.50s~", Val, e.what() );
      throw std::runtime_error( ErrorMsg );
    }

  if(ret < nI + nO)
    {
      delete [] Tab;
      sprintf( ErrorMsg, "~ErrorInFISFile~\n~ErrorInRuleValues~: %.50s~", Val );
      throw std::runtime_error( ErrorMsg );
    }

  try
    {
      Prem->SetAProps(Tab);
      Conclu->SetConcs(&Tab[nI]);
      if(ret > (nI+nO)) SetExpertWeight(Tab[nI+nO]);
      delete [] Tab;
      Active = true;
    }
  catch( std::exception &e )
    {
      delete [] Tab;
      sprintf( ErrorMsg, "~ErrorInFISFile~\n~ErrorInRuleValues~: %.50s\n%.100s", Val, e.what() );
      throw std::runtime_error( ErrorMsg );
    }
} 

RULE::RULE(RULE &regle, FISIN ** E, FISOUT ** S)
  //********************************************
{
  Init();
  Active = regle.Active;
  Weight = regle.Weight;
  ExpertWeight = regle.ExpertWeight;
  Prem = regle.Prem->Clone(E);
  Conclu = new CONCLUSION(*regle.Conclu, S);
}

RULE::RULE(RULE &regle, FISIN ** E)
  //*******************************
{
  Init();
  Active = regle.Active;
  Weight = regle.Weight;
  ExpertWeight = regle.ExpertWeight;
  Prem = regle.Prem->Clone(E);
}

int RULE::operator != (const RULE &regle)
  //*************************************
{
  if( (Active != regle.Active) || (ExpertWeight != regle.ExpertWeight) || (*Prem != *regle.Prem) || (*Conclu != *regle.Conclu) )
    return 1;
  return 0;
}


void RULE::SetPremise( int nI, FISIN ** E, char * cConj )
  //*****************************************************
{
  PREMISE *temp = NULL;

  if(!strcmp(cConj, RULE::PremiseProd()))
    temp = new PREMISE_PROD(nI, E);
  else if(!strcmp(cConj, RULE::PremiseMin()))
    temp = new PREMISE_MIN(nI, E);
  else if(!strcmp(cConj, RULE::PremiseLuka()))
    temp = new PREMISE_LUKA(nI, E);
  else
    {
      sprintf( ErrorMsg, "~UnknownConjunction~: %.50s~", cConj );
      throw std::runtime_error( ErrorMsg );
    }

  delete Prem;
  Prem = temp;
}


void RULE::SetConclusion( int nO, FISOUT ** S )
  //*******************************************
{
  CONCLUSION *temp = new CONCLUSION( nO, S );
  delete Conclu;
  Conclu = temp;
}


//**************************    RULE.CPP   *******************************



