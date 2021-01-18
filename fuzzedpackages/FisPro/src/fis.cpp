//***********************************************************************
//
//
//                              FIS.CPP
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : FIS class functions used by FISPRO, part of library fispro


//**********************************************************************

#include <errno.h>
#include "fis.h"



//*********************      Constructors and init functions **********
FIS::FIS()
  //******
{
  Init();
  SetName("");
  SetConjunction(RULE::PremiseMin());
  SetMissingValues(FIS::RandomMissingValues());
  SetErrorIndex(FIS::RmseErrorIndex());
}



FIS::FIS( const FIS & sif )
  //***********************
{
  int i;

  Init();
  SetName(sif.Name);
  SetConjunction(sif.TypeConj());
  SetMissingValues(sif.MissingValues());
  SetErrorIndex(sif.ErrorIndex());
  NbIn = sif.NbIn;
  NbOut = sif.NbOut;
  NbRules = sif.NbRules;
  NbActRules = sif.NbActRules;
  NbExceptions = sif.NbExceptions;


  if( NbIn != 0 )
    {
      In = new FISIN *[NbIn];
      for(i = 0 ; i < NbIn ; i++ )
	In[i] = NULL;
      for(i = 0 ; i < NbIn ; i++ )
	In[i] = new FISIN(*(sif.In[i]));
    }


  if( NbOut != 0 )
    {
      Out = new FISOUT *[NbOut];
      for(i = 0 ; i < NbOut ; i++ )
	Out[i] = NULL;
      for(i = 0 ; i < NbOut ; i++ )
	Out[i] = sif.Out[i]->Clone();
      OutValue = new double[NbOut];
      OutErr = new double[NbOut];
    }


  if( NbRules != 0 )
    {
      Rule = new RULE *[NbRules];
      for(i = 0 ; i < NbRules ; i++ )
	Rule[i] = NULL;
      for(i = 0 ; i < NbRules ; i++ )
	Rule[i] = new RULE(*sif.Rule[i], In, Out);
    }

  for(i = 0 ; i < NbOut ; i++ )  Out[i]->InitPossibles(Rule, NbRules, i);
}


void FIS::Init() //throw()
  //************
{
  Out = NULL;
  In = NULL;
  Rule = NULL;
  Name = NULL;
  OutValue = NULL;
  OutErr = NULL;
  cConjunction = NULL;
  strMissingValues = NULL;
  strErrorIndex = NULL;
  NbIn = 0;
  NbOut = 0;
  NbRules = 0;
  NbActRules = 0;
  PIn = RMSE = MAE = 0.;
}


void FIS::InitSystem(const char * fichier, int Cover) //throw(std::runtime_error)
  //*************************************************
{
  int i, j, bsize;
  ifstream f(fichier);

  if (! f)
    {
      sprintf( ErrorMsg, "~CannotOpenFISFile~: %.100s~", fichier);
      throw std::runtime_error( ErrorMsg );
    }

  bsize = MaxLineSize(f);
  ReadHdr(f, bsize);
  NbActRules = NbRules;

  if( NbIn != 0 )
    {
      In = new FISIN * [NbIn];
      for(i = 0 ; i < NbIn ; i++ ) In[i] = NULL;
    }

  if( NbOut != 0 )
    {
      Out = new FISOUT * [NbOut];
      for(i = 0 ; i < NbOut ; i++ ) Out[i] = NULL;
      OutValue = new double  [NbOut];
      OutErr = new double  [NbOut];
    }

  if( NbRules != 0 )
    {
      Rule = new RULE * [NbRules];
      for(i = 0 ; i < NbRules ; i++ ) Rule[i] = NULL;
    }

  for(i = 0; i < NbIn; i++)     ReadIn(f, bsize, i);
  for(i = 0; i < NbOut; i++)    ReadOut(f, bsize, i, Cover);
  ReadRules(f, bsize);
  NbActRules = NbRules;
  for(i = 0; i < NbOut; i++) // At least one implicative rule base
    if(!strcmp(Out[i]->Defuzzify(), OUT_FUZZY::ImpFuzzyDefuz()))
      for(j = 0; j < NbRules; j++) Rule[j]->SetExpertWeight(1.0);

  if(NbExceptions) ReadExcep(f, bsize);

  for(i = 0; i < NbOut; i++) Out[i]->InitPossibles(Rule, NbRules, i);

  SetErrorIndex(FIS::RmseErrorIndex());
}


FIS& FIS::operator=( const FIS & sif )
  //***********************
{
  // Standish version
  int i;

  //copied from destructor
  {
    for( i=0 ; i<NbIn ; i++) delete In[i];
    delete [] In;
    for( i=0 ; i<NbOut ; i++) delete Out[i];
    delete [] Out;
    for( i=0 ; i<NbRules ; i++) delete Rule[i];
    delete [] Rule;
    delete [] OutValue;
    delete [] OutErr;
    delete [] Name;
    delete [] cConjunction;
    delete [] strMissingValues;
    delete [] strErrorIndex;
  }

  Init();
  SetName(sif.Name);
  SetConjunction(sif.TypeConj());
  SetMissingValues(sif.MissingValues());
  SetErrorIndex(sif.ErrorIndex());
  NbIn = sif.NbIn;
  NbOut = sif.NbOut;
  NbRules = sif.NbRules;
  NbActRules = sif.NbActRules;
  NbExceptions = sif.NbExceptions;


  if( NbIn != 0 )
    {
      In = new FISIN *[NbIn];
      for(i = 0 ; i < NbIn ; i++ )
	In[i] = new FISIN(*(sif.In[i]));
    }


  if( NbOut != 0 )
    {
      Out = new FISOUT *[NbOut];
      for(i = 0 ; i < NbOut ; i++ )
	Out[i] = sif.Out[i]->Clone();
      OutValue = new double[NbOut];
      OutErr = new double[NbOut];
    }


  if( NbRules != 0 )
    {
      Rule = new RULE *[NbRules];
      for(i = 0 ; i < NbRules ; i++ )
	Rule[i] = new RULE(*sif.Rule[i], In, Out);
    }

  for(i = 0 ; i < NbOut ; i++ )  Out[i]->InitPossibles(Rule, NbRules, i);
  return *this;
}




int FIS::operator != (const FIS &sif) const
  //*********************************
  // returns 1 if the  2 FIS are different, 0 if identical
{
  int i;
  if( (strcmp( Name, sif.Name ) != 0) ||
      (strcmp( TypeConj(), sif.TypeConj() ) != 0) ||
      (strcmp( MissingValues(), sif.MissingValues() ) != 0) ||
      (strcmp( ErrorIndex(), sif.ErrorIndex() ) != 0) ||
      (NbIn != sif.NbIn) ||
      (NbOut != sif.NbOut) ||
      (NbRules != sif.NbRules) )
    return 1;


  for( i=0 ; i<NbIn ; i++ )
    if( *In[i] != *sif.In[i] )
      return 1;


  for( i=0 ; i<NbOut ; i++ )
    if( *Out[i] != *sif.Out[i] )
      return 1;


  for( i=0 ; i<NbRules ; i++ )
    if( *Rule[i] != *sif.Rule[i] )
      return 1;
  return 0;
}




int FIS::RulePos(RULE *R, int depart, int conc)
  //*******************************************
{
  int i;

  for(i = depart; i < NbRules; i ++)
    if(conc && ! R->Compare(Rule[i]) ) return i;
    else if(!conc && ! R->ComparePremises(Rule[i]) ) return i;

  return -1;
}


void FIS::RuleWeights (double * values, double * weights)
  //*****************************************************
{
  int i;

  //! Mfdeg array filling.
  for(i = 0; i < NbIn ; i++)
    {
      if(In[i]->IsActive() == false) continue;
      In[i]->GetDegsV(values[i]);
    }

  //! Mfdeg array filling.
  for(i = 0; i < NbIn ; i++)
    {
      if(In[i]->IsActive() == false) continue;
      In[i]->GetDegs(values[i]);
    }

  //! Computation of the matching degree of each of the active rules.
  for(i = 0; i < NbRules; i++)
    if ( Rule[i]->IsActive() )
      weights[i] = Rule[i]->MatchDeg();
    else weights[i] = -1;
}

//**********************    INFER     ******************************


double FIS::Infer(double * v, int out_number, FILE * fic, FILE *display, double deg) const
  //******************************************************************
  // v is the item description (input values and optionally output ones)
  // out_number is the output to infer, if -1 all outputs are infered
  // fic is the file in which the result is printed
  // display is the file in which intermediate result are displayed
{
  int i;
  double MaxWeight;


  if(NbRules == 0)
    {
      sprintf( ErrorMsg, "~NoRuleToInfer~");
      throw std::runtime_error( ErrorMsg );
    }

  if(! NbActRules)
    {
      for(i = 0; i < NbOut; i++)
	{
	  if(! Out[i]->IsActive() ) continue;
	  OutValue[i]  = Out[i]->DefaultValue();
	}
      return 0.;
    }

  MaxWeight = 0.0;

  if(display) fprintf(display, "\n");
#ifdef _OPENMP
#pragma omp parallel
#endif
  //! Step 1 : input fuzzification.
  //! Mfdeg array filling.
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for(i = 0; i < NbIn; i++)
    {
      if(In[i]->IsActive() == false) continue;

      if(FisIsnan(v[i]))            // missing value
	{
	  if(!strcmp(strMissingValues, FIS::RandomMissingValues()))
	    In[i]->GetRandDegs(v[i]);
	  else if(!strcmp(strMissingValues, FIS::MeanMissingValues()))
	    In[i]->SetEqDegs(v[i]);
	  else
	    {
	      sprintf( ErrorMsg, "~UnknownMissingValueStrategy~: %.50s", strMissingValues );
	      throw std::runtime_error( ErrorMsg );
	    }
	}
      In[i]->GetDegsV(v[i]);
      if(display) In[i]->PrintDeg(display);
    }


  //! Step 2 : Inference.
  //! Computation of the matching degree of each of the active rules.
#ifdef _OPENMP
  double lMaxWeight=MaxWeight; //thread local MaxWeight
  {
    double& MaxWeight=lMaxWeight;
#pragma omp for schedule(dynamic)
#endif
    for(int i = 0; i < NbRules; i++)
      {
	if ( Rule[i]->IsActive() )
	  {
	    Rule[i]->ExecRule();
	    if(Rule[i]->Weight > MaxWeight) MaxWeight = Rule[i]->Weight;
	  }
      }
#ifdef _OPENMP
  }
#pragma omp critical
  {if (lMaxWeight > MaxWeight) MaxWeight=lMaxWeight;} //reduce to shared variable
#endif

  //! Step 3 : Defuzzification.
  //! Rule conclusions are weighted and aggregated
  //! according to the  defuzzification strategy chosen for each output
  //! One numerical value per output is computed
  //#ifdef _OPENMP
  //#pragma omp for schedule(dynamic)
  //#endif
  for(i = 0; i < NbOut; i++)
    {
      if(out_number > -1 && i != out_number) continue;
      if(! Out[i]->IsActive() ) continue;
      OutValue[i] = Out[i]->InferOut(Rule, NbRules, fic, display, deg);
    }

  return MaxWeight;
} // End of Infer()


double FIS::Infer(MF ** v, int out_number, FILE *fic, FILE *display) const
  //***************************************************************
  // v is the item description (input values and optionally output ones)
  // out_number is the output to infer, if -1 all outputs are infered
  // fic is the file in which the result is printed
  // display is the file in which intermediate result are displayed
  // This version of Infere deals with items described by fuzzy sets
  // instead of crisp values.
  // It differs from the classical one mainly in the fuzzification step
{
  int i;
  double MaxWeight;

  if(NbRules == 0)
    {
      sprintf( ErrorMsg, "~NoRuleToInfer~");
      throw std::runtime_error( ErrorMsg );
    }

  if(! NbActRules)
    {
      for(i = 0; i < NbOut; i++)
	{
	  if(! Out[i]->IsActive() ) continue;
	  OutValue[i]  = Out[i]->DefaultValue();
	}
      return 0.;
    }

  MaxWeight = 0.0;

  if(display) fprintf(display, "\n");

  //! Step 1 : Fuzzification
  for(i = 0; i < NbIn; i++)
    {
      if(In[i]->IsActive() == false) continue;
      In[i]->MFMatchDegs(v[i]);
      if(display) In[i]->PrintDeg(display);
    }


  //! Step 2 : Inference.
  //! Computation of the matching degree of each of the active rules.
  for(i = 0; i < NbRules; i++)
    {
      if ( Rule[i]->IsActive() )
	{
	  Rule[i]->ExecRule();
	  if(Rule[i]->Weight > MaxWeight) MaxWeight = Rule[i]->Weight;
	}
    }

  //! Step 3 : Defuzzification.
  //! Rule conclusions are weighted and aggregated
  //! according to the  defuzzification strategy chosen for each output
  //! One numerical value per output is computed
  for(i = 0; i < NbOut; i++)
    {
      if(out_number > -1 && i != out_number) continue;
      if(! Out[i]->IsActive() ) continue;
      OutValue[i] = Out[i]->InferOut(Rule, NbRules, fic, display);
    }

  return MaxWeight;
} // End of Infer()

int FIS::CheckConsistency(void)
  //***************************
{
  int i, j;
  int min, max, v=0;

  // Check consistency between rules and inputs
  // Number of values in the premise part equal to number of inputs
  // CALLS InitPossibles for all outputs

  // case of no rules: exit
  if (NbRules<=0)
    return 0;

  if(Rule[0]->GetNbProp() != NbIn) return -100;

  // Linguistic labels ranging from 0 to nmf
  for(j = 0; j < NbIn; j++)
    {
      min = 10;
      max = -1;
      for(i = 0; i < NbRules; i++)
	{
	  Rule[i]->GetAProp(v, j);
	  if(v < min) min = v;
	  if(v > max) max = v;
	}
      if(min < 0 || max > In[j]->GetNbMf()) return -101 + j;
    }

  // Check consistency between rules and outputs
  if(Rule[0]->GetNbConc() != NbOut) return -200;

  for(j = 0; j < NbOut; j++)
    {
      // Linguistic labels ranging from 0 to nmf for fuzzy outputs
      if(Out[j]->GetNbMf()) // Fuzzy output
	{
	  min = 10;
	  max = -1;
	  for(i = 0; i < NbRules; i++)
	    {
	      v = (int) Rule[i]->GetAConc(j);
	      if(v < min) min = v;
	      if(v > max) max = v;
	    }
	  if(min < 0 || max > Out[j]->GetNbMf()) return -200 + j;
	}

      // Update (or set) array pointeurs
      Out[j]->InitPossibles(Rule, NbRules, j);

      // update NbActRules (member of fis)
      ComputeNbActRule();
    }
  return 0;
}

double FIS::InferCheck(double * v, double ** Val, int nb, int out_number, FILE * fic, FILE *display)
  //*********************************************************************************************
{
  int j;

  // add protection for no rule case - bch July 21, 2009
  if (NbRules<=0)
    {
      sprintf( ErrorMsg, "~No rule - inference is not possible~");
      throw std::runtime_error( ErrorMsg );
    }

  // CheckConsistency also allocates Possibles array
  if( (j = CheckConsistency()) ) return j;

  // Reallocate Classes array from data (if Val) or from rule conclusions
  InitClassLabels(Val, nb);

  // call to infer
  return Infer(v, out_number, fic, display);
}

double FIS::InferCheck(MF ** v, double ** Val, int nb, int out_number, FILE * fic, FILE *display)
  //*********************************************************************************************
{
  int j;
  // add protection for no rule case - bch July 21, 2009
  if (NbRules<=0)
    {
      sprintf( ErrorMsg, "~No rule - inference is not possible~");
      throw std::runtime_error( ErrorMsg );
    }

  // CheckConsistency also allocates Possibles array
  if( (j = CheckConsistency()) ) return j;

  // Reallocate Classes array from data (if Val) or from rule conclusions
  InitClassLabels(Val, nb);

  // call to infer
  return Infer(v, out_number, fic, display);
}

int FIS::InferFatiCheck(MFDPOSS ** v, int out_number, int nalphacut, double ** Val, int nb, FILE * fic, FILE *display)
  //*********************************************************************************************
{
  int j;
  MFDPOSS *dposs=NULL;

  // free MFConc array - it will be reallocated in InitPossibles, called by CheckConsistency
  DeleteMFConc(out_number);
  // CheckConsistency also allocates Possibles array  and MFConc array
  if( (j = CheckConsistency()) ) return j;

  // Reallocate Classes array from data (if Val) or from rule conclusions
  InitClassLabels(Val, nb);

  // call to infer
  dposs=InferFati(v, nalphacut, out_number, fic, display);

  // reset MFConc to empty dposs for correct display
  for (j=0;j<NbRules;j++)
    {
      delete Out[out_number]->MfConc[j];
      Out[out_number]->MfConc[j]= NULL;
      Out[out_number]->MfConc[j]= new MFDPOSS();
    }
  delete dposs;
  dposs=NULL;
  return 0;
}
void FIS::SetClassLabels(int num, double *Val, int nb)
  //**************************************************
{
  if(Val && Out[num]->Classification() &&
     (! strcmp(Out[num]->GetOutputType(), OUT_CRISP::OutputType())))
    {
      if (! strcmp(Out[num]->Defuzzify(), OUT_CRISP::SugenoDefuz()))
	{
	  DEFUZ_SugenoClassif * D = (DEFUZ_SugenoClassif *) Out[num]->Def;
	  if(Val) D->InitClasses(Val, nb);
	}

      if (! strcmp(Out[num]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))
	{
	  DEFUZ_MaxCrisp * D = (DEFUZ_MaxCrisp *) Out[num]->Def;
	  D->InitClasses(Val, nb);
	}
    }
}


void FIS::InitClassLabels(double **Val, int nb)
  //*******************************************
{
  int NumS;
  double * col;

  col = NULL;

  for(NumS = 0;  NumS < NbOut; NumS++)
    // for each output in the classification case
    {
      if(Out[NumS]->Classification() &&
	 (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType())))
	{
	  if(Val)
	    {
	      col = new double [nb];
	      GetColumn(Val, nb, NbIn + NumS, col);
	    }

	  if (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz()))
	    {
	      DEFUZ_SugenoClassif * D = (DEFUZ_SugenoClassif *) Out[NumS]->Def;
	      if(Val) D->InitClasses(col, nb);
	      else D->SetClasses(Out[NumS]->GetPossibles(), Out[NumS]->GetNbPossibles());
	    }

	  if (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))
	    {
	      DEFUZ_MaxCrisp * D = (DEFUZ_MaxCrisp *) Out[NumS]->Def;
	      if(Val) D->InitClasses(col, nb);
	      else D->SetClasses(Out[NumS]->GetPossibles(), Out[NumS]->GetNbPossibles());
	    }
	}
    }

  if(col) delete [] col;
}

int FIS::ClassCheck(int *& ResClassif, double *& Lab, double ** Data, int NbEx, int NumS)
  //*************************************************************************************
{
  //called once before calling Performance on data array
  // from java, fisopt or other module

  int j;

  // Check fis consistency and init possible values for rule conclusions
  // Delete array ptrs depnding on output type,   if necessary
  // and allocate new arrays
  //
  // Possibles Array in all cases
  if( (j = CheckConsistency()) ) return j;

  // calls ClassifCheck to update Classes and ResClassif arrays
  // (useful only in classification crisp case-does nothing otherwise)
  ClassifCheck(Data, NbEx, NumS);
  ResClassifAlloc(ResClassif,Lab,NumS);

  return 0;
}

int FIS::ClassCheckNoAllocResClassif(double ** Data, int NbEx, int NumS)
  //********************************************************************
{
  //called once before calling Performance on data array
  // from java, fisopt or other module
  //only difference with ClassCheck:  does not reallocate ResClassif

  int j;

  if( (j = CheckConsistency()) ) return j;

  // Delete array ptrs depnding on output type,   if necessary
  // and allocate new arrays
  //
  // Possibles Array in all cases
  // calls ClassifCheck to update Classes and ResClassif arrays
  // (useful only in classification crisp case-does nothing otherwise)

  //Out[NumS]->InitPossibles(Rule, NbRules);
  ClassifCheck(Data, NbEx, NumS);
  return 0;
}


int FIS::ClassifCheck(double ** Data, int NbEx, int NumS) const
  //*****************************************************
{
  //does something only in classification crisp case
  // deletes and recreates Classes arrays
  // casts defuz object according to defuz type
  //does not allocate resclassif and lab arrays (classification results)

  //does nothing otherwise

  int i;
  double * tmp;

  // casts defuz object
  // what follows does not depend on defuzzification type

  if(Out[NumS]->Classification() &&
     (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType())) )
    {
      if (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz()))
	{
	  DEFUZ_SugenoClassif * D = (DEFUZ_SugenoClassif *) Out[NumS]->Def;
	  if (D != NULL)
	    {
	      tmp = new double [NbEx];
	      for(i = 0; i < NbEx; i ++)  tmp[i] = Data[i][NbIn + NumS];
	      D->InitClasses(tmp, NbEx);
	      delete [] tmp;
	    }
	  else
		  throw std::runtime_error("error in ClassifCheck, Defuz object not initialized");
	}

      else if (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))
	{
	  DEFUZ_MaxCrisp * D = (DEFUZ_MaxCrisp *) Out[NumS]->Def;
	  if (D != NULL)
	    {
	      tmp = new double [NbEx];
	      for(i = 0; i < NbEx; i ++)  tmp[i] = Data[i][NbIn + NumS];
	      D->InitClasses(tmp, NbEx);
	      delete [] tmp;
	    }
	  else
		  throw std::runtime_error("error in ClassifCheck, Defuz object not initialized");
	}
    }
  return 0;
}
//*******************       ResClassifAlloc     *********************************
int FIS::ResClassifAlloc(int *& ResClassif, double *& Lab, int NumS) const
  //****************************************************************
{
  // allocates  ResClassif array
  // and updates Lab pointer to point on class labels from defuz object in output object
  // (for printing classification results)
  int i, NClas;

  if(Out[NumS]->Classification() &&
     (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType())) &&
     ( (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz()))	 ||
       (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))) )
    {
      // initialize array for misclassified items
      NClas = Out[NumS]->Def->NbClasses();
      if (NClas >0)
	{
	  // if resclassif already allocated delete it
	  if (ResClassif != NULL) delete [] ResClassif;
	  ResClassif=NULL;
	  Lab=NULL;
	  ResClassif = new int  [NClas];
	  for(i = 0; i < NClas; i ++) ResClassif[i] = 0;
	  // set Lab pointer on class labels
	  if (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz()))
	    {
	      if ((Out[NumS]->Def) != NULL)
		{
		  DEFUZ_SugenoClassif * D = (DEFUZ_SugenoClassif *) Out[NumS]->Def;
		  Lab=D->ClasLabel();
		}
	    }

	  else if (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))
	    {
	      if ((Out[NumS]->Def) != NULL)
		{
		  DEFUZ_MaxCrisp * D = (DEFUZ_MaxCrisp *) Out[NumS]->Def;
		  Lab=D->ClasLabel();
		}
	    }
	}
      else
    	  throw std::runtime_error("error in ResClassifAlloc:  classification case and no classes!");

    }
  return 0;
}

//*******************        PERFORMANCE    *********************************
double FIS::Performance(int NumS, const char * fdata, double &Couvert, double & MaxErr,
			double MuSeuil, const char * fres, FILE *display)
  //*****************************************************************************
  // Return : -2 when not enough input variables or not active output
  //          -1 when the observed output is not part of the sample file
  //          the number of misclassified items when classif
  //          the root mean square error
{
  if( (NumS < 0) || (NumS > NbOut - 1) || (! Out[NumS]->IsActive()) )
    {
      Couvert = 0.;
      sprintf( ErrorMsg, "~InvalidOutputNumber~: %d~", NumS);
      throw std::runtime_error( ErrorMsg );
    }

  double ** Data, err;
  int i, nbcol, classif, NbEx, ref;
  int * MisClassified=NULL;
  double * Lab=NULL;

  FILE * f;

  Data = NULL;
  nbcol = 0;
  MaxErr = 0.;
  err = -1.;
  Couvert = 0.;

  f = NULL;
  if(fres)
    if((f = fopen(fres, "wt")) == NULL)
      {
	sprintf( ErrorMsg, "~CannotOpenResultFile~: %.100s~", fres);
	throw std::runtime_error( ErrorMsg );
      }
  //read data -allocate Data array
  Data = ReadSampleFile(fdata, nbcol, NbEx);

  if(nbcol < NbIn)  return -2;

  ref = true;
  if(nbcol < NbIn + 1 + NumS) ref = false; // Observed output is not part of the sample file

  // Write a comment line in the result file
  WriteHeader(NumS, f, ref);

  // following two calls : always in the same order
  //allocate Classes array in crisp classification case
  ClassifCheck(Data, NbEx, NumS);
  // allocate  MisClassified and lab arrays in crisp classification case
  ResClassifAlloc(MisClassified, Lab, NumS);

  classif = Out[NumS]->Classification() &&
    (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType())) &&
    ( (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz()))	 ||
      (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz())));

  if(display && classif) fprintf(display, "\nThis is a classification case\n");


  //call to Perf function
  err = Perf(NumS, Data, NbEx, Couvert, MaxErr, MuSeuil, MisClassified, Lab, ref, f, display);

  if(f) fclose(f);

  if(display)
    {
      fprintf(display, "\n");
      if(ref)
	{
	  if(classif)
	    {
	      fprintf(display, "Number of misclassified items : %6d, in percentage %2d %%. \nDetail by classes : ",
		     (int) err, (int) (100. * err / NbEx));

	      for(i = 0; i < Out[NumS]->Def->NbClasses(); i ++) fprintf(display, "%6d ", MisClassified[i]);
	      fprintf(display, "\n");
	    }

	  else fprintf(display, "Mean square error: %11.2f\n", err);
	}
      int RuleWgt = false;
      for(i = 0; i < NbRules; i++)
	if(fabs(Rule[i]->GetExpertWeight() -1.0) > EPSILON) 
	  {
	    RuleWgt = true;
	    break;
	  }
      if(RuleWgt) fprintf(display, "\nWarning:  the rules are weighted.\n");
    }
  // delete 2D data array
  for (i = 0; i < NbEx; i++)
    {
      if (Data[i] != NULL) delete [] Data[i];
    }
  if (Data != NULL) delete [] Data;
  // delete MisClassified array
  if(MisClassified != NULL) delete [] MisClassified;

  return err;
} // End of Performance(char *)



  
double FIS::Perf(int NumS, double ** Data, int NbEx, double & Couvert, double & MaxErr,
		 double MuSeuil, int * MisClassified, double * Lab, int ref, FILE * f, FILE *display)
  //***********************************************************************************************
  // Measures the FIS performance du Sif for the NumS output for the data set fdata
  // Returns the number of misclassified for a class output
  // otherwise the error index
  // Also writes the MisClassified array which must be preallocated
  // if the output is a class
  // the Lab array contains labels for classes, must be preallocated
  // modif March 15 2009 - implicative inference case
  // flag impli=1 if implicative case, 0 otherwise
{
  if( (NumS < 0) || (NumS > NbOut - 1) || (! Out[NumS]->IsActive()) )
    {
      Couvert = 0.;
      sprintf( ErrorMsg, "~InvalidOutputNumber~: %d~", NumS);
      throw std::runtime_error( ErrorMsg );
    }

  int i, l, classif, diag, nPos, ValBlanc, NClas, NbBlancs, impli,nbOutMf;
  int countBl=0;
  double * muOut=NULL, * muOutInf=NULL;
  double diff, abserr, err, Mu, symatch=0.0,mulk=0.0,murk=0.0,lk=0.0,rk=0.0;
  double muThreshInf=MUTHIMPLI,muThreshObs=MUTHIMPLI,maxp=1.0;

  double * Pos, tmp;

  classif = (Out[NumS]->Classification() &&
	     (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType())) &&
	     ( (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz())) ||
	       (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))) );

  diag = ( (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType()))
	   && (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))
	   && (! Out[NumS]->Classification()) );

  impli= ! strcmp(Out[NumS]->Defuzzify(),OUT_FUZZY::ImpFuzzyDefuz());
  if (display)
    {
      if (impli) 
	fprintf(display, "\nimplicative inference");
      else
	fprintf(display, "\nconjunctive inference");
    }

  Pos = Out[NumS]->GetPossibles();
  nPos = Out[NumS]->GetNbPossibles();
  NClas = Out[NumS]->Def->NbClasses();

  diff = 0.;
  abserr = 0.;
  NbBlancs = 0;
  MaxErr = 0.;

  if(classif)
    for(i = 0; i < NClas; i++) MisClassified[i] = 0;

  if(! ref && (classif || diag))
    {
      if(! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz()))
	{
	  DEFUZ_SugenoClassif * D = (DEFUZ_SugenoClassif *) Out[NumS]->Def;
	  D->SetClasses(Pos, nPos);
	  // modif April 4, 2008 - assign Lab ptr to Classes + update Nclas
	  Lab=D->ClasLabel();
	  NClas=nPos;
	}
      else if(! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))
	{
	  DEFUZ_MaxCrisp * D = (DEFUZ_MaxCrisp *) Out[NumS]->Def;
	  D->SetClasses(Pos, nPos);
	  Lab=D->ClasLabel();
	  // modif April 4, 2008 - assign Lab ptr to Classes + update Nclas
	  NClas=nPos;
	}


    }
  // implicative case - write observed output value symb. labels
  nbOutMf=Out[NumS]->GetNbMf();
  if (impli)
    {
      muOut=new double[nbOutMf];
      muOutInf=new double[nbOutMf];
      for(int imf = 0; imf < nbOutMf; imf++)
	{
	  muOut[imf]=0.0;
	  muOutInf[imf]=0.0;
	}
    }
 

  // inference loop on data
  for(i = 0; i < NbEx; i ++)
    {
      Out[NumS]->MfGlob=NULL;
      OutErr[NumS] = 0.0;
      if(ref && f) fprintf(f, FORMAT_DOUBLE, Data[i][NbIn + NumS]);
      // implicative case
      if ( impli && ref )
	{
	  // get MF membership of observed output value and use Museuil for label membership

	  for(int imf = 0; imf < nbOutMf; imf++)
	    {
	      muOut[imf]=Out[NumS]->GetADeg(imf,Data[i][NbIn + NumS]);
	      if (f)
		fprintf(f, FORMAT_DOUBLE,muOut[imf]);
	    }
	}

      ValBlanc = 0;
      // same Infer function in all cases
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Mu = Infer(Data[i], NumS, f, display);
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(display) fprintf(display, "\nMax rule weight=%g",Mu);

      // blank calculation is more complex if implicative inference: blank example if >2 output MFs covered
      // plus inferred output sym. memberships are written in res. file

      if(Mu <= MuSeuil)
	{
	  NbBlancs++;
	  ValBlanc = 1;
	}
      else
	//implicative case: blank examples + sym. membership computation  for inferred d. poss
	if (impli)
	  {

	    // if MfGlob (inferred output poss. distrib. ) intersects more than 2 MFs, example is blank
	    countBl=0;
	    // compute MFGlob kernel
	    if (Out[NumS]->MfGlob !=NULL)
	      {
		Out[NumS]->MfGlob->AlphaKernel(lk,rk,maxp);
		for(int imf = 0; imf < nbOutMf; imf++)
		  {
		    mulk=Out[NumS]->GetADeg(imf,lk);
		    murk=Out[NumS]->GetADeg(imf,rk);
		    // consider maximum membership to each mf
		    muOutInf[imf]=mulk;
		    if (murk>mulk)
		      muOutInf[imf]=murk;

		    if (muOutInf[imf]>muThreshInf)
		      countBl++;
		  }
		if(countBl>2)
		  {
		    ValBlanc = 1;
		    NbBlancs++;
		  }
	      }// end case output dposs not empty
	    // case output dposs empty - what to do ????blank or not blank...
	    // end implicative case for blank cases
	  }
      //implicative case for symbolic match
      if (impli && ref)

	{
	  symatch=0.0;

	  if (Out[NumS]->MfGlob !=NULL)
	    {
	      symatch=((OUT_FUZZY *)Out[NumS])->SymbMatch(muOut,muOutInf,nbOutMf,muThreshObs,muThreshInf,display);
	      // for now symatch=-1->0
	      if (symatch<0) symatch=0;
	    }

	  if (f)
	    fprintf(f,FORMAT_DOUBLE,symatch);
	}// end implicative case for symbolic match

      // Special mode for OUT_CRISP, MaxCrisp defuz and classification
      if(diag)
	{
	  // Distance to the prototype (1 - MuInfer)
	  if(OutValue[NumS] == Data[i][NbIn + NumS])
	    {
	      for(l = 0; l < nPos; l++)
		if(OutValue[NumS] == Pos[l])
		  {
		    OutErr[NumS] = 1. - Out[NumS]->MuInfer[l];
		    break;
		  }
	    }
	  // Non ordered classes
	  else OutErr[NumS] = 1.;
	}

      
      else if(ref)
	OutErr[NumS] = OutValue[NumS] - Data[i][NbIn + NumS];

      if(display)
	{
	  fprintf(display, "\n Output  n %d  observed      inferred     delta    blank\n",  NumS+1);
	  if(ref)
	    fprintf(display, "             %8.3f      %8.3f     %8.3f", Data[i][NbIn + NumS], OutValue[NumS], OutErr[NumS]);
	  else
	    fprintf(display, "             %8.3f      %8.3f     %8.3f", FisMknan(),OutValue[NumS],FisMknan());

	  fprintf(display, "   %d", ValBlanc);
	}

      if((! ValBlanc) && ref)
	{
	  tmp = fabs(OutErr[NumS]);
	  if(tmp > fabs(MaxErr)) MaxErr = OutErr[NumS];

	  if(classif && tmp)
	    {
	      diff += 1.;
	      for(l = 0; l < NClas; l++)
		if(Data[i][NbIn + NumS] == Lab[l])
		  { MisClassified[l]++; break; }
	    }
	  else 
	    {
	      diff += (OutErr[NumS] * OutErr[NumS]);
	      abserr += tmp;
	    }
	}

      if(f)
	{
	  if(ref) fprintf(f, FORMAT_DOUBLE, OutErr[NumS]);
	  fprintf(f, "%5d  %8.3f\n", ValBlanc, diff);
	}
      // implicative case: delete mfconc for each rule (mfconc= new dposs for each item)
      if (impli)
	{
	  DeleteMFConc(NumS);
	  delete Out[NumS]->MfGlob;
	  Out[NumS]->MfGlob=NULL;
	}

    } // End of item loop

  if(classif)     err = diff;
  else if(diff) 
    {
      tmp = (double) (NbEx - NbBlancs);
      PIn = sqrt(diff) / tmp;
      RMSE = sqrt(diff / tmp);
      MAE = abserr / tmp;

      if(strErrorIndex != NULL)
	{
	  if(!strcmp(strErrorIndex, FIS::MaeErrorIndex())) err = MAE;
	  else if(!strcmp(strErrorIndex, FIS::PiErrorIndex())) err = PIn;
	  else err = RMSE;
	}
      else err = RMSE;
    }
  else err = 0.; // All blanks

  Couvert = 1. - ((double) NbBlancs / NbEx);

  delete [] muOut;
  delete [] muOutInf;
  if(!ref) return -1;
  return err;
} // End of Perf



void FIS::WriteHeader(int o, FILE *p, int obs) const
  //******************************************
{
  if(!p) return;
  if(obs)
    {
      fprintf(p, "    %s", OBSERVED);
      // implicative inference - result file has symbolic fields for observed values
      if (! strcmp(Out[o]->Defuzzify(),OUT_FUZZY::ImpFuzzyDefuz()))
	{
	  for(int i = 0; i < Out[o]->GetNbMf(); i++) fprintf(p, "      MF%d", i+1);
	}
    }
  // implicative inference - result file has symbolic fields for inferred values
  Out[o]->Def->WriteHeader(p, Out[o]);

  if(obs) fprintf(p, "    %s", FIS_ERROR);
  fprintf(p, "    %s", BLANK);
  fprintf(p, "    %s", CUMUL);
  fprintf(p, "\n");
}


//**********************         GETBREAKPOINTS          *********************

int FIS::GetBreakPoints(char * archive, int NbMax)
  //**********************************************
{
  int  *NbBp, *cBp;
  double ** Bp;
  int i, n;
  FILE * f;

  n = 1;
  for(i = 0; i < NbIn; i++)
    n *= (In[i]->GetNbMf() * 2 - 1);
  if(n > NbMax) return n;

  NbBp = new int [NbIn];       // Number of breakpoints for each input
  Bp = new double * [NbIn];    // Breakpoint values for each input
  cBp = new int [NbIn];        // Current index for each input

  for(i = 0; i < NbIn; i++)
    {
      In[i]->GetBreakPoints(Bp[i], NbBp[i]);
      cBp[i] = 0;
    }
  f = fopen(archive, "wt");
  if(f != NULL)
	  GenereCombi(0, f, NbBp, cBp, Bp);
  fclose(f);

  for(i = 0; i < NbIn; i++) delete [] Bp[i];
  delete [] Bp; delete [] cBp; delete [] NbBp;

  return 0;
}

void FIS::GenereCombi(int i, FILE *f, int  *NbBp, int *cBp, double ** Bp)
  //*********************************************************************
{
  int m;

  if(i == NbIn - 1)    // The last input
    {
      for(m = 0; m < NbBp[i]; m++)
	{
	  cBp[i] = m;
	  PrintBreakPoints(f, cBp, Bp);
	}
    }
  else
    {
      for(m = 0; m < NbBp[i]; m++)
	{
	  cBp[i] = m;
	  GenereCombi(i+1, f, NbBp, cBp, Bp);
	}
    }
} // End of GenereCombi()


void FIS::PrintBreakPoints(FILE *f, int *cBp, double ** Bp) const
  //*******************************************************
{
  int i;
  for(i = 0; i < NbIn; i++)
    {
      fprintf(f, FORMAT_DOUBLE, Bp[i][cBp[i]]);
      if(i == NbIn - 1) fprintf(f, "\n");
      else fprintf(f, "%c", SEPARE);
    }
}

//**********************         Componenent management  ********************


void FIS::SetName( const char *name )
  //*********************************
{
  delete [] Name;
  Name = new char[strlen(name)+1];
  sprintf( Name, "%s", name );
}



void FIS::SetConjunction( const char *conjunction )
  //***********************************************
{
  delete [] cConjunction;
  cConjunction = new char[strlen(conjunction)+1];
  sprintf( cConjunction, "%s", conjunction );

  if( Rule == NULL )    return;
  int *facteurs = NULL;
  try
    {
      facteurs = new int[NbIn];
      for( int i=0 ; i<NbRules ; i++ )
	{
	  Rule[i]->GetProps( facteurs );
	  Rule[i]->SetPremise( NbIn, In, cConjunction );
	  Rule[i]->SetAProps( facteurs );
	}
      delete [] facteurs;
    }
  catch( std::exception &e )
    {
      delete [] facteurs;
      throw;
    }
}


void FIS::SetMissingValues( const char *missing_values )
  //****************************************************
{
  delete [] strMissingValues;
  strMissingValues = new char[strlen(missing_values)+1];
  sprintf( strMissingValues, "%s", missing_values );
}

void FIS::SetErrorIndex( const char *index )
  //****************************************
{
  delete [] strErrorIndex;
  strErrorIndex = new char[strlen(index)+1];
  sprintf( strErrorIndex, "%s", index);
}


void FIS::AddInput( FISIN *entree )
  //*******************************
{
  FISIN **temp = new FISIN * [NbIn];
  int i;

  // Add input in In array
  for( i=0 ; i<NbIn ; i++ ) temp[i] = In[i];
  NbIn++;
  if (In !=NULL) delete [] In;

  In = new FISIN * [NbIn];
  for( i=0 ; i<NbIn - 1 ; i++ )
    In[i] = temp[i];
  In[NbIn-1] = entree;

  if (temp!=NULL) delete [] temp;

  // Add input in Rules
  int *facteurs = NULL;
  try
    {
      facteurs = new int[NbIn];
      for( i=0 ; i<NbRules ; i++ )
	{
	  Rule[i]->GetProps( facteurs );
	  facteurs[NbIn-1] = 0;
	  Rule[i]->SetPremise( NbIn, In, cConjunction );
	  Rule[i]->SetAProps( facteurs );
	}
      delete [] facteurs;
    }
  catch( std::exception &e )
    {
      delete [] facteurs;
      throw;
    }
}

void FIS::RemoveInput( int input_number )
  //*************************************
{
  // Remove input from In array
  FISIN **temp = new FISIN * [NbIn-1];
  delete In[input_number];

  for( int i=0, j=0 ; i<NbIn ; i++ )
    if( i != input_number )
      {
	temp[j] = In[i];
	j++;
      }
  NbIn--;
  delete [] In;
  In = temp;

  // remove input from rule descriptions
  int *facteurs = NULL;
  int *new_facteurs = NULL;
  try
    {
      facteurs = new int[NbIn+1];
      new_facteurs = new int[NbIn];
      for( int i=0 ; i<NbRules ; i++ )
	{
	  Rule[i]->GetProps( facteurs );
	  Rule[i]->SetPremise( NbIn, In, cConjunction );
	  for( int k=0, j=0 ; k<NbIn+1 ; k++ )
	    if( k != input_number )
	      {
		new_facteurs[j] = facteurs[k];
		j++;
	      }
	  Rule[i]->SetAProps( new_facteurs );
	}
      delete [] facteurs;
      delete [] new_facteurs;
    }
  catch( std::exception &e )
    {
      delete [] facteurs;
      delete [] new_facteurs;
      throw;
    }
}


void FIS::ReplaceInput( int input_number, FISIN *new_input )
  //********************************************************
{
  // removes the premise corresponding to a removed MF
  int prop;
  for( int i=0 ; i<NbRules ; i++ )
    if( Rule[i]->GetAProp( prop, input_number ) > new_input->GetNbMf() )
      Rule[i]->SetAProp( 0, input_number );


  // input replacement
  delete In[input_number];
  In[input_number] = new_input;
}

void FIS::DeleteMFConc(int output_number) const
  //********************************************************
{
  //delete MFConc array
  if  (Out[output_number]->MfConc !=NULL)
    {
      for (int irule=0;irule<NbRules;irule++)
	{
	  delete Out[output_number]->MfConc[irule];
	  Out[output_number]->MfConc[irule]=NULL;
	}
    }
}
void FIS::DeleteMFConcArray(int output_number)
  //********************************************************
{
  delete [] Out[output_number]->MfConc;
  Out[output_number]->MfConc=NULL;
}

void FIS::ReplaceOutput( int output_number, FISOUT *new_output )
  //************************************************************
{
  //PROTECTION
  if (output_number<0 || output_number>NbOut) return;
  //
  // check output MF types for implicative output-throw exception
  new_output->CheckImpliMFs();
  // the output type has changed
  if( strcmp(Out[output_number]->GetOutputType(), new_output->GetOutputType()) != 0 )
    for( int i=0 ; i<NbRules ; i++ )
      Rule[i]->SetAConc( output_number, 1 );
  else
    // if the new output is fuzzy
    if( strcmp(new_output->GetOutputType(), OUT_FUZZY::OutputType()) == 0 )
      // removes actions corresponding to a removed MF
      for( int i=0 ; i<NbRules ; i++ )
	if( (int)Rule[i]->GetAConc(output_number) > new_output->GetNbMf() )
	  Rule[i]->SetAConc( output_number, 1 );
  // free MFConc arrays
  DeleteMFConc(output_number);
  DeleteMFConcArray(output_number);
  // output replacement
  // OUTPUT is removed, together with its aggreg and defuz objects
  delete Out[output_number];
  Out[output_number] = new_output;
  // create Possibles, Classes, MFConc array for new output
  Out[output_number]->InitPossibles(Rule,NbRules,output_number);
}




void FIS::RemoveMFInInput( int input_number, int mf_number )
  //********************************************************
{
  //PROTECTION
  if ((input_number<0) || (input_number>=NbIn)) return;
  if ((mf_number<0) || (mf_number>In[input_number]->GetNbMf())) return;
  // retire le SEF de l'entr�e
  In[input_number]->RemoveMF( mf_number );
  // retire le SEF des r�gles et renum�rote les num�ros de SEF sup�rieur � mf_number
  mf_number++;
  for( int i=0 ; i<NbRules ; i++ )
    {
      int prop=0;
      Rule[i]->GetAProp( prop, input_number );
      if( prop == mf_number )
	Rule[i]->SetAProp( 0, input_number );
      if( prop > mf_number )
	Rule[i]->SetAProp( prop-1, input_number );
    }
}


void FIS::RemoveMFInOutput( int output_number, int mf_number )
  //**********************************************************
{
  //PROTECTION
  if ((output_number<0) || (output_number>=NbOut)) return;
  if ((mf_number<0) || (mf_number>Out[output_number]->GetNbMf())) return;
  // retire le SEF de la sortie
  Out[output_number]->RemoveMF( mf_number );
  // retire le SEF des r�gles et renum�rote les num�ros de SEF sup�rieur � mf_number
  mf_number++;
  for( int i=0 ; i<NbRules ; i++ )
    {
      int conc = (int)Rule[i]->GetAConc( output_number );
      if( conc == mf_number )
	Rule[i]->SetAConc( output_number, 1 );
      if( conc > mf_number )
	Rule[i]->SetAConc( output_number, conc-1 );
    }
  // free MFConc array
  DeleteMFConc(output_number);
  //DeleteMFConcArray(output_number);// done in InitPossibles
  // delete and recreate Possible and Classes arrays, MFConc array
  Out[output_number] -> InitPossibles(Rule, NbRules,output_number);

}


void FIS::AddOutput( FISOUT *sortie )
  //*********************************
{
  FISOUT **temp =NULL;
  // Add output in Out array
  if (NbOut>0)
    temp = new FISOUT * [NbOut];

  int i, j;

  for(i=0 ; i<NbOut ; i++)  temp[i] = Out[i];

  NbOut++;
  if (Out !=NULL) delete [] Out;
  if (OutValue != NULL) delete [] OutValue;
  if (OutErr != NULL) delete [] OutErr;

  Out=NULL;
  Out=new FISOUT * [NbOut];
  for(i=0 ; i<(NbOut-1) ; i++)
    Out[i]=temp[i];

  if (temp!=NULL) delete [] temp;

  // ckeck MF types-will throw exception if implicative output and MF types not allowed
  sortie->CheckImpliMFs();
  //
  Out[NbOut-1] = sortie;

  OutValue = new double [NbOut];
  OutErr = new double [NbOut];

  // Add output in rules
  double *actions = NULL;
  try
    {
      actions = new double[NbOut];
      for(i = 0 ; i < NbRules ; i++ )
	{
	  for(j = 0 ; j < NbOut-1 ; j++ )
	    actions[j] = Rule[i]->GetAConc( j );
	  actions[NbOut-1] = 1;
	  Rule[i]->SetConclusion( NbOut, Out );
	  Rule[i]->SetConcs( actions );
	}
      if (actions != NULL) delete [] actions;
      // delete and recreate Possible, Classes and MFConc arrays for output
      // which has just been added
      for(i = 0 ; i < NbOut ; i++)
	Out[i]->InitPossibles(Rule, NbRules, i);
    }
  catch( std::exception &e )
    {
      if (actions != NULL) delete [] actions;
      throw;
    }
}



void FIS::RemoveOutput( int output_number )
  //***************************************
{
  // normally nothing to be done for Possible and Classes arrays
  FISOUT **temp = NULL;
  int i,j,k;

  try
    {
      if ((output_number<0) || (output_number>NbOut)) return;
	temp = new FISOUT * [NbOut-1];

      for( i=0, j=0 ; i<NbOut ; i++ )
	if( i != output_number )
	  {
	    temp[j] = Out[i];
	    j++;
	  }
      // free MFConc array for output to be removed
      DeleteMFConc(output_number);
      DeleteMFConcArray(output_number);
      //
      if (Out[output_number]!=NULL) delete Out[output_number];

      NbOut--;
      if (Out != NULL) delete [] Out;
      if (OutValue != NULL) delete [] OutValue;
      OutValue = NULL;
      if (OutErr != NULL)  delete [] OutErr;
      OutErr = NULL;

      Out=NULL;
      if (NbOut>0)
	{
	  Out=new FISOUT * [NbOut];
	  for(i=0 ; i<NbOut ; i++)
	    Out[i]=temp[i];
	  OutValue = new double [NbOut];
	  OutErr = new double [NbOut];
	}
      if (temp!=NULL) delete [] temp;

    }
  catch( std::exception &e )
    {
      if (temp!=NULL) delete [] temp;
      throw;
    }

  // remove output from rules
  double *actions = NULL;
  double *new_actions = NULL;


  try
    {
      actions = new double[NbOut+1];
      if (NbOut>0)
	new_actions = new double[NbOut];
      for( i=0 ; i<NbRules ; i++ )
	{
	  for( j=0 ; j<NbOut+1 ; j++ )
	    actions[j] = Rule[i]->GetAConc( j );
	  Rule[i]->SetConclusion( NbOut, Out );
	  for( k=0, j=0 ; k<NbOut+1 ; k++ )
	    if( k != output_number )
	      {
		new_actions[j] = actions[k];
		j++;
	      }
	  Rule[i]->SetConcs( new_actions );
	}
      delete [] actions;
      delete [] new_actions;
      for(i = 0 ; i < NbOut ; i++)
	{
	  DeleteMFConc(i);
	  Out[i]->InitPossibles(Rule, NbRules, i);
	}
    }
  catch( std::exception &e )
    {
      if (actions != NULL) delete [] actions;
      if (new_actions != NULL) delete [] new_actions;
      throw;
    }

}



void FIS::AddRule( RULE *regle )
  //****************************
{
  int i;
  RULE **temp = new RULE * [NbRules+1];
  for( i=0 ; i<NbRules ; i++ )
    temp[i] = Rule[i];
  temp[NbRules] = regle;

  // free MFConc arrays for all outputs
  // bug correction - September 22, 2009 - first delete mfconcs with former rule #
  // need to call function DeleteMFConc in class FISOUT
  // and set the ptrs to NULL
  for( i=0 ; i<NbOut;i++)
    {
      DeleteMFConc(i);
      DeleteMFConcArray(i);
    }
  //
  NbRules++;
  //
  delete [] Rule;
  Rule = NULL;
  Rule = new RULE * [NbRules];
  // clone all rules
  for( i=0 ; i<NbRules ; i++ )
    Rule[i]= new RULE(*temp[i], In, Out);
  for( i=0 ; i<(NbRules-1) ; i++ )
    if (temp[i] != NULL) delete temp[i];
  if (temp != NULL) delete [] temp;
  
  // recreate arrays for rule conclusions
  for(i = 0 ; i < NbOut ; i++)
    Out[i]->InitPossibles(Rule, NbRules, i);
  //
  ComputeNbActRule();

}


void FIS::RemoveAllRules(void)
  //*************************
{
  int i;

  for( i=0 ; i<NbRules ; i++)
    if (Rule[i]!= NULL)  delete Rule[i];

  for( i=0 ; i<NbOut;i++)
    {
      DeleteMFConc(i);
      DeleteMFConcArray(i);
    }

  NbRules = 0;

}

void FIS::RemoveRule( int rule_number )
  //***********************************
{
  int i,j;
  // create temporary array for rules
  RULE **temp = NULL;

  //PROTECTION
  if ((rule_number<0) || (rule_number>NbRules)) return;
  if (NbRules>1)
    temp=new RULE * [NbRules-1];
  // transfer other rules into temp. array
  for( i=0, j=0 ; i<NbRules ; i++ )
    if( i != rule_number )
      {
	temp[j] = Rule[i];
	j++;
      }
  // delete mfonc array
  for( i=0 ; i<NbOut;i++)
    {
      DeleteMFConc(i);
      DeleteMFConcArray(i);
    }
  
  NbRules--;

  if (Rule != NULL) delete Rule[rule_number];
  if (Rule != NULL) delete [] Rule;
  Rule = NULL;
  Rule = new RULE * [NbRules];
  // clone all rules
  for( i=0 ; i<NbRules ; i++ )
    Rule[i]= new RULE(*temp[i], In, Out);
  //delete temp = also delete former rules
  for( i=0 ; i<NbRules ; i++ )
    if (temp[i] != NULL) delete temp[i];
  if (temp != NULL) delete [] temp;
  // update active rule number
  // delete removed rule
  ComputeNbActRule();

  // recreate arrays for rule conclusions
  for(i = 0 ; i < NbOut ; i++)
    Out[i]->InitPossibles(Rule, NbRules, i);

}

void FIS::Crisp2Fuz(int o, const char * DefuzType, double * c, int nc)
  //******************************************************************
{
  double *Val, *ValRange;
  int n, n_range, new_conc, i, j;
  FISOUT * newout;
  double Rmin, Rmax, vDef;

  Rmin = Out[o]->min();
  Rmax = Out[o]->max();
  vDef = Out[o]->DefaultValue();
  n = n_range = 0;

  if((o < 0) || (o > NbOut-1)) return;

  if(! strcmp(Out[o]->GetOutputType(), OUT_FUZZY::OutputType()))
    return;

  if(c == NULL)
    {
      Out[o]->InitPossibles(Rule, NbRules,o);
      Val = Out[o]->GetPossibles();
      n = Out[o]->GetNbPossibles();
    }
  else
    {
      if(NbRules>=1)
	{
	  sprintf(ErrorMsg, "~NbRules=~%d~in~Crisp2Fuz~function~incompatible~with~c~array\n~",NbRules);
	  throw std::runtime_error( ErrorMsg );
	}
      Val = c;
      n = nc;
    }
  // keeps the possible values within the output range

  if(n > MAX_MF)
    {
      sprintf(ErrorMsg, "~TooManyMFs~%d~ForOutput~%d~MaxAllowed~%d \n", n, o+1, MAX_MF);
      throw std::runtime_error( ErrorMsg );
    }

  ValRange = NULL;
  if (n >= 0) ValRange = new double[n];
  for(i = 0 ; i < n  ; i++ )
    if( (Val[i] >= Rmin) && (Val[i] <= Rmax) )
      {
	ValRange[n_range] = Val[i];
	n_range++;
      }
  newout = new OUT_FUZZY(ValRange, n_range, Rmin, Rmax, 1, DefuzType, OUT_FUZZY::DisjSum(), 0, vDef);

  newout->SetName(Out[o]->Name);
  newout->Classification(Out[o]->Classification());
  delete Out[o];
  Out[o] = newout;

  for(i = 0 ; i < NbRules ; i++ )
    {
      new_conc = 1;
      for(j = 0; j < n_range ; j++ )
	{
	  if( fabs(Out[o]->GetMF(j)->Corner() - Rule[i]->GetAConc(o)) < EPSILON )
	    new_conc = j+1;
	}
      Rule[i]->SetAConc(o, new_conc);
    }

  // delete and recreate Possible and Classes arrays
  Out[o] -> InitPossibles(Rule, NbRules,o);
  if (ValRange != NULL) delete [] ValRange;
}


void FIS::Fuz2Crisp(int o)
  //**********************
{
  if((o < 0) || (o > NbOut-1)) return;

  if(! strcmp(Out[o]->GetOutputType(), OUT_CRISP::OutputType()))
    return;

  double l, r, *newconc;
  int i, n;

  n = Out[o]->GetNbMf();
  newconc = new double [n];

  for(i = 0 ; i < n; i++ )
    newconc[i] = Out[o]->GetMF(i)->Corner();

  l = Out[o]->min();
  r = Out[o]->max();
  FISOUT *newout = new OUT_CRISP();
  newout->SetName(Out[o]->Name);
  newout->Classification(Out[o]->Classification());

  Out[o]->DeleteMFConc(NbRules);
  Out[o]->DeleteMFConcArray();
  delete Out[o]->MfGlob;
  delete Out[o];
  Out[o] = newout;
  Out[o]->SetRange(l, r);

  for(i = 0 ; i < NbRules ; i++ )
    {
      int index = (int) Rule[i]->GetAConc(o) - 1;
      Rule[i]->SetAConc(o, (index >= 0) ? newconc[index] : 0 );
    }
  // delete and recreate Possible and Classes arrays
  Out[o] -> InitPossibles(Rule, NbRules,o);

  delete [] newconc;
}


//*******************************************************************
//
//                      Config file reading functions
//
//*******************************************************************


void FIS::ReadHdr(ifstream & f, int bufsize) //throw(std::runtime_error)
  //***************************************
  // Returns 0 if everything OK otherwise 3
{
  char *tmp=NULL, *buf=NULL;

  try
    {
      tmp = new char[bufsize];
      buf = new char[bufsize];

      do{
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // skip empty lines and ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf( tmp, "[Interface]" );
      if( !strncmp(tmp, buf, strlen(tmp)) ) //Obsolete config file
	{
	  f.getline(buf, bufsize);
	  do{                           // skip empty lines
	    f.getline(buf, bufsize);
	  }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
	  // ignores the ^M of MS-DOS files read in Linux-Unix
	  // ignores comment lines beginning with # or %
	}

      sprintf( tmp, "[System]" );
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
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
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if(SearchStr(buf, tmp))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~~StringSeparatorNotFoundInString~: %.50s~", buf );
	  throw std::runtime_error( ErrorMsg );
	}
      SetName(tmp);

      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf( tmp, "Ninputs=" );
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      NbIn = atoi(buf + strlen(tmp));
      if(NbIn < 0)
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~~InvalidNumberOfInputs~: %-3d~", NbIn );
	  NbIn = 0;
	  throw std::runtime_error( ErrorMsg );
	}


      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf( tmp, "Noutputs=" );
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      NbOut = atoi(buf + strlen(tmp));
      if(NbOut < 0)
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~~InvalidNumberOfOutputs~: %-3d~", NbOut );
	  NbOut = 0;
	  throw std::runtime_error( ErrorMsg );
	}


      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %


      sprintf( tmp, "Nrules=" );
      if( strncmp(tmp, buf,  strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      NbRules = atoi(buf +  strlen(tmp));
      if(NbRules < 0)
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~~InvalidNumberOfRules~: %-3d~", NbRules );
	  NbRules = 0;
	  throw std::runtime_error( ErrorMsg );
	}


      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %
      sprintf( tmp, "Nexceptions=" );
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      NbExceptions = atoi(buf + strlen(tmp));




      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf( tmp, "Conjunction=" );
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if(SearchStr(buf, tmp))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~~StringSeparatorNotFoundInString~: %.50s~", buf );
	  throw std::runtime_error( ErrorMsg );
	}
      SetConjunction(tmp);


      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf(tmp,"MissingValues=");
      if( strncmp(tmp, buf,  strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if(SearchStr(buf , tmp))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~~StringSeparatorNotFoundInString~: %.50s~", buf );
	  throw std::runtime_error( ErrorMsg );
	}
      SetMissingValues(tmp);
      delete [] tmp;
      delete [] buf;
    }
  catch( std::exception &e )
    {
      delete [] tmp;
      delete [] buf;
      throw;
    }
} // end of ReadHdr()



void FIS::ReadIn(ifstream &f, int bufsize, int num) //throw(std::runtime_error)
  //***********************************************
{
  char *tmp=NULL, *buf=NULL;

  try
    {
      tmp = new char[bufsize];
      buf = new char[bufsize];

      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %


      sprintf(tmp,"[Input%d]", num + 1);    // Tag
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}

      In[num] = new FISIN(f, bufsize, num + 1);     // i + 1 for error messages
      delete [] tmp;
      delete [] buf;
    }
  catch( std::exception &e )
    {
      delete [] tmp;
      delete [] buf;
      throw;
    }
}


void FIS::ReadOut(ifstream &f, int bufsize, int num, int Cover)
  //***********************************************************
  // Returns 0 if everything OK otherwise 5
{
  double vDefault;
  char *tmp=NULL, *buf=NULL, *nature=NULL, *cDefuz=NULL, *cDisj=NULL;
  int cclas;

  try
    {
      tmp = new char[bufsize];
      buf = new char[bufsize];
      nature = new char[bufsize];
      cDefuz = new char[bufsize];
      cDisj = new char[bufsize];

      do{                               // Skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %


      sprintf(tmp,"[Output%d]", num + 1);    // Tag
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}

      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf(tmp,"Nature=");
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", num+1, tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if(SearchStr(buf, nature))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~StringSeparatorNotFoundInString~: %.50s~", num+1, buf );
	  throw std::runtime_error( ErrorMsg );
	}


      // read the defuzzification type
      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23) || (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf(tmp,"Defuzzification=");
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", num+1, tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if(SearchStr(buf, cDefuz))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~StringSeparatorNotFoundInString~: %.50s~", num+1, buf );
	  throw std::runtime_error( ErrorMsg );
	}

      // read the disjunction type
      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf(tmp,"Disjunction=");
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", num+1, tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if(SearchStr(buf, cDisj))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~StringSeparatorNotFoundInString~: %.50s~", num+1, buf );
	  throw std::runtime_error( ErrorMsg );
	}

      // read the default value
      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf(tmp,"DefaultValue=");
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", num+1, tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      vDefault = strtod(buf + strlen(tmp), NULL);

      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %

      sprintf(tmp,"Classif=");
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", num+1, tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if(SearchStr(buf, tmp))
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~StringSeparatorNotFoundInString~: %.50s~", num+1, buf );
	  throw std::runtime_error( ErrorMsg );
	}
      if( ! strncmp(tmp, "no", 4) )  cclas = false;
      else  if( ! strncmp(tmp, "yes", 4) ) cclas = true;
      else
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~ExpectedString~: Classif=yes or no\n~ReadString~: %.50s~", num+1, tmp );
	  throw std::runtime_error( ErrorMsg );
	}

      if(!strcmp(nature, OUT_CRISP::OutputType()))
	Out[num] = new OUT_CRISP(f, bufsize, num + 1, cDefuz, cDisj, cclas, vDefault);
      // i+1 for error messages
      else if(!strcmp(nature, OUT_FUZZY::OutputType()))
	Out[num] = new OUT_FUZZY(f, bufsize, num + 1, cDefuz, cDisj, cclas, vDefault, Cover);
      else
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n~UnknownNature~:~%.50s~", num+1, nature );
	  throw std::runtime_error( ErrorMsg );
	}
      if (tmp !=NULL) delete [] tmp;
      if (buf != NULL) delete [] buf;
      if (nature != NULL) delete [] nature;
      if (cDefuz != NULL) delete [] cDefuz;
      if (cDisj != NULL) delete [] cDisj;
    }
  catch(std::runtime_error &e)
    {
      if (tmp !=NULL) delete [] tmp;
      if (buf != NULL) delete [] buf;
      if (nature != NULL) delete [] nature;
      if (cDefuz != NULL) delete [] cDefuz;
      if (cDisj != NULL) delete [] cDisj;
      throw;
    }
}


void FIS::ReadRules(ifstream &f, int bufsize) //throw(std::runtime_error)
  //*****************************************
  // Returns 0 if everything OK otherwise 6
{
  int i;
  char *tmp=NULL, *buf=NULL;

  try
    {
      tmp = new char[bufsize];
      buf = new char[bufsize];

      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %


      sprintf(tmp, "[Rules]");
      if( strncmp(tmp, buf, strlen(tmp)) )
	{
	  sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	  throw std::runtime_error( ErrorMsg );
	}

      if(! NbRules)
	{
	  delete [] tmp;
	  delete [] buf;
	  return;
	}

      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %
      if(SearchStr(buf, tmp))   // It is not a file name but values.
	{
	  Rule[0] = new RULE( NbIn, In, NbOut, Out, cConjunction, buf);
	  for(i = 1; i < NbRules; i++)
	    {
	      do{                           // skip empty lines
		f.getline(buf, bufsize);
	      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
	      // ignores the ^M of MS-DOS files read in Linux-Unix
	      // ignores comment lines beginning with # or %
	      Rule[i] = new RULE( NbIn, In, NbOut, Out, cConjunction, buf);
	    }
	}
      else
	{
	  ifstream freg(tmp);
	  if(!freg)
	    {
	      sprintf( ErrorMsg, "~ErrorInFISFile~\n~CannotOpenRulesFile~: %.100s~", tmp );
	      throw std::runtime_error( ErrorMsg );
	    }

	  bufsize = MaxLineSize( freg );
	  delete [] buf;
	  buf = new char[bufsize];

	  for(i = 0; i < NbRules; i++)
	    {
	      freg.getline(buf, bufsize);
	      Rule[i] = new RULE( NbIn, In, NbOut, Out, cConjunction, buf);
	    }
	}
      delete [] tmp;
      delete [] buf;
    }
  catch( std::exception &e )
    {
      delete [] tmp;
      delete [] buf;
      throw;
    }
}



void FIS::ReadExcep(ifstream &f, int bufsize)
  //*****************************************
{
  RULE * RegleTmp;
  int i, pos, depuis;
  char *tmp, *buf;


  tmp = new char[bufsize];
  buf = new char[bufsize];

  do{                           // skip empty lines
    f.getline(buf, bufsize);
  }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
  // ignores the ^M of MS-DOS files read in Linux-Unix
  // ignores comment lines beginning with # or %


  sprintf(tmp, "[Exceptions]");
  if( strncmp(tmp, buf, strlen(tmp)) )
    if( strncmp(tmp, buf, strlen(tmp)) )
      {
	sprintf( ErrorMsg, "~ErrorInFISFile~\n~ExpectedString~: %.50s\n~ReadString~: %.50s~", tmp, buf );
	throw std::runtime_error( ErrorMsg );
      }

  for(i = 0; i < NbExceptions; i++)
    {
      do{                           // skip empty lines
	f.getline(buf, bufsize);
      }while((strlen(buf) == 0) || (buf[0]==0x0D) || (buf[0]==0x23)|| (buf[0]==0x25));
      // ignores the ^M of MS-DOS files read in Linux-Unix
      // ignores comment lines beginning with # or %
      RegleTmp = new RULE( NbIn, In, NbOut, Out, cConjunction, buf);
      pos = depuis = 0;
      while( (pos = RulePos(RegleTmp, depuis)) != -1)
	{
	  Rule[pos]->Deactivate();
	  depuis = pos + 1;
	};
      delete RegleTmp;
    }

  delete [] tmp;
  delete [] buf;
}


//*******************************************************************
//
//                      Print functions
//
//*******************************************************************
void FIS::Print(FILE *f) const
  //********************
{
  int i;

  fprintf(f, "\nSystem : %s", Name);
  fprintf(f, "\nNumber of Inputs: %d\tNumber of outputs : %d\n", NbIn, NbOut);

  fprintf(f, "\nNumber of rules : %d\tNumber of exceptions : %d", NbRules, NbExceptions) ;

  fprintf(f, "\nConjunction : %s", cConjunction);
  fprintf(f, "\nMissing values handling, membership : %s\n", strMissingValues);

  for(i = 0; i < NbIn; i++)   In[i]->Print(f);
  for(i = 0; i < NbOut; i++)  Out[i]->Print(f);

  fprintf(f, "\nRules : \n");
  if(NbRules < 30)
    for(i = 0; i < NbRules; i++)   Rule[i]->Print(f);

  else
    {
      char * fich;
      fich = new char[strlen(Name) + 10];
      sprintf(fich, "%s.rules", Name);
      fprintf(f, "\nsee file %s\n", fich);
      FILE * g;
      g = fopen(fich, "wt");
      if( ! g)
	{
      sprintf(ErrorMsg, "\nFile opening failed: %s\n", fich);
      throw std::runtime_error( ErrorMsg );
	}
      for(i = 0; i < NbRules; i++)  Rule[i]->Print(g);

      delete [] fich;
    }
}


void FIS::PrintCfg(FILE *f, const char *fd) const
  //***************************************
{
  int i, j;
  bool ExpWgt = false;

  j = 0;
  for(i = 0; i < NbRules; i++)
    {
      if(Rule[i]->IsActive()) j++;
      if(fabs(Rule[i]->GetExpertWeight() -1.0) > EPSILON) ExpWgt = true;
    }
  fprintf(f, "[System]\n");
  fprintf(f, "Name=%c%s%c\n", STRING_SEP, Name, STRING_SEP);
  fprintf(f, "Ninputs=%d\n", NbIn);
  fprintf(f, "Noutputs=%d\n", NbOut);
  fprintf(f, "Nrules=%d\n", j);
  fprintf(f, "Nexceptions=0\n");
  fprintf(f, "Conjunction=%c%s%c\n", STRING_SEP, cConjunction, STRING_SEP);
  fprintf(f, "MissingValues=%c%s%c\n", STRING_SEP, strMissingValues, STRING_SEP);

  for(i = 0; i < NbIn; i++)   In[i]->PrintCfg(i + 1, f, fd);

  for(i = 0; i < NbOut; i++)  Out[i]->PrintCfg(i + 1, f, fd);

  fprintf(f, "\n[Rules]\n");
  for(i = 0; i < NbRules; i++)  
    if(Rule[i]->IsActive()) Rule[i]->PrintCfg(f, fd, ExpWgt);

  fprintf(f, "\n[Exceptions]\n");
}


//*******************************************************************
//
//                      Rule base analysis
//
//*******************************************************************
int FIS::WriteHeaderPerfRB(int nout, FILE * f)
  //******************************************
{
  InfoRB t;
  int ret;
  ret = AnalyzeRB(t, nout);

  /*ret = AnalyzeRB(t, nout);

  ret = AnalyzeRB(t, nout);
  */

  if(ret) return ret;

  fprintf(f, " Name  &  PI  &   CI  &   maxE  & ");
  t.WriteHeader(f);

  return 0;
}


int FIS::PerfRB(double pi, double ci, double errmax, int nout, FILE * f)
  //********************************************************************
{
  int ret;
  InfoRB t;
  ret = AnalyzeRB(t, nout);
  if(ret) return ret;

  fprintf(f, "%s & %f & %f & %f & ", Name, pi, ci, errmax);
  t.Print(f);



  return 0;
}

int FIS::AnalyzeRB(InfoRB & i, int n, double ** Varray, int nb)
  //***********************************************************
{
  if(n > NbOut) return n;

  int p, j, k, classif, tmp, *val, *use;
  double conc;
  double * loclabels=NULL;

  InitClassLabels(Varray, nb);

  loclabels = NULL;
  val = use = NULL;
  classif = false;

  i.out = n;
  i.nIn = NbIn;
  i.nOut = NbOut;

  i.nMf = new int [NbIn + NbOut];

  // loclabels is set to the Def->Classes pointer without any memory allocation
  // but val is allocated within ResClassifAlloc
  ResClassifAlloc(val, loclabels, n);

  if(val)   // Non NULL when crisp output and classif
    {
      i.nClass = Out[n]->Def->NbClasses();
      i.labels = new double [i.nClass];
      for(p = 0; p < i.nClass; p++) i.labels[p] = loclabels[p];
      classif = true;
    }

  else if(! (strcmp(Out[n]->GetOutputType(), OUT_FUZZY::OutputType())))
    {
      i.nClass = Out[n]->GetNbMf();
      i.labels = new double [i.nClass];
      for(p = 0; p < i.nClass; p++) i.labels[p] = p+1;
      classif = true;
    }

  if(classif)
    {
      i.nRc = new int [i.nClass];
      for(p = 0; p < i.nClass; p++) i.nRc[p] = 0;
    }

  if (val) delete [] val;

  val = new int [NbIn];
  use = new int [NbIn];

  i.nR = NbRules;
  tmp = 0;
  for(p = 0; p < NbIn; p++) val[p] = use[p] =0;

  for(p = 0; p < NbRules; p++)
    {
      k = 0;
      Rule[p]->GetProps(val);

      for (j = 0; j < NbIn; j++)
	if(val[j])
	  {
	    use[j] = 1;  // variable i is used in the rule.
	    k ++;        // increase the number of variable in use in the rule
	  }
      if (k > i.maxVr)  i.maxVr = k;
      tmp += k;         // increase the overall number of variable in use

      if(classif)
	{
	  conc = Rule[p]->GetAConc(n);
	  for (j = 0; j < i.nClass; j++)
	    if(conc == i.labels[j])
	      { (i.nRc[j])++; break; }
	}
    }

  if(NbRules) i.meanVr = (double) tmp / NbRules;
  i.maxR = 1;
  i.meanMF = 0;
  tmp = k = 0;
  for(p = 0; p < NbIn; p++)
    {
      j = In[p]->GetNbMf();
      if(In[p]->IsActive())
	{
	  i.maxR *= j;
	  i.nMf[p] = j;
	}
      else i.nMf[p] = -j;
      if(use[p])  { k++; i.meanMF += (double) j;}
    }
  if(k) i.meanMF /= k;
  i.nVar = k;

  for(p = 0; p < NbOut; p++)
    {
      i.nMf[p+NbIn] = Out[p]->GetNbMf();
      if(p == i.out && i.nMf[p+NbIn] == 0) // crisp
	{
	  if(classif) i.nMf[p+NbIn] = i.nClass;
	  else
	    i.nMf[p+NbIn] = Out[p]->GetNbPossibles();
	}
    }

  delete [] val;
  delete [] use;
  return 0;
}


void FIS::Normalize(double **SampleData,int nbrow)
  //**********************************************
{
  int i;
  for(i = 0; i < NbIn; i++)
    {
      if ( SampleData != NULL ) ::Normalize( SampleData , i ,nbrow, In[i]->min(),In[i]->max());
      In[i]->Normalize();
    }
  for(i = 0 ; i < NbOut ; i++ )
    {
      if ( SampleData != NULL ) ::Normalize( SampleData , i+NbIn ,nbrow, Out[i]->min(),Out[i]->max());
      if ( strcmp(Out[i]->GetOutputType(), OUT_FUZZY::OutputType() ) != 0 )
	{
	  // normalize rule conclusions
	  for ( int varnrule = 0 ; varnrule < GetNbRule() ; varnrule ++ )
	    Rule[varnrule]->Normalize(i,Out[i]->min(),Out[i]->max());
	}
      Out[i]->Normalize();
    }

  return;
}


void FIS::UnNormalize(double **SampleData,int nbrow)
  //************************************************
{
  for(int i = 0; i < NbIn; i++)
    {
      In[i]->UnNormalize();
      // modif bch August 2012 - for optimization unnormalize input data
      if ( SampleData != NULL ) ::UnNormalize( SampleData , i ,nbrow, In[i]->OLower,In[i]->OUpper);
    }
  for(int i = 0 ; i < NbOut ; i++ )
    {
      if ( SampleData != NULL ) ::UnNormalize( SampleData , NbIn+i ,nbrow, Out[i]->OLower,Out[i]->OUpper);
      if ( strcmp(Out[i]->GetOutputType(), OUT_FUZZY::OutputType() ) != 0 )
	{
	  for ( int varnrule = 0 ; varnrule < GetNbRule() ; varnrule ++ )
	    Rule[varnrule]->UnNormalize(i,Out[i]->OLower,Out[i]->OUpper);
	}
      Out[i]->UnNormalize();
    }
  return;
}

int FIS::ComputeNbActRule(void)
  //***************************
{
  int i;

  NbActRules = 0;

  for(i = 0; i < NbRules; i++)
    if( Rule[i]->IsActive() ) NbActRules ++;

  return NbActRules;
}



//****************************************************************
//
//
//             Global variable and Fonction  (Sorting)
//
//
//****************************************************************

double * CumG;  // Used for rule sorting

int CmpCumDec(const void * a, const void *b)
  //****************************************
  // Out of class function called by  qsort()
{
  if(CumG[*(unsigned int *)a] > CumG[*(unsigned int *)b]) return -1;
  if(CumG[*(unsigned int *)a] < CumG[*(unsigned int *)b]) return  1;
  return 0;
}

int CmpCumInc(const void * a, const void *b)
  //****************************************
  // Out of class function called by  qsort()
{
  if(CumG[*(unsigned int *)a] < CumG[*(unsigned int *)b]) return -1;
  if(CumG[*(unsigned int *)a] > CumG[*(unsigned int *)b]) return  1;
  return 0;
}

//  RULE SORT

void FIS::SortRules(double **dat, int n, int order)
  //***********************************************
  // order following cumulated weight : >0 decreasing, <0 increasing, 0 no change
{

  if(order == 0) return;

  int i, j, k;
  unsigned int * sorted;

  CumG = new double [NbRules];

  for(k = 0; k < NbRules; k ++) CumG[k] = 0.;

  for(j = 0; j < n; j ++)
    {
      for(i = 0; i  < NbIn; i++)
	{
	  if(In[i]->IsActive() == false) continue;

	  if(FisIsnan(dat[j][i]))            // missing value
	    {
	      if(!strcmp(strMissingValues, FIS::RandomMissingValues()))
		In[i]->GetRandDegs(dat[j][i]);
	      else if(!strcmp(strMissingValues, FIS::MeanMissingValues()))
		In[i]->SetEqDegs(dat[j][i]);
	      else
		{
		  sprintf( ErrorMsg, "~UnknownMissingValueStrategy~: %.50s", strMissingValues );
		  throw std::runtime_error( ErrorMsg );
		}
	    }
	  In[i]->GetDegs(dat[j][i]);
	}

      for(k = 0; k < NbRules; k ++)
	{
	  if ( Rule[k]->IsActive() )	  Rule[k]->ExecRule();
	  CumG[k] += Rule[k]->Weight;
	}
    }

  sorted = new unsigned int [NbRules];
  for(i = 0; i < NbRules; i++)   sorted[i] = i;

  if(order > 0)
    qsort(sorted, NbRules, sizeof(unsigned int), CmpCumDec);
  else qsort(sorted, NbRules, sizeof(unsigned int), CmpCumInc);

  RULE ** tmp;
  tmp = new RULE *[NbRules];

  for(i = 0; i < NbRules; i++)
    tmp[i] = new RULE(*(Rule[sorted[i]]), In, Out);

  for(i = 0; i < NbRules; i++)     delete Rule[i];
  delete [] Rule;

  Rule = tmp;

  for(i = 0; i < NbOut; i++) Out[i]->InitPossibles(Rule, NbRules,i);

  delete [] CumG;
  delete [] sorted;
}

/*****************************************************************
 *                                                               *
 *                                                               *
 *                    Vocabulary Reduction                       *
 *                                                               *
 *                                                               *
 *****************************************************************/

int FIS::VocReduc(int NumS, double **Data, int NbEx, double MuSeuil, double PerfLoss, int NConc, int ExtVoc)
  //*********************************************************************
{
  if( (NumS < 0) || (NumS > NbOut - 1) || (! Out[NumS]->IsActive()) )
    {
      sprintf( ErrorMsg, "~InvalidOutputNumber~: %d~", NumS);
      throw std::runtime_error( ErrorMsg );
    }
  // Nclas:  number of classes
  // nData:  number of data to be clustered (NbEx or NPos when ExtVoc= true)
  int i, j, n, nc, classif, NClas, nPos, nData;
  // Conc:  array containing modified Conclusions
  // SaveConc:  array od the initial conclusions
  // Centres:  cluster centres (result of the k-means)
  double *Conc, *OutData, *SaveConc, *Centres, *Pos;
  // PerfRef:  Initial performance
  double cov,err, PerfRef, pi, mini, maxi;
  int *MisClassified=NULL;
  double *Lab=NULL;

  Conc = OutData = SaveConc = Centres = Pos= NULL;
  Conc = new double[NbRules];

  // following two calls : always in the same order
  // allocate Classes array in crisp classification case and
  // allocate  MisClassified and lab arrays in crisp classification case
  ClassifCheck(Data, NbEx, NumS);
  ResClassifAlloc(MisClassified, Lab, NumS);

  PerfRef = Perf(NumS,Data,NbEx,cov,err,MuSeuil,MisClassified,Lab);

  classif = (Out[NumS]->Classification() &&
	     (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType())) &&
	     ( (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz())) ||
	       (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))) );

  if(!classif && NConc<2 && NConc!=0)
    {
      sprintf(ErrorMsg, "~Must Be At Least Two Conclusions (min-max)~");
      throw std::runtime_error( ErrorMsg );
    }

  for(i=0;i<NbRules;i++) Conc[i]=Rule[i]->GetAConc(NumS);
  Out[NumS]->InitPossibles(Rule,NbRules,NumS);
  NClas = 0;

  // Init OutData array and nData value
  // ExtVoc=false: data values
  // ExtVoc=true:  possible (rule conclusion) values
  if (!ExtVoc)
    {
      OutData = new double[NbEx];
      for(i=0;i<NbEx;i++) OutData[i]=Data[i][NbIn + NumS];
      nData = NbEx;
      mini=Out[NumS]->min();
      maxi=Out[NumS]->max();
    }
  else
    {
      Pos = Out[NumS]->GetPossibles();
      nPos = Out[NumS]->GetNbPossibles();
      OutData = new double[nPos];
      for(i=0;i<nPos;i++) OutData[i]=Pos[i];
      nData = nPos;
      mini=Pos[0];
      maxi=Pos[nPos-1];
    }

  if(NConc && (!classif))         // number of conclusions pre-defined
    {
      InitCentres(Centres, NConc, mini, maxi);
      Kmeans(OutData,nData,Centres,NConc,false);
      nc = NConc;
      // Keep only the nc centers which correspond to non empty groups
      i = KmeansNE(OutData,nData,Centres,nc);
      Centres[0]=mini;
      Centres[nc-1]=maxi;
      NewConc(Conc,Centres,nc); // Compute new conclusions
      delete [] Centres; // allocated within InitCentres

      for(j=0;j<NbRules;j++) Rule[j]->SetAConc(NumS, Conc[j]); // Conclusions update
      Out[NumS]->InitPossibles(Rule,NbRules,NumS);
    }
  else if(classif) // classif case
    {
      NClas = Out[NumS]->Def->NbClasses();
      for(i = 0; i < NClas; i++) MisClassified[i] = 0;
      NewConc(Conc,Lab,NClas);
      for(j=0;j<NbRules;j++) Rule[j]->SetAConc(NumS, Conc[j]);
      Out[NumS]->InitPossibles(Rule,NbRules,NumS);
    }
  else  // Find the minimum number of distinct conclusions within a given performance loss
    {
      SaveConc = new double[NbRules];
      for(i=0;i<NbRules;i++) SaveConc[i]=Conc[i];

      if(NbRules < NbEx) n = NbRules;
      else n = NbEx;
      for(i=2;i<n;i++)
	{
	  InitCentres(Centres, i, mini, maxi);
	  Kmeans(OutData,nData,Centres,i,false);
	  nc = i;
	  KmeansNE(OutData,nData,Centres,nc);
	  Centres[0]=mini;
	  Centres[nc-1]=maxi;
	  NewConc(Conc,Centres,nc);
	  delete [] Centres;

	  for(j=0;j<NbRules;j++) Rule[j]->SetAConc(NumS, Conc[j]);
	  Out[NumS]->InitPossibles(Rule,NbRules,NumS);

	  pi=(Perf(NumS,Data,NbEx,cov,err,MuSeuil,MisClassified,Lab)-PerfRef)/PerfRef;
	  if(pi < PerfLoss) break;   // check the performance loss
	  else if(i!=(n-1))
	    for(j=0;j<NbRules;j++)
	      {
		Conc[j]=SaveConc[j];
		Rule[j]->SetAConc(NumS, Conc[j]);
	      }
	}
      delete [] SaveConc;
    }

  delete [] OutData;
  delete [] Conc;
  if(MisClassified != NULL) delete [] MisClassified;

  //reset ptr
  MisClassified=NULL;

  return Out[NumS]->GetNbPossibles();
}

void FIS::NewConc(double *&Conc, double *Centres, int nconc)
  //********************************************************
{
  int index, j, k;
  double min;
  // replace old conclusions by the nearest new ones
  for(j=0;j<NbRules;j++)
    {
      index=0;
      min=fabs(Conc[j]-Centres[0]);
      for(k=1;k<nconc;k++)
	{
	  if(min>fabs(Conc[j]-Centres[k]))
	    {
	      min=fabs(Conc[j]-Centres[k]);
	      index=k;
	    }
	}
      Conc[j]=Centres[index];
    }
}

/*****************************************************************
 *                                                               *
 *                                                               *
 *                    Weighted Performance                       *
 *                                                               *
 *                                                               *
 *****************************************************************/


double FIS::WeightedPerf(int NumS, char * fdata, int NPart, char *DomBreakpoints, char *PartWeight, double &WeightedCov, double &MaxError, double MuSeuil, int ErrorType, char *fres, FILE *display)
  //***************************************************************
{
  double *Result;
  double *Cov;
  double *MaxErr;
  double *NSam;
  double WeightedPerf;

  if( (ErrorType < 1) || (ErrorType > 5) )
    {
      sprintf( ErrorMsg, "~ErrorType must be 1,2,3,4 or 5~");
      throw std::runtime_error( ErrorMsg );
    }
  // compute the performance of each part and the global performance, initialize all arrays
  NPart = Performance(NumS,fdata,NPart,DomBreakpoints,Result,Cov,MaxErr,NSam,MuSeuil,ErrorType,fres,display);

  MaxError = MaxErr[NPart];
  // compute the Weighted perf
  WeightedPerf = ComputeWeightedPerf(PartWeight,NPart,Result,Cov,WeightedCov);

  delete [] Result;
  delete [] Cov;
  delete [] MaxErr;
  delete [] NSam;

  return WeightedPerf;

}

double FIS::ComputeWeightedPerf(char *PartWeight, int NPart, double *&ResultTab, double *&Couvert, double &WeightedCov)
  //***********************************************************
{
  int i;
  double WeightedPerf = 0;
  double * Weights=NULL;
  // Init the weights
  Weights = new double[NPart];
  InitWeights(NPart,PartWeight,Weights);

  WeightedPerf = ResultTab[NPart];
  // compute the weighted performance
  for(i=0;i<NPart;i++)
    {
      WeightedPerf+=(Weights[i]*ResultTab[i]);
    }
  WeightedPerf=WeightedPerf/2;

  WeightedCov = Couvert[NPart];
  // compute the wieghted coverage
  for(i=0;i<NPart;i++)
    {
      WeightedCov+=(Weights[i]*Couvert[i]);
    }
  WeightedCov=WeightedCov/2;

  delete [] Weights;

  return WeightedPerf;
}




int FIS::Performance(int NumS, char * fdata, int NPart, char *DomBreakPoints,
                     double *&ResultTab, double *&Couvert, double *&MaxErr, double *&NItem,
                     double MuSeuil, int ErrorType, char * fres, FILE *display)
  //*****************************************************************************
  // Return : -2 when not enough input variables or not active output
  //          -1 when the observed output is not part of the sample file
  //          the number of misclassified items when classif
  //          the root mean square error
{
  if( (NumS < 0) || (NumS > NbOut - 1) || (! Out[NumS]->IsActive()) )
    {
      Couvert[NPart] = 0.;
      sprintf( ErrorMsg, "~InvalidOutputNumber~: %d~", NumS);
      throw std::runtime_error( ErrorMsg );
    }


  double ** Data;
  int i,j, nbcol, classif, NbEx, ref;
  int * MisClassified=NULL;
  double * Lab=NULL;
  double * BreakPoints=NULL;

  FILE * f;

  Data = NULL;
  nbcol = 0;

  f = NULL;
  if(fres)
    if((f = fopen(fres, "wt")) == NULL)
      {
	sprintf( ErrorMsg, "~CannotOpenResultFile~: %.100s~", fres);
	throw std::runtime_error( ErrorMsg );
      }

  //read data -allocate Data array
  Data = ReadSampleFile(fdata, nbcol, NbEx);

  if(nbcol < NbIn)  return -2;

  ref = true;
  if(nbcol < NbIn + 1 + NumS) ref = false; // Observed output is not part of the sample file

  // Write a comment line in the result file
  WriteHeader(NumS, f, ref);

  // following two calls : always in the same order
  //allocate Classes array in crisp classification case
  ClassifCheck(Data, NbEx, NumS);
  // allocate  MisClassified and lab arrays in crisp classification case
  ResClassifAlloc(MisClassified, Lab, NumS);

  classif = Out[NumS]->Classification() &&
    (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType())) &&
    ( (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz()))	 ||
      (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz())));

  if(display && classif) fprintf(display, "\nThis is a classification case\n");


  if((!classif) && (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType())))
    {
      BreakPoints = new double[NPart-1];
      for(j=0;j<NPart-1;j++) BreakPoints[j]=0;
      InitBreakPoints(NumS,NPart,DomBreakPoints,BreakPoints);
    }

  //call to Perf function
  NPart = Perf(NumS, Data, NbEx, NPart, ResultTab, Couvert, MaxErr, NItem, MuSeuil, ErrorType, BreakPoints ,MisClassified, Lab, ref, f, display);

  if(f) fclose(f);

  if(display)
    {
      fprintf(display, "\n");
      if(ref)
	{
	  if(classif)
	    {
	      fprintf(display, "Number of misclassified items : %6d, in percentage %2d %%. \nDetail by classes : ",
		     (int) ResultTab[NPart], (int) (100. * ResultTab[NPart]/ NbEx));

	      for(i = 0; i < Out[NumS]->Def->NbClasses(); i ++)
		fprintf(display, "Number of misclassified items : %6d, in percentage %2d %%. ",
		       MisClassified[i], (int) (100. * ResultTab[i]/NItem[i]));
	      fprintf(display, "\n");
	    }

	  else
	    {
	      fprintf(display, "Mean Square error : %11.6f\n", ResultTab[NPart]);
	      if(NPart!=1)
		{
		  for(i=0;i<NPart;i++) fprintf (display, "Mean Square error of part %i : %11.6f\n", i, ResultTab[i]);
		}
	    }
	}
      int RuleWgt = false;
      for(i = 0; i < NbRules; i++)
	if(fabs(Rule[i]->GetExpertWeight() -1.0) > EPSILON) 
	  {
	    RuleWgt = true;
	    break;
	  }
      if(RuleWgt) fprintf(display, "\nWarning:  the rules are weighted.\n");
    }

  // delete 2D data array
  for (i = 0; i < NbEx; i++)
    {
      if (Data[i] != NULL) delete [] Data[i];
    }
  if (Data != NULL) delete [] Data;
  // delete MisClassified array
  if(MisClassified != NULL) delete [] MisClassified;

  if(BreakPoints!=NULL) delete [] BreakPoints;

  if(Lab!=NULL) delete [] Lab;

  return NPart;
} // End of Performance(char *)


void FIS::InitWeights(int NPart, char *PartWeights, double *&Weights)
  //*****************************************************************
{
  int i;
  double weightsum=0;
  // init weights array with given values
  if(PartWeights!=NULL)
    {
      SearchNb(PartWeights,Weights,NPart,',','[',']');
      for(i=0;i<NPart;i++) weightsum+=Weights[i];
      for(i=0;i<NPart;i++) Weights[i]=Weights[i]/weightsum;
    }
  // if no given values, each part has an equal weight
  else
    {
      for(i=0;i<NPart;i++) Weights[i]=(1.)/NPart;
    }
}

void FIS::InitBreakPoints(int NumS, int Npart, char *DomBreakPoints, double *&BreakPoints)
  //**************************************************************************************
{
  int i;
  double incfront;

  //Reading and filling array of breakpoints
  if(DomBreakPoints!=NULL)
    {
      SearchNb(DomBreakPoints,BreakPoints,(Npart-1),',','[',']');
      for(i=0;i<(Npart-1);i++)
	{
	  if(BreakPoints[i]<=(Out[NumS]->min()))
	    {
	      sprintf(ErrorMsg, "~Invalid BreakPoint : %f Less Than or Equal To Inferior Bound %f~\n",
		      BreakPoints[i], Out[NumS]->min());
	      throw std::runtime_error(ErrorMsg);
	    }
	  else if(BreakPoints[i]>=(Out[NumS]->max()))
	    {
	      sprintf(ErrorMsg, "~Invalid BreakPoint : %f Higher Than or Equal To Superior Bound %f~\n",
		      BreakPoints[i], Out[NumS]->max());
	      throw std::runtime_error(ErrorMsg);
	    }
	}
    }
  //if no given breakpoints, variation domain is cut into equal parts
  else
    {
      incfront=(Out[NumS]->max()-Out[NumS]->min())/Npart;
      for(i=1;i<Npart;i++)
	{
	  BreakPoints[i-1]=(Out[NumS]->min())+i*incfront;
	}
    }
}

int FIS::Perf(int NumS, double ** Data, int NbEx, int NPart, double *&ResultTab,  double *&Couvert,
	      double *&MaxErr, double *&NItem, double MuSeuil, int ErrorType, double *BreakPoints,
	      int * MisClassified, double * Lab, int ref, FILE * f, FILE *display)
  //***********************************************************************************************
  // Measures the FIS performance du Sif for the NumS output for the data set fdata
  // Returns the number of misclassified for a class output
  // otherwise the mean error
  // Also writes the MisClassified array which must be preallocated
  // if the output is a class
  // the Lab array contains labels for classes, must be preallocated
{
  if( (NumS < 0) || (NumS > NbOut - 1) || (! Out[NumS]->IsActive()) )
    {
      Couvert[NPart] = 0.;
      sprintf( ErrorMsg, "~InvalidOutputNumber~: %d~", NumS);
      throw std::runtime_error( ErrorMsg );
    }

  int i, j, l, classif, diag, nPos, NClas, Blank, ExPart;
  double Mu;
  double * Pos,*NbBlank=NULL;

  ExPart = 0;
  //check for classif case
  classif = (Out[NumS]->Classification() &&
	     (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType())) &&
	     ( (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz())) ||
	       (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))) );
  //Npart and arrays initialization
  if(classif)
    {
      NPart=Out[NumS]->Def->NbClasses();
    }
  else if(! strcmp(Out[NumS]->GetOutputType(), OUT_FUZZY::OutputType()))
    {
      NPart=Out[NumS]->GetNbMf();
    }

  ResultTab = new double[NPart+1];
  Couvert = new double[NPart+1];
  MaxErr = new double[NPart+1];
  NItem = new double[NPart+1];
  for(j=0;j<NPart+1;j++)
    {
      ResultTab[j]=0;
      Couvert[j]=0;
      MaxErr[j]=0;
      NItem[j]=0;
    }

  NbBlank = new double [NPart+1];
  for(i=0;i<(NPart+1);i++) NbBlank[i]=0;

  NItem[NPart]=NbEx;

  diag = ( (! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType()))
	   && (! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))
	   && (! Out[NumS]->Classification()) );

  Pos = Out[NumS]->GetPossibles();
  nPos = Out[NumS]->GetNbPossibles();
  NClas = Out[NumS]->Def->NbClasses();

  Blank = 0;

  if(classif)
    for(i = 0; i < NClas; i++) MisClassified[i] = 0;

  if(! ref && (classif || diag))
    {
      if(! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::SugenoDefuz()))
	((DEFUZ_SugenoClassif *) (Out[NumS]->Def))->SetClasses(Pos, nPos);
      else if(! strcmp(Out[NumS]->Defuzzify(), OUT_CRISP::MaxCrispDefuz()))
	((DEFUZ_MaxCrisp *) (Out[NumS]->Def))->SetClasses(Pos, nPos);
    }

  for(i = 0; i < NbEx; i ++)
    {
      if(ref && f) fprintf(f, FORMAT_DOUBLE, Data[i][NbIn + NumS]);
      //check the part in which the data is (depending on the case) and increment NItem
      if(classif)
	{
	  for(j=0;j<NClas;j++)
	    {
	      if(Data[i][NbIn + NumS]==Lab[j]) ExPart=j;
	      NItem[j]+=1;
	    }
	}
      //fuzzy output case (membership degree weighting)
      else if(! strcmp(Out[NumS]->GetOutputType(), OUT_FUZZY::OutputType()))
	{
	  for(j=0;j<(Out[NumS]->GetNbMf());j++)
	    NItem[j]+=(1*(Out[NumS]->GetADeg(j,Data[i][NbIn + NumS])));
	}
      //crisp regression output case
      else if(! strcmp(Out[NumS]->GetOutputType(), OUT_CRISP::OutputType()))
	{
	  for(j=0;j<NPart;j++)
	    {
	      //for each part, lower bound is part of the domain, upper is not (except for last part)
	      if(j!=(NPart-1))
		{
		  if((Data[i][NbIn + NumS])<BreakPoints[j])
		    {
		      ExPart=j;
		      NItem[j]+=1;
		      break;
		    }
		}
	      else
		{
		  ExPart=j;
		  NItem[j]+=1;
		}
	    }
	}

      Blank=0;
      Mu = Infer(Data[i], NumS, f, display);
      if(display) fprintf(display, "\n Output  n %d  observed      inferred     delta    blank",  NumS+1);
      //check if sample is active
      if(Mu <= MuSeuil)
	{
	  Blank=1;
	  NbBlank[NPart]=(NbBlank[NPart])+1;
	  //fuzzy output case (membership degree weighting)
	  if(! strcmp(Out[NumS]->GetOutputType(), OUT_FUZZY::OutputType()))
	    {
	      for(j=0;j<(Out[NumS]->GetNbMf());j++)
		NbBlank[j]+=(1*(Out[NumS]->GetADeg(j,Data[i][NbIn + NumS])));
	    }
	  else
	    NbBlank[ExPart]+=1;
	}

      // Special mode for OUT_CRISP, MaxCrisp defuz and classification
      if(diag)
	{
	  // Distance to the prototype (1 - MuInfer)
	  if(OutValue[NumS] == Data[i][NbIn + NumS])
	    {
	      for(l = 0; l < nPos; l++)
		if(OutValue[NumS] == Pos[l])
		  {
		    OutErr[NumS] = 1. - Out[NumS]->MuInfer[l];
		    break;
		  }
	    }
	  // Non ordered classes
	  else OutErr[NumS] = 1.;
	}
      //output part of the testing file
      else if(ref)
	{
	  OutErr[NumS] = OutValue[NumS] - Data[i][NbIn + NumS];
	}

      if(display)
	{
	  if(ref) fprintf(display, "  %8.3f      %8.3f     %8.3f", Data[i][NbIn + NumS], OutValue[NumS], OutErr[NumS]);
	  fprintf(display, "   %d", Blank);
	}
      //error incrementation and maximum error check (depend on the chosen error type)
      if(! Blank && ref)
	{
	  //fuzzy output case
	  if(! strcmp(Out[NumS]->GetOutputType(), OUT_FUZZY::OutputType()))
	    {
	      for(j=0;j<(Out[NumS]->GetNbMf());j++)
		{
		  if(fabs(OutErr[NumS]*(Out[NumS]->GetADeg(j,Data[i][NbIn + NumS]))) > fabs(MaxErr[j]))
		    {
		      MaxErr[j] = OutErr[NumS]*(Out[NumS]->GetADeg(j,Data[i][NbIn + NumS]));
		      if(fabs(OutErr[NumS]) > fabs(MaxErr[NPart])) MaxErr[NPart] = OutErr[NumS];
		    }
		  if(ErrorType==1 || ErrorType==2 || ErrorType==3)
		    {
		      ResultTab[j] += (OutErr[NumS] * OutErr[NumS] * (Out[NumS]->GetADeg(j,Data[i][NbIn + NumS])) * (Out[NumS]->GetADeg(j,Data[i][NbIn + NumS])));
		      ResultTab[NPart] += (OutErr[NumS] * OutErr[NumS]);
		    }
		  else if (ErrorType==4 || ErrorType==5)
		    {
		      ResultTab[j] += fabs(OutErr[NumS] * (Out[NumS]->GetADeg(j,Data[i][NbIn + NumS])));
		      ResultTab[NPart] += fabs(OutErr[NumS]);
		    }
		  else if (ErrorType==4 || ErrorType==5)
		    {
		      ResultTab[j] += fabs((OutErr[NumS]/OutValue[NumS]) * (Out[NumS]->GetADeg(j,Data[i][NbIn + NumS])));
		      ResultTab[NPart] += fabs(OutErr[NumS]/OutValue[NumS]);
		    }
		}
	    }
	  else
	    {
	      if(fabs(OutErr[NumS]) > fabs(MaxErr[ExPart]))
		{
		  MaxErr[ExPart] = OutErr[NumS];
		  if(fabs(OutErr[NumS]) > fabs(MaxErr[NPart])) MaxErr[NPart] = OutErr[NumS];
		}
	      //classification case
	      if(classif && fabs(OutErr[NumS]) )
		{
		  ResultTab[ExPart] += 1.;
		  ResultTab[NPart] +=1.;
		  for(l = 0; l < NClas; l++)
		    if(Data[i][NbIn + NumS] == Lab[l])
		      {MisClassified[l]++; break; }
		}
	      //crisp regression case
	      else
		{
		  if(ErrorType==1 || ErrorType==2 || ErrorType==3)
		    {
		      ResultTab[ExPart] += (OutErr[NumS] * OutErr[NumS]);
		      ResultTab[NPart] += (OutErr[NumS] * OutErr[NumS]);
		    }
		  else if (ErrorType==4)
		    {
		      ResultTab[ExPart] += fabs(OutErr[NumS]);
		      ResultTab[NPart] += fabs(OutErr[NumS]);
		    }
		  else if (ErrorType==5)
		    {
		      ResultTab[ExPart] += fabs(OutErr[NumS]/OutValue[NumS]);
		      ResultTab[NPart] += fabs(OutErr[NumS]/OutValue[NumS]);
		    }
		}
	    }

	  if(f)
	    {
	      if(ref) fprintf(f, FORMAT_DOUBLE, OutErr[NumS]);
	      fprintf(f, "%5d  %8.3f\n", Blank, ResultTab[NPart]);
	    }

	}
    } //End of Item loop
  //performance indexes computation
  if(ResultTab[NPart] && (!classif))
    {
      for(j=0;j<NPart+1;j++)
	{
	  if(ResultTab[j])
	    {
	      if(ErrorType==1) ResultTab[j] = sqrt(ResultTab[j])/(NItem[j] - NbBlank[j]);
	      else if(ErrorType==2 || ErrorType==4 || ErrorType==5) ResultTab[j] = (ResultTab[j])/(NItem[j] - NbBlank[j]);
	      else if(ErrorType==3) ResultTab[j] = sqrt((ResultTab[j])/(NItem[j] - NbBlank[j]));
	    }
	}
    }

  for(j=0;j<NPart+1;j++) Couvert[j] = 1. - ((double) NbBlank[j] / NItem[j]);

  delete [] NbBlank;
  return NPart;
} // End of Perf


// **********    FATI INFERENCE FOR IMPLICATIVE RULES AND FUZZY INPUTS **********

void FIS::InferFatiPrep(int numout)
  //*******************************
{
  int impli;
  impli= ! strcmp(Out[numout]->Defuzzify(),OUT_FUZZY::ImpFuzzyDefuz());
  if (impli)
    {
      int i;
      std::list<double> ** dpl = new list<double>*[NbIn];
      // check inputs - add universal MF to any empty input
      for (i = 0; i < NbIn; i++)
	{
	  if ((In[i]->GetNbMf())==0)
	    {
	      In[i]->AddMF(new MFUNIV(In[i]->min(),In[i]->max()));
	    }
	}

      for(i = 0; i < NbIn; i++) dpl[i]= new list<double>;
      if(NbIn == 2) KinkPoints(dpl, numout);
      for(i = 0; i < NbIn; i++) In[i]->DecomposePart(dpl[i]);
      for(i = 0; i < NbIn; i++) delete dpl[i];
      delete [] dpl;
    }
}

void FIS::KinkPoints(std::list<double> ** dl, int nout)
  //*****************************************************
{
  MFDPOSS *mf1, *mf2, *mfi;
  int i, j, n;

  n = Out[nout]->GetNbMf();
  for(i = 0; i < n; i++)
    {
      mf1 = new MFDPOSS(Out[nout]->GetMF(i),0);
      for(j = i+1; j < n; j++)
	{
	  mf2 = new MFDPOSS(Out[nout]->GetMF(j),0);
	  mfi = mf1->Inter(mf2);
	  if (mfi != NULL)
	    {
	      UpdatePartList(nout, dl, mfi->getMaxposs(), i, j);
	      delete mfi;
	    }
	  delete mf2;
	}
      delete mf1;
    }
}

void FIS::UpdatePartList(int iout, std::list<double> **dL, double mposs, int m1, int m2)
  //**************************************************************************************
{
  if(DBL_INF_EQUAL(mposs, 0.5)) return;
  double l, r, conc;
  int i, j, id = 0;
  for (j = 0; j < NbRules; j++)
    {
      conc = Rule[j]->GetAConc(iout);
      if (DBL_EQUAL(conc-1, m1))
	{

	  for (i = 0; i < NbIn; i++)
	    {
	      Rule[j]->GetAProp(id,i);//updates id - returns id=0 if input #i has no MFs
	      // protection for empty premise
	      if (id>0)
		{
		  In[i]->GetMF(id-1)->AlphaKernel(l,r, mposs);
		}
	      else
		{
		  l=In[i]->min();
		  r=In[i]->max();
		}
	      dL[i]->push_back(l);
	      dL[i]->push_back(r);
	    }
	}

      if (DBL_EQUAL(conc-1, m2))
	{
	  for (i = 0; i < NbIn; i++)
	    {
	      Rule[j]->GetAProp(id,i);//updates id - returns id=0 if input #i has no MFs
	      // protection for input with no MFs
	      if (id>0)
		{
		  In[i]->GetMF(id-1)->AlphaKernel(l,r, mposs);
		}
	      else
		{
		  l=In[i]->min();
		  r=In[i]->max();
		}
	      dL[i]->push_back(l);
	      dL[i]->push_back(r);
	    }
	}
    }
}

void FIS::BuildFuzIn(double * val, MFDPOSS ** tpl, MFDPOSS ** fuzval)
  //*******************************************************************
{
  int i;
  for(i = 0; i < NbIn; i++)
    fuzval[i] = tpl[i]->translate(val[i], In[i]->min(), In[i]->max());
}

void FIS::BuildFuzIn(double * crispin, double * KW, double * SW,  MFDPOSS ** &fuzval, double maxposs)
  //**********************************************************************************************************
{
  // get vmin,vmax=input range as boundaries for trapezoid fuzzy number
  // KW: kernel width, SW: support width, cripsIn: kernel center
  double left,right,topleft,topright,vmin,vmax;
  double KWD2,SWD2;
  int i;
  LIST* lst =NULL;

  if(fuzval != NULL)
    {
      for(i = 0; i < NbIn; i++)
	delete fuzval[i];
      delete [] fuzval;
    }

  fuzval = new MFDPOSS*[NbIn];

  for(i = 0; i < NbIn; i++)
    {
      KWD2=KW[i]*0.5;
      SWD2=SW[i]*0.5;
      left=crispin[i]-SWD2;
      right=crispin[i]+SWD2;
      topleft=crispin[i]-KWD2;
      topright=crispin[i]+KWD2;
      vmin=In[i]->min();
      vmax=In[i]->max();
      if (left<vmin) left=vmin;
      if (right>vmax) right=vmax;
      if (topleft<vmin) topleft=vmin;
      if (topright>vmax) topright=vmax;
      lst= new LIST;
      lst->home();
      lst->add(left,0.0);
      lst->add(topleft,maxposs);
      lst->add(topright,maxposs);
      lst->add(right,0.0);
      fuzval[i] =new MFDPOSS(lst);// simplification of dposs is included in constructor
      if (lst!=NULL)
	{
	  delete lst;
	  lst=NULL;
	}
    }//loop on inputs
}

MFDPOSS * FIS::InferAcut(double *binf, double *bsup, int iout, FILE *fg, double alpha, FILE *display)
  //************************************************************************************
{
  MFDPOSS *mfinf, *mfsup, *mfjoin;



  Infer(binf, iout, fg, 0, alpha);
  if (Out[iout]->MfGlob != NULL)  mfinf = Out[iout]->MfGlob->Clone();
  else mfinf =NULL;
  if(display)
    {
      fprintf(display, "\nin InferAcut after infer with binf mfinf=\n");
      if (mfinf)
	mfinf->Print(display);
      else
	fprintf(display, "\nmfinf is NULL");
    }
  //save MFConc for binf

  Infer(bsup, iout, fg, 0, alpha);

  if (Out[iout]->MfGlob != NULL)  mfsup = Out[iout]->MfGlob->Clone();
  else mfsup = NULL;
  if(display)
    {
      fprintf(display, "\nin InferAcut after infer with bsup mfsup=\n");
      if (mfsup)
	mfsup->Print(display);
      else
	fprintf(display, "\nmfsup is NULL");
    }
  // join MFConc for binf and MFConc for sup
  if (mfinf && mfsup)
    {
      mfjoin = mfinf->Join(mfsup);
      if (display)
	{
	  fprintf(display, "\nin InferAcut join dposs=\n");
	  mfjoin->Print(display);
	}
      delete mfinf;
      delete mfsup;
      return mfjoin;
    }

  if (mfinf) delete mfinf;
  if (mfsup) delete mfsup;
  return NULL;
}

MFDPOSS * FIS::InferFatiAlpha(MFDPOSS ** v, int a, int nout, FILE * f, FILE *display)
  //*****************************************************************************
{
  int i, j, *nb;
  double ** inf, ** sup, *vinf, *vsup, alpha;
  MFDPOSS * res=NULL;   // Inferred result of a given alpha-cut
  // Inferred result of the union of the alpha-cuts
  std::list<MFDPOSS> * ualf = NULL, *tmp = NULL;

  inf = new double *[NbIn];
  sup = new double *[NbIn];
  vinf = new double [NbIn];
  vsup = new double [NbIn];
  nb = new int [NbIn];

  for(i = 0; i < NbIn; i++)
    {
      // bug correction  August 4, 2009
      //inf[i] = new double[2*In[i]->GetNbMf()-1];
      //sup[i] = new double[2*In[i]->GetNbMf()-1];
      inf[i] = new double[In[i]->nPart];
      sup[i] = new double[In[i]->nPart];

      nb[i] = In[i]->getIntersect(v[i]->acuts[a], inf[i], sup[i]);
      if (display)
	{
	  fprintf(display, "\nin InferFatiAlpha input=%d, nb[%d]=%d",i+1,i,nb[i]);
	  for (int k=0;k<nb[i];k++)
	    fprintf(display, "\ninf[%d][%d]=%g,sup[%d][%d]=%g",i,k,inf[i][k],i,k,sup[i][k]);
	}

    }

  alpha = v[0]->acuts[a].alpha;

  // This loop is only valid for one or two non empty inputs
  // changed May 7 - added forgotten union in loop for two input case
  for(i = 0; i < nb[0]; i++)
    {
      vinf[0] =inf[0][i];
      vsup[0]= sup[0][i];
      //two input case
      for(j = 0; NbIn == 2 && j < nb[1]; j++)
	{
	  vinf[1] =inf[1][j];
	  vsup[1]= sup[1][j];
	  res = InferAcut(vinf, vsup, nout, f, alpha, display);
	  if(res)
	    {
	      tmp  = res->Union(ualf);
	      if(ualf) { ualf->clear(); delete ualf;}
	      ualf = tmp;
	      delete res;
	      res=NULL;
	    }
	  if (display)
	    {
	      fprintf(display, "\nin InferFatiAlpha after InferAcut in loop i=%d j=%d\n",i,j);
	      (*ualf->begin()).Print(display);
	    }

	}
      if(NbIn == 1) res = InferAcut(vinf, vsup, nout, f, alpha, display);

      //single input case
      if(res)
	{
	  tmp  = res->Union(ualf);
	  if(ualf) { ualf->clear(); delete ualf;}
	  ualf = tmp;
	  delete res;
	  res=NULL;
	}
      if (display)
	{
	  fprintf(display, "\nin InferFatiAlpha after InferAcut in loop i=%d\n",i);
	  (*ualf->begin()).Print(display);
	}

    }

  delete [] nb; delete [] vinf; delete [] vsup;
  for(i = 0; i < NbIn; i++)
    {
      delete [] inf[i];
      delete [] sup[i];
    }
  delete [] inf; delete [] sup;

  if(!ualf || (ualf->size() == 0)) return NULL;
  if (ualf->size()>1)  fprintf(display, "WARNING separated unions for alpha union\n");
  res = (*ualf->begin()).Clone();
  ualf->clear();   delete ualf;
  return res;
}

MFDPOSS * FIS::InferFati(MFDPOSS ** v, int nalf, int numout, FILE * f, FILE *display)
  //****************************************************************************
  // For inference of fuzzy inputs with implicative rules
{
  if(NbIn < 0 || NbIn > 2)
    {
      sprintf( ErrorMsg, "~Invalid#InputsInferFatiLimitedTo2~:  %d", NbIn);
      throw std::runtime_error( ErrorMsg );
    }

  if(NbRules == 0)
    {
      sprintf( ErrorMsg, "~NoRuleToInfer~");
      throw std::runtime_error( ErrorMsg );
    }
  // implicative case
  if (strcmp(Out[numout]->Defuzzify(),OUT_FUZZY::ImpFuzzyDefuz()))
    {
      // conjunctive system - call InferCheck
      sprintf( ErrorMsg, "~OUTPUT~MUST~BE~IMPLICATIVE~FOR~FUZZY~INPUT~INFERENCE");
      throw std::runtime_error( ErrorMsg );
    }

  int i, j;
  MFDPOSS * resalf, * res;
  std::list<MFDPOSS> *resfin = NULL, *tmp = NULL;

  for(i = 0; i < NbIn; i++)
    {

      v[i]->DecompAcut(nalf);
    }

  for  (j = nalf -1; j >=0; j--) // alpha-cut loop
    {
      resalf = InferFatiAlpha(v, j, numout, f, display);
      if(resalf)
	{
	  tmp  = resalf->Union(resfin);
	  if(resfin) {resfin->clear(); delete resfin;}
	  resfin = tmp;
	  delete resalf;
	}
      else break;
    }

  if(!resfin || (resfin->size() == 0)) return NULL;
  if (resfin->size()>1)  fprintf(display, "WARNING separated unions for alpha union\n");
  res = (*resfin->begin()).Clone();
  resfin->clear(); delete resfin;
  //clone inf. dposs result to keep in MfGlob
  if (Out[numout]->MfGlob != NULL)
    delete Out[numout]->MfGlob;
  Out[numout]->MfGlob=res->Clone();
  // defuzzified value=middle of kernel of output dposs
  OutValue[numout] = Out[numout]->Def->EvalOut(Rule, NbRules, Out[numout], f, display);
  return res;
} // End of InferFati()

int FIS::Conj2Imp(int numoutput, const char * DisjType, bool transfPart)
  //**********************************************************************
{
  //**********************************************************************
  // Convert a conjuncttive output into an implicative one
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
  // if transfPart=false returns 1 if partition is already Qsp, 0 otherwise
  //**********************************************************************

  int rettrans=0;
  bool ret=false;
  int nbmf;
  int allowedimpliMFs=0;

  if((numoutput < 0) || (numoutput > NbOut-1)) return -1;
  // An impli output: nothing to do
  if(!strcmp(Out[numoutput]->Defuzzify(), OUT_FUZZY::ImpFuzzyDefuz())) return -2;
  // A crisp output: nothing to do
  if(strcmp(Out[numoutput]->GetOutputType(), OUT_FUZZY::OutputType())) return -3;
  //Fuzzy output but no MFs
  if ((Out[numoutput]->GetNbMf())<=0) return -4;

  // fuzzy output with forbidden MF for implicative output
  nbmf=  Out[numoutput]->GetNbMf();

  for (int imf=0;imf<nbmf;imf++)
    {
      allowedimpliMFs=(!strcmp(Out[numoutput]->GetMF(imf)->GetType(),"trapezoidal"))||(!strcmp(Out[numoutput]->GetMF(imf)->GetType(),"triangular"));
      allowedimpliMFs=allowedimpliMFs||(!strcmp(Out[numoutput]->GetMF(imf)->GetType(),"SemiTrapezoidalSup"))||(!strcmp(Out[numoutput]->GetMF(imf)->GetType(),"SemiTrapezoidalInf"));
      allowedimpliMFs=allowedimpliMFs||(!strcmp(Out[numoutput]->GetMF(imf)->GetType(),"universal"));
      allowedimpliMFs=allowedimpliMFs||(!strcmp(Out[numoutput]->GetMF(imf)->GetType(),"door"));
      if (!allowedimpliMFs)
	{
	  sprintf( ErrorMsg, "ForbiddenMFshape~in~implicative~Systems");
	  throw std::runtime_error( ErrorMsg );
	  //return -5;
	}
    }
  //all other cases: fuzzy output, cases=no MF,1MF or more
  // in all cases update inference operators
  Out[numoutput]->SetOpDefuz(OUT_FUZZY::ImpFuzzyDefuz());
  if(DisjType)  Out[numoutput]->SetOpDisj(DisjType);
  else Out[numoutput]->SetOpDisj(OUT_FUZZY::DisjIrg());

  // Reset ExpertWeight for implicative rules
  for(int i = 0; i < NbRules; i++) Rule[i]->SetExpertWeight(1.0);

  if (transfPart)
    {
      rettrans=FIS2Qsp(numoutput,DisjType);//InitPossibles called from there
      return rettrans;
    }
  else
    {
      // delete and recreate MFConc  arrays because of implicative output
      Out[numoutput]->InitPossibles(Rule, NbRules, numoutput);
      ret=((OUT_FUZZY*)Out[numoutput])->IsQsp();
      return ret;
    }

}

int FIS::FIS2Qsp(int numoutput , const char * DisjType)
  //*****************************************************
{
  //********************************************************
  //turns the quasi strong fuzzy partition (QSP) into a strong fuzzy partition (SFP)
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


  int *sorted, i, ret, retsfp,retqsp;
  ret = 0;
  retsfp=0;
  sorted = NULL;

  if((numoutput < 0) || (numoutput > NbOut-1)) return -1;

  // A crisp output: nothing to do
  if(strcmp(Out[numoutput]->GetOutputType(), OUT_FUZZY::OutputType())) return -3;
  //Fuzzy output but no MFs
  if ((Out[numoutput]->GetNbMf())<=0) return -4;

  // test if output partition is already Qsp
  retqsp=((OUT_FUZZY *) Out[numoutput])->IsQsp();
  if (retqsp) return 2;

  // Turn Sfp into Qsp
  // Sfp2Qsp returns -1 if Nmf=0, -2 if initial output partition is not SFP
  // it modifies its argument (sorted) if necessary + it returns 0 if partition was transformed
  retsfp=((OUT_FUZZY *) Out[numoutput])->Sfp2Qsp(sorted);

  if(retsfp<0)
    {
      // delete and recreate MFConc  arrays because of implicative output
      Out[numoutput]->InitPossibles(Rule, NbRules, numoutput);
      return retsfp-3; //-1,-2 ->-4,-5
    }

  // else sort rules if necessary
  if(sorted)
    // The MF have been re-ordered, update rule conclusions
    // according the new order within the sfp
    {
      for(i = 0 ; i < NbRules ; i++ )
	Rule[i]->SetAConc(numoutput, sorted[(int)Rule[i]->GetAConc(numoutput) - 1] + 1);
      delete [] sorted;
      ret=1;
    }

  // Update rule conclusions according to the Qsp
  for(i = 0 ; i < NbRules ; i++ )
    Rule[i]->SetAConc(numoutput, (int) Rule[i]->GetAConc(numoutput) *2 -1);

  // delete and recreate Possible, MFConc, Classes arrays
  Out[numoutput]->InitPossibles(Rule, NbRules, numoutput);

  return ret;
}




int FIS::Imp2Conj(int outputnumber, const char * DefuzType, const char * DisjType, bool transfPart)
  //*************************************************************************************************
{
  //**********************************************************************************************
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
  //             -5 : Problem in function Qsp2Sfp-initial output partition is not Qsp
  //              0 if Qsp2Sfp has been done OK
  //              1 if Qsp2Sfp has been done OK and if MFs have been re-sorted
  //              2 if output partition is already an Sfp partition
  // if transfPart=false returns 1 if partition is already Sfp, 0 otherwise
  //**********************************************************************************************
  int ret = 0;

  // invalid output number-nothing to do
  if((outputnumber < 0) || (outputnumber > NbOut-1)) return -1;

  // Already a conjunctive (not implicative) output-nothing to do
  if(strcmp(Out[outputnumber]->Defuzzify(), OUT_FUZZY::ImpFuzzyDefuz())) return -2;

  // A crisp output: nothing to do
  if(strcmp(Out[outputnumber]->GetOutputType(), OUT_FUZZY::OutputType())) return -3;

  //Fuzzy output but no MFs
  if ((Out[outputnumber]->GetNbMf())<=0) return -4;



  // Update inference operators
  if(DefuzType) Out[outputnumber]->SetOpDefuz(DefuzType);
  else Out[outputnumber]->SetOpDefuz(OUT_FUZZY::AreaDefuz());
  if(DisjType)  Out[outputnumber]->SetOpDisj(DisjType);
  else Out[outputnumber]->SetOpDisj(OUT_FUZZY::DisjMax());


  // transfPart=true: try turning Qsp into Sfp
  // else:            only test if output partition is Qsp
  if(transfPart)
    //Imp2ConjSfp returns
    // -1 if invalid outputnumber
    // -2 if already implicative output
    // -3 the output is not fuzzy
    // -4 if  Nmf=0
    // -5 if transf. into Sfp does not work
    // 0 if  transf. into Sfp has been done OK
    // 1 if transf. into SfpQsp has been done OK and if MFs have been re-sorted
    // 2 if output partition is already an Sfp partition

    ret=FIS2Sfp(outputnumber,DefuzType,DisjType);//InitPossibles called from there
  else
    {
      //test mode only - partition is not transformed- IsSfp returns 1 if output partition is Qsp, 0 if not.
      ret=((OUT_FUZZY *) Out[outputnumber])->IsQsp();
      // delete and recreate Possible, MFConcs and Classes arrays
      Out[outputnumber]->InitPossibles(Rule, NbRules,outputnumber);
    }



  return ret;
}

int FIS::FIS2Sfp(int outputnumber, const char * DefuzType, const char * DisjType)
  //*******************************************************************************
  //Imp2ConjSfp returns
  // -1 if invalid outputnumber
  // -3 the output is not fuzzy
  // -4 if  Nmf=0
  // -5 if transf. into Sfp does not work
  // 0 if  transf. into Sfp has been done OK
  // 1 if transf. into SfpQsp has been done OK and if MFs have been re-sorted
  // 2 if output partition is already an Sfp partition
{
  int *sorted, i, ret,testsfp;
  ret = 0;
  sorted = NULL;

  // invalid output number-nothing to do
  if((outputnumber < 0) || (outputnumber > NbOut-1)) return -1;

  // A crisp output: nothing to do
  if(strcmp(Out[outputnumber]->GetOutputType(), OUT_FUZZY::OutputType())) return -3;

  //Fuzzy output but no MFs
  if ((Out[outputnumber]->GetNbMf())<=0) return -4;

  // test if partition is already sfp
  testsfp=((FISIN*) Out[outputnumber])->IsSfp(sorted);
  if(testsfp) return 2;
  delete [] sorted;
  sorted = NULL;

  // Turn Qsp into Sfp
  if(! ((OUT_FUZZY *) Out[outputnumber])->Qsp2Sfp(sorted)) return -5;

  // Update rule conclusions according to the Sfp
  for(i = 0 ; i < NbRules ; i++ )
    {
      if((int) Rule[i]->GetAConc(outputnumber) % 2)
	Rule[i]->SetAConc(outputnumber, (Rule[i]->GetAConc(outputnumber) +1) / 2);
      else Rule[i]->SetAConc(outputnumber, 1);
    }

  if(sorted)
    // The MF have been re-ordered, update rule conclusions
    // according the new order within the sfp
    {
      for(i = 0 ; i < NbRules ; i++ )
	Rule[i]->SetAConc(outputnumber, sorted[(int)Rule[i]->GetAConc(outputnumber) - 1] + 1);
      ret = 1;
      delete [] sorted;
    }

  // Update inference operators
  if(DefuzType) Out[outputnumber]->SetOpDefuz(DefuzType);
  else Out[outputnumber]->SetOpDefuz(OUT_FUZZY::AreaDefuz());
  if(DisjType)  Out[outputnumber]->SetOpDisj(DisjType);
  else Out[outputnumber]->SetOpDisj(OUT_FUZZY::DisjMax());

  // delete and recreate Possible and Classes arrays
  Out[outputnumber]->InitPossibles(Rule, NbRules,outputnumber);

  // Ensure the range limits be inferred.
  ((OUT_FUZZY *) Out[outputnumber])->OutCoverage();

  return ret;
}

//! *******************************************stable rule module******************************************

//! Global variable used by qsort
int * OccurG;

//! Global function used by qsort
int CmpOccur(const void * a, const void *b)
  //***************************************
{
  if(OccurG[*(unsigned int *)a] > OccurG[*(unsigned int *)b]) return -1;
  if(OccurG[*(unsigned int *)a] < OccurG[*(unsigned int *)b]) return  1;
  return 0;
}


int MergeRules(const char * Fis1, const char *Fis2, const char * Merge, const char *Occur, double **&tabconc, int conc)
  //*******************************************************************************************
{
  FIS R(Fis1);
  FIS M(Fis2);
  int i,j,k, ret;
  FILE * FisMerged;
  FILE * fOccur;
  char buf[15];
  int *TabOccur;
  double *tmprule=NULL;
  double **tmpmat=NULL;

  if(R.GetNbIn() != M.GetNbIn())
    {
      sprintf( ErrorMsg, "~InMergeRules~, ~NbInMustBeTheSameInBothSystems~\n~Values~: %d %d\n",
	       R.GetNbIn(), M.GetNbIn());
      throw std::runtime_error( ErrorMsg );
    }

  if(R.GetNbOut() != M.GetNbOut())
    {
      sprintf( ErrorMsg, "~InMergeRules~, ~NbOutMustBeTheSameInBothSystems~\n~Values~: %d %d\n",
	       R.GetNbOut(), M.GetNbOut());
      throw std::runtime_error( ErrorMsg );
    }

  for(i = 0; i < R.GetNbIn(); i++)
    if(R.In[i]->GetNbMf() != M.In[i]->GetNbMf())
      {
	sprintf( ErrorMsg, "~InMergeRules~, ~NbMfMustBeTheSameInBothSystems~\n~ValuesForInput~ %d: %d %d\n",
		 i+1, R.In[i]->GetNbMf(), M.In[i]->GetNbMf());
	throw std::runtime_error( ErrorMsg );
      }

  for(i = 0; i < R.GetNbOut(); i++)
    if(R.Out[i]->GetNbMf() != M.Out[i]->GetNbMf())
      {
	sprintf( ErrorMsg, "~InMergeRules~, ~NbMfMustBeTheSameInBothSystems~\n~ValuesForOutput~ %d: %d %d\n",
		 i+1, R.Out[i]->GetNbMf(), M.Out[i]->GetNbMf());
	throw std::runtime_error( ErrorMsg );
      }

  TabOccur = new int [R.GetNbRule() + M.GetNbRule()];
  // if conclusions are not taken into account for the rule comparison, initialisation of the conclusions array
  if(tabconc==NULL && conc==false)
    {
      tabconc = new double *[R.GetNbRule()];
      for(i=0; i < R.GetNbRule(); i++) tabconc[i]= new double [1];
    }

  if(! (fOccur = fopen(Occur, "rt")))
    {
      for(i = 0; i < R.GetNbRule(); i++)
	{
	  TabOccur[i] = 1;
	  if(conc == false)
	    tabconc[i][0]=R.Rule[i]->GetAConc(0);
	}
    }
  else
    for(i = 0; i < R.GetNbRule(); i++)
      {
	char *tmp = fgets(buf, 15, fOccur);
	if(tmp != NULL) TabOccur[i] = atoi(buf);
      }
  if(fOccur) fclose(fOccur);

  for(i = 0; i < M.GetNbRule(); i++)
    {
      ret = R.RulePos(M.Rule[i], 0, conc);
      if(ret == -1)
	{
	  R.AddRule(M.Rule[i]);
	  TabOccur[R.GetNbRule()-1] = 1;
	  //increasing conclusion array size
	  if(conc==false)
	    {
	      tmpmat = new double *[R.GetNbRule()];
	      for(j=0;j<R.GetNbRule();j++)
		{
		  tmpmat[j]=new double [TabOccur[j]];
		  if(j!=(R.GetNbRule()-1))
		    for(k=0;k<TabOccur[j];k++) tmpmat[j][k]=tabconc[j][k];
		  else
		    tmpmat[j][0]=M.Rule[i]->GetAConc(0);
		}
	      for(j=0;j<R.GetNbRule()-1;j++) delete [] tabconc[j];
	      delete [] tabconc;
	      tabconc = new double *[R.GetNbRule()];
	      for(j=0;j<R.GetNbRule();j++)
		{
		  tabconc[j]=new double [TabOccur[j]];
		  for(k=0;k<TabOccur[j];k++) tabconc[j][k]=tmpmat[j][k];
		}
	      for(j=0;j<R.GetNbRule();j++) delete [] tmpmat[j];
	      delete [] tmpmat;
	    }
	}
      else
	{
	  TabOccur[ret] ++;
	  if(conc==false)
	    {
	      tmprule = new double [TabOccur[ret]];
	      for(j=0;j<TabOccur[ret]-1;j++) tmprule[j]=tabconc[ret][j];
	      delete [] tabconc[ret];
	      tmprule[TabOccur[ret]-1]=M.Rule[i]->GetAConc(0);
	      tabconc[ret]=new double [TabOccur[ret]];
	      for(j=0;j<TabOccur[ret];j++) tabconc[ret][j]=tmprule[j];
	      delete [] tmprule;
	    }
	}

    }

  // write temporary files
  FisMerged = fopen(Merge, "wt");
  R.PrintCfg(FisMerged);
  fclose(FisMerged);

  fOccur = fopen(Occur, "wt");
  FILE * fres = fopen("merge.res", "wt");
  for(i = 0; i < R.GetNbRule(); i++)
    fprintf(fOccur, "%d\n", TabOccur[i]);

  fclose(fOccur);
  fclose(fres);

  delete [] TabOccur;
  return 0;
}

int StableRules(char * firstfispart, int n, char * lastfispart, char * res, int & nR, double & meanO, double & stdO, int conc)
  //***********************************************************************************
{
  if(n < 2)
    {
      sprintf( ErrorMsg, "~InStableRules~, ~NbOfFisToBeAnalyzedLessThan2~: %d\n", n);
      throw std::runtime_error( ErrorMsg );
    }

  int i, j, nfis;
  unsigned int * sorted;
  char * f1, * f2;
  char buf[15];
  double median, min, max;
  double * tmp;
  double **concmat=NULL;
  FIS *s;
  FILE *f, *o;

  i = strlen(firstfispart);
  if(lastfispart) i += strlen(lastfispart);
  f1 = new char [i + 4];
  f2 = new char [i + 4];

  j = nfis = 0;

  // name of the two first fis files
  for(i = 0; i < n; i++)
    {
      if(lastfispart!=NULL)
	sprintf(f1, "%s.%d.%s", firstfispart, i, lastfispart);
      else
	sprintf(f1, "%s.%d", firstfispart, i);

      if((f = fopen(f1, "rt")) != NULL)
	{ fclose(f); break;}
    }
  j = ++i;

  for(i = j; i < n; i++)
    {
      if(lastfispart!=NULL)
	sprintf(f2, "%s.%d.%s", firstfispart, i, lastfispart);
      else
	sprintf(f2, "%s.%d", firstfispart, i);

      if((f = fopen(f2, "rt")) != NULL)
	{ fclose(f); break;}
    }
  if(i == n)
    {
      sprintf( ErrorMsg, "~InStableRules~, ~NbOfValidFisLessThan2~: %d\n", n);
      throw std::runtime_error( ErrorMsg );
    }
  j = ++i;

  remove("occur.tab");

  MergeRules(f1, f2, "merge.tmp", "occur.tab", concmat, conc);
  nfis = 2;

  for(i = j; i < n; i++)
    {
      if(lastfispart!=NULL)
	sprintf(f2, "%s.%d.%s", firstfispart, i, lastfispart);
      else
	sprintf(f2, "%s.%d", firstfispart, i);
      if((f = fopen(f2, "rt")) != NULL) fclose(f);
      else continue;
      MergeRules( "merge.tmp", f2, "merge.tmp",  "occur.tab", concmat, conc);
      nfis ++;
    }

  s = new FIS("merge.tmp");
  f = fopen(res, "wt");
  o = fopen("occur.tab", "rt");

  nR = s->GetNbRule();
  OccurG = new int [nR];
  sorted = new unsigned int [nR];
  // sorting rules by their decreasing number of occurence
  for(i = 0; i < nR; i++)
    {
      sorted[i] = i;
      char *tmp = fgets(buf, 15, o);
      if(tmp != NULL) {
    	  OccurG[i] = atoi(buf);
    	  s->Rule[i]->NbOccur = OccurG[i];
      }
    }

  qsort(sorted, s->GetNbRule(), sizeof(unsigned int), CmpOccur);

  // Print results : if distinct conclusions are taken into account,
  // each rule is followed by its conclusion
  // If distinct conclusions are not taken into account, each rule is followed
  // by the mean of the conclusions and the standard deviation
  for(i = 0; i < nR; i++)
    {
      fprintf(f, "%d, ", s->Rule[sorted[i]]->NbOccur);
      if(conc) s->Rule[sorted[i]]->Print(f);
      else
	{
	  s->Rule[sorted[i]]->PrintPrems(f);
	  StatArray(concmat[sorted[i]], s->Rule[sorted[i]]->NbOccur, 0, median, meanO, stdO, min, max, 0);
	  fprintf(f, "%f, %f \n", meanO, stdO);
	}
    }
  fprintf(f, "Number of valid fis %d \n", nfis);

  fclose(o); fclose(f);

  tmp = new double [nR];
  for(i = 0; i < nR; i++) tmp[i] = OccurG[i];
  // final result file
  StatArray(tmp, nR, 0, median, meanO, stdO, min, max, 0);

  if(concmat!=NULL)
    {
      for(i = 0; i < nR; i++) if(concmat[i]!=NULL) delete [] concmat[i];
      delete [] concmat;
    }

  delete [] tmp;
  delete [] OccurG;
  OccurG = NULL;
  delete [] sorted;
  delete s;
  delete [] f1;  delete [] f2;

  return 0;
}
//*********************************************** END FIS.CPP*********************************
