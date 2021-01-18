//***********************************************************************

//
//
//                              OUT.CPP
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC 
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : MF class functions used by FISPRO, part of library fispro

//**********************************************************************

#include "fis.h"

void FISOUT::Init( ifstream &f, int bufsize, int num, const char *OpDefuz, const char *OpDisj, int ccl, double VDef )
  //*****************************************************************************************************   
{
  try 
    { 
      Init( OpDefuz, OpDisj, ccl, VDef ); 
    }
  catch( std::exception &e )
    {
      sprintf( ErrorMsg, "~ErrorInFISFile~\n~Output~: %-3d\n%.100s", num, e.what() );
      throw std::runtime_error( ErrorMsg );
    }

  int impli=0;
 
  FISIN::Init( f, bufsize, num );
  // MFs allowed for implicative output : triangular, (semi)trapezoidal, universal, door
  impli= ! strcmp(Defuzzify(),OUT_FUZZY::ImpFuzzyDefuz());
  if (impli)
    CheckImpliMFs();//throws exception if forbidden MF types are used
  // end of implicative ouput case
}

void  FISOUT::CheckImpliMFs()
{
  int impli=0;
  // MFs allowed for implicative output : triangular, (semi)trapezoidal, universal, door
  impli= ! strcmp(Defuzzify(),OUT_FUZZY::ImpFuzzyDefuz());
  if (impli)
      for (int imf=0;imf<Nmf;imf++)
	CheckImpliMF(Fp[imf]); //throws exception if MF type not allowed for implicative output
}

//check candidate MF type
void  FISOUT::CheckImpliMF(MF* mfcandidate)
{
  int impli=0,allowedimpliMFs=0;
  // MFs allowed for implicative output : triangular, (semi)trapezoidal, universal, door
  impli= ! strcmp(Defuzzify(),OUT_FUZZY::ImpFuzzyDefuz());
  if (impli)
    {
      allowedimpliMFs=(!strcmp(mfcandidate->GetType(),"trapezoidal"))||(!strcmp(mfcandidate->GetType(),"triangular"));
      allowedimpliMFs=allowedimpliMFs||(!strcmp(mfcandidate->GetType(),"SemiTrapezoidalSup"))||(!strcmp(mfcandidate->GetType(),"SemiTrapezoidalInf"));
      allowedimpliMFs=allowedimpliMFs||(!strcmp(mfcandidate->GetType(),"universal"));
      allowedimpliMFs=allowedimpliMFs||(!strcmp(mfcandidate->GetType(),"door"));
      if (!allowedimpliMFs)
	{
	  sprintf( ErrorMsg, "ForbiddenMFshape~in~implicative~Systems");
	  throw std::runtime_error( ErrorMsg );
	}
    }
}

int FISOUT::operator != (const FISOUT & o)
  //**************************************
{
  if( (FISIN::operator!=(o)) ||
      (strcmp(GetOutputType(), o.GetOutputType()) !=0) ||
      (strcmp(Defuz, o.Defuz) != 0) ||
      (strcmp(Disj, o.Disj) != 0) ||
      (Default != o.Default) ||
      (Classif != o.Classif) )
    return 1;
  return 0;
}

void FISOUT::SetOpDefuz( const char *op_defuz )
  //*******************************************
{
  delete [] Defuz;
  Defuz = new char[strlen(op_defuz)+1];
  sprintf( Defuz, "%s", op_defuz );

  if(Def) delete Def;
  Def = NULL;

}

void FISOUT::Classification( int classif )
  //**************************************
{
  Classif = classif;
  char *tmp_defuz = new char[strlen(Defuzzify())+1];
  strcpy( tmp_defuz, Defuzzify());
  SetOpDefuz( tmp_defuz );
  delete []tmp_defuz;
}

void FISOUT::ReplaceMF( int sef_number, MF *new_sef )
  //************************************************
{
  // check that output is fuzzy
  if(! strcmp(GetOutputType(), OUT_CRISP::OutputType()))
    return;
  // check that new MF is allowed
  CheckImpliMF(new_sef); // throws exception if MF type is not allowed
  // now call FISIN method if MF type is allowed
  FISIN::ReplaceMF(sef_number,new_sef );
}

void FISOUT::AddMF( MF *newsef )
  //**************************
{
   // check that new MF is allowed
  CheckImpliMF(newsef); // throws exception if MF type is not allowed
  // now call FISIN method if MF type is allowed
  FISIN::AddMF(newsef);
}

void OUT_FUZZY::SetOpDefuz( const char *op_defuz )
  //**********************************************
{
  if( ! strcmp( op_defuz, OUT_FUZZY::AreaDefuz() ) || ! strcmp( op_defuz, OUT_FUZZY::MeanMaxDefuz())  || 
      ! strcmp( op_defuz, OUT_FUZZY::SugenoDefuz())|| ! strcmp( op_defuz, OUT_FUZZY::ImpFuzzyDefuz())  )
    {
      FISOUT::SetOpDefuz( op_defuz );
      if(!strcmp(Defuz, OUT_FUZZY::SugenoDefuz()))
	Def = new DEFUZ_SugenoFuzzy();
      else if(!strcmp(Defuz, OUT_FUZZY::AreaDefuz()))
	Def = new DEFUZ_WeArea();
      else if(!strcmp(Defuz, OUT_FUZZY::MeanMaxDefuz()))
	Def = new DEFUZ_MeanMax();
	  else if(!strcmp(Defuz, OUT_FUZZY::ImpFuzzyDefuz()))
	Def = new DEFUZ_ImpFuzzy();
	
    }  
  else
    {
      sprintf( ErrorMsg, "~Output~%.50s~:~Defuzzification~%.50s~NotAllowed~", GetOutputType(), op_defuz );
      throw std::runtime_error( ErrorMsg );
    }
}

void OUT_CRISP::SetOpDefuz( const char *op_defuz )
  //**********************************************
{
  if( !strcmp( op_defuz, OUT_CRISP::SugenoDefuz()) || ! strcmp( op_defuz, OUT_CRISP::MaxCrispDefuz()) )
    {
      FISOUT::SetOpDefuz( op_defuz );
      if(!strcmp(Defuz, OUT_CRISP::SugenoDefuz()))
	{
	  if(Classif) Def = new DEFUZ_SugenoClassif();
	  else Def = new DEFUZ_Sugeno();
	}
      else if(!strcmp(Defuz, OUT_CRISP::MaxCrispDefuz()))
	Def = new DEFUZ_MaxCrisp();
    }
  else
    {
      sprintf( ErrorMsg, "~Output~%.50s~:~Defuzzification~%.50s~NotAllowed~", GetOutputType(), op_defuz ); 
      throw std::runtime_error( ErrorMsg );
    }
}

void FISOUT::SetOpDisj( const char *op_disj )
  //*****************************************
{
  if (Disj != NULL) delete [] Disj;
  Disj = new char[strlen(op_disj)+1];
  sprintf( Disj, "%s", op_disj );

  if(Ag!=NULL) delete Ag;
  Ag = NULL;
}

void OUT_FUZZY::SetOpDisj( const char *op_disj ) //throw( std::runtime_error )
  //********************************************
{
  if( ! strcmp( op_disj, OUT_FUZZY::DisjSum() ) || ! strcmp( op_disj, OUT_FUZZY::DisjMax() ) ||
      ! strcmp( op_disj, OUT_FUZZY::DisjIgg() ) || ! strcmp( op_disj, OUT_FUZZY::DisjIgd() ) ||
      ! strcmp( op_disj, OUT_FUZZY::DisjIrg() ) )
    {
      FISOUT::SetOpDisj( op_disj );
      if(!strcmp(Disj, OUT_FUZZY::DisjSum()))
	Ag = new AGGREGSUM();
      if(!strcmp(Disj, OUT_FUZZY::DisjMax()))
	Ag = new AGGREGMAX();
	  if(!strcmp(Disj, OUT_FUZZY::DisjIgd()))
	Ag = new AGGREGIMP(new IMPLIGD());
	  if(!strcmp(Disj, OUT_FUZZY::DisjIrg()))
	Ag = new AGGREGIMP(new IMPLIRG());
	  if(!strcmp(Disj, OUT_FUZZY::DisjIgg()))
	Ag = new AGGREGIMP(new IMPLIGG());	
		
	}
  else
    {
      sprintf( ErrorMsg, "~Output~%.50s~:~Disjunction~%.50s~NotAllowed~", GetOutputType(), op_disj );
      throw std::runtime_error( ErrorMsg );
    }
}

void OUT_CRISP::SetOpDisj( const char *op_disj ) //throw( std::runtime_error )
  //********************************************
{
  if( ! strcmp( op_disj, OUT_CRISP::DisjSum() ) || ! strcmp( op_disj, OUT_CRISP::DisjMax() ) )
    {
      FISOUT::SetOpDisj( op_disj );
      if(!strcmp(Disj, OUT_CRISP::DisjSum()))
	Ag = new AGGREGSUM();
      else if(!strcmp(Disj, OUT_CRISP::DisjMax()))
	Ag = new AGGREGMAX();
    }
  else
    {
      sprintf( ErrorMsg, "~Output~%.50s~:~Disjunction~%.50s~NotAllowed~", GetOutputType(), op_disj );
      throw std::runtime_error( ErrorMsg );
    }
}
void FISOUT::DeleteMFConc(int numberofrules) const
  //********************************************************
{
  //delete MFConc array
  if  (MfConc !=NULL)
    {
      for (int irule=0;irule<numberofrules;irule++)
	{
	  delete MfConc[irule];
	  MfConc[irule]=NULL;
	}
    }
}
void FISOUT::DeleteMFConcArray()
  //********************************************************
{
  delete [] MfConc;
  MfConc=NULL;
}

void FISOUT::InitPossibles(RULE ** r, int nr, int nb)
  //*************************************************
{                                                   
#define IMPOSSIBLE (-INFINI-0.0005)  //  conclusion label = impossible
  double *val;
  int tmp, j, k;
  bool impli=false;

  //PROTECTION in case of no (active) rules
  if(! active || nr <= 0) return;  
  
  DeletePossibles(nr);
  tmp = 0;
  
  val = new double [nr];
  for(j = 0; j < nr; j ++) val[j] = IMPOSSIBLE;
  
  for(j = 0; j < nr; j ++)
    {
      for(k = 0; k < tmp; k++)
	if(fabs(r[j]->GetAConc(nb) - val[k]) < EQUAL_CONC) break;
      if(k == tmp)
	{
	  val[tmp] = r[j]->GetAConc(nb);
	  tmp++;
	}
    }  

  qsort(val, tmp, sizeof(double), CmpDblAsc);

  NbPossibles = tmp;  
  Possibles = new double [NbPossibles];
  for(k = 0; k < NbPossibles; k++)        Possibles[k] = val[k];
  delete [] val;
  
  MuInfer = new double [NbPossibles];
  RuleInfer = new int [NbPossibles];
  ConcInfer = new int [nr];

  DeleteMFConcArray();//delete array of ptrs-assume that MFConc have already been deleted
  impli= ! strcmp(Defuzzify(),OUT_FUZZY::ImpFuzzyDefuz());
  if (impli)
    {
      MfConc = new MFDPOSS*[nr]; 
      for(k = 0; k < nr; k++) MfConc[k] =NULL;
    }
  else
    MfConc=NULL;
    
  MfGlob =NULL;
  
  InitTabRes();

  for(k = 0; k < nr; k++)
    { 
  
      if(! r[k]->IsActive() ) continue; 
      for(j = 0; j < NbPossibles; j++)
	if(fabs(r[k]->GetAConc(nb) - Possibles[j]) < EQUAL_CONC)
	  {
	    ConcInfer[k] = j;
	    break;
	  }
      if(j == NbPossibles) 
	{
	  sprintf( ErrorMsg, "~ErrorInInitPossibles~\n~Output~: %50s\n", Name);
	  throw std::runtime_error( ErrorMsg );
	}
    }
  
}
void FISOUT::UpdatePossibles(RULE ** r, int nr, int changedrulenum, int nb)
  //***********************************************************************
{                                                   
  int j;
  int concexist=0;

  // inactive output or no rules
  if(! active || nr <= 0) return;  
  // PROTECTION (changedrulenum)
  if (changedrulenum<0 || changedrulenum>=nr) return;
  // examine all possible values to see 
  // if new conclusion is already implemented
  for(j = 0; j < NbPossibles; j++)
    if(fabs(r[changedrulenum]->GetAConc(nb) - Possibles[j]) < EQUAL_CONC)
      {
	ConcInfer[changedrulenum] = j;
	concexist=1;
	break;
      }
  // it is a new conclusion, all arrays must be deleted and reallocated
  if (concexist==0)   InitPossibles(r,nr,nb);   
}

void FISOUT::DeletePossibles()
  //*************************************
{
  if( Possibles != NULL ) delete [] Possibles;
  if( MuInfer != NULL ) delete [] MuInfer;
  if( RuleInfer != NULL ) delete [] RuleInfer;
  if( ConcInfer != NULL ) delete [] ConcInfer;
  
  Possibles = NULL;
  MuInfer = NULL;
  RuleInfer = NULL;
  ConcInfer = NULL;
  NbPossibles = 0;  
}

void FISOUT::DeletePossibles(int nbrules)
  //*************************************
{
  if( Possibles != NULL ) delete [] Possibles;
  if( MuInfer != NULL ) delete [] MuInfer;
  if( RuleInfer != NULL ) delete [] RuleInfer;
  if( ConcInfer != NULL ) delete [] ConcInfer;
  
  Possibles = NULL;
  MuInfer = NULL;
  RuleInfer = NULL;
  ConcInfer = NULL;
  NbPossibles = 0;  
  // for implicative output ONLY
  if  (MfConc !=NULL)
    {
      for (int irule=0;irule<nbrules;irule++)
	{
	  delete MfConc[irule];
	  MfConc[irule]=NULL;
	}
    }
  delete [] MfConc;
  MfConc=NULL;
  
  if( MfGlob != NULL ) delete MfGlob;
  MfGlob=NULL;
  
}
  

void OUT_FUZZY::InitDiscrete(double *t, int n, double min, double max)
  //******************************************************************
{
  int i;

  SetRange( min, max );
  Nmf = n;
  if(Nmf == 0)
    return;
  Fp = new MF *[Nmf];
  for( int i=0 ; i<Nmf ; i++ )
    Fp[i] = NULL;
  //Mfdeg = new double [Nmf];
  for(i = 0; i < Nmf; i++)
    Fp[i] = new MFDISCRETE(t[i]);
    
}

void OUT_FUZZY::OutCoverage(void)
  //*****************************
{
  double *p, *pe;

  if(Nmf < 2) return;
  if(strcmp(Fp[0]->GetType(), MFTRAPINF::Type()) ||
     strcmp(Fp[Nmf-1]->GetType(), MFTRAPSUP::Type()) )
    {
      sprintf(ErrorMsg, "~ErrorInOutCoverage~~InOutput~%50s\n~PartitionEndShouldBeSemitrapezoidalShaped~", Name);
      throw std::runtime_error( ErrorMsg );
    }
  p = new double [3];
  pe = new double [3];

  Fp[0]->GetParams(p);
  Fp[Nmf-1]->GetParams(pe);

  if(ValInf > p[1] || ValSup < pe[1])
    {
      sprintf(ErrorMsg, "~ErrorInOutCoverage~~InOutput~%50s\n~UnreachableTarget~,~BothValinfAndValsupMustBelongToTheKernels", Name);
      throw std::runtime_error( ErrorMsg );
    }

  if(! strcmp(Defuz, OUT_FUZZY::MeanMaxDefuz()) || ! strcmp(Defuz, OUT_FUZZY::SugenoDefuz()) )
    {
      // min
      MF *new_mf = new MFTRAPINF(ValInf * 2 - p[1], p[1], p[2]);
      new_mf->SetName(Fp[0]->Name);
      ReplaceMF(0, new_mf);
      // max
      new_mf = new MFTRAPSUP(pe[0], pe[1], ValSup * 2 - pe[1]);
      new_mf->SetName(Fp[Nmf - 1]->Name);
      ReplaceMF(Nmf - 1, new_mf);
    }

  else if(! strcmp(Defuz, OUT_FUZZY::AreaDefuz()) )  
    {
      // min
      // ValInf = [(b-a)(a+(b-a)/2) + ((c-b)/2)(b+(c-b)/3)]/[(b-a) + (c-b)/2]
      double k, root;

      k = (p[2] - p[1]) * (ValInf - p[1] - (p[2]-p[1])/3) + ValInf*p[1]*2 - p[1]*p[1];
      // a^2 + -2*ValInf*a + k = 0
      root = sqrt(ValInf*ValInf - k);
      MF *new_mf = new MFTRAPINF(ValInf-root, p[1], p[2]);
      new_mf->SetName(Fp[0]->Name);
      ReplaceMF(0, new_mf);

      // max
      // ValSup = [((b-a)/2)*(a+(b-a)*2/3) + (c-b)*(b+(c-b)/2)]/[(b-a)/2 + (c-b)]
      k = (ValSup - pe[0] - 2.*(pe[1]-pe[0])/3) * (pe[1]-pe[0]) - ValSup*pe[1]*2 + pe[1]*pe[1];
      // c^2 - ValSup*2*c - k = 0
      root = sqrt(ValSup*ValSup + k);
      new_mf = new MFTRAPSUP(pe[0], pe[1], ValSup+root);
      new_mf->SetName(Fp[Nmf - 1]->Name);
      ReplaceMF(Nmf - 1, new_mf);
    }

  delete [] p; delete [] pe;
}
double OUT_FUZZY::SymbMatch(double * muObs, double * muInf, int nOutMF, double mutObs,double mutInf,FILE *display)
{
  // nOutMF number of output MFs
  // muInf=maximum membership of inferred dposs to each output MF
  // comment calculer ?????????????????????????????????????????????????????????????
  // muObs=membership of observed crisp value to each output MF
  // non symmetrical match
  // first examine k= number of intersected kernels by Inferred dposs, 
  // then m=number of covered MFs by Inferred dposs (beyond a threshold)
  // returns -1 if result is too imprecise (k>=2)
  //         0<val<1 otherwise
  // if m=1, val=1 if match otherwise val=0
  // if m=2 
  double val=0.0,lk,rk;
  int imf,k=0,m=0, count;
  int * infMF=NULL, *obsMF=NULL;
  MFDPOSS * kerneldposs=NULL, *inters =NULL;;
  ACUT * acut=NULL;
  
  
  if (Nmf<=0) return val;//case no MFs
  
  infMF=new int[Nmf];
  obsMF=new int[Nmf];
  
  if (MfGlob == NULL) return val;// case of null dposs output
  if (display)
    {
      fprintf(display, "\nIn symbmatch mutObs=%g\tmutInf=%g\tMFGlob:\n",mutObs,mutInf);
      MfGlob->Print(display);
      for (imf=0;imf<Nmf;imf++)
	{
	  fprintf(display, "muObs[%d]=%g\t",imf,muObs[imf]);
	  fprintf(display, "muInf[%d]=%g\t",imf,muInf[imf]);
	}
    }
  

  for (imf=0;imf<Nmf;imf++)
    {
      if (muObs[imf]>=mutObs)
	obsMF[imf]=1;
      else
	obsMF[imf]=0;
    }
  for (imf=0;imf<Nmf;imf++)
    {
      // compute k=number of intersected kernels
      // for each output MF compute its kernel
      Fp[imf]->Kernel(lk,rk);
      acut=new ACUT(lk,rk,1.0);
      kerneldposs=new MFDPOSS(acut);
      inters=MfGlob->Inter(kerneldposs);
      if (inters !=NULL)
	k=k+1;
      delete kerneldposs;
      delete acut;
      delete inters;
      kerneldposs=NULL;
      acut=NULL;
      inters=NULL;
      //compute m=number of intersected MFs with membership > threshold
      if (muInf[imf]>=mutInf)
	{
	  m=m+1;
	  infMF[imf]=1;
	}
      else
	infMF[imf]=0;
    }
  if (display)
    {
      fprintf(display, "\nIn symbmatch #intersected output MFs (with threshold %g) m=%d #intersected output MF kernels k=%d",mutInf,m,k);
      for (imf=0;imf<Nmf;imf++)
	fprintf(display, "\nobsMF[%d]=%d\tinfMF[%d]=%d",imf,obsMF[imf],imf,infMF[imf]);
    }
  
  
  val=0.0;
  count=0;

  if (k>=2) 
    val=-1.0;//imprecise inferred output
  else
    switch (m) 
      {
      case 1:
	{
	  // single MF intersected by MFGlob-match/not match with observed value labels
	  
	  for (imf=0;imf<Nmf;imf++)
	    {
	      if ((obsMF[imf]==1) && (infMF[imf]==1))
		  val=1.0;
	    }
	  break;
	}
      case 2:
	{
	  // 2 MFs intersected by MFGlob - count number of matches with observed value labels
	   for (imf=0;imf<Nmf;imf++)
	     {
	       if ((obsMF[imf]==1) && (infMF[imf]==1))
		 count++;
	     }
	   if (count==2)
	     val=1.0;
	   if (count==1)
	     val=0.75;
	   
	  break; 
	}
      case 3:
	{
	  // 3 MFs intersected by MFGlob - count number of matches with observed value labels
	  for (imf=0;imf<Nmf;imf++)
	    if ((obsMF[imf]==1) && (infMF[imf]==1))
	      count++;
	  if (count==3)
	     val=1.0;
	  if (count==2)
	    val=0.75;
	  if (count==1)
	    val=0.5;
	   break;
	}
      default:
	{
	  val=0.0;
	  break;
	}
     
      }
if (display) fprintf(display, "\nEnd of Symbmatch val=%g\n",val);
 delete [] infMF;
 delete [] obsMF;
return val;	  
}

int OUT_FUZZY::Sfp2Qsp(int *& s)
//*******************************
// Sfp2Qsp returns -1 if Nmf=0 or 1, -2 if output partition is not SFP
  // and it modifies its argument (sorted) if necessary + if returns 0
{
 
  int i, j;
  double *par; 
  char* nameMF=NULL;

  if(Nmf == 0) return -1; 
  if(Nmf == 1) return -1; 
  if(! IsSfp(s)) return -2;
 
  MF **tmp = new MF * [Nmf*2-1];
  //double *TmpDeg = new double [2*Nmf-1];
  nameMF=new char[15];

  // NewNum = OldNum * 2
  par = new double[4];
  j = 0;

  for(i = 0; i < Nmf-1; i++)
    {
      Fp[i]->GetParams(par);
      if(i) tmp[j++] = Fp[i]->Clone();
      // Set semi-trapezoidal limits to range
      else tmp[j++] = new MFTRAPINF(ValInf, par[1], par[2]);

      if(! strcmp(Fp[i]->GetType(), MFTRAP::Type())) 
	tmp[j++] = new MFTRI(par[2], (par[2]+par[3])/2, par[3]);
      else tmp[j++] = new MFTRI(par[1], (par[1]+par[2])/2, par[2]);
    }

  // Set semi-trapezoidal limits to range
  Fp[i]->GetParams(par);
  tmp[j++] = new MFTRAPSUP(par[0], par[1], ValSup );

  //for (i=0;i < j; i++) tmp[i]->PrintCfg(i);

  delete [] par;  
  for (i=0;i < Nmf; i++)
    {
      delete Fp[i];
      Fp[i]=NULL;
    }
  delete [] Fp;
  Fp=NULL;
  //delete [] Mfdeg;
  //Mfdeg=NULL;
  
  Nmf = j;
  Fp = tmp;
  //Mfdeg = TmpDeg;  
  Mfdeg().resize(Nmf);
  for (i=0;i < Nmf; i++) 
    {
      // MFs names
      if (i<1000)
	{
	  sprintf(nameMF,"MF%d",i+1);
	  Fp[i]->SetName(nameMF);
	}
      else
	Fp[i]->SetName("MF");
      
    }
  delete [] nameMF;
  return 0;
}

bool OUT_FUZZY::IsQsp()
//********************************************
//TEST IF PARTITION IS qsp
{
  int * s=NULL;
  if(Nmf == 1) return true;
  return Qsp2Sfp(s,true);
}
  



bool OUT_FUZZY::Qsp2Sfp(int *& s, bool onlytestQSP)
  //***********************************************
{

  if(Nmf == 0) return false; 
  if(Nmf == 1) return false; 
  if(! (Nmf%2)) return false; // Nmf must be odd

  FISIN TmpIn(*this);
  int i, j;
  MF **tmp = new MF * [(Nmf+1)/2];
  //double *TmpDeg = new double [(Nmf+1)/2];
  bool testSFP;
  j = 0;

  // Only even MF are kept
  for(i = 0; i < Nmf-1; i+=2) tmp[j++] = Fp[i]->Clone();
  tmp[j++] = Fp[Nmf-1]->Clone();
  for (i=0;i < Nmf; i++)
    {
      delete Fp[i];
      Fp[i]=NULL;
    }
  delete [] Fp;
  Fp=NULL;
   
  //delete [] Mfdeg;
  //Mfdeg=NULL;

  Nmf = j;
  Fp = tmp;
  //Mfdeg = TmpDeg;  
  Mfdeg().resize(Nmf);
  //for (i=0;i < Nmf; i++) Fp[i]->SetName("");

  //
  testSFP=IsSfp(s);
  //
  if (! testSFP||onlytestQSP) // restore if transformed partition is not SFP or if testQSP argument
    {
      for (i=0;i < Nmf; i++)
	{
	  delete Fp[i];
	  Fp[i]=NULL;
	}
      delete [] Fp;
      Fp=NULL;
   
      //delete [] Mfdeg;
      //Mfdeg=NULL;

      Nmf = TmpIn.GetNbMf();
      Fp = new MF * [Nmf];
      //Mfdeg = new double [Nmf];  
      Mfdeg().resize(Nmf);
      for(i = 0; i < Nmf; i++) Fp[i] = (TmpIn.GetMF(i))->Clone();
     
    }
  if (testSFP)
    // the transformed partition is SFP so the initial one was QSP  
    return true;
  else
    return false;
}

//**************************    OUT.CPP   *******************************

