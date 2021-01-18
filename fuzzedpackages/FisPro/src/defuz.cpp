//***********************************************************************
//
//
//                              DEFUZ.CPP
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : DEFUZ class functions used by FISPRO, part of library fispro

//**********************************************************************

#include "fis.h"

void DEFUZ::GetMax(FISOUT * O, double & max1, double & max2, int & i_max1, int & i_max2) const
//************************************************************************************
{
	int i;

	max1 = max2 = -1.;
	i_max1 = i_max2 = -1;

	for(i = 0; i < O->GetNbPossibles(); i++)
	{
		if(O->MuInfer[i] && (O->MuInfer[i] > max1 - Thres) )
		{
			if(max1 == -1)
			{
				max1 = O->MuInfer[i];
				i_max1 = i;
			}
			else if(O->MuInfer[i] > max1)
			{
				max2 = max1;
				i_max2 = i_max1;
				max1 = O->MuInfer[i];
				i_max1 = i;
			}
			else if(O->MuInfer[i] <= max1)
			{
				max2 = O->MuInfer[i];
				i_max2 = i;
			}
		}
	}

	if(max1 - max2 > Thres)
	{ max2 = -1; i_max2 = -1; }
}

double DEFUZ_Sugeno::EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display)
//*******************************************************************************
{
	double ret, wgtsum, pondsum;
	int i;
	double * ValConc = O->GetPossibles();

	Alarm = NOTHING;
	wgtsum = 0.0;
	pondsum = 0.0;
	for(i = 0; i < O->GetNbPossibles(); i++)
	{
		wgtsum += O->MuInfer[i];
		pondsum += O->MuInfer[i] * ValConc[i];
	}

	ret = O->DefaultValue();
	if(wgtsum) ret = pondsum / wgtsum;
	else Alarm = NO_ACTIVE_RULE;

	if(display)  fprintf(display, "Inferred output:  %f Alarm: %d\n", ret, Alarm);
	if(fa)
	{
		fprintf(fa, FORMAT_DOUBLE, ret);
		fprintf(fa, "%5d", Alarm);
	}
	return ret;
}  // End of DEFUZ_Sugeno::EvalOut()


DEFUZ_SugenoClassif::DEFUZ_SugenoClassif()
//**************************************
: DEFUZ_Sugeno()
{
	//! The out inferred by the DEFUZ_SugenoClassif is a class label (double).
	//! The possible class labels can be specified using the function SetClasses.
	//! The function FIS::Performance initializes the class labels with the distinct
	//! values of the data file using the function InitClasses(double *T, int n).
	//! The attribute Classes should not be mixed up with FISOUT::Possibles.
	//! initialized with the rule conclusions.
	Classes = NULL;

	Thres = AMBIGU_THRES;
}

double DEFUZ_SugenoClassif::EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display)
//**************************************************************************************
{
	double val, min, max, *tmp;
	int i, i_min;

	val = DEFUZ_Sugeno::EvalOut(TabR, NbR, O, fa, display);
	if(Classes == NULL)
		throw std::runtime_error("Classes non initialized in object DEFUZ_SugenoClassif");

	if(Alarm == NO_ACTIVE_RULE)
	{
		if(fa)
		{
			fprintf(fa, FORMAT_DOUBLE, val);
			fprintf(fa, "%5d", Alarm);
		}
		return val;
	}
	tmp = new double [NbClas];

	min = INFINI;
	max = -INFINI;
	i_min = -1;

	for(i = 0; i < NbClas; i++)
	{
		tmp[i] = fabs(val - Classes[i]);
		if(tmp[i] < min )
		{
			min = tmp[i];
			i_min = i;
		}
		if(tmp[i] > max ) max = tmp[i];
	}

	if(i_min == -1) val = O->DefaultValue();
	else  val = Classes[i_min];

	if(i_min > -1)
	{
		min = INFINI;
		for(i = 0; i < NbClas; i++)
		{
			if(i == i_min) continue;
			if(tmp[i] <= min ) min = tmp[i];
		}
		if(((min - tmp[i_min]) / (max -  tmp[i_min])) <= Thres) Alarm = AMBIGUITY;
	}

	if(display)  fprintf(display, "Inferred class label %f Alarm: %d \n", val, Alarm);
	if(fa)
	{
		fprintf(fa, FORMAT_DOUBLE, val);
		fprintf(fa, "%5d", Alarm);
	}

	delete [] tmp;
	return val;
}

DEFUZ_SugenoFuzzy::DEFUZ_SugenoFuzzy() : DEFUZ()
//********************************************
{
	Consequences = NULL;
	Thres = 0.;       // Not used
}

void DEFUZ_SugenoFuzzy::InitConsequences(FISOUT * O)
//************************************************
{
	int i, nmf;
	double b, e;

	delete [] Consequences;

	nmf = O->GetNbMf();
	Consequences = new double [nmf];

	for(i = 0; i < nmf; i++)
		Consequences[i] = O->Kernel(i, b, e);
}

double DEFUZ_SugenoFuzzy::EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display)
//************************************************************************************
{
	double ret, wgtsum, pondsum;
	int i;
	double * ValConc = O->GetPossibles();

	Alarm = NOTHING;
	if(Consequences == NULL) InitConsequences(O);

	wgtsum = 0.0;
	pondsum = 0.0;
	for(i = 0; i < O->GetNbPossibles(); i++)
	{
		wgtsum += O->MuInfer[i];
		pondsum += O->MuInfer[i] * Consequences[(int) ValConc[i]-1];
	}

	ret = O->DefaultValue();
	if(wgtsum) ret = pondsum / wgtsum;
	else Alarm = NO_ACTIVE_RULE;

	if(display)  fprintf(display, "Inferred output %f Alarm %d\n", ret, Alarm);
	if(fa)
	{
		fprintf(fa, FORMAT_DOUBLE, ret);
		fprintf(fa, "%5d", Alarm);
	}

	if(O->Classification())
	{
		O->GetDegsV(ret);
		if(fa) 	for(i = 0; i < O->GetNbMf(); i++) fprintf(fa, FORMAT_DOUBLE, O->Mfdeg()[i]);
	}

	return ret;
}

void DEFUZ_SugenoFuzzy::WriteHeader(FILE *p, FISOUT *O) const
//***************************************************
{
	int i;
	fprintf(p, "     %s", INFERRED);
	fprintf(p, "     %s", ALARM);

	if(O->Classification())
		for(i = 0; i < O->GetNbMf(); i++) fprintf(p, "      MF%d", i+1);
}

DEFUZ_WeArea::DEFUZ_WeArea() : DEFUZ()
//**********************************
{
	Thres = MIN_THRES;
}

double DEFUZ_WeArea::EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display)
//*******************************************************************************
{
	Trapeze *ct;
	double ret, center, card, massum, cog;
	double begi, endi, begj, endj;
	int i, j;

	ct = new Trapeze;
	ct->lk = ct->rk = ct->ls = ct->rs = -1.;
	double * ValConc = O->GetPossibles();

	Alarm = NOTHING;
	massum = 0.0;
	cog = 0.0;

	for(i = 0; i < O->GetNbPossibles(); i++)
	{
		O->Centroid((int) ValConc[i] - 1,  O->MuInfer[i], center, card, ct) ;
		massum += card;
		cog += card * center;
		if(display) fprintf(display, "MF %d  : Weight %f Mass %f cog %f  Trapeze Kernel : %f %f Support : %f %f \n", i+1, O->MuInfer[i], card, center, ct->lk, ct->rk, ct->ls, ct->rs);
	}

	if(massum) ret = cog / massum;
	else
	{
		Alarm = NO_ACTIVE_RULE;
		ret = O->DefaultValue();
	}

	double npos = O->GetNbPossibles();
	double nomf = O->GetNbMf();
	int k;
	for(i = 0; i < nomf - 1; i++)
	{
		for(k = 0; k < npos - 1; k++)
			if(i == (int) ValConc[k] - 1) break;
		if(k == npos -1 || O->MuInfer[k] < Thres) continue;
		O->Support(i, begi, endi);
		for(j = i+1; j < nomf; j++)
		{
			for(k = 0; k < npos; k++)
				if(j == (int) ValConc[k] - 1) break;
			if(k == npos || O->MuInfer[k] < Thres) continue;
			O->Support(j, begj, endj);
			if((endi - begj) < EPSILON)
				Alarm = NON_CONNEX_AREA;
			else break;
			// The area is connex from MF #i. Increment i.
		}
	}

	if(display)  fprintf(display, "Inferred output %f Alarm %d\n", ret, Alarm);
	if(fa)
	{
		fprintf(fa, FORMAT_DOUBLE, ret);
		fprintf(fa, "%5d", Alarm);
	}

	if(O->Classification())
	{
		O->GetDegsV(ret);
		if(fa) 	for(i = 0; i < O->GetNbMf(); i++) fprintf(fa, FORMAT_DOUBLE, O->Mfdeg()[i]);
	}

	delete ct;
	return ret;
}  // End of DEFUZ_WeArea::EvalOut()

void DEFUZ_WeArea::WriteHeader(FILE *p, FISOUT *O) const
//**********************************************
{
	int i;
	fprintf(p, "     %s", INFERRED);
	fprintf(p, "     %s", ALARM);

	if(O->Classification())
		for(i = 0; i < O->GetNbMf(); i++) fprintf(p, "      MF%d", i+1);
}

DEFUZ_MeanMax::DEFUZ_MeanMax(void)
//******************************
: DEFUZ()
{
	Thres = EQUALITY_THRES;
}

double DEFUZ_MeanMax::EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display)
//********************************************************************************
{
	double max_1, max_2, d1, d2, ker_r, ker_l, val;
	int i, i_max_1, i_max_2;

	Trapeze * ct;
	double * ValConc = O->GetPossibles();

	ct = new Trapeze;
	ct->lk = ct->rk = ct->ls = ct->rs = -1.;

	Alarm = NOTHING;
	val = O->DefaultValue();

	GetMax(O, max_1, max_2, i_max_1, i_max_2);

	// No rule is activated
	if(max_1 == -1)
		Alarm = NO_ACTIVE_RULE;

	// There is no ambiguity
	else if(max_2 == -1)
	{
		O->Centroid((int) ValConc[i_max_1] - 1, max_1, d1, d2, ct) ;
		val =  (ct->rk -  ct->lk) /2. + ct->lk;
	}

	// Is the segment connex?
	else
	{
		O->Centroid((int) ValConc[i_max_1] - 1, max_1, d1, d2, ct) ;
		ker_r = ct->rk;
		ker_l = ct->lk;

		O->Centroid((int) ValConc[i_max_2] - 1, max_2, d1, d2, ct) ;
		if((ct->lk - ker_r > EPSILON)  || (ker_l - ct->rk > EPSILON))
		{
			val =  (ker_r -  ker_l) /2. + ker_l;
			Alarm = NON_CONNEX_SEGMENT;
		}
		else
		{
			if(ct->lk < ker_r) val = (ct->rk - ker_l) / 2. + ker_l;
			else // ker_l < ct->rk
				val = (ker_r - ct->lk) / 2. + ct->lk;
		}
	}

	delete ct;

	if(display)  fprintf(display, "Inferred output %f Alarm %d\n", val, Alarm);
	if(fa)
	{
		fprintf(fa, FORMAT_DOUBLE, val);
		fprintf(fa, "%5d", Alarm);
	}

	if(O->Classification())
	{
		//      O->GetDegsV(val);
		if(fa)
			for(i = 0; i < O->GetNbMf(); i++)
				fprintf(fa, FORMAT_DOUBLE, O->MuInfer[i]);
	}

	return val;
}

void DEFUZ_MeanMax::WriteHeader(FILE *p, FISOUT *O) const
//***********************************************
{
	int i;
	fprintf(p, "     %s", INFERRED);
	fprintf(p, "     %s", ALARM);

	if(O->Classification())
		for(i = 0; i < O->GetNbMf(); i++) fprintf(p, "      MF%d", i+1);
}

double DEFUZ_MaxCrisp::EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display)
//*********************************************************************************
{
	double max_1, max_2, val;
	int i_max_1, i_max_2;
	double *p;

	p = O->GetPossibles();

	Alarm = NOTHING;

	GetMax(O, max_1, max_2, i_max_1, i_max_2);

	// No rule is activated
	if(max_1 == -1)
	{
		Alarm = NO_ACTIVE_RULE;
		val = O->DefaultValue();
	}

	else
	{
		val =  p[i_max_1];
		if(max_2 != -1 && i_max_2 != i_max_1)  Alarm = AMBIGUITY;
	}

	if(display)  fprintf(display, "Inferred output %f Alarm %d\n", val, Alarm);
	if(fa)
	{
		fprintf(fa, FORMAT_DOUBLE, val);
		fprintf(fa, "%5d", Alarm);
	}

	if(O->Classification())
	{
		int i;
		if(fa)
			for(i = 0; i < O->GetNbPossibles(); i++)
				fprintf(fa, FORMAT_DOUBLE, O->MuInfer[i]);
	}

	return val;
}

void DEFUZ_MaxCrisp::WriteHeader(FILE *p, FISOUT *O) const
//***********************************************
{
	int i;
	if(p != NULL)
	{
		fprintf(p, "     %s", INFERRED);
		fprintf(p, "    %s", ALARM);

		if(O->Classification())
			for(i = 0; i < O->GetNbPossibles(); i++) fprintf(p, "      MF%d", i+1);
	}
}


double DEFUZ_ImpFuzzy::EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display)
  //*********************************************************************************
{
  double val = O->DefaultValue();
  double mulk=0.0,murk=0.0;
  double ls, rs, lk, rk;
  double maxp=1.0;
  Alarm = NOTHING;
  // writes implicative inf. results in perf.res file
  // mean(kernel), MF# from output partition which are covered for mu>museuil, symbolic match result
  // minkernel, maxkernel, minsupport, maxsupport
  if(O->MfGlob != NULL)
    {
      //fprintf(f, FORMAT_DOUBLE,
      O->MfGlob->AlphaKernel(lk,rk,maxp);
      val = 0.5*(lk+rk);//middle of kernel=defuzzified value
      O->MfGlob->Support(ls,rs);
      if (fa)
	{
	  // defuzzified output
	  fprintf(fa,FORMAT_DOUBLE,val);
	  fprintf(fa, "%5d", Alarm);

	  // covered output MFs : number, or -1 if MF not covered above threshold
	  //mfi = mf1->Inter(mf2);

	  for(int i = 0; i < O->GetNbMf(); i++)
	    {
	      mulk=O->GetADeg(i,lk);
	      murk=O->GetADeg(i,rk);
	      if (murk > mulk)
		fprintf(fa,FORMAT_DOUBLE, murk);
	      else
		fprintf(fa,FORMAT_DOUBLE, mulk);
	    }

	  //kernel
	  fprintf(fa,FORMAT_DOUBLE,lk);
	  fprintf(fa,FORMAT_DOUBLE,rk);
	  //support
	  fprintf(fa,FORMAT_DOUBLE,ls);
	  fprintf(fa,FORMAT_DOUBLE,rs);
	}
    }//end dposs not null case
  else
    if(fa) 
      {
	// null dposs- default output & alarm
	fprintf(fa,FORMAT_DOUBLE,val);
	fprintf(fa, "%5d", Alarm);
	//write NaN string because format %8.3f writes nan not NaN
	for(int i = 0; i < O->GetNbMf(); i++)
	  fprintf(fa,"   NaN  ");
	fprintf(fa,"   NaN  ");//lk
	fprintf(fa,"   NaN  ");//rk
	fprintf(fa,"   NaN  ");//ls
	fprintf(fa,"   NaN  ");//rs
      }

  if(display)  fprintf(display, "Inferred output %f Alarm %d\n", val, Alarm);
  return val;
}

void DEFUZ_ImpFuzzy::WriteHeader(FILE *p, FISOUT *O) const
//***********************************************
{
	if( p != NULL)
	{
		fprintf(p, "     %s", INFERRED);
		fprintf(p, "     %s", ALARM);
		for(int i = 0; i < O->GetNbMf(); i++) fprintf(p, "      MF%d", i+1);
		fprintf(p, "     %s", MINKERNEL);
		fprintf(p, "     %s", MAXKERNEL);
		fprintf(p, "     %s", MINSUPPORT);
		fprintf(p, "     %s", MAXSUPPORT);
		fprintf(p, "     %s", INFERREDSYMBMATCH);
	}

}


//**************************    DEFUZ.CPP   *******************************

