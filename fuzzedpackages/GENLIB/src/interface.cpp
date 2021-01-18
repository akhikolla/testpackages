/*! \file interface.cc
\brief Implementation toutes les fonctions d'interface de Splus

	C'est dans ce fichier que sont placer toutes les interfaces pour Splus.
	Seulement les fonctions d'ici sont caller à partir des wrappers des .ssc

\author Sébastien Leclerc
\contributor Jean-François Lefebvre
*/

#include "base.h"
#include "userInterface.h"
#include "probsimul.h"
#include "apparentement.h"
#include "congen.h"
#include "consanguinite.h"
#include "fondateur.h"
#include "outils.h"
#include "genphi.h"
#include "statanal.h"
#include "interface.h"

#include <string>
#define ALLOWPRINTPROGRESS
#include <R.h>
//#include <Rdefines.h>
#include <Rmath.h> 
#include <Rcpp.h>
#include <Rcpp/as.h>
#include <RcppCommon.h>
//#define EXPORTTYPE extern "C"  -> remplacer par RcppExport

#define R_NO_REMAP
#define R_NO_REMAP_RMATH

/// Valeur textuel de l'enumeration typenoeud_t
/** \sa typenoeud_t*/
const char *stype[]=
{
	"NON-EXPLORER",
	"INUTILE",
	"NOEUD",
	"DEPART",
	"PROPOSANT",
	"GENPROPOSANTINUTILE"
};

//********************************************
//     FONCTION DE CONTROLE
//******************************************** 

RcppExport SEXP SPLUSFlushCacheGenealogie()
{
	//Flush la cache de la genealogie et autre
	FlushGenealogie();
	return R_NilValue;
}

RcppExport SEXP SPLUSGetTimer(SEXP sTimeInSec)
{
	//Flush la cache de la genealogie et autre
	int * TimeInSec;
	//std::vector<int> dat = Rcpp::as<std::vector<int>> (sTimeInSec);
	Rcpp::IntegerVector dat(sTimeInSec);
	TimeInSec = &dat[0];
	*TimeInSec = getLastTimer();
	return R_NilValue;
}

RcppExport SEXP SPLUSValidateGenealogie(SEXP RGenealogie, SEXP RisValid)
{
	STARTTIMER;
	int * Genealogie, * tmp;
	Rcpp::IntegerVector gen(RGenealogie);
//	std::vector<int> genT = Rcpp::as<std::vector<int>>(gen);
	Genealogie = INTEGER(gen); //&gen[0]; // genT
	
	//std::vector<int> isV = Rcpp::as<std::vector<int>>(RisValid);
	Rcpp::IntegerVector isV(RisValid);
	tmp = INTEGER(isV);//&isV[0];

	*tmp = ValidateGenealogie(Genealogie);
	//for(int i=0; i<gen.size(); i++) gen[i] = Genealogie[i];
	STOPTIMER;
     return (Rcpp::List::create( Rcpp::Named("Data")    = gen, //Rcpp::wrap(tmp1),
						   Rcpp::Named("isValid") = RisValid ) );
}

/// Fonction d'interface Splus pour change le temps maximum des fonctions longue
/** \sa setCurrentMaxTime() getCurrentMaxTime() */
RcppExport SEXP SPLUSChangeMaxProcessingTime(SEXP snewMaximum,SEXP soldMaximum)
{	
	double * newMaximum, * oldMaximum;
	//newMaximum = NUMERIC_POINTER(snewMaximum);
	//Rcpp::NumericVector v_snewMaximum (snewMaximum);
	newMaximum = REAL(snewMaximum); //&v_snewMaximum[0];

	//Rcpp::NumericVector v_soldMaximum(soldMaximum);
	//oldMaximum = REAL(soldMaximum) //&v_soldMaximum[0];
	oldMaximum = REAL(soldMaximum); //NUMERIC_POINTER(soldMaximum);
	getCurrentMaxTime(oldMaximum);
	if (*newMaximum>=0)
		setCurrentMaxTime(*newMaximum);
	return R_NilValue;
}

/// **********
//	APPARENTEMENT
// *********
#define USEALTERNATEEXCEPTION

/// Fonction d'interface Splus pour PhiMatrix
/** \sa PhiMatrix()*/
RcppExport SEXP SPLUSPhiMatrix(SEXP sGenealogie, SEXP sProposant, SEXP sNProposant, SEXP sNiveau, SEXP sPDRetour, SEXP sPrintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * Niveau, * printit ;
	double * pdRetour;

	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sProposant);
	Rcpp::NumericVector lpdRetour(sPDRetour);

	Genealogie = INTEGER(lGenealogie); //&lGenealogie[0];
	proposant = INTEGER(lproposant); //&lproposant[0];
	pdRetour	= REAL(lpdRetour);
	
	//int NProposant = Rcpp::as<int>(sNProposant);
	NProposant = INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	//int Niveau	= Rcpp::as<int>(sNiveau);
	Niveau	= INTEGER(sNiveau); //INTEGER_POINTER(sNiveau);
	//int printit	= Rcpp::as<int>(sPrintit);
	printit	= INTEGER(sPrintit); //INTEGER_POINTER(sPrintit);
	
	PhiMatrix(Genealogie, proposant,*NProposant,*Niveau, pdRetour,*printit);
//	PhiMatrix(Genealogie, proposant, NProposant, Niveau, pdR, printit);
	STOPTIMER;
	return R_NilValue;
}

/// Fonction d'interface Splus pour PhiMatrix
/** \sa PhiMatrix()*/
RcppExport SEXP SPLUSPhiMatrixMT(SEXP sGenealogie, SEXP sproposant, SEXP sNProposant, SEXP sNiveau, SEXP spdRetour, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * Niveau, * printit;
	double * pdRetour;

	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	Rcpp::NumericVector lpdRetour(spdRetour);
	
	Genealogie= INTEGER(lGenealogie); //&lGenealogie[0];
	proposant = INTEGER(lproposant); //&lproposant[0];
	pdRetour	= REAL(lpdRetour);
	
	//int NProposant = Rcpp::as<int>(sNProposant);
	NProposant = INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	//int Niveau	= Rcpp::as<int>(sNiveau);
	Niveau	= INTEGER(sNiveau); //INTEGER_POINTER(sNiveau);
	///int printit	= Rcpp::as<int>(sprintit);
	printit	= INTEGER(sprintit); //INTEGER_POINTER(sprintit);
	
	PhiMatrixMT(Genealogie, proposant, *NProposant, *Niveau, pdRetour,*printit);
	//NProposant, Niveau, pdR,printit);
	
	STOPTIMER;
	return R_NilValue;
}

/// Fonction d'interface Splus pour Phis2
/** \sa Phis2() Phis()*/
RcppExport SEXP SPLUSPhis(SEXP sGenealogie, SEXP sproposant, SEXP sNProposant,SEXP sNiveauMin, SEXP sNiveauMax, SEXP spdRetour, 
					 SEXP sMatrixArray, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * NiveauMin, * NiveauMax, * printit;
	double * pdRetour, * MatrixArray;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	Rcpp::NumericVector lpdRetour(spdRetour);
	Rcpp::NumericVector lMatrixArray(sMatrixArray);
	
	Genealogie = INTEGER(lGenealogie); //&lGenealogie[0];
	proposant  = INTEGER(lproposant); //&lproposant[0];
	pdRetour	 = REAL(spdRetour);
	MatrixArray= REAL(sMatrixArray);

	//int NProposant		= Rcpp::as<int>(sNProposant); 
	NProposant	= INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	//int NiveauMin		= Rcpp::as<int>(sNiveauMin);  
	NiveauMin		= INTEGER(sNiveauMin); //INTEGER_POINTER(sNiveauMin);
	//int NiveauMax		= Rcpp::as<int>(sNiveauMax);  
	NiveauMax		= INTEGER(sNiveauMax); //INTEGER_POINTER(sNiveauMax);
	//int printit		= Rcpp::as<int>(sprintit); 	
	printit		= INTEGER(sprintit); //INTEGER_POINTER(sprintit);
	
	//Peut utilise Phis en guise de comparaison  
	Phis(Genealogie, proposant, *NProposant,*NiveauMin,*NiveauMax,pdRetour,MatrixArray,*printit);
	//NProposant,NiveauMin,NiveauMax,pdR,MArray,printit);
	STOPTIMER;

	return R_NilValue;
}

RcppExport SEXP SPLUSPhisMT(	SEXP sGenealogie , SEXP sproposant, SEXP sNProposant, SEXP sNiveauMin, SEXP sNiveauMax, SEXP spdRetour, 
						SEXP sMatrixArray, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * NiveauMin, * NiveauMax, * printit;
	double * pdRetour, * MatrixArray;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	Rcpp::NumericVector lpdRetour(spdRetour);
	Rcpp::NumericVector lMatrixArray(sMatrixArray);
	
	Genealogie = INTEGER(lGenealogie); //&lGenealogie[0];
	proposant  = INTEGER(lproposant); //&lproposant[0];
	pdRetour	 = REAL(spdRetour);
	MatrixArray= REAL(sMatrixArray);

	//int NProposant	= Rcpp::as<int>(sNProposant);
	NProposant  = INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	//int NiveauMin	= Rcpp::as<int>(sNiveauMin); 
	NiveauMin	  = INTEGER(sNiveauMin); //INTEGER_POINTER(sNiveauMin);
	//int NiveauMax	= Rcpp::as<int>(sNiveauMax);
	NiveauMax	  = INTEGER(sNiveauMax); //INTEGER_POINTER(sNiveauMax);
	//int printit	= Rcpp::as<int>(sprintit);
	printit	  = INTEGER(sprintit); //INTEGER_POINTER(sprintit);
	

	//Peut utilise Phis en guise de comparaison  
	PhisMT(Genealogie, proposant,*NProposant,*NiveauMin,*NiveauMax, pdRetour,MatrixArray,*printit);
	//NProposant,NiveauMin,NiveauMax, pdR,MArray,printit);
	STOPTIMER;

	return R_NilValue;
}


/// **********
//	CONSANGUINITE
// *********

/// Fonction d'interface Splus pour consan (F)
/** \sa consan()*/
RcppExport SEXP SPLUSF(SEXP sGenealogie, SEXP sproposant, SEXP sNProposant, SEXP sNiveau, SEXP spdRetour, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * printit;
	double * pdRetour, * Niveau;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	Rcpp::NumericVector lpdRetour(spdRetour);
	
	Genealogie = INTEGER(lGenealogie); //&lGenealogie[0];
	proposant  = INTEGER(lproposant); //&lproposant[0];
	pdRetour   = REAL(lpdRetour); //&lpdRetour[0];

	//int NProposant	= Rcpp::as<int>(sNProposant);
	NProposant = INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	//int Niveau	= Rcpp::as<int>(sNiveau);
	Niveau	 = REAL(sNiveau); //NUMERIC_POINTER(sNiveau);
	//int printit	= Rcpp::as<int>(sprintit);
	printit	 = INTEGER(sprintit); //INTEGER_POINTER(sprintit);
	
	consan(Genealogie, proposant,*NProposant,*Niveau, pdRetour,*printit);
	//,NProposant,Niveau, pdRetour,printit);
	STOPTIMER;
	return R_NilValue;
}

/// Fonction d'interface Splus pour consanFs 
/** \sa consanFs()*/
RcppExport SEXP SPLUSFS(SEXP sGenealogie, SEXP sproposant, SEXP sNProposant, SEXP sNiveauMin, SEXP sNiveauMax, SEXP spdRetour, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie, * proposant, * NProposant, * printit;
	double * pdRetour, * NiveauMin, * NiveauMax;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sproposant);
	Rcpp::NumericVector lpdRetour(spdRetour);
	
	Genealogie = INTEGER(lGenealogie); //&lGenealogie[0];
	proposant = INTEGER(lproposant); //&lproposant[0];
	pdRetour = REAL(lpdRetour); //&lpdRetour[0];

	//int NProposant	= Rcpp::as<int>(sNProposant);
	NProposant  = INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	//int NiveauMin	= Rcpp::as<int>(sNiveauMin);
	NiveauMin	  = REAL(sNiveauMin); //NUMERIC_POINTER(sNiveauMin);
	//int NiveauMax	= Rcpp::as<int>(sNiveauMax);
	NiveauMax	  = REAL(sNiveauMax); //NUMERIC_POINTER(sNiveauMax);
	//int printit	= Rcpp::as<int>(sprintit);
	printit	  = INTEGER(sprintit); //INTEGER_POINTER(sprintit);
	
	consanFs(Genealogie, proposant,*NProposant,*NiveauMin,*NiveauMax, pdRetour,*printit);
	//NProposant,NiveauMin,NiveauMax, pdRetour,printit);
	STOPTIMER;
	return R_NilValue;
}

/// **********
//	OUTILS
// *********

/// Fonction d'interface Splus pour CountChild
/** \sa CountChild()*/
RcppExport SEXP SPLUSChild(SEXP sGenealogie, SEXP splProposant,SEXP slNProposant, SEXP sretour)
{
	STARTTIMER;
	int * Genealogie, * plProposant, * retour, * lNProposant;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lplProposant(splProposant);
	Rcpp::IntegerVector lretour(sretour);
	
	Genealogie  = INTEGER(lGenealogie); //&lGenealogie[0];
	plProposant = INTEGER(lplProposant); //&lplProposant[0];
	retour	  = INTEGER(lretour);

	//int lNProposant= Rcpp::as<int>(slNProposant);
	lNProposant = INTEGER(slNProposant); //INTEGER_POINTER(slNProposant);
	
	CountChild(Genealogie, plProposant, *lNProposant, retour);
	//, lNProposant, retour);
	STOPTIMER;
	return R_NilValue;
} 


/// Fonction d'interface Splus pour ebranche
/** \sa ebranche()*/
RcppExport  SEXP  SPLUSebranche(SEXP sGenealogie, SEXP sProposant, SEXP sNProposant, SEXP sAncetre, SEXP sNAncetre, 
						  SEXP sRetour, SEXP sTaille)
{
	STARTTIMER;
	int * Genealogie, * proposant, * ancetre, * retour, * taille, * nproposant, * nancetre ;

	Rcpp::IntegerVector lGenealogie(sGenealogie);
	Rcpp::IntegerVector lproposant(sProposant);
	Rcpp::IntegerVector lancetre(sAncetre);
	Rcpp::IntegerVector lretour(sRetour);
	
	Genealogie = INTEGER(lGenealogie);
	proposant  = INTEGER(lproposant);
	ancetre    = INTEGER(lancetre);
	retour     = INTEGER(lretour);

	nproposant= INTEGER(sNProposant); // INTEGER_POINTER(sNProposant);
	nancetre	= INTEGER(sNAncetre); // INTEGER_POINTER(sNAncetre);
	taille    = INTEGER(sTaille);// INTEGER_POINTER(sTaille);

	ebranche(Genealogie, proposant, *nproposant, ancetre, *nancetre, retour, taille);
	//nproposant, ancetre, nancetre, retour, taille);
	STOPTIMER;
	return R_NilValue;
}


/// Fonction d'interface Splus pour numeroGen
/** \sa compareGen()*/
RcppExport  SEXP  SPLUSnumeroGen(SEXP sGenealogie, SEXP splProposant, SEXP sNProposant, SEXP sretour)		   
{

	STARTTIMER;
	int * Genealogie, * plProposant, * retour, * NProposant;

	Rcpp::IntegerVector lGenealogie (sGenealogie);
	Rcpp::IntegerVector lplProposant(splProposant);
	Rcpp::IntegerVector lretour     (sretour);

	Genealogie  = INTEGER(lGenealogie); //&lGenealogie[0];
	plProposant = INTEGER(lplProposant); //&lplProposant[0];
	retour 	  = INTEGER(lretour); //&lretour[0];

	//int NProposant	= Rcpp::as<int>(sNProposant);
	NProposant  = INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	
	numeroGen(Genealogie, plProposant,*NProposant, retour);
	//NProposant, retour);
	STOPTIMER;
	return R_NilValue;
}
/// Fonction d'interface Splus pour numeroGenMin
/** \sa compareGen()*/
RcppExport  SEXP  SPLUSnumGenMin(SEXP sGenealogie, SEXP splProposant,SEXP sNProposant, SEXP sretour)		   
{

	STARTTIMER;
	int * Genealogie, * plProposant, * retour, * NProposant;

	Rcpp::IntegerVector lGenealogie (sGenealogie);
	Rcpp::IntegerVector lplProposant(splProposant);
	Rcpp::IntegerVector lretour     (sretour);

	Genealogie  = INTEGER(lGenealogie); //&lGenealogie[0];
	plProposant = INTEGER(lplProposant); //&lplProposant[0];
	retour 	  = INTEGER(lretour); //&lretour[0];

	//int NProposant	= Rcpp::as<int>(sNProposant);
	NProposant  = INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	
	numeroGenMin(Genealogie, plProposant,*NProposant, retour);
	//NProposant, retour);
	STOPTIMER;

	return R_NilValue;
}

/// Fonction d'interface Splus pour numeroGenMoy
/** \sa compareGen()*/
RcppExport  SEXP  SPLUSnumGenMoy(SEXP sGenealogie, SEXP splProposant,SEXP sNProposant, SEXP sretour)		   
{

	STARTTIMER;
	int * Genealogie, * plProposant, * NProposant;
	double * retour;

	Rcpp::IntegerVector lGenealogie (sGenealogie);
	Rcpp::IntegerVector lplProposant(splProposant);
	Rcpp::NumericVector lretour     (sretour);

	Genealogie  = INTEGER(lGenealogie); //&lGenealogie[0];
	plProposant = INTEGER(lplProposant); //&lplProposant[0];
	retour	  = REAL(sretour);

	//int NProposant	= Rcpp::as<int>(sNProposant);
	NProposant  = INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	
	numeroGenMoy(Genealogie, plProposant,*NProposant, retour);
	//NProposant, retour);
	STOPTIMER;

	return R_NilValue;
}
/// **********
//	CONTRIBUTION GENETIQUE
// ***********

/// Fonction d'interface Splus pour Congen
/** \sa Congen()*/
RcppExport SEXP SPLUSConGen(SEXP sGenealogie, SEXP slProposant, SEXP sNProposant, SEXP slAncetre, SEXP sNAncetre, SEXP sdRetour, 
					   SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie , * plProposant, * plAncetre, * NProposant, * NAncetre, * printit ;
	double * dRetour;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie); // conversion automatique avec Rcpp
	Rcpp::IntegerVector lProposant(slProposant);
	Rcpp::IntegerVector lAncetre(slAncetre);
	Rcpp::NumericVector ldRetour(sdRetour);
	
	Genealogie  = INTEGER(lGenealogie); //&lGenealogie[0];
	plProposant = INTEGER(lProposant); //&lProposant[0];
	plAncetre   = INTEGER(lAncetre); //&lAncetre[0];
	dRetour	  = REAL(ldRetour);
	
	//int NProposant	= Rcpp::as<int>(sNProposant); //	
	NProposant	= INTEGER(sNProposant); //INTEGER_POINTER(sNProposant);
	//int NAncetre	= Rcpp::as<int>(sNAncetre); //	
	NAncetre		= INTEGER(sNAncetre); //INTEGER_POINTER(sNAncetre);
	//int printit	= Rcpp::as<int>(sprintit); //	
	printit		= INTEGER(sprintit); //INTEGER_POINTER(sprintit);
	
	Congen(Genealogie, plProposant , *NProposant, plAncetre, *NAncetre, dRetour, *printit);
	//NProposant, plAncetre, NAncetre, dRetour, printit);
	STOPTIMER;
	return R_NilValue;
}

/// Fonction d'interface Splus pour Congen
/** \sa Congen()*/
RcppExport SEXP SPLUSConGenPLUS(SEXP sGenealogie, SEXP splProposant,SEXP slNProposant, SEXP splAncetre, SEXP slNAncetre, SEXP spdSexe, 
						  SEXP spdRetour, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie , * plProposant, * plAncetre, * lNProposant, * lNAncetre, * printit ;
	double * pdRetour, * pdSexe;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie); // conversion automatique avec Rcpp
	Rcpp::IntegerVector lProposant (splProposant);
	Rcpp::IntegerVector lAncetre   (splAncetre);
	Rcpp::NumericVector lpdSexe	 (spdSexe);
	Rcpp::NumericVector lpdRetour  (spdRetour);
	
	Genealogie	= INTEGER(lGenealogie); //&lGenealogie[0]; // INTEGER_POINTER(
	plProposant	= INTEGER(lProposant); //&lProposant[0]; // INTEGER_POINTER(
	plAncetre		= INTEGER(lAncetre); //&lAncetre[0]; // INTEGER_POINTER(
	pdRetour	 	= REAL(lpdRetour);
	pdSexe		= REAL(lpdSexe);
	
	//int lNProposant= Rcpp::as<int>(slNProposant);
	lNProposant	= INTEGER(slNProposant); //INTEGER_POINTER(slNProposant);
	//int lNAncetre	= Rcpp::as<int>(slNAncetre);
	lNAncetre		= INTEGER(slNAncetre); //INTEGER_POINTER(slNAncetre);
	//int printit	= Rcpp::as<int>(sprintit);
	printit		= INTEGER(sprintit); //INTEGER_POINTER(sprintit);
	
	CongenPLUS(Genealogie, plProposant, *lNProposant, plAncetre, *lNAncetre, pdSexe, pdRetour, *printit);
	// lNProposant, plAncetre, lNAncetre, pdSexe, pdRetour, printit);
	STOPTIMER;
	return R_NilValue;
}

/// Fonction d'interface Splus pour CongenCumul
/** \sa CongenCumul()*/
RcppExport SEXP SPLUSCGCumul(	SEXP sGenealogie, SEXP splProposant,SEXP slNProposant, SEXP splAncetre, SEXP slNAncetre, SEXP sAncRet,
						SEXP spdRetour, SEXP spdRetourCumul, SEXP sprintit)
{
	STARTTIMER;
	int * Genealogie , * plProposant, * plAncetre, * AncRet, * lNProposant, * lNAncetre, * printit ;
	double * pdRetour, * pdRetourCumul;
	
	Rcpp::IntegerVector lGenealogie(sGenealogie); // conversion automatique avec Rcpp
	Rcpp::IntegerVector lProposant(splProposant);
	Rcpp::IntegerVector lAncetre(splAncetre);
	Rcpp::NumericVector lpdRetour(spdRetour);
	Rcpp::NumericVector lpdRetourCumul(spdRetourCumul);
	
	Genealogie  = INTEGER(lGenealogie); //&lGenealogie[0];
	plProposant = INTEGER(lProposant); //&lProposant[0];
	plAncetre	  = INTEGER(lAncetre); //&lAncetre[0];
	pdRetour      = REAL(lpdRetour);
	pdRetourCumul = REAL(lpdRetourCumul);

	AncRet		= NULL;
	//int lNProposant= Rcpp::as<int>(slNProposant);
	lNProposant= INTEGER(slNProposant); //INTEGER_POINTER(slNProposant);
	//int lNAncetre	= Rcpp::as<int>(slNAncetre);
	lNAncetre	 = INTEGER(slNAncetre); //INTEGER_POINTER(slNAncetre);
	//int printit	= Rcpp::as<int>(sprintit);
	printit	 = INTEGER(sprintit); //INTEGER_POINTER(sprintit);

	CongenCumul(Genealogie, plProposant, *lNProposant, plAncetre, *lNAncetre, AncRet, pdRetour, pdRetourCumul, *printit);
	//lNProposant, plAncetre, lNAncetre, AncRet, pdRetour, pdRetourCumul, printit);
	STOPTIMER
	return R_NilValue;
}

/// Fonction d'interface Splus pour CongenCumul
/** \sa CongenCumuldirect()*/
RcppExport SEXP SPLUSCGCumuldirect(SEXP smatriceCG, SEXP slNProposant, SEXP splAncetre, SEXP slNAncetre, SEXP sAncRet, 
							SEXP spdSomAnc, SEXP spdSomCumul)
{
	STARTTIMER;
	int * matriceCG, * plAncetre, * AncRet, * lNProposant, * lNAncetre;
	double * pdSomAnc, * pdSomCumul;
	Rcpp::IntegerVector lmatriceCG ( smatriceCG );
	Rcpp::IntegerVector lplAncetre ( splAncetre );
	Rcpp::IntegerVector lAncRet    ( sAncRet );
	
	matriceCG		= INTEGER(lmatriceCG); //&lmatriceCG[0]; // INTEGER_POINTER( lmatriceCG );
	plAncetre		= INTEGER(lplAncetre); //&lplAncetre[0]; // INTEGER_POINTER( lplAncetre );
	AncRet		= INTEGER(lAncRet); //&lAncRet[0]; // INTEGER_POINTER( lAncRet );
	
	//int lNProposant= Rcpp::as<int>(slNProposant);
	lNProposant	= INTEGER(slNProposant); //INTEGER_POINTER( slNProposant );
	//int lNAncetre= Rcpp::as<int>(slNAncetre);
	lNAncetre		= INTEGER(slNAncetre); //INTEGER_POINTER( slNAncetre );
	
	pdSomAnc   = REAL(spdSomAnc);
	pdSomCumul = REAL(spdSomCumul);
	
	CongenCumuldirect(matriceCG, *lNProposant, plAncetre, *lNAncetre, AncRet,pdSomAnc,pdSomCumul);
	//lNProposant, plAncetre, lNAncetre, AncRet, pdSomAnc, pdSomCumul);
	STOPTIMER
	return R_NilValue;
}

/// **********
//	DIVERS
// *********

/*FONCTION D'INTERFACE POUR SPLUS*/

/// Fonction d'interface Splus pour simul
/** \sa simul()*/
//RcppExport SEXP SPLUSSimul(SEXP sGenealogie, SEXP sproposant, SEXP setatproposant, SEXP snproposant, SEXP sancetre, SEXP setatancetre,
RcppExport SEXP SPLUSSimul(SEXP sGenealogie, SEXP sproposant, SEXP setatproposant, SEXP snproposant, SEXP sancetre, SEXP setatancetre,
					  SEXP snancetre, SEXP snSimul, SEXP spdRetConj, SEXP spdRetSimul, SEXP spdRetProp, SEXP sprobRecomb,
					  SEXP sprobSurvieHomo, SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie, * proposant, * etatproposant, * ancetre, * etatancetre, * nproposant, * nancetre, * nSimul, * PrintProgress;
	double * pdRetConj, * pdRetSimul, * pdRetProp, * probRecomb;
	
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector letatproposant	( setatproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::IntegerVector letatancetre	( setatancetre );
	Rcpp::NumericVector lpdRetSimul	( spdRetSimul );
	Rcpp::NumericVector lpdRetProp	( spdRetProp );
	Rcpp::NumericVector lprobRecomb	( sprobRecomb );
	
	Genealogie	= INTEGER(lGenealogie); //&lGenealogie[0]; // INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER(lproposant); //&lproposant[0]; // INTEGER_POINTER( lproposant );
	etatproposant	= INTEGER(letatproposant); //&letatproposant[0]; // INTEGER_POINTER( letatproposant );
	ancetre		= INTEGER(lancetre); //&lancetre[0]; // INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER(letatancetre); //&letatancetre[0]; // INTEGER_POINTER( letatancetre );
	pdRetSimul = REAL(lpdRetSimul); //&lpdRetSimul[0]; // NUMERIC_POINTER( lpdRetSimul );
	pdRetProp	 = REAL(lpdRetProp); //&lpdRetProp[0]; // NUMERIC_POINTER( lpdRetProp );
	probRecomb = REAL(lprobRecomb); //&lprobRecomb[0]; // NUMERIC_POINTER( lprobRecomb );

	//int nproposant	= Rcpp::as<int>( snproposant ); 
	nproposant	= INTEGER(snproposant);// INTEGER_POINTER( snproposant );
	//int nancetre		= Rcpp::as<int>( snancetre );
	nancetre		= INTEGER(snancetre);// INTEGER_POINTER( snancetre );
	//int nSimul		= Rcpp::as<int>( snSimul ); 
	nSimul		= INTEGER(snSimul);// INTEGER_POINTER( snSimul );
	//int PrintProgress	= Rcpp::as<int>( sPrintProgress ); 
	PrintProgress	= INTEGER(sPrintProgress);// INTEGER_POINTER( sPrintProgress );

	pdRetConj = REAL(spdRetConj);
	
	double probSurvieHomo  = Rcpp::as<double>(sprobSurvieHomo);
	
	simul(Genealogie, proposant, etatproposant, *nproposant, ancetre, etatancetre, *nancetre, *nSimul, pdRetConj, pdRetSimul,
		 pdRetProp,  probRecomb, probSurvieHomo, *PrintProgress);
	//nproposant, ancetre, etatancetre, nancetre, nSimul, pdRetConj, pdRetSimul, pdRetProp,  probRecomb, probSurvieHomo, PrintProgress);
	STOPTIMER;
	//return;
	return Rcpp::wrap(getLastTimer());
}

/// Fonction d'interface Splus pour simulsingle
/** \sa simulsingle()*/
RcppExport SEXP SPLUSSimulSingle(SEXP sGenealogie, SEXP sproposant, SEXP snproposant, SEXP sancetre, SEXP setatancetre, SEXP snancetre, 
						   SEXP sNSimul, SEXP spdRetour, SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie, * proposant, * ancetre, * etatancetre, * nproposant, * nancetre, * NSimul, * PrintProgress;
	double * pdRetour;
	Rcpp::IntegerVector lGenealogie ( sGenealogie );
	Rcpp::IntegerVector lproposant  ( sproposant );
	Rcpp::IntegerVector lancetre	  ( sancetre );
	Rcpp::IntegerVector letatancetre( setatancetre );
	
	Genealogie	= INTEGER(lGenealogie); //&lGenealogie[0]; // INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER(lproposant); //&lproposant[0];  // INTEGER_POINTER( lproposant );
	ancetre		= INTEGER(lancetre); //&lancetre[0];    // INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER(letatancetre); //&letatancetre[0]; // INTEGER_POINTER( letatancetre );

	//int nproposant	= Rcpp::as<int>( snproposant );
	nproposant	= INTEGER(snproposant); //INTEGER_POINTER( snproposant );
	//int nancetre		= Rcpp::as<int>( snancetre );
	nancetre		= INTEGER(snancetre); //INTEGER_POINTER( snancetre );
	//int NSimul		= Rcpp::as<int>( sNSimul );
	NSimul		= INTEGER(sNSimul); //INTEGER_POINTER( sNSimul );
	//int PrintProgress	= Rcpp::as<int>( sPrintProgress );
	PrintProgress	= INTEGER(sPrintProgress); //INTEGER_POINTER( sPrintProgress );
	
	pdRetour		= REAL(spdRetour);
	
	simulsingle(Genealogie, proposant, *nproposant, ancetre, etatancetre, *nancetre, *NSimul, pdRetour, *PrintProgress);
	//nproposant, ancetre, etatancetre, nancetre, NSimul, pdRetour, PrintProgress);
	STOPTIMER;
	return R_NilValue;
}

/// Fonction d'interface Splus pour simulsingleFreq
/** \sa simulsingleFreq()*/
RcppExport SEXP SPLUSSimulSingleFreq(SEXP sGenealogie, SEXP sproposant, SEXP snproposant, SEXP sancetre, SEXP setatancetre, SEXP snancetre,
							  SEXP sNSimul, SEXP spdRetour,SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie, * proposant, * ancetre, * etatancetre, * nproposant , * nancetre, * NSimul, * PrintProgress;
	double * pdRetour;
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::IntegerVector letatancetre	( setatancetre );
	Rcpp::NumericVector lpdRetour	( spdRetour );
	
	Genealogie	= INTEGER(lGenealogie); //&lGenealogie[0]; // INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER(lproposant); //&lproposant[0]; // INTEGER_POINTER( lproposant );
	ancetre		= INTEGER(lancetre); //&lancetre[0]; // INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER(letatancetre); //&letatancetre[0]; // INTEGER_POINTER( letatancetre );
	pdRetour		= REAL(lpdRetour);

	//int nproposant	= Rcpp::as<int>( snproposant );
	nproposant	= INTEGER(snproposant); //INTEGER_POINTER( snproposant );
	//int nancetre		= Rcpp::as<int>( snancetre );
	nancetre		= INTEGER(snancetre); //INTEGER_POINTER( snancetre );
	//int NSimul		= Rcpp::as<int>( sNSimul );
	NSimul		= INTEGER(sNSimul); //INTEGER_POINTER( sNSimul );
	//int PrintProgress	= Rcpp::as<int>( sPrintProgress );
	PrintProgress	= INTEGER(sPrintProgress); //INTEGER_POINTER( sPrintProgress );

	
	simulsingleFreq(Genealogie, proposant, *nproposant, ancetre, etatancetre, *nancetre, *NSimul, pdRetour, *PrintProgress);
	//nproposant, ancetre, etatancetre, nancetre, NSimul, pdRetour, PrintProgress);
	STOPTIMER;
	return R_NilValue;
}

/// Fonction d'interface Splus pour simulsingleProb
/** \sa simulsingleProb()*/
RcppExport SEXP SPLUSSimulSingleProb( SEXP sGenealogie,SEXP sproposant, SEXP snproposant, SEXP sancetre,SEXP snancetre,
							   SEXP setatancetre,SEXP smtProb, SEXP sNSimul, SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie, * proposant, * ancetre, * etatancetre, * nproposant, * nancetre, * NSimul, * PrintProgress; //* mtProb, 
	
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::IntegerVector letatancetre	( setatancetre );
	
	Genealogie	= INTEGER(lGenealogie); //&lGenealogie[0]; // INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER(lproposant); //&lproposant[0]; // INTEGER_POINTER( lproposant );
	ancetre		= INTEGER(lancetre); //&lancetre[0]; // INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER(letatancetre); //&letatancetre[0]; // INTEGER_POINTER( letatancetre );

	//int nproposant	= Rcpp::as<int>( snproposant );
	nproposant	= INTEGER(snproposant); //INTEGER_POINTER( snproposant );
	//int nancetre		= Rcpp::as<int>( snancetre );
	nancetre		= INTEGER(snancetre); //INTEGER_POINTER( snancetre );
	//int NSimul		= Rcpp::as<int>( sNSimul );
	NSimul		= INTEGER(sNSimul); //INTEGER_POINTER( sNSimul );
	//int PrintProgress	= Rcpp::as<int>( sPrintProgress );
	PrintProgress	= INTEGER(sPrintProgress); //INTEGER_POINTER( sPrintProgress );
	
	SEXP ret = simulsingleProb(Genealogie, proposant, *nproposant, ancetre, *nancetre, etatancetre, smtProb, *NSimul, *PrintProgress );
	//nproposant, ancetre, nancetre, etatancetre, smtProb, NSimul, PrintProgress );
	STOPTIMER; 
	return ret;
}

/// Fonction d'interface Splus pour simulsingleFct
/** \sa simulsingleFct()*/
RcppExport SEXP SPLUSSimulSingleFct(SEXP sGenealogie, SEXP sproposant, SEXP sancetre, SEXP sancEtatAll1, SEXP sancEtatAll2, 
							 SEXP snancetre, SEXP sNSimul, SEXP sfctSousGrp, SEXP sPrintProgress)
{
	STARTTIMER;
	int * Genealogie,  * proposant, * ancetre, * ancEtatAll1, * ancEtatAll2, * nancetre, * NSimul, * PrintProgress;

	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	
	Genealogie	= INTEGER(lGenealogie); //&lGenealogie[0]; // INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER(lproposant); //&lproposant[0]; // INTEGER_POINTER( lproposant );
	ancetre		= INTEGER(lancetre); //&lancetre[0]; // INTEGER_POINTER( lancetre );

	//int nancetre		= Rcpp::as<int>( snancetre );
	nancetre		= INTEGER(snancetre); //INTEGER_POINTER( snancetre );
	//int NSimul		= Rcpp::as<int>( sNSimul );
	NSimul		= INTEGER(sNSimul); //INTEGER_POINTER( sNSimul );
	//int PrintProgress	= Rcpp::as<int>( sPrintProgress );
	PrintProgress	= INTEGER(sPrintProgress); //INTEGER_POINTER( sPrintProgress );

	ancEtatAll1	= INTEGER(sancEtatAll1);// INTEGER_POINTER( sancEtatAll1 );
	ancEtatAll2	= INTEGER(sancEtatAll2);// INTEGER_POINTER( sancEtatAll2 );
	int nproposant = lproposant.size();
	
	SEXP ret = simulsingleFct(Genealogie, proposant, nproposant, ancetre, ancEtatAll1, ancEtatAll2, *nancetre, *NSimul,
						 sfctSousGrp, *PrintProgress);
	//ancEtatAll1, ancEtatAll2, nancetre, NSimul, sfctSousGrp, PrintProgress);
	STOPTIMER;
	return ret;
}

/// Fonction d'interface Splus pour prob
/** \sa prob()*/
RcppExport SEXP SPLUSProb(SEXP sGenealogie, SEXP sproposant, SEXP setatproposant,SEXP snproposant, SEXP sancetre, SEXP setatancetre, SEXP snancetre,
					 SEXP spdRetConj, SEXP spdRetSimul,SEXP sPrintProgress,SEXP sonlyConj)
{
	STARTTIMER;
	int * Genealogie, * proposant, * etatproposant, * ancetre, * etatancetre, * nproposant, * nancetre, * PrintProgress, * onlyConj;
	double * pdRetConj, * pdRetSimul;
	
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector letatproposant	( setatproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::IntegerVector letatancetre	( setatancetre );
	Rcpp::NumericVector lpdRetSimul	( spdRetSimul );
	
	Genealogie	= INTEGER(lGenealogie); //&lGenealogie[0]; // INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER(lproposant); //&lproposant[0]; // INTEGER_POINTER( lproposant );
	etatproposant	= INTEGER(letatproposant); //&letatproposant[0]; // INTEGER_POINTER( letatproposant );
	ancetre		= INTEGER(lancetre); //&lancetre[0]; // INTEGER_POINTER( lancetre );
	etatancetre	= INTEGER(letatancetre); //&letatancetre[0]; // INTEGER_POINTER( letatancetre );
	pdRetSimul	= REAL(lpdRetSimul); //&lpdRetSimul[0]; // NUMERIC_POINTER( lpdRetSimul );
	
	//int nproposant		= Rcpp::as<int>( snproposant );
	nproposant	= INTEGER(snproposant); //INTEGER_POINTER( snproposant );
	//int nancetre		= Rcpp::as<int>( snancetre );
	nancetre		= INTEGER(snancetre); //INTEGER_POINTER( snancetre );
	//int PrintProgress	= Rcpp::as<int>( sPrintProgress );
	PrintProgress	= INTEGER(sPrintProgress); //INTEGER_POINTER( sPrintProgress );
	//int onlyConj		= Rcpp::as<int>( sonlyConj );
	onlyConj		= INTEGER(sonlyConj); //INTEGER_POINTER( sonlyConj );

	pdRetConj = REAL(spdRetConj);
	
	SEXP ret = prob(Genealogie, proposant, etatproposant, *nproposant, ancetre, etatancetre, *nancetre, pdRetConj, 
				 pdRetSimul,*PrintProgress, *onlyConj);
	//nproposant, ancetre, etatancetre, nancetre, pdRetConj, pdRetSimul, PrintProgress, onlyConj);
	STOPTIMER;	
	return ret;
}


/// Fonction d'interface Splus pour CoefApparentement
/** \sa CoefApparentement()*/

RcppExport SEXP SPLUSCoeffApparentement(SEXP sGenealogie, SEXP sproposant, SEXP snproposant, SEXP sancetre, SEXP sretour, 
								SEXP sDuppDetection, SEXP sprintprogress)
{
	STARTTIMER
	int * Genealogie, * proposant, * ancetre, * nproposant, * DuppDetection, * printprogress;
	double * retour;
	
	Rcpp::IntegerVector lGenealogie	( sGenealogie );
	Rcpp::IntegerVector lproposant	( sproposant );
	Rcpp::IntegerVector lancetre		( sancetre );
	Rcpp::NumericVector lretour		( sretour );
	
	Genealogie	= INTEGER(lGenealogie); //&lGenealogie[0]; // INTEGER_POINTER( lGenealogie );
	proposant		= INTEGER(lproposant); //&lproposant[0]; // INTEGER_POINTER( lproposant );
	ancetre		= INTEGER(lancetre); //&lancetre[0]; // INTEGER_POINTER( lancetre );
	retour		= REAL(lretour); //&lretour[0]; // NUMERIC_POINTER( lretour );
	
	//int nproposant	= Rcpp::as<int>( snproposant );
	nproposant	= INTEGER(snproposant); //INTEGER_POINTER( snproposant );
	//int DuppDetection	= Rcpp::as<int>( sDuppDetection );
	DuppDetection	= INTEGER(sDuppDetection); //INTEGER_POINTER( sDuppDetection );
	//int printprogress	= Rcpp::as<int>( sprintprogress );
	printprogress	= INTEGER(sprintprogress); //INTEGER_POINTER( sprintprogress );
	
	CoefApparentement(Genealogie, proposant, *nproposant, ancetre, retour, *DuppDetection, *printprogress);
	//nproposant, ancetre, retour, DuppDetection, printprogress);
	STOPTIMER;
	return R_NilValue; 
} 



//********************************************
//     FONCTION D'INTERFACE SPLUS CALL
//********************************************

//#ifndef MODETEST


/*! 
	\brief SPLUSCALL: Creer une genealogie binaire à l'aide d'une genealogie de forme classique No Ind, No Pere, No Mere

     Cette fonction est une interface SPLUS de la fonction CreerGenealogie()
	
	\param  SIndividu	[in] SEXP de type vecteur int contenant no individu

	\param  SPere		[in] SEXP de type vecteur int contenant no pere
								
	\param  SMere		[in] SEXP de type vecteur int contenant no mere
	
	\return Un pointeur vers un SEXP qui contient la genealogie en entier sous forme d'un long vecteur de int
	
	\remarque La genealogie composee de ind,pere,mere doit etre valide pour que la fonction s'execute correctement

	\sa CreerGenealogie()
*/
RcppExport SEXP SPLUSCALLCreerObjetGenealogie(SEXP SIndividu,SEXP SPere, SEXP SMere, SEXP SSexe)
{
	STARTTIMER
	//VARIABLE OPERATIONNEL		   conversion automatique avec Rcpp
	int * plIndividu, * plPere, * plMere, * plSexe;

	Rcpp::IntegerVector lIndividu(SIndividu);
	Rcpp::IntegerVector lPere(SPere);
	Rcpp::IntegerVector lMere(SMere);
	Rcpp::IntegerVector lSexe(SSexe);
//	try{ lSexe = Rcpp::as<Rcpp::IntegerVector>(SSexe); } catch(std::exception &ex) { }
	
	plIndividu = INTEGER(lIndividu); //&lIndividu[0]; // INTEGER_POINTER(lIndividu);
	plPere     = INTEGER(lPere); //&lPere[0]; // INTEGER_POINTER(lPere);
	plMere     = INTEGER(lMere); //&lMere[0]; // INTEGER_POINTER(lMere);
	plSexe     = INTEGER(lSexe); //&lSexe[0]; // INTEGER_POINTER(lSexe);
	int lNIndividu = lIndividu.size();
	
	//Trois vecteur de meme taille
	if (lNIndividu != lPere.size() || lNIndividu != lMere.size())
	{
		ErrorHandler();
		//PROBLEM "%s","LES TROIS VECTEURS (ind,pere & mere) DOIVENT ETRE DE MEME TAILLE";
		//RECOVER(NULL_ENTRY);
	}
	
	//Verifie si le sexe est une information utile
	if (lSexe.size() != lNIndividu)	plSexe=NULL;
	
	int NombreEnfant=0; 
	int i;  

	//COMPLETE LA GENEALOGIE AVEC LES ELEMENTS MANQUANTS
	INITGESTIONMEMOIRE;	
	int* fInd	=(int*)  memalloc(lNIndividu*3,sizeof(int));	
	int* fpere	=(int*)  memalloc(lNIndividu*3,sizeof(int));	
	int* fmere	=(int*)  memalloc(lNIndividu*3,sizeof(int));
	int* fsexe 	= NULL;
	if (plSexe)
		fsexe	=(int*)  memalloc(lNIndividu*3,sizeof(int));
	CompleteGenealogie(plIndividu,plPere,plMere,plSexe,fInd,fpere,fmere,fsexe,&lNIndividu);

	//COMPTER LE NOMBRE D'ENFANT
	for(i=0;i<lNIndividu;i++) 	
	{
		if (fpere[i]!=0)	++NombreEnfant;
		if (fmere[i]!=0)	++NombreEnfant;
	}  

	//CREER UN OBJET XPLUS QUI CONTIENT TOUTE LES DONNEES
	const int TAILLESAUVEGARDE=TAILLEGENVERSION7(lNIndividu,NombreEnfant);
//	int	* saveobj = new int[TAILLESAUVEGARDE];		  /*chgt JFL*/
//	int * saveptr = saveobj;		  					/*chgt JFL*/
	int * saveptr = new int[TAILLESAUVEGARDE];		  /*chgt JFL*/
//return (Rcpp::wrap( TAILLESAUVEGARDE ));

	//Creation de la genealogie
	CreerGenealogie(fInd,fpere,fmere,fsexe,lNIndividu,saveptr);
	Rcpp::IntegerVector retour(TAILLESAUVEGARDE);
//	for(int i=0;i<TAILLESAUVEGARDE;i++) retour[i] = saveobj[i];		  /*chgt JFL*/
	for(int i=0;i<TAILLESAUVEGARDE;i++) retour[i] = saveptr[i];		  /*chgt JFL*/
	STOPTIMER;
	for(int i=0;i<lNIndividu;i++) 
	{
		plIndividu[i] = fInd[i];
		plPere[i] = fpere[i];
		plMere[i] = fmere[i];
	}
	if(saveptr) { delete [] saveptr; }		  /*chgt JFL*/
	return (Rcpp::wrap( retour )); //saveobj;
} 

//#endif


//********************************************
//     FONCTION EXTRACTION DE GENEALOGIE
//******************************************** 

/*! 
	\brief Extraction de la genealogie sous la forme: No individu, No Pere, No mere

	La fonction Extrait la genealogie et retourne 3 vecteur soit le No Individu, le No du pere et le No de la Mere.
	
	\param genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\retval plRetIndividu [out] Vecteur de taille LengthGenealogie
								En cas de succes contient les No d'individus
	
	\retval plRetPere	[out] Vecteur de taille LengthGenealogie
								En cas de succes contient les No des peres

  	\retval plRetMere	[out] Vecteur de taille LengthGenealogie
								En cas de succes contient les No des meres

	\remark Il faut utiliser ( .C LengthGenealogie() ou SPLUS. gen.length) pour connaitre le nombre d'individu de la genealogie

	\sa LengthGenealogie()
*/
RcppExport SEXP SPLUSOutgen(SEXP Rgenealogie, SEXP RplRetIndividu, SEXP RplRetPere, SEXP RplRetMere, SEXP RplRetSexe, SEXP Rmustsort)
{
	STARTTIMER;
	int * genealogie, * plRetIndividu, * plRetPere, * plRetMere, * plRetSexe, * mustsort;

	Rcpp::IntegerVector a(Rgenealogie);
	Rcpp::IntegerVector b(RplRetIndividu);
	Rcpp::IntegerVector c(RplRetPere);
	Rcpp::IntegerVector d(RplRetMere);
	Rcpp::IntegerVector e(RplRetSexe);
	
	genealogie    = INTEGER(a); //&a[0]; // INTEGER_POINTER(a); //par["Data"]);
	plRetIndividu = INTEGER(b); //&b[0]; // INTEGER_POINTER(b); //par["ind"]);
	plRetPere     = INTEGER(c); //&c[0]; // INTEGER_POINTER(c); //par["pere"]);
	plRetMere     = INTEGER(d); //&d[0]; // INTEGER_POINTER(d); //par["mere"]);
	plRetSexe     = INTEGER(e); //&e[0]; // INTEGER_POINTER(e); //par["sexe"]);

	mustsort = INTEGER(Rmustsort);
	
	//CREATION DU TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(genealogie, GFALSE, &lNIndividu, &Noeud); //Pas d'enfant  int test = 
	
	int havesex = (LoadNIndMasc()!=-1);
	//creation vecteur de sortie
	for(int i=0;i<lNIndividu;i++)
	{
		plRetIndividu[i]=Noeud[i].nom;
		if (Noeud[i].pere!=NULL)	plRetPere[i]=Noeud[i].pere->nom;		
		else						plRetPere[i]=0;
			
		if (Noeud[i].mere!=NULL)	plRetMere[i]=Noeud[i].mere->nom;			
		else						plRetMere[i]=0;
		//Recuperation du sexe
		if (havesex)				plRetSexe[i]=Noeud[i].sex;
		else						plRetSexe[i]=-1; //Pas utilisable
	}
	//Trie si nécessaire
	if (*mustsort)
		SortGenealogie3Vecteur(plRetIndividu,plRetPere,plRetMere,plRetSexe,lNIndividu);
	STOPTIMER;
	return (Rcpp::List::create(	Rcpp::Named("Data") = a,
							Rcpp::Named("ind")  = b,
							Rcpp::Named("father") = c,
							Rcpp::Named("mother") = d,
							Rcpp::Named("sex") = e
							));
}

/*! 
	\brief Extraction de la genealogie sous la forme: No individu, Indice Pere, Indice mere

	La fonction Extrait la genealogie et retourne 3 vecteur soit le No Individu, le Indice du pere et le Indice de la Mere.
	
	\param genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\retval plRetIndividu [out] Vecteur de taille LengthGenealogie
								En cas de succes contient les No d'individus
	
	\retval plRetPere	[out] Vecteur de taille LengthGenealogie
								En cas de succes contient les Indice des peres

  	\retval plRetMere	[out] Vecteur de taille LengthGenealogie
								En cas de succes contient les Indice des meres

	\remark Il faut utiliser ( .C LengthGenealogie() ou SPLUS. gen.length) pour connaitre le nombre d'individu de la genealogie

	\sa LengthGenealogie()
*/
RcppExport SEXP SPLUSOutIndice(SEXP sgenealogie, SEXP splRetIndividu, SEXP splRetPere, SEXP splRetMere, SEXP splRetSexe, SEXP smustsort)
{

	STARTTIMER;			
	int * genealogie, * plRetIndividu, * plRetPere, * plRetMere, * plRetSexe, * mustsort;
	
	Rcpp::IntegerVector lgenealogie	( sgenealogie );
	Rcpp::IntegerVector lplRetIndividu	( splRetIndividu );
	Rcpp::IntegerVector lplRetPere	( splRetPere );
	Rcpp::IntegerVector lplRetMere	( splRetMere );
	Rcpp::IntegerVector lplRetSexe	( splRetSexe );
	
	genealogie	= INTEGER(lgenealogie); //&lgenealogie[0]; // INTEGER_POINTER( lgenealogie );
	plRetIndividu	= INTEGER(lplRetIndividu); //&lplRetIndividu[0]; // INTEGER_POINTER( lplRetIndividu );
	plRetPere 	= INTEGER(lplRetPere); //&lplRetPere[0]; // INTEGER_POINTER( lplRetPere );
	plRetMere 	= INTEGER(lplRetMere); //&lplRetMere[0]; // INTEGER_POINTER( lplRetMere );
	plRetSexe 	= INTEGER(lplRetSexe); //&lplRetSexe[0]; // INTEGER_POINTER( lplRetSexe );

	mustsort = INTEGER(smustsort);
	
	//CREATION DU TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(genealogie,GFALSE,&lNIndividu,&Noeud); //Pas d'enfant
	
	int havesex = (LoadNIndMasc()!=-1);
	//creation vecteur de sortie
	for(int i=0;i<lNIndividu;i++)
	{
		plRetIndividu[i]=Noeud[i].nom;
		if (Noeud[i].pere!=NULL)	plRetPere[i]=Noeud[i].pere->noind+1;		
		else						plRetPere[i]=0;
			
		if (Noeud[i].mere!=NULL)	plRetMere[i]=Noeud[i].mere->noind+1;			
		else						plRetMere[i]=0;
		//Recuperation du sexe
		if (havesex)				plRetSexe[i]=Noeud[i].sex;
		else						plRetSexe[i]=-1; //Pas utilisable
	}

	//Trie si nécessaire
	if (*mustsort)
		SortGenealogie3Vecteur(plRetIndividu,plRetPere,plRetMere,plRetSexe,lNIndividu);
	STOPTIMER;

	return R_NilValue;
}



//********************************************
//     FONCTIONS POUR L ANALYSE STATISTIQUE
//******************************************** 


/// Fonction d'interface Splus pour FondParGen
/** \sa initImplexe()*/
RcppExport SEXP SPLUSFondParGen(SEXP sgenealogie, SEXP sprop, SEXP snbProp, SEXP sretour)
{
	STARTTIMER;
	int * genealogie, * prop, * retour, *nbProp;
	
	Rcpp::IntegerVector lgenealogie	( sgenealogie );
	Rcpp::IntegerVector lprop		( sprop );
	Rcpp::IntegerVector lretour		( sretour );
	
	genealogie = INTEGER(lgenealogie); //&lgenealogie[0]; // INTEGER_POINTER( lgenealogie );
	prop		 = INTEGER(lprop); //&lprop[0]; // INTEGER_POINTER( lprop );
	retour 	 = INTEGER(lretour); //&lretour[0]; // INTEGER_POINTER( lretour );

	//int nbProp = Rcpp::as<int>( snbProp );
	nbProp = INTEGER(snbProp); //INTEGER_POINTER( snbProp );

	FondParGen(genealogie, prop, *nbProp, retour); //nbProp, 
	STOPTIMER;

	return R_NilValue;
}

