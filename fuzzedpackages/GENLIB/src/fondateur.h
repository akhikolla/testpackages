/*! \file fondateur.h
\brief Interface des fonctions de simulation calcul de probabilite

Interface de toutes les fonctions en rapport avec le gene fondateur

\author Sébastien Leclerc 
\contributor Jean-François Lefebvre

*/

#ifndef GENFOND
#define GENFOND

#include <RcppCommon.h>

int simul(int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
		int lSimul, double* pdRetConj,double* pdRetSimul,double* pdRetProp,double* probRecomb,double probSurvieHomo,int printprogress);

int simulsingle(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
			 int lSimul, double* pdRetour,int printprogress);

int simulsingleFreq(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
				int lSimul, double* pdRetour,int printprogress);

SEXP simulsingleFct(int* SGenealogie, int * proposant, int lproposant, int* SplAncetre, int* SplAncEtatAll1, int* SplAncEtatAll2, int SlNAncetre,
				int SlSimul, SEXP SfctSousGrp, int Sprintprogress);

SEXP simulsingleProb(int* SGenealogie, int* SplProposant, int SlNProposant, int* SplAncetre,int SlNAncetres, int* SplAncEtat,SEXP mtProb,
				 int SlSimul, int Sprintprogress);

SEXP prob(int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
	    double* pdRetConj,double* pdRetSimul,int printprogress,int onlyConj);

int CoefApparentement(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre,	double* pdRetour,int DuppDetection, int printprogress);

#endif



