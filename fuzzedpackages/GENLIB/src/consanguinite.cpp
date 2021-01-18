/*! \file consanguinite.cc
\brief Implementation des fonctions de calcul de F

Calcul et Analyse de diverse valeur dérivé de F et Fmoyen

\author Sébastien Leclerc
\contributor Jean-François Lefebvre

*/


/// Authorise l'affichage du niveau de progression sur la sortie standard
/** 
	Si ALLOWPRINTPROGRESS est defini, les fonctions suivantes peuvent
	indiquer leur niveau de progression sur la sortie standard stdout
*/
#define ALLOWPRINTPROGRESS

#include "base.h"
#include "outils.h"
#include "consanguinite.h"
#include "apparentement.h"
#include "userInterface.h"
#include <limits.h>
#include <math.h>
#include <cstdlib>
#include <Rcpp.h>
#define R_NO_REMAP
// ********************************************************************
//
//			PUBLIC 
//
// ********************************************************************

/*! 
	\brief Calcul de la consanguinité (F)

	Calcul le coefficient de consanguinité (F) de chaque proposant.

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param proposant	[in] Vecteur des no de proposant à étudier
	\param NProposant	[in] Nombre d'élément du vecteur proposant
  
	\param Niveau		[in] Profondeur maximal de la recherche

	\retval pdConsan	[out] Pointeur vers un vecteur de intueur NProposant
						En cas de succes, la valeur de la consanguinite de chaque proposant
						est ecrite dans le vecteur

	\param printprogress [in] Imprime un message indiquant les progress accomplies

	\return 0 si la fonction est executé avec succès

	\remark lors de l'appel de kinship sur les 2 parents, le Niveau est reduit de 1 pour montrer la remonte
			necessaire pour ce rendre au parent. (valide?)
*/
int consan(int* Genealogie, int* proposant, int NProposant,int Niveau, double* pdConsan, int printprogress)
{
	try{
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GFALSE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(proposant,NProposant,&NoeudPro);
	
	//Remise à zero des valeurs utilisé
	int i=0;
	for(i=0;i<lNIndividu;i++)
		Noeud[i].prob[1]=-1;	//CONTRIBUTION DES 2 PARENTS	
	
	//Mise en place du niveau maximal
	if (Niveau<=0)
		Niveau=SHRT_MAX;

	if (Niveau>SHRT_MAX) {
//		GENError("Niveau must be smaller than %d",SHRT_MAX);
		char erreur[TAILLEDESCRIPTION];
		sprintf(erreur, "Niveau must be smaller than %d",SHRT_MAX);
		throw std::range_error(erreur);
		//GENError("Le niveau doit-être inférieur à %d",SHRT_MAX);
	}
	const short niveauMax = short(Niveau);

	//Pour chaque proposant
	CREATE_PROGRESS_BAR(NProposant,printprogress)
	for(i=0;i<NProposant;i++)
	{
		if (NoeudPro[i]->pere && NoeudPro[i]->mere)	
			pdConsan[i] = Kinship(NoeudPro[i]->pere,NoeudPro[i]->mere,niveauMax-1,niveauMax-1);			
		else
			pdConsan[i]	=	0;

		//Affichage des progress
		INCREMENT_PROGRESS_BAR();
	}
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			}
 			  catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


/*! 
	\brief Liste de F moyens pour différentes profondeurs

  Donne une liste de F moyen pour différentes profondeurs

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param proposant	[in] Vecteur des no de proposant à étudier
	\param NProposant	[in] Nombre d'élément du vecteur proposant
  
	\param NiveauMin	[in] Nombre de génération à prendre en compte
	\param NiveauMax	[in] Nombre de génération à prendre en compte
  
	\retval pdDeepConsan [out] Un pointeur vers un vecteur de intueur NProposant*(NiveauMax-NiveauMin+1)
						  En cas de success, contient la consanguinite de chaque proposant et ce pour chaque niveau 
						  de NiveauMin a NiveauMax.
						  <br>Ex:
						  <br>&nbsp; &nbsp;IndPro: Indice du proposant
						  <br>&nbsp; &nbsp;NiveauMin : 2
						  <br>&nbsp; &nbsp;NiveauMax : 4
						  <br>&nbsp; &nbsp;NiveauVoulu :3   -> donc IndVoulu= 3-NiveauMin=1
						  <br>&nbsp; &nbsp;pdDeepConsan[IndVoulu*NProposant+IndPro]=0;

	\param printprogress [in] Imprime un message indiquant les progress accomplies

    \return 0 si la fonction est executé avec succès

*/
int consanFs(int* Genealogie, int* proposant, int NProposant,int NiveauMin,int NiveauMax,
		   double* pdDeepConsan, int printprogress)
{
	try{
	//TEST D'ERREUR DE BASE
	if (NProposant<1){
//		GENError("At least two probands are required for this function");
		throw std::range_error("At least one proband is required for this function");
		//GENError("Il faut au minimum 2 proposant pour utilise cette fonction");
	}
	if (NiveauMin<1){
//		GENError("depthmax and depthmin must be greater than one.");
		throw std::range_error("depthmax and depthmin must be greater than one.");
		//GENError("Le niveau minimum et le niveau maximum doivent-être supérieur à un");
	}
	if (NiveauMax<NiveauMin){
//		GENError("depthmax must be greater or equal to depthmin");
		throw std::range_error("depthmax must be greater or equal to depthmin");
		//GENError("Le niveau maximum doit-être supérieur ou égal au niveau minimum");
	}
	//Mise en place du niveau maximal
	if (NiveauMax>SHRT_MAX){
//		GENError("depthmax must be smaller than %d",SHRT_MAX);
		char erreur[TAILLEDESCRIPTION];
		sprintf(erreur, "depthmax must be smaller than %d",SHRT_MAX);
		throw std::range_error(erreur);
		//GENError("Le NiveauMax doit-être inférieur à %d",SHRT_MAX);
	}
	const short trueNiveauMax = NiveauMax-1;
	const int trueNiveauMin = NiveauMin-1;

	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GFALSE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(proposant,NProposant,&NoeudPro);
	
	//Creation du tableau Noeud,proposant et du tableau de resultat pour FS
	const int LongEcart		= NiveauMax-NiveauMin+1;

	INITGESTIONMEMOIRE;	
	double *Phideep		       	= (double*)	memalloc(NiveauMax,sizeof(double));
	int i;
		
	//Initialisation 
	//Initialisation des noeuds
	for(i=0;i<lNIndividu;i++)	
		Noeud[i].pGen=NULL;			

	
	//CALCUL LUI MEME		
	Kinship4Struct Element(trueNiveauMax,Phideep);	
	
	//Pour chaque proposant
//	CREATE_PROGRESS_BAR(NProposant,printprogress)
	for(int cPro=0;cPro<NProposant;++cPro)
    {
			
		if (NoeudPro[cPro]->pere && NoeudPro[cPro]->mere)	
		{	//Ce proposant a un pere et une mere		

			//remise à zero de phi deep
			for(int a=0;a<=trueNiveauMax;a++)
				Phideep[a]=0.0;

			Kinship4(NoeudPro[cPro]->pere,NoeudPro[cPro]->mere,trueNiveauMax,trueNiveauMax,Element);
			
			for(int a=0;a<LongEcart;a++)			
				pdDeepConsan[a*NProposant+cPro]=Phideep[a+trueNiveauMin];							
		}
		else
		{   //Ce proposant n'as pas de parent (facon de parler)
			for(int a=0;a<LongEcart;a++)			
				pdDeepConsan[a*NProposant+cPro]=0.0;							
		}

		//Affichage des progress
//		INCREMENT_PROGRESS_BAR();

	}//Fin chaque proposant

	return 0;
 	} catch(std::exception &ex) {
 		forward_exception_to_r(ex);
 	}
 	  catch(...){
 		::Rf_error("c++ exception (unknown reason)"); 
 	} 
 	return 0;
}

