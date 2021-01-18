/*! \file fondateur.cc
\brief Implementation des fonctions de simulation calcul de probabilite

Calcul et Analyse de diverse valeur dérivé du gene fondateur


\author Sébastien Leclerc 
\contributor Jean-Francois Lefebvre
*/


/// Authorise l'affichage du niveau de progression sur la sortie standard
/** 
	Si ALLOWPRINTPROGRESS est defini, les fonctions suivantes peuvent
	indiquer leur niveau de progression sur la sortie standard stdout
*/
#define ALLOWPRINTPROGRESS
#include "base.h"
#include "outils.h" 
#include "cbignum.h"
#include "hashtable.h"
#include "userInterface.h"
#include "fondateur.h"
#include <iostream>
#include <time.h>
#include <string.h> 
#include <math.h> 
#include <limits.h>
#include <vector>
#include <random>
#include <cstdlib>
#include <RcppCommon.h>
#include <R.h>
//#include <Rdefines.h>

#include <Rcpp.h>
#include <Rcpp/as.h>
#include <Rcpp/Function.h>
//#include <boost/random/random_device.hpp>

#define R_NO_REMAP
//using namespace std;

#ifdef NEG
	#undef NEG
#endif
extern "C"
{
	#include "mpi.h"
	#include "mplogic.h"
}

// ******************************************************************** 
//
//			CONSTANTE & STRUCTURE
//
// ********************************************************************

const int MEGAOCTET = 1048576;//octet

///Liste de tableau de unsigned char.... Utilisé pour representer une liste de nombre binaire de longueur variable
struct CApPath
{
	///Grand nombre binaire en format mpi représentant le chemin....
	mp_int num;
	///Pointeur vers l'element suivant de la liste
	CApPath *next;
};

//CONSTRUCTION ET DESTRUCTION DES PATHS
static CIndSimul** g_ExpCoeff_CheminParcouru=NULL; 
///Utilise par ExploreCoeff comme étant le dernier chemin remplis (ou le premiers)
static CApPath ** g_ExpCoeff_Path=NULL;
///Utilise par ExploreCoeff comme étant la cible de l'exploration
static CIndSimul* g_ExpCoeff_Cible=NULL;
static void FASTCALL ExploreCoeff(CIndSimul* Noeud);
static void PathDestruction(CApPath **Path,int npath);

// ********************************************************************
//
//			PUBLIC
//
// ********************************************************************
/*! 
	\brief Execute une simulation pour determiner les probabilites du passage d'un allele a une serie de proposant

	Calcule les probabilite qu'un proposants recoivent 0,1-2,2 allele a partir d'un groupe d'ancetre
	Calcule la probabilite conjointe que chaque proposant soit atteint
	Calcule la probatilite que de un..n proposant soit atteint
	Calcule la probabilie que chaque proposant soit atteint

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant à étudier
	\param plProEtat    [in] vecteur de taille lNProposant et representant l'etat a considerer pour chaque proposant
			<br>&nbsp; &nbsp;&nbsp; &nbsp;0: La condition est remplie si se proposant n'est pas malade 
			<br>&nbsp; &nbsp;&nbsp; &nbsp;1: La condition est remplie si le proposant recois 1-2 allele 
			<br>&nbsp; &nbsp;&nbsp; &nbsp;2: La condition est remplie si le proposant recois 2 allele 
	\param lNProposant	[in] Nombre d'élément du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant à étudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre

	\param lSimul		[in] Nombre de simulation a effectuer

	\retval pdRetConj	[out] Un pointeur vers un double. 
					En cas de succes, le double represente la probabilite conjointe que la condition de chaque proposant soit remplis. 
	\retval pdRetSimul	[out] Un pointeur vers une vecteur de taille lNProposant.. 
					En cas de succes, Probabilite que la condition de chaque proposant soit remplis
	\retval pdRetProp	[out] Un pointeur vers une vecteur de taille lNProposant+1.
					En cas de succes represente la probabilite que 0,1,2..n condition soit remplis
									
	\param printprogress [in] imprime un message indiquant les progress accomplies

	\return 0 si la fonction est executé avec succès 
*/ 
int simul(int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
		int lSimul, double* pdRetConj, double* pdRetSimul, double* pdRetProp, double* probRecomb, double probSurvieHomo,
		int printprogress)
{	
	try{
	//Validation genealogie
	if (lSimul<=0){
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulation doit-être supérieur à zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);

	//Creation des tableau
	INITGESTIONMEMOIRE;
	CIndSimul** Ordre		=(CIndSimul**) memalloc(lNIndividu,sizeof(CIndSimul*));	
	
	//Pour le sort spécial		
	int*		OrdreSaut	=(int*) memalloc(lNIndividu,sizeof(int*));				
	int NOrdre;
	
	//RECEUIL DES STATISTIQUES
	int *ProCompteur		=(int*)memalloc(lNProposant,sizeof(int*));
	int *NCompteur			=(int*)memalloc(lNProposant+1,sizeof(int*));
	int bConj;
	
	//initialisation
	//initrand();
	int i,j;
	int ap,am;
	//int *VecteurPosition=NULL;

	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(i=0;i<lNIndividu;i++)
	{
		Noeud[i].allele=0;	
		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;	   				
	}
	
	//identifier et etiqueter les proposant
	for(i=0;i<lNProposant;i++)
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
					
	//identifier et etiqueter les points de departs
	for(i=0;i<lNAncetre;i++)
	{		
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);		
	}

	//identifier et marque les noeuds utile et ceux inutile a la recherche
	for(i=0;i<lNAncetre;i++)
		ExploreArbre(NoeudAnc[i]);

	
	//creation d'un ordre d'execution et calcul des sauts
	PrepareSortPrioriteArbre(Noeud,lNIndividu);	
	NOrdre=0;
	memset(OrdreSaut,0,sizeof(int)*lNIndividu);
	for(i=0;i<lNAncetre;i++)		
		StartSortPrioriteArbre(NoeudAnc[i],Ordre,&NOrdre,OrdreSaut);
	
//	unsigned seed2 = time(0);
//	std::mt19937 gen(seed2);

//	boost::random_device rd;
//	std::random_device rd;		// **chgt IGES**
//	std::mt19937 gen(rd());
	
	#if defined _WIN32 || defined _WIN64
	  std::mt19937 gen(time(0));
	#else
//	  boost::random_device gen;
	  std::random_device gen;
	#endif
 	int nbannulee = 0;		// **chgt IGES**
 	int nbCasHomo = 0;		// **chgt IGES**

	//Simulation
	memset(ProCompteur,0,sizeof(int)*lNProposant);
	memset(NCompteur,0,sizeof(int)*(lNProposant+1));
	
	//CREATE_PROGRESS_BAR(lSimul,printprogress)
	for(int csimul=0;csimul<lSimul; csimul++)
	{
		//Par ordre du parent -> enfant
		//les 2 sorts n'ordonne pas dans le meme sens
		bool simAnnulee = false;
		for(int i=0;i<NOrdre;i++)
		{
			//croisement p/r au parent
			//double iRandom =urand();
			double iRandom = (double)gen()/(double)gen.max();
			//double iRandom = (double)rd()/(double)rd.max();		// pour le passage de l'allele.  **chgt IGES**
			
			if (Ordre[i]->pere!=NULL)	ap=	Ordre[i]->pere->allele;	// Ordre=vecteur de ptr vers les elts tries. allelePere (ap)
			else						ap=0;					// si pas de ptr, ap = 0.
			
			if (Ordre[i]->mere!=NULL)	am=	Ordre[i]->mere->allele;	// alleleMere (am) existe et on l'attribe
			else						am=0;					// existe pas donc am = 0.
			
			if (iRandom<TransGenCum[ap][am][0])
			{
				Ordre[i]->allele=0;
				//Saute un certain nombre de noeud
				j=i+OrdreSaut[i];
				while (i!=j)	Ordre[++i]->allele=0;
			}
			else
			{
				double alea    = (double)gen()/(double)gen.max();//(double)rd()/(double)rd.max(); // pour la recombinaison.**chgt IGES**
				int sex = Ordre[i]->sex;
				if(alea < probRecomb[1]) // tx femme (plus élevé)
				{
				// si femme ou si prob < au taux male et inconnu (le plus petit)
				// si le morceau recombine (selon tx male et aussi sexe inconnu) OU si le morceau recombine (selon tx femelle)
				if( sex == GEN_FEM || alea < probRecomb[0] )
				{
					Ordre[i]->allele=0;
					//Saute un certain nombre de noeud
					j=i+OrdreSaut[i];
					while (i!=j)	Ordre[++i]->allele=0;
				}
				}
				else if (iRandom<TransGenCum[ap][am][1])
					Ordre[i]->allele=1;			
				else
				{
					Ordre[i]->allele=2;
					nbCasHomo++;
					//Descendance conditionnelle aux nombre d'alleles recus
					double alea2 = (double)gen()/(double)gen.max(); // 3e aleatoire pour determiner les cas homo
					if(alea2 > probSurvieHomo){ simAnnulee = true; break; }
				}
			}
		}
		if(!simAnnulee) {
			//Modification des compteurs
			bConj=0;
			for(int i=0;i<lNProposant;i++)
			{
//				if ( (plProEtat[i]==0 && NoeudPro[i]->allele==0) || (plProEtat[i]!=0 && NoeudPro[i]->allele>=plProEtat[i]))
				if ( (plProEtat[i]==0 && NoeudPro[i]->allele==0) || 
					(plProEtat[i]==1 && NoeudPro[i]->allele==1) ||
					(plProEtat[i]==2 && NoeudPro[i]->allele==2) ||
					(plProEtat[i]==3 && NoeudPro[i]->allele>=1) )
				{
					++ProCompteur[i];
					++bConj;
				}
			}
			++NCompteur[bConj];
			//Barre de progress
			//INCREMENT_PROGRESS_BAR()
		}
		else { // la simulation est annulee pcq la sim ne concorde pas avec la genealogie
			csimul--;
			nbannulee++;
		}
	}

	//ECRITURE DE LA VALEUR DE RETOUR
	//double* pdRetConj,double* pdRetSimul,double* pdRetProp
	for(int i=0;i<lNProposant;i++)
	{
		//printf("%f\n", double(ProCompteur[i]));
		pdRetSimul[i]=double(ProCompteur[i])/double(lSimul);
		pdRetProp[i]=double(NCompteur[i])/double(lSimul);
	}
	pdRetProp[lNProposant]=double(NCompteur[lNProposant])/double(lSimul);
	*pdRetConj=double(NCompteur[lNProposant])/double(lSimul);
	
	//FIN
	//outrand();
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}




/*! 
	\brief Execute une ou plusieurs simulation et retourne le nombre d'allele transmit a chaque proposant (pour chaque simulation)

	Calcule un etat possible pour chaque proposant en tenant compte de chaque ancetre et son etat

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant à étudier
	\param lNProposant	[in] Nombre d'élément du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant à étudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre
	
	\param lSimul		[in] Nombre de simulation à effectuer

	\retval pdRetour	[out] Pointeur vers un vecteur de NProposant x lSimul
								En cas de succes, ce vecteur le nombre d'allele assigne a chaque proposant pour la simulation
	
	\param printprogress imprime un message indiquant les progress accomplies

	\return 0 si la fonction est executé avec succès
*/
int simulsingle(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre, int lSimul,
			 double* pdRetour,int printprogress)
{
	try{
	//VALIDATION GENEALOGIE
	if (lSimul<=0){
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulation doit-être supérieur à zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	/**D***/
	//Creation des tableau
	INITGESTIONMEMOIRE;
	CIndSimul** Ordre		=(CIndSimul**) memalloc(lNIndividu,sizeof(CIndSimul*));	
	
	//Pour le sort spécial		
	int*		OrdreSaut	=(int*) memalloc(lNIndividu,sizeof(int*));				
	int NOrdre;
	/**F***/
	//CREATION DES TABLEAU       
	int i;
	//int *VecteurPosition=NULL;

	/**D***/
	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(i=0;i<lNIndividu;i++)
	{
		Noeud[i].allele=0;	
		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;	   				
	}
	//INITIALISATION DE LA STRUCTURE DE NOEUD
//	for(i=0;i<lNIndividu;i++)
//		Noeud[i].allele=0;			
	/**F***/
	
	//IDENTIFIER ET ETIQUETER LES PROPOSANT
	for(i=0;i<lNProposant;i++)
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;					

	//IDENTIFIER ET ETIQUETER LES POINTS DE DEPARTS
	for(i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);		
	}

	/**D***/
	//identifier et marque les noeuds utile et ceux inutile a la recherche
	for(i=0;i<lNAncetre;i++)
		ExploreArbre(NoeudAnc[i]);

	
	//creation d'un ordre d'execution et calcul des sauts
	PrepareSortPrioriteArbre(Noeud,lNIndividu);	
	NOrdre=0;
	memset(OrdreSaut,0,sizeof(int)*lNIndividu);
	for(i=0;i<lNAncetre;i++)		
		StartSortPrioriteArbre(NoeudAnc[i],Ordre,&NOrdre,OrdreSaut);
	/**F***/

	//INITIALISATION
	//initrand();
//	unsigned seed2 = time(0);
//	std::mt19937 gen(seed2);
//	boost::random_device rd;
//	std::random_device rd;		// **chgt IGES**
//	std::mt19937 gen(rd());
	#if defined _WIN32 || defined _WIN64
	  std::mt19937 gen(time(0));
	#else
	std::random_device gen;		// **chgt IGES**
//	  boost::random_device gen;
	#endif
	int ap,am;
		
	//Partie 3: SIMULATION
	//CREATE_PROGRESS_BAR(lSimul,printprogress)
	for(int csimul=0;csimul<lSimul; csimul++)
	{
/* Utilisant la facon de faire de la fonction simul */
		for(int i=0;i<NOrdre;i++)
		{
			if (Ordre[i]->pere!=NULL) ap=	Ordre[i]->pere->allele;
       		else					ap=0;
				
			if (Ordre[i]->mere!=NULL) am=	Ordre[i]->mere->allele;
       		else					am=0;
						
			if (ap==0 && am==0)		Ordre[i]->allele=0;
			else {
				double iRandom = (double)gen()/(double)gen.max();
				if (iRandom<TransGenCum[ap][am][0])		Ordre[i]->allele=0;
				else
					if (iRandom<TransGenCum[ap][am][1]) Ordre[i]->allele=1;			
					else							 Ordre[i]->allele=2;
			}
		}
	
/* Fin de la nouvelle facon de faire */ 

/* Code remplace par la facon de faire de la fonction simul
		for(i=0;i<lNIndividu;i++)
		{	  
			if (Noeud[i].etat!=GENDEPART)
			{
				if (Noeud[i].pere!=NULL) ap=	Noeud[i].pere->allele;
	       		else					 ap=0;
				
				if (Noeud[i].mere!=NULL) am=	Noeud[i].mere->allele;
	       		else					 am=0;
						
				if (ap==0 && am==0)		Noeud[i].allele=0;
				else {
					//double iRandom =urand();
					//double iRandom = (double)gen()/(double)rd.max();
					double iRandom = (double)gen()/(double)gen.max();//(double)rd()/(double)rd.max();		// pour le passage de l'allele.  **chgt IGES**
					if (iRandom<TransGenCum[ap][am][0])		Noeud[i].allele=0;
					else
						if (iRandom<TransGenCum[ap][am][1]) Noeud[i].allele=1;			
						else								Noeud[i].allele=2;
				} //fin ap==am==0
			}//fin if noeud depart
		}// fin for pour chaque noeud
** Fin du code remplace */		

        //Modification des compteurs
		const int tmp = csimul*lNProposant;
		for(i=0;i<lNProposant;i++)
			pdRetour[tmp+i]=NoeudPro[i]->allele;		

		//Barre de progress
		//INCREMENT_PROGRESS_BAR();

	}// Fin du pour chaque simulation
	///nettoyer
	//outrand();
	return 0;
 } catch(std::exception &ex) {
 	forward_exception_to_r(ex);
 } catch(...){
 	::Rf_error("c++ exception (unknown reason)"); 
 } 
 return 0;
}
/*! 
	\brief Execute une plusieur simulation et retourne le nombre d'allele transmit a chaque proposant (pour chaque simulation)

	Calcule un etat possible pour chaque proposant en tenant en comple chaque ancetre et son etat

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant à étudier
	\param lNProposant	[in] Nombre d'élément du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant à étudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre
	
	\param lSimul		[in] Nombre de simulation à effectuer
	\param mtProb		[in] Matrice S-PLUS: tableau des probabilités

	\retval pdRetour	[out] Pointeur vers un vecteur de NProposant x lSimul
								En cas de succes, ce vecteur le nombre d'allele assigne a chaque proposant pour la simulation
	
	\param printprogress imprime un message indiquant les progress accomplies

	\return 0 si la fonction est executé avec succès
*/

SEXP simulsingleProb(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre, int* plAncEtat, SEXP mtProb,
				 int lSimul, int printprogress)
{
	try{
	//Conversion des paramètres
	Rcpp::NumericMatrix matprob(mtProb);

	//VALIDATION GENEALOGIE
	if (lSimul<=0){
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulations doit-être supérieur à zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	//CREATION DES TABLEAUX       
	int i;
	//long *VecteurPosition=NULL;

	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(i=0;i<lNIndividu;i++)
		Noeud[i].allele=0;			
	
	//IDENTIFIER ET ETIQUETER LES PROPOSANT
	for(i=0;i<lNProposant;i++)
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;					

	//IDENTIFIER ET ETIQUETER LES POINTS DE DEPARTS
	for(i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);		
	}

	//INITIALISATION
	//initrand();
//	unsigned seed2 = time(0);
//	std::mt19937 gen(seed2);
//	boost::random_device rd;
//	std::random_device rd;		// **chgt IGES**
//	std::mt19937 gen(rd());
	#if defined _WIN32 || defined _WIN64
	  std::mt19937 gen(time(0));
	#else
//	  boost::random_device gen;
	  std::random_device gen;
	#endif
	int ap,am;
		
	//Pointeur de retour
	//CSPnumeric tmp(lNProposant*lSimul);
	Rcpp::IntegerVector tmp(lNProposant*lSimul);
	
	//Partie 3: SIMULATION
	//CREATE_PROGRESS_BAR(lSimul,printprogress)
	for(int csimul=0;csimul<lSimul; csimul++)
	{
		for(i=0;i<lNIndividu;i++)
		{	  
			if (Noeud[i].etat!=GENDEPART)
			{	
				if (Noeud[i].pere!=NULL)	ap=	Noeud[i].pere->allele;
	       		else						ap=0;
				
				if (Noeud[i].mere!=NULL)	am=Noeud[i].mere->allele;
	       		else						am=0;
						
				if (ap==0 && am==0)								Noeud[i].allele=0;
				else {
//					am += 1;
//					ap += 1;
					if (Noeud[i].sex == GEN_FEM) am += 6;
//					double iRandom =urand(); 											
					//double iRandom = (double)gen()/(double)rd.max();
					double iRandom = (double)gen()/(double)gen.max();//(double)rd()/(double)rd.max();		// pour le passage de l'allele.  **chgt IGES**
					if (iRandom<(double)matprob(ap,am))			Noeud[i].allele=0;		       
					else
						if (iRandom<(double)matprob(ap,am+=3))	Noeud[i].allele=1;			
						else									Noeud[i].allele=2;
				} //fin ap==am==0
			}//fin if noeud depart
		}// fin for pour chaque noeud
	
        //Modification des compteurs
		for(i=0;i<lNProposant;i++)
			tmp(csimul*lNProposant+i)=(int)NoeudPro[i]->allele;	 //csimul*lNProposant+i+1
				
		//Barre de progress
		//INCREMENT_PROGRESS_BAR();

	}// Fin du pour chaque simulation
	///nettoyer
	//outrand();
	return tmp; //.Detach();
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}

/*! 
	\brief Execute une plusieur simulation et retourne le nombre d'allele transmit a chaque proposant en tableau
			de fréquences pour toutes les simulations.

	Calcule un etat possible pour chaque proposant en tenant en compte chaque ancetre et son etat

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant à étudier
	\param lNProposant	[in] Nombre d'élément du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant à étudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre
	
	\param lSimul		[in] Nombre de simulations à effectuer

	\retval pdRetour	[out] Pointeur vers un vecteur de NProposant x 0-1-2 pour la fréquence d'allèles tramises pour toutes les
							  simulations. 
	
	\param printprogress imprime un message indiquant les progress accomplies

	\return 0 si la fonction est executé avec succès
*/
int simulsingleFreq(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
				int lSimul, double* pdRetour,int printprogress)
{
	try{
	//VALIDATION GENEALOGIE
	if (lSimul<=0) {
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulation doit-être supérieur à zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	//CREATION DES TABLEAU    
	int i;
	//int *VecteurPosition=NULL;

	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(i=0;i<lNIndividu;i++)
		Noeud[i].allele=0;		
	
	//IDENTIFIER ET ETIQUETER LES PROPOSANT
	for(i=0;i<lNProposant;i++)
		NoeudPro[i]->etat = GENPROPOSANTINUTILE;					

	//IDENTIFIER ET ETIQUETER LES POINTS DE DEPARTS
	for(i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);		
	}

	//Indice de tableau pour la variable de sortie
	const int tmp0 = 0*lNProposant; //Fréquence de 0 allèle
	const int tmp1 = 1*lNProposant; //Fréquence de 1 allèle
	const int tmp2 = 2*lNProposant; //Fréquence de 2 allèle

	//INITIALISATION
	//initrand(); 
//	unsigned seed2 = time(0);
//	std::mt19937 gen(seed2);
//	boost::random_device rd;
//	std::random_device rd;		// **chgt IGES**
//	std::mt19937 gen(rd());
	#if defined _WIN32 || defined _WIN64
	  std::mt19937 gen(time(0));
	#else
//	  boost::random_device gen;
	  std::random_device gen;
	#endif
	int ap,am;
		
	//Partie 3: SIMULATION
	//CREATE_PROGRESS_BAR(lSimul,printprogress)
	for(int csimul=0;csimul<lSimul; csimul++)
	{
		for(i =0;i<lNIndividu;i++)
		{	  
			if (Noeud[i].etat!=GENDEPART)
			{
				if (Noeud[i].pere!=NULL) ap = Noeud[i].pere->allele;
	       		else					 ap=0;
				
				if (Noeud[i].mere!=NULL) am = Noeud[i].mere->allele;
	       		else					 am=0;
						
				if (ap==0 && am==0)							Noeud[i].allele=0;
				else {
					//double iRandom = urand();
					//double iRandom = (double)gen()/(double)rd.max();
					double iRandom = (double)gen()/(double)gen.max();//(double)rd()/(double)rd.max();		// pour le passage de l'allele.  **chgt IGES**
					if (iRandom<TransGenCum[ap][am][0])		Noeud[i].allele=0;		       
					else
						if (iRandom<TransGenCum[ap][am][1])	Noeud[i].allele=1;			
						else								Noeud[i].allele=2;
				} //fin ap==am==0
			}//fin if noeud depart
		}// fin for pour chaque noeud

		for(i=0;i<lNProposant;i++)
		{
			if (NoeudPro[i]->allele == 0)		pdRetour[tmp0+i]+= 1;
			else if (NoeudPro[i]->allele == 1)	pdRetour[tmp1+i]+= 1;
			else								pdRetour[tmp2+i]+= 1;
		}
		//Barre de progress
		//INCREMENT_PROGRESS_BAR();
	}// Fin du pour chaque simulation
	///nettoyer
	//outrand();
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}
/*! 
	\brief Execute une plusieur simulation et retourne le nombre d'allele transmit a chaque proposant (pour chaque simulation)

	Calcule un etat possible pour chaque proposant en tenant en comple chaque ancetre et son etat

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant à étudier
	\param lNProposant	[in] Nombre d'élément du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant à étudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre
	
	\param lSimul		[in] Nombre de simulation à effectuer

	\retval pdRetour	[out] Pointeur vers un vecteur de NProposant x lSimul
								En cas de succes, ce vecteur le nombre d'allele assigne a chaque proposant pour la simulation
	
	\param printprogress imprime un message indiquant les progress accomplies

	\return 0 si la fonction est executé avec succès
*/
SEXP simulsingleFct(int * Genealogie, int * proposant, int lproposant, int * plAncetre, int * plAncEtatAll1, int * plAncEtatAll2, int lNAncetre, int lSimul, SEXP SfctSousGrp, int printprogress)
{	
	try{
	//VALIDATION GENEALOGIE
	if (lSimul<=0){
//		GENError("Number of simulation must be greater than zero");
		throw std::range_error("Number of simulation must be greater than zero");
		//GENError("Le nombre de simulation doit-être supérieur à zero");
	}
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul ** NoeudPro=NULL; //[i]
	LoadProposant(proposant,lproposant,&NoeudPro); //[i]

	//CREATION D'UN VECTEUR D'ANCETRES
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);

	//CREATION DES TABLEAU       
	//long *VecteurPosition=NULL;

	//INITIALISATION DE LA STRUCTURE DE NOEUD
	for(int i=0;i<lNIndividu;i++)
	{
		Noeud[i].allele2Pos[0] = 0;
		Noeud[i].allele2Pos[1] = 0;
	}

	//IDENTIFIER ET ETIQUETER LES POINTS DE DEPARTS
	for(int i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele2Pos[0]=interval(plAncEtatAll1[i],0,5);	
		NoeudAnc[i]->allele2Pos[1]=interval(plAncEtatAll2[i],0,5);		
	}	

	//INITIALISATION
//	initrand();
	int lptrAp[2];
	int lptrAm[2];

	//Partie 3: SIMULATION	
	//Déclaration de la liste de résultats de la fonction de l'utilisation pour chaque simulation
	//CSPlist resultFct;
	Rcpp::List resultFct; //(lSimul);
	Rcpp::Function f(SfctSousGrp);
	//CREATE_PROGRESS_BAR(lSimul,printprogress)

	for(int csimul=0;csimul<lSimul; csimul++)
	{
		for(int i=0;i<lNIndividu;i++)
		{	  	
			if (Noeud[i].etat!=GENDEPART)
			{	
				if (Noeud[i].pere!=NULL)
				{
					lptrAp[0] = Noeud[i].pere->allele2Pos[0];
					lptrAp[1] = Noeud[i].pere->allele2Pos[1];
				}
	       		else
					lptrAp[0] = lptrAp[1] = 0;
				if (Noeud[i].mere!=NULL)
				{
					 lptrAm[0] = Noeud[i].mere->allele2Pos[0];
					 lptrAm[1] = Noeud[i].mere->allele2Pos[1];
				}
	       		else
					lptrAm[0] = lptrAm[1] = 0;
						
				int iRandom = irand(0,1);  // utilise random_device a travers irand()
				Noeud[i].allele2Pos[0]=lptrAp[iRandom];
				iRandom = irand(0,1);
				Noeud[i].allele2Pos[1]=lptrAm[iRandom];

			}//fin if noeud depart

		}// fin for pour chaque noeud	
		/*
		//Déclaration de la liste de résultats
		//CSPlist resultGrp;
		Rcpp::List resultGrp(igrp);
		//Pour chaque groupe
		for (int i=0;i<igrp;i++)
		{	
			SEXP grp = indSelect[i];//extraction des individus du groupe

			//CSPnumeric lind(grp);
			Rcpp::IntegerVector lind(grp);
			
			//CSPnumericMatrix ans(lind.length(), 2); //output matrix.
			Rcpp::IntegerMatrix ans(lind.size(), 2);
			
			//ans.SetRowNames(CSPcharacter(grp));
			ans.attr("rownames") = Rcpp::as<Rcpp::CharacterVector>(grp);
			
			//CSPnumeric tmp(lind.length());
			Rcpp::IntegerVector tmp (lind.size());
			
			//CSPnumeric tmp2(lind.length());
			Rcpp::IntegerVector tmp2 (lind.size());
			
			for(int j=0;j<lind.length();j++)//Pour chaque individus du groupe
			{
				 // **tmp(j+1) = NoeudPro[i][j]->allele2Pos[0];	//Résultat dans un tableau	
				 ans(j+1, 0) = NoeudPro[i][j]->allele2Pos[0];
				 // **tmp2(j+1) = NoeudPro[i][j]->allele2Pos[1];	//Résultat dans un tableau	
				 ans(j+1, 1) = NoeudPro[i][j]->allele2Pos[1];
			}
			// **ans.SetJthColumnDirect(0, tmp.Detach());
			
			// **ans.SetJthColumnDirect(1, tmp2.Detach());
			//En liste
			// **resultGrp.Add(ans);
			resultGrp.push_back(ans);
		}
		// **CSPfunction f(SfctSousGrp); //Call the S function.
		//Rcpp::Function f(SfctSousGrp);
		//En liste de résultats de la fct de l'utilisateur
		// **resultFct.Add(f.Eval(resultGrp));
		resultFct.push_back(f(resultGrp));
		*/
		
		Rcpp::IntegerMatrix ans(lproposant, 2);

		Rcpp::CharacterVector rowNames(lproposant);
		for(int i=0; i<lproposant; i++) { char nomLigne [10]; /*int n = */sprintf(nomLigne, "%d", proposant[i]); rowNames[i] = nomLigne; }
		//for(int i=0; i<lproposant; i++) sprintf(rowNames[i], "%d", proposant[i]); // itoa(proposant[i], &rowNames[i], 10)

		Rcpp::List dimnms = Rcpp::List::create( rowNames, //Rcpp::as<Rcpp::CharacterVector>(indSelect[i]),
										Rcpp::CharacterVector::create("1", "2"));
		ans.attr("dimnames") = dimnms;
		
		for(int j=0;j<lproposant;j++)//Pour chaque individus du groupe
		{
			 ans(j, 0) = NoeudPro[j]->allele2Pos[0]; //[i] j+1
			 ans(j, 1) = NoeudPro[j]->allele2Pos[1]; //[i] j+1
		}
		resultFct.push_back(f(ans));
		//Barre de progress
		//INCREMENT_PROGRESS_BAR();

	}// Fin du pour chaque simulation

	///nettoyer
	//outrand();
	//Retourne la liste de résultat pour chaque simulation de la fct de l'utilisateur
	return Rcpp::wrap(resultFct); //.Detach();
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}

// *************************************************************************
//	CALCUL DE LA PROBABILITE EXACTE
// ************************************************************************* 
//const int PBARINTERVAL_PROB=9;
//100 devrais être plus que suffisant pour occupé un ordinateur moderne un bon bout de temps
//const int PROB_MAXIMUM_NORDRE = 100; //646 Taille maximale d'un double (utiliser numeric_limits<double>::max() ? )

/*! 
	\brief Evalue la probabilite exacte de transfert d'un allele a une serie de proposant

	Cette fonction détermine la probabilité conjointe qu'un ou plusieurs proposant reçoivent  un certain nombre d'allèle malade.
	En prenant pour acquis qu'un ou plusieurs ancêtre possède un ou deux d'allèle malade. 
	Pour ce faire, cette fonction calcule la valeur de toute les branches de manière à obtenir la probabilité EXACT. 
	Le temps de calcul peut-être très prohibitif.

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant à étudier
	\param plProEtat    [in] vecteur de taille lNProposant et representant l'etat a considerer pour chaque proposant
			<br>&nbsp; &nbsp;&nbsp; &nbsp;0: La condition est remplie si se proposant n'est pas malade 
			<br>&nbsp; &nbsp;&nbsp; &nbsp;1: La condition est remplie si le proposant recois 1-2 allele 
			<br>&nbsp; &nbsp;&nbsp; &nbsp;2: La condition est remplie si le proposant recois 2 allele 
	\param lNProposant	[in] Nombre d'élément du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant à étudier
	\param plAncEtat	[in] Vecteur de taille plAncetre representant le nombre d'allele atteint pour chaque Ancetre (0,1,2) 
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre

	\param OrdreMaximum [in] Le nombre de noeud touche au maximum, si le nombre de noeud touche est trop grand alors la procedure s'interromp automatiquement

	\retval pdRetConj	[out] Un pointeur vers un double. 
								En cas de succes, le double represente la probabilite conjointe que la condition de chaque proposant soit remplis. 
	\retval pdRetSimul	[out] Un pointeur vers une vecteur de taille lNProposant.. 
								En cas de succes, Probabilite que la condition de chaque proposant soit remplis
							
	\return 0 si la fonction est executé avec succès
*/
SEXP prob( int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
		double* pdRetConj,double* pdRetSimul,int printprogress,int onlyConj)
{
	try{
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	if(LoadProposant(plProposant,lNProposant,&NoeudPro) == -1) return Rcpp::wrap(-1);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	//Creation du tableau de Noeud, d'ancetre et proposant
	INITGESTIONMEMOIRE;			
	int NOrdre;
	CIndSimul** Ordre	=(CIndSimul**) memalloc(lNIndividu,sizeof(CIndSimul*));		
	int   *OrdreSaut	= (int*)	memalloc(lNIndividu,sizeof(int));	
	
	//receuil des statistiques
	int i;
	int ap,am;

	//creation et initialisation de la structure noeud et l'ordre
	for(i=0;i<lNIndividu;i++)
	{
		Noeud[i].allele=0;
		
		Noeud[i].prob[0]=0.;
		Noeud[i].prob[1]=0.;
		Noeud[i].prob[2]=0.;			
		Noeud[i].iind=-1; //Pas un proposant

		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;	
	}
	
	//identifier et etiqueter les proposant
	for(i=0;i<lNProposant;i++)
	{
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
		NoeudPro[i]->iind=interval(plProEtat[i],0,2);
	}

	//identifier et etiqueter les points de departs
	for(i=0;i<lNAncetre;i++)
	{
		NoeudAnc[i]->etat=GENDEPART;
		NoeudAnc[i]->allele=interval(plAncEtat[i],0,2);
	}

	//identifier et marque les noeuds utile et ceux inutile a la recherche
	for(i=0;i<lNAncetre;i++)
		ExploreArbre(NoeudAnc[i]);
	
	//creation d'un ordre d'execution et calcul des sauts
	PrepareSortPrioriteArbre(Noeud,lNIndividu);
	NOrdre=0;
	memset(OrdreSaut,0,sizeof(int)*lNIndividu);
	for(i=0;i<lNAncetre;i++)		
		StartSortPrioriteArbre(NoeudAnc[i],Ordre,&NOrdre,OrdreSaut);

	if (NOrdre==-1){		
//		GENError("There is no link between any ancetres and any probands");
		throw std::range_error("There is no link between any ancetres and any probands");
		//GENError("Il n'y a pas de lien entre aucun des ancetres et aucun proposant");
	}
	//Liste de tableau dépendante de l'ordre
	double *Cumul		= (double*) memalloc(NOrdre+1,sizeof(double));

	// Valeur cumulative courante
	//double *PrCumul	  = (double*) memalloc(NOrdre+1,sizeof(double));	
	//Valeur actuelle pour chaque position pour chaque allele
	//double (*PrValue)[3] = (double (*)[3]) memalloc((NOrdre+1),sizeof(double[3]));	
	
	// Debut du calcul
	int n=0;
	Cumul[0]=1;
	int iteration=0;
	const int Lastindex=NOrdre-1;
	double ConjProb=0;
	
	//Initialisation de début de simulation
	for(i=0;i<NOrdre;i++)	
	{
		Ordre[i]->allele=-1;		
	}
	
	//On réutilise la variable bFlagSort pour indique si la condition d'un proposant
	//est acceptable ou non (1 oui, non =0)
	int nbCritereValid=0;
	for(i=0;i<lNProposant;i++)
	{
		NoeudPro[i]->bFlagSort=0; //0: non satisfait
		if (NoeudPro[i]->etat==GENPROPOSANTINUTILE)
		{
			NoeudPro[i]->prob[0]=1.0; //Si c'est un proposant inutile fait une correction
			if (NoeudPro[i]->iind==0) 
			{
				NoeudPro[i]->bFlagSort=1;
				++nbCritereValid;
			}
		}
	}
	
	//Verifie que la durée d'execution ne sera pas trop intue
//	if (NOrdre > PROB_MAXIMUM_NORDRE){
//		GENError("Execution time is too great to launch this function call");
//		throw std::exception();
		//GENError("La durée d'exécution de cet appel de fonction serait surement comparable à l'âge de l'univers");
//	}
	//PROGRESS BAR	
	//Construction du vecteur de valeur
/*	const int PBar_position = MAX(NOrdre - PBARINTERVAL_PROB,0);	
	const double maximumiteration = pow(3,PBar_position+1);

	PrCumul[0]=0.0;
	PrCumul+=1; //Pour simplifier
	for(i=0;i<=PBar_position;i++)	
	{			
		PrCumul[i]=0;
		//Calcul de valeur
		const double val= pow(3,PBar_position-i);
		PrValue[i][0]=0.0;
		PrValue[i][1]=val;
		PrValue[i][2]=val+val;
	}
	CREATE_FLOAT_PROGRESS_BAR(maximumiteration,&PrCumul[PBar_position],printprogress)
*/

#ifdef USESDEBUG
	//Informatin de débuggage
	SXLONG iter =0;
	//printf("\nThe order is = %d\n",NOrdre);
	//printf("Maximum number of iterations used:%.15G\n",pow(3,NOrdre)/2);
	//printf("\nTaskbar order  = %d\n",PBar_position);
	//printf("Progress bar iteration number:%.15G\n",maximumiteration);
#endif

	while(n>=0)
	{
			//n niveau de l'ordre qui est en traitement
			//a nombre d'allele du noeud courant

			//pour économiser du code
			CIndSimul& nd = *Ordre[n];
			int &a = nd.allele;

			//test
			iteration++;
			++a;			
			if (a==3)
			{
				a=-1;
				n--;
				/*if (n==PBar_position){ INCREMENT_FLOAT_PROGRESS_BAR();	}*/
				#ifdef USESDEBUG
					++iter;
				#endif
			}
			else {
				//Calcul pour la progressbar
				//if (n<=PBar_position) PrCumul[n]=PrCumul[n-1]+PrValue[n][a];					

				//Trouve le nombre d'allele des parents
				if (nd.pere!=NULL)	ap = nd.pere->allele;
				else				ap = 0;
				
				if (nd.mere!=NULL)	am=	nd.mere->allele;
				else				am=0;
				
				//ebranchage (C'est ce qui empeche de calcule atteint)
				if (TransGen[ap][am][a]!=0.0) {
					//Probabilité d'être dans l'état courant					
					Cumul[n+1]=Cumul[n]*TransGen[ap][am][a];
					
					//Si on désirer la probabilié individuel décoché en dessous		
					if (nd.etat==GENPROPOSANT){
						nd.prob[a]+=Cumul[n+1];							

						//Evaluation des criteres pour probabilité conjointe
						if (nd.bFlagSort==1) {
							--nbCritereValid;
							nd.bFlagSort=0;
						}
						if ( ((nd.iind==0 && a==0) || (nd.iind!=0 && a>=nd.iind)))
						{
							++nbCritereValid;
							nd.bFlagSort=1;							
						}						

						//Si tous les proposants sont considéré comme valide alors...
						//Compet dans la propabilité conjointe
						if (n==Lastindex && nbCritereValid==lNProposant)												
							ConjProb+=Cumul[n+1];

						//Passe à la sous-étape suivante	
						if (n!=Lastindex && (!onlyConj || (onlyConj && nd.bFlagSort==1) ) )	++n;
					}//fin si proposant
					else
					{
						//Passe à la sous-étape suivante	
						if (n!=Lastindex)	++n;
					}//fin else proposant
				} //Fin ébranchage
			}//Fin de la boucle pour iteration valide
	}//Fin de la boucle principale
	//END_FLOAT_PROGRESS_BAR();

#ifdef USESDEBUG	
	//printf("\nprogressbar value: %.15G\n",PrCumul[PBar_position]);
	//printf("Number of iterations executed:%I64d\n",iter);
#endif
	
	//CALCUL DE LA PROBABILITE CONJOINTE
	*pdRetConj=ConjProb;
		
	//CALCUL DE LA PROBABILITE INDIVIDUEL (OBSOLETE) OPTIONNELLE
	for(i=0;i<lNProposant;i++)
	{
		switch(plProEtat[i])
		{
		case 0:
			pdRetSimul[i]=NoeudPro[i]->prob[0];break; 
		case 1:
			pdRetSimul[i]=NoeudPro[i]->prob[1]+NoeudPro[i]->prob[2];break;
		case 2:
			pdRetSimul[i]=NoeudPro[i]->prob[2];break;
		}
	}	
	return Rcpp::wrap(0);
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			}catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


// *************************************************************************
//	APPARENTEMENT MODIFIER 
// ************************************************************************* 

//Facteur de charge maximale de la table de hashage utilise par CoefApparentement
/**
	si ex= 0.8 alors la table de hashage sera 20% plus grande que le nombre maximale de possibilite
	Avec ce facteur on peut se permettre d'echange memoire/performance..

  Mais dans tous les cas COAPPMOD_HASHAGECHARGEMAXIMAL ne devrais pas etre egale a 1 sinon sa vas etre extremement lent
 */
const double COAPPMOD_HASHAGECHARGEMAXIMAL=.8;

///Mémoire utilisable au grand maximum en octet;
const long double COAPPMOD_MAXMEMORY_USABLE=(long double) 4000*MEGAOCTET; //800*MEGAOCTET;
///Nombre maximal de possibilité autorisé par l'algorithme...
//COAPPMOD_NBPOSSMAX_DUPP doit-être plus petit qu'un unsigned long (d'une bonne marge)
const long double COAPPMOD_NBPOSSMAX_DUPP=(long double)ULONG_MAX*COAPPMOD_HASHAGECHARGEMAXIMAL*.7;  //.7 = Facteur sécurité
//COAPPMOD_NBPOSSMAX_DUPP doit-être plus petit qu'un XLONG (pour la progresse bar et autre)
const long double COAPPMOD_NBPOSSMAX_NODUPP=(long double)1E16;  //SERIEUSEMENT LONG SUR QUOI QUE CE SOIT  ( avant 1E12 ) **JFL**


/*! 
	\brief Calcul de l'apparentement modifier entre un ancetre et quelques proposants

	Cette fonction calcule l'apparentement modifier d'entre un ancêtre et n proposant. 
	Plusieurs calcul à partir d'ancêtre différent vers le même nombre n de proposant peuvent-être fait lors du même appel de fonction.

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plInput	 [in] Vecteur representant une matrice d'entier (Selon Formule: col * NLigne + ligne)
						Matrice d'entier de n colonne,  sur chaque ligne, on retrouve n proposant à étudier p/r à un ancêtre. 
						<br>Chaque ligne représente un calcul complètement distinct
						<br>Le nombre de ligne de la matrice doit être égale au nombre d'élément du vecteur ancêtre pour faire la correspondance.
						<br>Ex :  
							<br>ancetre =  Vecteur 10, 20
							<br>
							<table>
								<tr>
								<td><B>10</B></td>
								<td>1</td>
								<td>2</td>
								<td>3</td>
								</tr>
								<tr>
								<td><B>20</B></td>
								<td>4</td>
								<td>5</td>
								<td>6</td>
								</tr>
							</table>
						<br>Dans ce cas, la fonction calculera.
							<br>&nbsp;&nbsp;1- 	L'apparentement modifié des proposant 1, 2, 3 avec l'ancêtre 10
							<br>&nbsp;&nbsp;2- 	Et l'apparentement modifié des proposant 4, 5, 6 avec l'ancêtre 20
						<br><br>La valeur de retour sera un vecteur avec les deux résultats précédents.

	\param lNColonne [in] Nombre de Colone de la matrice (Incluant la 1e colonne des ancetres)
	\param lNLigne	 [in] Nombre de ligne de la matrice
  

	\retval pdRetVecAppMod	[out] Un pointeur vers un vecteur de double de taillel lNLigne
							En cas de success, contient la valeur de l'apparentement modifier pour chaque ligne de la matrice
		 		
	\param DuppDetection	[in] Si !=0 la detection de dupplicata dans les chemins sera activer

	\param MaxMoTableHash	[in] Nombre max de Mo que la table de hash peut utiliser avant de declenche une erreur
	  
	\param Maxcombinaison	[in] Nombre maximum de combinaision a tester avant d'imprimer un message d'erreur

	\param printprogress	[in] imprime un message indiquant les progress accomplies

	\return 0 si la fonction est executé avec succès 
	
	 \sa CoPro
*/ 
int CoefApparentement(int* Genealogie,	int* plProposant, int NProposant, int* plAncetre, double* pdRetour,int DuppDetection, int printprogress)
{
	try{
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,NProposant,&NoeudPro);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **TmpNoeudAnc=NULL;	
	LoadAncetre(plAncetre,1,&TmpNoeudAnc);
	CIndSimul *NoeudAnc=*TmpNoeudAnc;	

	//Creation tableau
    INITGESTIONMEMOIRE;

	//OPTIONNEL
	CIndSimul** Ordre		=(CIndSimul**) 	memalloc(lNIndividu,sizeof(CIndSimul*));		
	CApPath**  Path			=(CApPath**) 	memalloc(NProposant,sizeof(CApPath*));
	CApPath **Current 		=(CApPath**) 	memalloc(NProposant,sizeof(CApPath*));	
	int i;	
	
	//initialisation de la structure noeud
	for(i=0;i<lNIndividu;i++)
	{		
		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;	   		
		Noeud[i].iind=0;
		
		Ordre[i]=NULL;
	}
			
	//#2 : CREATION DU BITFIELD OU iind EST LE NO DU BITS	
	//initialisation
	for(i=0;i<NProposant;i++)
	{
 		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
		Current[i]=NULL;
		Path[i]=NULL;		
	}
	NoeudAnc->etat=GENDEPART;

	//Trouve les noeuds utile
	ExploreArbre(NoeudAnc);

	int cbit=-1;
	for(i=0;i<lNIndividu;i++)
	{
		if (Noeud[i].etat!=GENNONEXPLORER && Noeud[i].etat!=GENINUTILE)
		{
			Noeud[i].iind=++cbit;
		}
	}
	//Nombre de noeud touché, taille de la chaine de char 8 bit
	//const long NOrdre=cbit;
	const int taille=(int) ceil(double(cbit)/MP_DIGIT_BIT); //MPI
	
	//AJUSTE LA PRECISION DU MPI
	mp_set_prec(taille); //Indique a mpi de faire des nombres de la bonne taille....
	
	//RECHERCHE DE TOUS LES CHEMINS POSSIBLES ENTRE PROPOSANT ET ANCETRE
	g_ExpCoeff_CheminParcouru=Ordre; 		
	for(int cCol=0;cCol<NProposant;++cCol)
	{						
		Path[cCol]=NULL;
		g_ExpCoeff_Path=&Path[cCol];  
		g_ExpCoeff_Cible=NoeudPro[cCol];

#ifdef USESDEBUG
		//printf("\nPath finder %d -> %d\n",NoeudAnc->nom,NoeudPro[cCol]->nom);
#endif
		ExploreCoeff(NoeudAnc);				
	}

	//EVALUATION DU NOMBRE DE POSSIBILITE
	//long double dnbposs=1;
	double dnbposs=1;
	//unsigned long NbIteration = 1;
	unsigned int NbIteration = 1;
	for(int cCol=0;cCol<NProposant;++cCol)
	{						
		int tmpcountpath=0;			
		CApPath* tmp=Path[cCol];
		while (tmp!=NULL)
		{
			++tmpcountpath;
			tmp=tmp->next;
		}
		dnbposs *= tmpcountpath;
		NbIteration *= tmpcountpath; //Peut débordé
	}
	
	//VALIDE SI LE NOMBRE DE COMBINAISON EST TROP GRAND
	const long double MaxMemoryUsed = dnbposs*(taille*sizeof(mp_digit)+sizeof(mp_int))/COAPPMOD_HASHAGECHARGEMAXIMAL; //En octet
	if (DuppDetection)
	{		
		if (dnbposs>COAPPMOD_NBPOSSMAX_DUPP)
		{
			PathDestruction(Path,NProposant);
//			GENError("Number of combination to evaluate is too great for duplicata detection\nDeactivate it if you want to continue.");
			throw std::range_error("Number of combination to evaluate is too great for duplicata detection\nDeactivate it if you want to continue.");
			//GENError("Le nombre de combinaison a évaluer est trop grand pour pouvoir utilisé la détection de dupplicata\n Déactivé la détection de dupplication si vous désiré continuer"); 
		}

		if (MaxMemoryUsed>COAPPMOD_MAXMEMORY_USABLE)
		{
			PathDestruction(Path,NProposant);
//			GENError("Memory usage is too great for duplicata detection\nDeactivate it if you want to continue\n"
//				    "Maximum memory allowed: %lG Mo  memory needed: %lG Mo\n\n", 
//					COAPPMOD_MAXMEMORY_USABLE/MEGAOCTET,MaxMemoryUsed/MEGAOCTET);
			char erreur[TAILLEDESCRIPTION];
			
			double maxMemAllowed = COAPPMOD_MAXMEMORY_USABLE/MEGAOCTET;
			double memUsed = MaxMemoryUsed/MEGAOCTET;
			sprintf(erreur, "Memory usage is too great for duplicata detection\nDeactivate it if you want to continue\nMaximum memory allowed: %f Mo  memory needed: %f Mo\n\n",
				maxMemAllowed, memUsed);
			//sprintf(erreur, "Memory usage is too great for duplicata detection\nDeactivate it if you want to continue\nMaximum memory allowed: %Lf Mo  memory needed: %Lf Mo\n\n",
				// (COAPPMOD_MAXMEMORY_USABLE/MEGAOCTET),(MaxMemoryUsed/MEGAOCTET));
			throw std::range_error(erreur);
			//GENError("La quantite de mémoire utilisé par la détection de dupplicata est trop importante\n Déactivé la détection de dupplication si vous désiré continuer"
					 //"\nTaille maximal permise : %lG Mo   Mémoire demandé : %lG Mo\n\n",	 COAPPMOD_MAXMEMORY_USABLE/MEGAOCTET,MaxMemoryUsed/MEGAOCTET);
		}
	}
	else
	{
		if (dnbposs>COAPPMOD_NBPOSSMAX_NODUPP)
		{
			PathDestruction(Path,NProposant);
//			GENError("Execution time of this function call is too great");
			throw std::range_error("Execution time of this function call is too great");
			//GENError("La durée d'exécution de cet appel de fonction serait surement comparable à l'âge de l'univers");
		}
	}

	//validation de solution triviale
	*pdRetour=0.0;//Initialise a zero la valeur de retour
	if (NbIteration==0)
	{
		PathDestruction(Path,NProposant);
		return 0;  //Operation réussi mais le résultat est simplement zero
	}

	//CREATION DE LA TABLE ANTI DUPPLICATA	
	//initrand();
	HashDouble<mp_int>* ptrHashTable=NULL;

	if (DuppDetection)
		ptrHashTable= new HashDouble<mp_int>(NbIteration,COAPPMOD_HASHAGECHARGEMAXIMAL);
	if (ptrHashTable==NULL || ptrHashTable->Error())
	{
		//PathDestruction(Path,NProposant);
//		GENError( "Insufficient memory to create an anti-duplicata table\n Deactivate duplicata detection if you want to continue\n"
//				"Memory needed: %lG Mo\n", MaxMemoryUsed/MEGAOCTET);
		char erreur[TAILLEDESCRIPTION];
		double memUsed = MaxMemoryUsed/MEGAOCTET;
		sprintf(erreur, "Insufficient memory to create an anti-duplicata table\n Deactivate duplicata detection if you want to continue\nMemory needed: %f Mo\n",
			memUsed);
		//sprintf(erreur, "Insufficient memory to create an anti-duplicata table\n Deactivate duplicata detection if you want to continue\nMemory needed: %Lf Mo\n",
		//	 (MaxMemoryUsed/MEGAOCTET));
		throw std::range_error(erreur);
		//GENError("Mémoire insuffisante pour utilisé créer une table anti-dupplicata\n Déactivé la détection de dupplication si vous désiré continuer"
				 //"\nMémoire nécessaire : %lG Mo\n",MaxMemoryUsed/MEGAOCTET);
	}

	//SOMMATION DE LA PROBABILITE CONJOINTE DE CHAQUE CHEMIN POSSIBLE
	//l'équivalent d'un ensemble de nbproposant boucle imbriqué.
	//Chaque boucle représente l'ensemble des chemins partant de l'ancetre
	//et ce rendant au proposant de cette boucle
	
	//Initialisation
	for(i=0;i<NProposant;++i)
		Current[i]=NULL;; //INITIALISE LE COMPTEUR DE BOUCLE

	int n=0; //Curseur qui indique sur que chemin on est actuellement
		
	mp_int tmpResult; //La solution courante à étudier
	mp_int tmp; 
	mp_init(&tmpResult);
	mp_init(&tmp);
	
	const int LastProposant = NProposant-1;	
	//CREATE_PROGRESS_BAR(NbIteration,printprogress);	
	while(n>=0)
	{
		//AVANCE Aincremente
		if (Current[n]==NULL) Current[n]=Path[n];
		else				  Current[n]=Current[n]->next;
						
		if (Current[n]!=NULL) {
			//EST-CE LA DERNIERE BOUCLE
			if (n==LastProposant) {
				//OUI 
				
				//INCREMENT_PROGRESS_BAR();

				//genere la solution & la distance																						
				//Pour chaque proposant
				mp_zero(&tmpResult); //Remise a zero
				for(int cCol=0;cCol<NProposant;++cCol) 
				{						
					//OR pour combiner le tout
					mpl_or(&tmpResult, &(Current[cCol]->num), &tmp);
					//Remet le résultat dans tmpResult					
					mp_exch(&tmpResult, &tmp);									
				}
				
				//Validation Anti-dupplicata											
				if (!DuppDetection ||
					 (DuppDetection && ptrHashTable->Add(tmpResult)==GTRUE) )
				{
					//Calcul et ajout de la distance au résultat		
					int distance=0;
					mpl_num_set(&tmpResult, &distance);  
					*pdRetour += pow2(distance-1); //Me rappelle plus pourquoi -1... 
				}
			}//Derniere boucle: NON : if (n==(NProposant-1))
			else
				++n;
		} //if (Current[n]!=NULL)
		else
			--n;		
	}//fin while(n>=0)	

	//VALEUR DE RETOUR
	//*pdRetour

	//NETTOYER LES ASSIGNEMENTS DE MÉMOIRES			
	//outrand();
	PathDestruction(Path,NProposant);

	mp_clear(&tmpResult);
	mp_clear(&tmp);
	if (ptrHashTable)
		delete ptrHashTable;	
	
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			}catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


// ********************************************************************
//
//			FONCTION PRIVE
//
// ********************************************************************


/*! 
	\brief (INTERNE) Detruit un tableau de liste CApPath

		Pour chaque element du vecteur Path. detruit la liste (et tout le element contenu dans celle-ci)

	\param Path		[in] Ptr vers un vecteur de pointeur vers des liste CApPath
	\param npath	[in] Taille du vecteur Path
	
*/
static void PathDestruction(CApPath **Path,int npath)
{
	//DESTRUCTION DES PATHS (A VERIFIER)
	CApPath *Current=NULL;	
	CApPath *tmp	=NULL;
	for(int j=0;j<npath;j++)
	{
		Current=Path[j];
		while(Current!=NULL)
		{
			
			tmp=Current->next;
			mp_clear(&Current->num); 
			memfreeIN(Current); 
			Current=tmp;								
		}
	}
}


///Utilise par ExploreCoeff pour représenté le chemin parcouru
/*
static CIndSimul** g_ExpCoeff_CheminParcouru=NULL; 
///Utilise par ExploreCoeff comme étant le dernier chemin remplis (ou le premiers)
static CApPath ** g_ExpCoeff_Path=NULL;
///Utilise par ExploreCoeff comme étant la cible de l'exploration
static CIndSimul* g_ExpCoeff_Cible=NULL;*/
/*! 
	\brief (INTERNE) Explore une genealogie et construit une liste CapPath
		
		Trouve tous les chemins entre le Noeud et la cible et sa les inscrit dans une liste CApPath
		sous forme binaire mpi_int.

		La forme binaire correspond a la somme de tout les 2<<iind des noeuds implique et peut s'etendre sur un nombre
		infini d'octet

	\attention Il faut que les noeuds soit correctement etiquette (etat)
				De plus, il faut qu'un iind unique soit assigne a chaque noeud qui sont Utile
				???? important avant mp_set_prec(taille);
				g_ExpCoeff_CheminParcouru
				g_ExpCoeff_InputPath
				g_ExpCoeff_Cible
	
	\param Noeud		[in] Noeud courant ou l'exploration est rendu (Pts de depart la majorite du temps)
	\param profondeur	[in] Profondeur actuel (au depart:0 normalement)
	\param Cible		[in] Noeud Cible 
	\param Array2		[in] Pointeur vers un vecteur de cindsimul* de taille nombre de noeud max implique 
							 (utilisé pour mémorisé le chemin utilisé)
	\param InputPath	[in] Taille du vecteur Path
	\param taille2		[in] Taille du vecteur Path
	
	  \remark Cette fonction est recursive
*/
static void FASTCALL ExploreCoeff(CIndSimul* Noeud)
{
	//Explore l'arbre et retourne tous les path utilisé		
	static int profondeur=0;

	//Trace le chemin parcouru
	g_ExpCoeff_CheminParcouru[profondeur]=Noeud;

	//on vient d'atteindre la cible?
	if (Noeud==g_ExpCoeff_Cible)
	{		
		//Creation d'un nouveau chemin
		CApPath *tmp=(CApPath*) memallocIN(1,sizeof(CApPath));
		
		//Initialise le nombre contenu....		
		mp_init( &tmp->num );		
		tmp->next=NULL;

		//Complete la série
		*g_ExpCoeff_Path=tmp;
		g_ExpCoeff_Path=&(tmp->next);	//Avance le curseur
		
		//Ajuste le nombre en conséquence....
		for(int i=0;i<=profondeur;i++)
		{
#ifdef USESDEBUG
			//printf("%d,",g_ExpCoeff_CheminParcouru[i]->nom);
#endif
			mpl_bit_set(&tmp->num, g_ExpCoeff_CheminParcouru[i]->iind);
		}
#ifdef USESDEBUG
		unsigned char buffer[100];
		mp_toradix(&tmp->num, buffer, 16);
		//printf(" == %s\n",buffer);
#endif

		return;	
	}	
	else
	{
		//Non.. on continue a cherché
		Clist *current=Noeud->fils;
		if (current!=NULL)
		{
			do
			{	
				CIndSimul *tmp=current->noeud;				
				if (tmp->etat!=GENNONEXPLORER && tmp->etat!=GENINUTILE)
				{
					++profondeur;
					ExploreCoeff(tmp);
					--profondeur;
				}
				//Fils suivant
				current=current->next;				
			}
			while(current!=NULL);
		}
	}	
	return;	
}

