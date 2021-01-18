//#define MODETEST

//contributor Jean-François Lefebvre

/// Authorise l'affichage du niveau de progression sur la sortie standard
/** 
	Si ALLOWPRINTPROGRESS est defini, les fonctions suivantes peuvent 
	indiquer leur niveau de progression sur la sortie standard stdout 
*/ 
#define ALLOWPRINTPROGRESS

#include "base.h"
#include "outils.h"
#include "apparentement.h"
#include "basemt.h"
#include "userInterface.h"

//Temporaire
#include <fstream>
// #include "windows.h"
#include <cstdlib>
#include <math.h>

#include <limits.h>
#include <Rcpp.h>

using namespace std;

#define R_NO_REMAP

#define DIV 1024
string divisor = "K";
#define WIDTH 7
/*! 
	\brief Donne la matrice Phi pour une profondeur fixe

	La fonction calcule la valeur de l'apparentement entre chaque proposant qui lui est fourni. 
	Le resultat retourne sous la forme d'une matrice NProposant x NProposant   

	\param Genealogie		[in] Une genealogie construite à l'aide de gen.genealogie 

	\param proposant		[in] Vecteur avec les numeros des proposants à étudier 
	\param NProposant		[in] Nombre d'éléments du vecteur proposant 
  
	\param Niveau			[in] Profondeur fixe 

	\retval pdMatricePhi	[out] Adresse d'un tableau de NProposant x NProposant. 
							En cas de Succes, pdMatricePhi contient la matrice Phi et 
		  					les differents colonnes et lignes corresponde au proposant fourni et ce dans le meme ordre 
							<br>Ex: <br> 
								soit:<br> 
									Ind1 : Indice du proposant 1 dans la liste proposant<br> 
									Ind2 : Indice du proposant 2 dans la liste proposant<br> 
							  Apparentement entre proposant 1 et proposant 2 <br> 
								= pdMatricePhi[ Ind1 x NProposant + Ind2 ]<br> 
									ou<br> 
								  pdMatricePhi[ Ind2 x NProposant + Ind1 ]<br> 
						
	\param printprogress	[in] Imprime une serie de message indiquant la progression du calcul 
 
	\return 0 si la fonction est executé avec succès 
*/ 
int PhiMatrix(int* Genealogie, int* proposant, int NProposant,int Niveau, double* pdMatricePhi, int printprogress)
{
	try{
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GFALSE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(proposant,NProposant,&NoeudPro);

	//DEBUT DU CALCUL
	if (Niveau==0)
		Niveau=SHRT_MAX;

	if (Niveau>SHRT_MAX){
//		GENError("Niveau must be smaller than %d", SHRT_MAX);
		char erreur[TAILLEDESCRIPTION];
		sprintf(erreur, "depthmin must be smaller than %d",SHRT_MAX);
		throw std::range_error(erreur);
		//GENError("Le niveau doit-être inférieur à %d",SHRT_MAX);
	}
	const short niveauMax = short(Niveau);

	//Barre de progression
	//CREATE_PROGRESS_BAR_MATRIX(NProposant,printprogress)
	for(int cPro1=0;cPro1<NProposant;++cPro1)
	{
		for(int cPro2=cPro1;cPro2<NProposant;++cPro2)
		{
			/*
			if (cPro1==cPro2)
			{
				// L'apparentement de qqun avec lui-meme = 0.5
				pdMatricePhi[cPro1*NProposant+cPro2]=0.5;
				pdMatricePhi[cPro2*NProposant+cPro1]=0.5;
			}
			else
			{*/						
				const double tmp=Kinship(NoeudPro[cPro1],NoeudPro[cPro2],niveauMax,niveauMax);								
				pdMatricePhi[cPro1*NProposant+cPro2]=tmp;
				pdMatricePhi[cPro2*NProposant+cPro1]=tmp;
				
				//Affichage des progress
				//INCREMENT_PROGRESS_BAR();
			//}//Fin IF
		}//Fin itérateur proposant 2
	}//Fin itérateur proposant 1
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


/* VERSION MT */
struct CBASEMTPhiMatrixMT : public CMtGlobalMessage
{
	CIndSimul* ind1;
	CIndSimul* ind2;
	double danswer;
	int indice1;			//Indice ou placer la reponse dans le tableau de reponse
	int indice2;			
	short niv;				//Nombre de remonte de la matrice
};

BASEMT_CREATE_GLOBALMESSAGE(CBASEMTPhiMatrixMT,1)

BASEMT_DEBUT_HELPERFCT(CBASEMTPhiMatrixMT,1)
		//LANCEMENT DU CALCUL ET RECUPERE LE RESULTAT
		BASEMT_HLPMES.danswer  = 
			Kinship(BASEMT_HLPMES.ind1,BASEMT_HLPMES.ind2,BASEMT_HLPMES.niv,BASEMT_HLPMES.niv); 
BASEMT_FIN_HELPERFCT() 


int PhiMatrixMT(int* Genealogie, int* proposant, int NProposant,int Niveau, double* pdMatricePhi, int printprogress)
{
	try{
	//CREATION DE TABLEAU D'INDIVIDU
	INITGESTIONMEMOIRE
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GFALSE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(proposant,NProposant,&NoeudPro);

	//DEBUT DU CALCUL
	if (Niveau==0)
		Niveau=SHRT_MAX;	

	if (Niveau>SHRT_MAX){
//		GENError("Niveau must be smaller than %d", SHRT_MAX);
		char erreur[TAILLEDESCRIPTION];
		sprintf(erreur, "depthmin must be smaller than %d",SHRT_MAX);
		throw std::range_error(erreur);
		//GENError("Le niveau doit-être inférieur à %d",SHRT_MAX);
	}
	const short niveauMax = short(Niveau);

	//Initialisation d'une structure multithread
	BASEMT_DEBUTBOUCLE_INITMESSAGE(1) 		
		//initialisation des données Phi de la structure
		BASEMT_MESSAGE(1).niv=niveauMax;	 	  
		BASEMT_MESSAGE(1).indice1 = -1;
		BASEMT_MESSAGE(1).indice2 = -1;	 
	BASEMT_FINBOUCLE_INITMESSAGE(1)

	//Demarre le calcul de phi
	//Barre de progression
	CREATE_PROGRESS_BAR_MATRIX(NProposant,printprogress);		
	for(int cPro1=0;cPro1<NProposant;++cPro1)
	{
		for(int cPro2=cPro1;cPro2<NProposant;++cPro2)
		{
			/*
			if (cPro1==cPro2)
			{
				// L'apparentement de qqun avec lui-meme = 0.5 
				pdMatricePhi[cPro1*NProposant+cPro2]=0.5;
			}
			else
			{*/
				BASEMT_DEBUT_REQUETE(1) 
					//Résupere le résultat de phi
					if (BASEMT_MESSAGE(1).indice1!=-1)
					{
						pdMatricePhi[BASEMT_MESSAGE(1).indice1*NProposant+BASEMT_MESSAGE(1).indice2]=BASEMT_MESSAGE(1).danswer;
						pdMatricePhi[BASEMT_MESSAGE(1).indice2*NProposant+BASEMT_MESSAGE(1).indice1]=BASEMT_MESSAGE(1).danswer;
					}
					//Parametre pour un nouveau calcul de phi
					BASEMT_MESSAGE(1).indice1 = cPro1;
					BASEMT_MESSAGE(1).indice2 = cPro2;
					BASEMT_MESSAGE(1).ind1=NoeudPro[cPro1];
					BASEMT_MESSAGE(1).ind2=NoeudPro[cPro2];				
				BASEMT_FIN_REQUETE(1)

				//Affichage des progress
				INCREMENT_PROGRESS_BAR();
			//}//fin if pour chaque case de la matrice
		}//Fin itérateur proposant 2
	}//Fin itérateur proposant 1
	
	//Fermeture des threads
	BASEMT_DEBUT_FERMETURE(1)
		//RECUPERE LES DERNIERS RESULTATS DE PHI S'IL SONT VALIDE
		if (BASEMT_MESSAGE(1).indice1!=-1)
		{
			pdMatricePhi[BASEMT_MESSAGE(1).indice1*NProposant+BASEMT_MESSAGE(1).indice2]=BASEMT_MESSAGE(1).danswer;
			pdMatricePhi[BASEMT_MESSAGE(1).indice2*NProposant+BASEMT_MESSAGE(1).indice1]=BASEMT_MESSAGE(1).danswer;
		}
	BASEMT_FIN_FERMETURE(1)

	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}

