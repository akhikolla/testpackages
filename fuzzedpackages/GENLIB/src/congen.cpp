/*! \file congen.cc
\brief Implementation des fonctions de calcul de la Contribution Genetique

Calcul et Analyse de diverse valeur dérivé de la contribution génétique

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
#include "apparentement.h"
#include "userInterface.h"

//#include <iostream>
#include <vector>

#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>

#include <cstdlib> //jfl
#include <ctime>   //jfl
#include <stdio.h> //jfl

using namespace std;

// ********************************************************************
//
//			CONSTANTE
//
// ********************************************************************


/// Paire Valeur(double) , name (int)  utilise pour trier en ordre croissant dans CongenCumul
/**
	Utilise pour permettre a qsort de classer des noms selons un ordre numerique
*/
struct PairAncValue
{
	///Etiquette pour cette paire
    int name;
	///Valeur numerique correspondant a ce nom
    double value;
};


static int WINCDECL PairCompare(const void *p1,const void *p2);
static void ExploreConGenProposant(CIndSimul* Noeud, int profondeur);
static void ExploreConGenProposantPLUS(CIndSimul* Noeud, int profondeur, double* pdSexe, vector<double>& vProb);

// ********************************************************************
//
//			PUBLIC
//
// ********************************************************************

/*! 
	\brief Calcul la contribution genetique entre des ancetres et des proposants

	 Calcule la contribution genetique entre une serie d'ancetres et une serie de proposants.
	 Le resultat se presente sous la forme d'une matrice de Nancetre*Nproposant elements.

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant à étudier
	\param lNProposant	[in] Nombre d'élément du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant à étudier
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre

	\retval pdCongen	[out] Un pointeur vers une vecteur de taille  lNAncetre*NProposant pour recevoir résultat
						En cas de succes, on peut allez cherche la contribution genetique de la maniere suivante.
						<br>Ex:
						<br>&nbsp; &nbsp; anc : Indice de l'ancetre
						<br>&nbsp; &nbsp; pro : Indice du proposant
						<br>alors:
						Congen = pdCongen[pro*lNAncetre+anc]

	\param printprogress [in] Imprime un message indiquant les progress accomplies

	\return 0 si la fonction est executé avec succès

	\remark les no d'individu ne sont que des étiquettes, ne sont utilisé qu'en référence au père et mère
*/
int Congen(int* Genealogie, int* plProposant,int lNProposant, int* plAncetre, int lNAncetre, double* pdCongen, int printprogress) 
{
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);
	
	//CREATION D'UN VECTEUR D'ANCETRE
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	//INITIALISE LES DONNES	
	//int i;
	for(int i=0;i<lNProposant;i++)
	{	
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
		NoeudPro[i]->prob[0]=0.0;
	}
	
	//POUR CHAQUE ANCETRE FAIRE LE CALCUL SUIVANT
	CREATE_PROGRESS_BAR(lNAncetre,printprogress)
	for(int cIndAnc=0;cIndAnc<lNAncetre;++cIndAnc) 
	{	     	

		//CALCUL CONTRIBUTION GENETIQUE
		ExploreConGenProposant(NoeudAnc[cIndAnc],0);
		for(int i=0;i<lNProposant;i++)
		{
			pdCongen[cIndAnc*lNProposant+i]=NoeudPro[i]->prob[0];
			NoeudPro[i]->prob[0]=0;
		}

		//BARRE DE PROGRESSION
		INCREMENT_PROGRESS_BAR();

	}//Fin pour chaque ligne

	return 0;
}
// ********************************************************************
//
//			PUBLIC
//
// ********************************************************************

/*! 
	\brief Calcul la contribution genetique entre des ancetres et des proposants

	 Calcule la contribution genetique entre une serie d'ancetres et une serie de proposants.
	 Le resultat se presente sous la forme d'une matrice de Nancetre*Nproposant elements.

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] Vecteur des no de proposant à étudier
	\param lNProposant	[in] Nombre d'élément du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant proposant à étudier
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre
	\param pdSexe		[in] Probabilité de transmission selon le sexe

	\retval pdCongen	[out] Un pointeur vers une vecteur de taille  lNAncetre*NProposant pour recevoir résultat
						En cas de succes, on peut allez cherche la contribution genetique de la maniere suivante.
						<br>Ex:
						<br>&nbsp; &nbsp; anc : Indice de l'ancetre
						<br>&nbsp; &nbsp; pro : Indice du proposant
						<br>alors:
						Congen = pdCongen[pro*lNAncetre+anc]

	\param printprogress [in] Imprime un message indiquant les progress accomplies

	\return 0 si la fonction est executé avec succès

	\remark les no d'individu ne sont que des étiquettes, ne sont utilisé qu'en référence au père et mère
*/
int CongenPLUS(int* Genealogie, int* plProposant,int lNProposant, int* plAncetre, int lNAncetre, double* pdSexe, 
			double* pdCongen,int printprogress) 
{
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);
	
	//CREATION D'UN VECTEUR D'ANCETRE
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	//INITIALISE LES DONNES	
	int i;
	for(i=0;i<lNProposant;i++)
	{	
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
		NoeudPro[i]->prob[0]=0.0;
	}
	
	//POUR CHAQUE ANCETRE FAIRE LE CALCUL SUIVANT
	CREATE_PROGRESS_BAR(lNAncetre,printprogress)
	for(int cIndAnc=0;cIndAnc<lNAncetre;++cIndAnc) 
	{	 
		vector<double> vProb(lNIndividu, 0);
		//CALCUL CONTRIBUTION GENETIQUE
		ExploreConGenProposantPLUS(NoeudAnc[cIndAnc],0, pdSexe, vProb);
		for(i=0;i<lNProposant;i++)
		{
			pdCongen[cIndAnc*lNProposant+i]=NoeudPro[i]->prob[0];
			NoeudPro[i]->prob[0]=0;
		}

		//BARRE DE PROGRESSION
		INCREMENT_PROGRESS_BAR();

	}//Fin pour chaque ligne

	return 0;
}


/*! 
	\brief Liste/sommation cumulative des contributions genetique par ancetre

  Produit une matrice de Nancetre*Nproposant ou la contribution genetique entre chaque ancetre/proposant est calculer   

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] vecteur des no de proposant à étudier 
	\param lNProposant	[in] Nombre d'élément du vecteur proposant 
  
	\param plAncetre	[in] Vecteur des no d'ancetre à étudier
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre 

	\retval AncRet		[out] Un pointeur vers un vecteur de int de taille lNAncetre
						 En cas de success, il contient le no de ancetre une par ordre decroissant de contribution genetique total (voir pdSomAnc) 						
	
	\retval pdSomAnc	[out] Un pointeur vers une vecteur de taille lNAncetre qui contient 
						la somme des contributions genetique de tous les proposants p/r a cet ancetre.
						Celle-ci sont classe en ordre decroissant (No Ancetre en ordre decroissant = AncRet)

	\retval pdSomCumul  [out] Somme cummulative des contributions genetique de pdSomAnc

	\param printprogress imprime un message indiquant les progress accomplies 

	\return 0 si la fonction est executé avec succès 	
*/
int CongenCumul(int* Genealogie,
	int* plProposant,int lNProposant,
	int* plAncetre, int lNAncetre,
    int* AncRet,double* pdSomAnc,double* pdSomCumul, int printprogress)
{
	//VARIABLE OPERATIONNEL
	
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,lNProposant,&NoeudPro);
	
	//CREATION D'UN VECTEUR D'ANCETRE
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);
	
	//Niveau nombre de ligne dans le tableau
	INITGESTIONMEMOIRE;	
	PairAncValue* AncPair	=(PairAncValue*) memalloc(lNAncetre,sizeof(PairAncValue));

	//REMISE EN PLACE DES DONNES DES NOEUDS	
	int i;
	for(i=0;i<lNProposant;i++)
	{
		NoeudPro[i]->etat=GENPROPOSANTINUTILE;
		NoeudPro[i]->prob[0]=0.0;
	}
	
	//POUR CHAQUE ANCETRE
	CREATE_PROGRESS_BAR(lNAncetre,printprogress)
	for(int cIndAnc=0;cIndAnc<lNAncetre;++cIndAnc) 
	{

		//CALCUL CONTRIBUTION GENETIQUE
		ExploreConGenProposant(NoeudAnc[cIndAnc],0);
			
		//CALCUL DE LA VALEUR CUMULATIVE
		AncPair[cIndAnc].name=plAncetre[cIndAnc];
		AncPair[cIndAnc].value=0;
		for(i=0;i<lNProposant;i++)
			AncPair[cIndAnc].value+=NoeudPro[i]->prob[0];

		//BARRE DE PROGRESSION
		INCREMENT_PROGRESS_BAR();

		//REMISE A ZERO
		for(i=0;i<lNProposant;i++)
			NoeudPro[i]->prob[0]=0;
	}//Fin pour chaque ligne

	//Maintenant on trie en ordre ascendant
	qsort(AncPair,lNAncetre,sizeof(PairAncValue),PairCompare);

    //Retourne resultat
	for(int cIndAnc=0;cIndAnc<lNAncetre;++cIndAnc) //pour chaque ancetre
	{
	       AncRet[cIndAnc]=AncPair[cIndAnc].name;
	       pdSomAnc[cIndAnc]=AncPair[cIndAnc].value;
	       if (cIndAnc!=0)
			pdSomCumul[cIndAnc]=pdSomCumul[cIndAnc-1]+AncPair[cIndAnc].value;
	       else
			pdSomCumul[cIndAnc]=AncPair[cIndAnc].value;
	}

	return 0;
}

/*! 
	\brief Liste/sommation cumulative des contributions genetique par ancetre

  Produit une matrice de Nancetre*Nproposant ou la contribution genetique entre chaque ancetre/proposant est calculer   

	\param Genealogie	[in] Une genealogie construite à l'aide de gen.genealogie 

	\param plProposant	[in] vecteur des no de proposant à étudier 
	\param lNProposant	[in] Nombre d'élément du vecteur proposant 
  
	\param plAncetre	[in] Vecteur des no d'ancetre à étudier
	\param lNAncetre	[in] Nombre d'élément du vecteur ancetre 

	\retval AncRet		[out] Un pointeur vers un vecteur de int de taille lNAncetre
						 En cas de success, il contient le no de ancetre une par ordre decroissant de contribution genetique total (voir pdSomAnc) 						
	
	\retval pdSomAnc	[out] Un pointeur vers une vecteur de taille lNAncetre qui contient 
						la somme des contributions genetique de tous les proposants p/r a cet ancetre.
						Celle-ci sont classe en ordre decroissant (No Ancetre en ordre decroissant = AncRet)

	\retval pdSomCumul  [out] Somme cummulative des contributions genetique de pdSomAnc

	\param printprogress imprime un message indiquant les progress accomplies 
	
	\remark marche mal

	\return 0 si la fonction est executé avec succès 	
*/
int CongenCumuldirect(int* matriceCG,
	int lNProposant,
	int* plAncetre, int lNAncetre,
    int* AncRet,double* pdSomAnc,double* pdSomCumul)
{
	//Niveau nombre de ligne dans le tableau
	INITGESTIONMEMOIRE;	
	PairAncValue* AncPair	=(PairAncValue*) memalloc(lNAncetre,sizeof(PairAncValue));
	
	//POUR CHAQUE ANCETRE
	int i;
	for(int cIndAnc=0;cIndAnc<lNAncetre;++cIndAnc) 
	{			
		//CALCUL DE LA VALEUR CUMULATIVE
		AncPair[cIndAnc].name=plAncetre[cIndAnc];
		AncPair[cIndAnc].value=0;
		for(i=0;i<lNProposant;i++)
		{
			AncPair[cIndAnc].value+=matriceCG[cIndAnc*lNProposant+i];					
		}
	}//Fin pour chaque ligne

	//Maintenant on trie en ordre ascendant
	qsort(AncPair,lNAncetre,sizeof(PairAncValue),PairCompare);

    //Retourne resultat
	for(int cIndAnc=0;cIndAnc<lNAncetre;++cIndAnc) //pour chaque ancetre
	{
	       AncRet[cIndAnc]=AncPair[cIndAnc].name;
	       pdSomAnc[cIndAnc]=AncPair[cIndAnc].value;
	       if (cIndAnc!=0)
			pdSomCumul[cIndAnc]=pdSomCumul[cIndAnc-1]+AncPair[cIndAnc].value;
	       else
			pdSomCumul[cIndAnc]=AncPair[cIndAnc].value;
	}
	return 0;
}
//******************************************************************************************
//	Contribution genetique par groupe 
//*****************************************************************************************/

//CongenGroupe existait dans le cvs version 1.1


// ********************************************************************
//
//			FONCTION PRIVE
//
// ********************************************************************

///Fonction de comparaison pour QSORT: Compare des PairAncValue par ordre de Value
static int WINCDECL PairCompare(const void *p1,const void *p2)
{
  
	const PairAncValue *e1=(PairAncValue*) p1;
	const PairAncValue *e2=(PairAncValue*) p2;

	if (e1->value > e2->value)
		return -1;
	else
		if (e1->value < e2->value)
			return 1;

	return 0;
}


/*! 
	\brief  Exploration et calcul contribution genetique

	Fonction recursive utilise par congen pour le calcul de la matrice de contribution genetique
	Sa valeur de la contribution genetique de chaque proposant est stocker dans sa propriete
	Prob[0]

	\attention	Tous les prob[0] devrais etre initialise a zero avant de lancer cette fonction

	\param Noeud		[in] Noeud Fesant partie d'une genealogie valide
	\param profondeur	[in] Profondeur actuel

	\remark Au lancement profondeur devrais etre egale a zero
*/
static void ExploreConGenProposant(CIndSimul* Noeud, int profondeur)
{
	//Calcule la contribution génétique de tout les descendants d'un noeud	
	if (Noeud->etat==GENPROPOSANTINUTILE)
	{
		Noeud->prob[0]+=pow2(profondeur);	
	}

	Clist *current=Noeud->fils;		
	while(current!=NULL)
	{	
		ExploreConGenProposant(current->noeud,profondeur+1);
		current=current->next;				
	}

	return;	
}
/*! 
	\brief  Exploration et calcul contribution genetique

	Fonction recursive utilise par congen pour le calcul de la matrice de contribution genetique
	Sa valeur de la contribution genetique de chaque proposant est stocker dans sa propriete
	Prob[0]

	\attention	Tous les prob[0] devrais etre initialise a zero avant de lancer cette fonction

	\param Noeud		[in] Noeud Fesant partie d'une genealogie valide
	\param profondeur	[in] Profondeur actuel

	\remark Au lancement profondeur devrais etre egale a zero
*/
static void ExploreConGenProposantPLUS(CIndSimul* Noeud, int profondeur, double* pdSexe, vector<double>& vProb)
{
	//Calcule la contribution génétique de tout les descendants d'un noeud	
	if (Noeud->etat==GENPROPOSANTINUTILE)
	{
		double dprob = 0;
		for (int i=0;i<profondeur;i++)
		{
			if (i != 0)
				dprob*=vProb[i];
			else
				dprob=vProb[i];
		}
		if (profondeur == 0)
			Noeud->prob[0]+=pow2(profondeur);
		else
			Noeud->prob[0]+=dprob;
	}

	Clist *current=Noeud->fils;	

	while(current!=NULL)
	{	
		int iprob = 0;

		if (Noeud->sex == GEN_MASC && current->noeud->sex == GEN_MASC) //homme-homme
			iprob = 0;
		else if (Noeud->sex == GEN_MASC && current->noeud->sex == GEN_FEM)//homme-femme
			iprob = 1;
		else if (Noeud->sex == GEN_FEM && current->noeud->sex == GEN_MASC) //femme-homme
			iprob = 2;
		else if (Noeud->sex == GEN_FEM && current->noeud->sex == GEN_FEM) //femme-femme
			iprob = 3;
		
		if (pdSexe[iprob] != 0)
		{
			vProb[profondeur] = pdSexe[iprob];
			ExploreConGenProposantPLUS(current->noeud,profondeur+1, pdSexe, vProb);
		}
		current=current->next;				
	}

	return;	
}
