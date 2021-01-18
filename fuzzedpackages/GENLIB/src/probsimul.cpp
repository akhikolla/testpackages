/*! \file probsimul.cc
\brief Library Genlib: fonction diverse


\author Sébastien Leclerc
\contributor Jean-Francois Lefebvre

*/

/// Authorise l'affichage du niveau de progression sur la sortie standard
/** 
	Si ALLOWPRINTPROGRESS est defini, les fonctions suivantes peuvent
	indiquer leur niveau de progression sur la sortie standard stdout
*/
#define ALLOWPRINTPROGRESS


/* POUR UTILISER LE PROGRAMME EN MODE CONSOLE SUR WINDOWS, IL FAUT DEFINIR MODETEST*/
#include "base.h"
#include "outils.h"
#include "outilanal.h"
#include "userInterface.h"

#include <time.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <R.h>
//#include <Rdefines.h>
#include <Rcpp.h>
#include <Rcpp/as.h>
#include <RcppCommon.h>

#define R_NO_REMAP

typedef std::vector<Tuple> vecTuple_t;

// *********************************************************** ///
//						STRUCTURE ET CLASSE
// *********************************************************** ///

/*FONCTION LOCALE*/
//static int numeroteInd(int nbT, vecTuple_t& vTuple);
//static int genereTuple(int nogen, vecTuple_t& vTuple, CIndSimul** tabgen);

//******************************************
//	OUTILS
//****************************************** 

/*! 
	\brief Compte le nombre d'enfant d'une serie de parent

	\param Genealogie		[in] Une genealogie construite à l'aide de gen.genealogie

	\param plProposant		[in] Vecteur avec les numeros des proposants à étudier
	\param NProposant		[in] Nombre d'éléments du vecteur proposant
  	\retval retour			[out] Adressse d'un vecteur de taille NProposant.
							En cas de Succes, contient le nombre d'enfant associer a chaque proposant
	\return 0 si la fonction est executé avec succès
*/
int CountChild(int* Genealogie, int* plProposant,int NProposant, int* retour)
{
	
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud); //Avec enfants

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,NProposant,&NoeudPro);	

	for(int i=0;i<NProposant;i++)
	{		
		//Ce Noeud est un proposant ou un depart il est forcément utile mais il faut explorer les enfants
		Clist* current=NoeudPro[i]->fils;
		int compteur=0;
		if (current!=NULL)
		{
			do
			{
				++compteur;							
				current=current->next;				
			}
			while(current!=NULL);
		}
		retour[i]=compteur;		
	}

	return 0;
} 
/*! 
	\brief Creer une nouvelle genealogie ne contenant que les noeuds ayant un entre au moins un proposant et un ancetre
	
	<br>Si le 1e proposant !=0.
	<br> La nouvelle genealogie ne contient que le individu qui font parti d'au moins un chemin
		  entre un des individu et un des proposant.
	<br><br>Si le 1e proposant ==0 
	<br> alors la nouvelle genealogie ne contient que les descendants de tous les ancetres


	\param Genealogie		[in] Une genealogie construite à l'aide de gen.genealogie

	\param plProposant	[in] Vecteur des no de proposant à étudier
							 Si le premier proposant est 0, alors tous les enfants
							 d'un des ancetre seront conserver
	\param lNProposant	[in] Nombre d'élément du vecteur proposant
  
	\param plAncetre	[in] Vecteur des no des ancetres correspondant au proposant à étudier
	\param lNAncetre	[in] Nombre d'élément du vecteur plAncetre


	\retval NouvelGenealogie	[out] Adresse d'un vecteur de int de taille de Genealogie
							En cas de Succes, contient la nouvelle genealogie

	\retval tailleNouvelGenealogie	[out] Adresse d'un int
							En cas de Succes, Contient la taille de la nouvelle genealogie
							ce qui permet de tronquer NouvelGenealogie

	\return 0 si la fonction est executé avec succès
*/

int ebranche(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre, int* NouvelGenealogie, 
		   int* tailleNouvelGenealogie)
{ 
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud); //Avec enfants

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	int allProposant=0;
	if (plProposant[0]==0) //Conserve tous les descendants des ancetres
		allProposant=1;
	else
		LoadProposant(plProposant,lNProposant,&NoeudPro);	
	//CREATION D'UN VECTEUR D'ANCETRE
	CIndSimul **NoeudAnc=NULL;
	LoadAncetre(plAncetre,lNAncetre,&NoeudAnc);

	//VARIABLE OPERATIONNEL 	
	INITGESTIONMEMOIRE;	
	int* plRetIndividu	=(int*)memalloc(lNIndividu,sizeof(int));	
	int* plRetPere		=(int*)memalloc(lNIndividu,sizeof(int));
	int* plRetMere		=(int*)memalloc(lNIndividu,sizeof(int));
	int* plRetSexe		= NULL;		
	if (LoadNIndMasc()>=0)
		plRetSexe		=(int*)memalloc(lNIndividu,sizeof(int));

	//RÉINITIALISE LES CHAMPS UTILE DE LA GENEALOGIE
	for(int i=0;i<lNIndividu;i++)
	{		
		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].bFlagSort=0;	   		
	}

	//ETIQUETER LES PROPOSANT 
	if (!allProposant) //Pas tous les proposants qui sont valide alors
	{
		for(int i=0;i<lNProposant;i++) 
			NoeudPro[i]->etat=GENPROPOSANTINUTILE;					 
	}
	
	//IDENTIFIER ET ETIQUETER LES POINTS DE DEPARTS
	for(int i=0;i<lNAncetre;i++) 	
		NoeudAnc[i]->etat=GENDEPART; 

	for(int i=0;i<lNAncetre;i++) 				
		if (allProposant)
			ExploreArbreTousDescendant(NoeudAnc[i]);				
		else
			ExploreArbre(NoeudAnc[i]);			

	//CREATION D'UN ORDRE D'EXECUTION ET CALCUL DES SAUTS 
	int countInd=0;
	int NombreEnfant=0;
	for(int i=0;i<lNIndividu;i++)  
	{ 
		const typenoeud_t et=Noeud[i].etat;
		if(et==GENNOEUD || et==GENDEPART || et==GENPROPOSANT || et==GENPROPOSANTINUTILE) 
		{
			plRetIndividu[countInd]=Noeud[i].nom; 
			if(plRetSexe)
				plRetSexe[countInd]=Noeud[i].sex; 
			if(Noeud[i].pere!=NULL) 
			{ 
				const typenoeud_t petat = Noeud[i].pere->etat; 
				if(petat==GENNOEUD || petat==GENDEPART || petat==GENPROPOSANT || petat==GENPROPOSANTINUTILE) 
				{
					plRetPere[countInd]=Noeud[i].pere->nom;++NombreEnfant;
				}
				else 
					plRetPere[countInd]=0; 
			}
			else 
				plRetPere[countInd]=0;

			if(Noeud[i].mere!=NULL) 
			{ 
				const typenoeud_t metat = Noeud[i].mere->etat; 
				if (metat==GENNOEUD || metat==GENDEPART || metat==GENPROPOSANT || metat==GENPROPOSANTINUTILE) 
				{
					plRetMere[countInd]=Noeud[i].mere->nom;++NombreEnfant;
				}
				else 
					plRetMere[countInd]=0; 
			}
			else 
				plRetMere[countInd]=0; 
			++countInd;				
		}
	}

	//CREATION DE LA GENEALOGIE
	*tailleNouvelGenealogie = TAILLEGENVERSION7(countInd,NombreEnfant);

	//#Utilise CompleteGenealogie... je me pose sérieusement la question...
	CreerGenealogie(plRetIndividu,plRetPere,plRetMere,plRetSexe,countInd,NouvelGenealogie);
	
	//FIN
	return 0; 
}

/*! 
	\brief Traite une partie de l'algorithme

	Pour chaque génération, traite l'ensemble des individus afin de modifier
	le vecteur de Tuple associe à la generation suivante correctement

	\param nogen	[in]	numero de la generation consideree

	\param vTuple	[in, out]	Vecteur STL des tuples

	\param tabgen	[in, out]	vecteur des individus tries par generations

	\return 0 si la fonction est execute avec succes

	\remark ATTENTION : Cette fonction modifie le vecteur de Tuples, mais également
	ptr[0], le tuple des parents des individus de cette génération
*/
/*  PAS Utilisée ..??
static int genereTuple(int nogen, vecTuple_t& vTuple, CIndSimul** tabgen)
{
	CIndSimul* indcour= tabgen[nogen];
	int taille = 0; //nombre réel de tuples dans vTuple

	while (indcour!=NULL) //traiter chaque individu de la liste
	{
		int decale = 1;
		CIndSimul* suivant = (CIndSimul*) indcour->pGen;

		//pas oublier de mettre un flag pour évite de remonte
		//indcour de deux générations

		if(indcour->pere != NULL)
		{
			// si le pere n'est pas dans la prochaine génération, 
			// reporter l'individu à la prochaine generation
			if(indcour->pere->bFlagSort > nogen+1)
			{
				indcour->pGen = (double*) tabgen[nogen+1];
				tabgen[nogen+1] = indcour;
				decale = 0;
				//ICI
				//printf("gen %d, décalage de : %d\n", nogen, indcour->nom);
			}
			else if(indcour->pere->bFlagSort == nogen+1)//sinon traiter le pere
			{
				//si le pere n'a pas de tuple associer,
				// lui en attribuer un
				if(indcour->pere->ptr[0]==NULL)
				{
					vTuple[taille].clear(indcour->pere);
					indcour->pere->ptr[0]=(void*)&vTuple[taille];
					taille++;
				}
				//ajouter la valeur de l'individu au pere
				((Tuple*)(indcour->pere->ptr[0]))->addtab(indcour->allele);
			}
		}
		if(indcour->mere != NULL)
		{
			//si mere n'est pas de generation nogen+1
			if(indcour->mere->bFlagSort>nogen+1)
			{
				//Si l individu n a pas deja ete ajoute a la
				// generation nogen+1
				//ICI
				if(decale)  //indcour!=tabgen[nogen+1]
				{
					indcour->pGen = (double*) tabgen[nogen+1];
					tabgen[nogen+1] = indcour;
					//ICI
					//printf("gen %d, décalage de : %d\n", nogen, indcour->nom);
				}
			}
			else if(indcour->mere->bFlagSort == nogen+1)//sinon traiter le mere
			{
				//si la mere n'a pas de tuple associer,
				// lui en attribuer un
				if(indcour->mere->ptr[0]==NULL)
				{
					vTuple[taille].clear(indcour->mere);
					indcour->mere->ptr[0]=(void*)&vTuple[taille];
					taille++;
				}
				//ajouter la valeur de l individu a la mere
				((Tuple*)(indcour->mere->ptr[0]))->addtab(indcour->allele);
			}
		}
		
		//passer au suivant:
		indcour = suivant;
	}
	
	return 0;
}
*/
/*! 
	\brief Traite une partie de l'algorithme

	Pour chaque individu associe aux nbT premiers Tuples, 
	calcule et attribut la valeur de allele.

	\param nbT		[in]	nombre de tuples a considerer 

	\param vTuple	[in, out]	Vecteur STL des tuples

	\return 0 si la fonction est execute avec succes

	\remark ATTENTION : Cette fonction modifie l attribut allele
	des individus associés aux tuples
*/
/* PAS UTILISÉE ..??
static int numeroteInd(int nbT, vecTuple_t& vTuple)
{
	int k= 1;
	Tuple* t = &vTuple[0];

	for(int i=0; i<nbT; i++)
	{
		if(!(vTuple[i]==*t))
		{
			t=&vTuple[i];
			k++;
		}
		vTuple[i].getNoeud()->allele= k;
		
	}

	return 0;
}
*/
/*! 
	\brief Calcul du numéro de génération

	Calcule et enregistre dans un premier temps le numéro de génération de tous les individus de la 
	généalogie dans l'attribut bFlagSort. Enregistre dans retour le numéro de génération
	de chaque proposant de plProposant.

	\param plProposant	[in, out]	Vecteur des proposants à considérer

	\param Gen			[in, out]	Une genealogie construite à l'aide de gen.genealogie

	\param retour		[out]	Adressse d'un vecteur de taille NProposant.
							En cas de Succes, contient le numéro de génération de chaque proposant

	\return 0 si la fonction est execute avec succes

	\remark ATTENTION : Cette fonction modifie l attribut bFlagSort de tous les individus de la généalogie

*/
int numeroGen(int* Genealogie, int* plProposant,int NProposant, int* retour)
{
	//CREATION DU TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,NProposant,&NoeudPro);

	for(int i=0;i<lNIndividu;i++)
		Noeud[i].bFlagSort=0;

	//numérotation des générations
	classeGen(Noeud, lNIndividu, NULL , NULL);

	// compléter le vecteur de retour
	for(int i = 0; i< NProposant; i++)
	{
		retour[i]=NoeudPro[i]->bFlagSort;
	}

	return 0;
}
/*! 
	\brief Calcul du numéro de génération minimum

	Calcule et enregistre dans un premier temps le numéro de génération de tous les individus de la 
	généalogie dans l'attribut bFlagSort. Enregistre dans retour le numéro de génération
	de chaque proposant de plProposant.

	\param plProposant	[in, out]	Vecteur des proposants à considérer

	\param Gen			[in, out]	Une genealogie construite à l'aide de gen.genealogie

	\param retour		[out]	Adressse d'un vecteur de taille NProposant.
							En cas de Succes, contient le numéro de génération de chaque proposant

	\return 0 si la fonction est execute avec succes

	\remark ATTENTION : Cette fonction modifie l attribut bFlagSort de tous les individus de la généalogie

*/
int numeroGenMin(int* Genealogie, int* plProposant,int NProposant, int* retour)
{
	//CREATION DU TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,NProposant,&NoeudPro);

	for(int i=0;i<lNIndividu;i++)
		Noeud[i].bFlagSort=0;

	//numérotation des générations
	classeGenMin(Noeud, lNIndividu, NULL, NULL);

	// compléter le vecteur de retour
	for(int i = 0; i< NProposant; i++)
	{
		retour[i]=NoeudPro[i]->bFlagSort;
	}

	return 0;
}
/*! 
	\brief Calcul du numéro de génération moyen

	Calcule et enregistre dans un premier temps le numéro de génération de tous les individus de la 
	généalogie dans l'attribut bFlagSort. Enregistre dans retour le numéro de génération
	de chaque proposant de plProposant.

	\param plProposant	[in, out]	Vecteur des proposants à considérer

	\param Gen			[in, out]	Une genealogie construite à l'aide de gen.genealogie

	\param retour		[out]	Adressse d'un vecteur de taille NProposant.
							En cas de Succes, contient le numéro de génération de chaque proposant

	\return 0 si la fonction est execute avec succes

	\remark ATTENTION : Cette fonction modifie l attribut bFlagSort de tous les individus de la généalogie

*/
int numeroGenMoy(int* Genealogie, int* plProposant,int NProposant, double* retour)
{
	//CREATION DU TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GTRUE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(plProposant,NProposant,&NoeudPro);

	for(int i=0;i<lNIndividu;i++)
		Noeud[i].dFlagSort=(double)0;

	//numérotation des générations
	classeGenMoy(Noeud, lNIndividu);

	// compléter le vecteur de retour
	for(int i = 0; i< NProposant; i++)
	{
		retour[i]=NoeudPro[i]->dFlagSort;
	}

	return 0;
}
