/*! \file outils.cc
\brief Outils d'utilisation generale pour toute la library

Groupes de fonctions d'utilité générale. comme la creation des classe INDsimul,
reordonnancement, creation d'ordre, ordre saut etc... 

  plus les petites fonction mathématiques

	\author Sébastien Leclerc
\contributor Jean-Francois Lefebvre

*/
#include "base.h"
#include "outils.h"
#include "stdlib.h"
#include "limits.h"
#include "md5.h"

#include <string.h>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <random>
#include <Rcpp.h>

using namespace std;

#define R_NO_REMAP

// TransGenCum[p][m][a] : chances d'avoir a allèles si le père a p
// allèles mutants et la mère m allèles mutants. (Cumulatif)
double TransGenCum[3][3][3] =
{
	{ {  1,  1,  1}, { .5,  1,  1}, {  0,  1,  1} },
	{ { .5,  1,  1}, {.25,.75,  1}, {  0,0.5,  1} },
	{ {  0,  1,  1}, {  0, .5,  1}, {  0,  0,  1} }
};

	
// TransGenCum[p][m][a] : chances d'avoir a allèles si le père a p
// allèles mutants et la mère m allèles mutants. (non-Cumulatif)
double TransGen[3][3][3] =
{
	{ {  1,  0,  0}, { .5, .5,  0}, {  0,  1,  0} },
	{ { .5, .5,  0}, {.25, .5,.25}, {  0, .5, .5} },
	{ {  0,  1,  0}, {  0, .5, .5}, {  0,  0,  1} }
};
	
//CONSTANTE POUR LA GENEALOGIE BINAIRE
const int GENBIN_NBIND_OFFSET = 12;

const int GENBIN_NIND_OFFSET = 8;
const int GENBIN_NENFANT_OFFSET = 9;
const int GENBIN_PROFMAX = 10;
const int GENBIN_NINDMASC = 11;

const int GENBIN_MD5_1 = 4;
const int GENBIN_MD5_2 = 5;
const int GENBIN_MD5_3 = 6;
const int GENBIN_MD5_4 = 7;



// *************************************************************************** //
//		STRUCTURE
// *************************************************************************** //

/// Utiliser par qsort pour trier est valeur d'un tableau mais memoriser leur indice(position) d'origine
/**
	Utilise dans plusieur fonction pour permettre de trier un vecteur par ordre numerique mais
	memoriser la position initiale de chaque element pour etablir des liens.
*/
struct CDuoPair
{
	///Etiquette pour cette paire
	int nom;	
	///Indice de cette paire
	int pos;
};

struct CIndPereMere
{
	int ind;
	int pere;
	int mere;
	int sex;
};

static int WINCDECL QsortDuoPaircompare(const void *p1, const void *p2);
static int WINCDECL QsortCIndPereMereCmp(const void *p1, const void *p2);
static int ReTrouverIndice(int nom, CDuoPair* Pair, int iNind , int* resultat);
//static void Exchange(void *a, void *b, int size);
static int temoin(unsigned int a, unsigned int n);
enum ENUMBANQUE {PROPOSANT,ANCETRE};
static int LoadVec(ENUMBANQUE banque,int* vec, int nb,CIndSimul*** NproAnc);
static int FlushProposantAncetre(ENUMBANQUE banque);
static int LoadVecGroupe(ENUMBANQUE banque,int* BorneGr,int nbGroupe, CIndSimul**** GrProAnc,int** nIndGr);
static int FlushGroupeProposantAncetre(ENUMBANQUE banque);

// *************************************************************************** //
//		UTILITAIRE SPLUS
// *************************************************************************** //
#ifndef MODETEST
	const char* DescIEEEValue(int* val)
	{
		char c[10];
		sprintf(c, "%d", *val);
		double valTmp = atof(c);
		const char* buffer[]={"NA","Not a Number","Infinite"};
		//if (is_na(val,S_MODE_INT))
			//return buffer[0];
		//else 
			//if (is_nan(val,S_MODE_INT))
			if (isnan(valTmp))
				return buffer[1];
			else
				//if (is_inf(val,S_MODE_INT))
				if (isinf(valTmp))
					return buffer[2];
				else
					//return NULL;
					return buffer[0];
	}
#else
	char* DescIEEEValue(int* val) {return NULL;}
#endif
// *************************************************************************** //
//		GESTION DE LA MEMOIRE
// *************************************************************************** //
///Taille de chaque tableau successif pour la gestion memoire
const int MAXALLOCINGESTIONMEMOIRE=100;

///Structure utiliser par GestionMemoire pour gere la "garbage collection"
/**
	A l'interne GestionMemoire travaille avec un liste de tableau de void*.
	A l'initialisation, GestionMemoire cree un GestionMemoireBlock qui contient un tableau de pointeur 
	de taille MAXALLOCINGESTIONMEMOIRE. Par la suite, chaque fois qu'une fonction effectue une allocation memoire memalloc
	l'adresse du bloc memoire nouvellement allouer est stocker dans le tableau. Si le nombre d'allocation
	depasse la capacite du tableau alors un nouveau GestionMemoireBlock est creer a la suite du precedent et l'operation
	continu.
  */
struct GestionMemoireBlock
{
	//Tableau de pointeur void* pour stocker l'adresse des bloc memoire allouer par GestionMemoire::alloc
	void**  tableau;							
	///GestionMemoireBlock suivant pour creer une liste de GestionMemoireBlock
	GestionMemoireBlock *next;  
};

/*! 
	\brief Constructeur de GestionMemoire

	Construit le premier GestionMemoireBlock et initialise la GestionMemoire
	GestionMemoire utilise le type d'allocation specifier par la definition de
	USESPLUSALLOC.

	\param UseStdMalloc	[in] Si UseStdMalloc >= 1 alors la valeur de USESPLUSALLOC est ignorer
							 et alloc fera toutes les allocations memoire a l'aide de malloc.
   \sa GestionMemoire::alloc
*/
GestionMemoire::GestionMemoire(char UseStdMalloc)
{
	try{
	n=-1;
	UseMalloc		= (UseStdMalloc==1);	
	startblock	= (GestionMemoireBlock*) malloc( sizeof(GestionMemoireBlock) );	
	if (!startblock){
//	 GENError("Insufficient memory"); 
	 throw std::range_error("Insufficient memory");
	}
	startblock->tableau	= (void**) malloc(  MAXALLOCINGESTIONMEMOIRE * sizeof(void*) );	
	startblock->next	= NULL;
	tableaublock		= startblock;
 	} catch(std::exception &ex) {
 		forward_exception_to_r(ex);
	} catch(...){
		::Rf_error("c++ exception (unknown reason)"); 
	}
}

/*! 
	\brief Alloue de la memoire suivant la methode specifier par USESPLUSALLOC
	
	Alloue de la memoire et retourne un pointeur vers l'adresse du nouveau bloc  
	
	Si USESPLUSALLOC est defini alors  alloc alloue de la memoire a l'aide de S_alloc
	dans le cas contraire, utilise malloc.
	
	\param nelement	[in] Nombre d'element du vecteur a allouer

	\param size	[in] Taille (en byte) de chaque element du vecteur.
	
	\param NLigne [in] Utiliser au débuggage, no de la ligne qui fait l'appel

	\return Un pointeur vers un bloc d'adresse de taille nelement*size
*/
void* GestionMemoire::alloc(int nelement, size_t size)
{
	try{
	void* tmp;
	//if (UseMalloc)
		tmp= (void*) malloc((nelement)*(size));		
	//else
	//	tmp= (void*) memallocIN(nelement,size);
	if (tmp!=NULL)
	{
		if ((++n)==MAXALLOCINGESTIONMEMOIRE)
		{
			tableaublock->next		= (GestionMemoireBlock*) malloc( sizeof(GestionMemoireBlock) );	
			tableaublock->next->tableau	= (void**) malloc(  MAXALLOCINGESTIONMEMOIRE * sizeof(void*) );			
			tableaublock->next->next	= NULL;
			n				= 0;
			tableaublock			= tableaublock -> next;			
		}
		tableaublock->tableau[n]	=tmp;
	}
	else{		
//		GENError("Insufficient memory"); //GENError("Memoire insuffisante");
		throw std::range_error("Insufficient memory");
	}
	return tmp;
	} catch(std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...){
		::Rf_error("c++ exception (unknown reason)"); 
	} 
	return 0;
}

/*! 
	\brief Ajoute la valeur d'un pointeur a la table de GestionMemoire
	
	Cette fonction peut-être utilise si on desire rajoute un pointeur dans la table de GestionMemoire.
	Car a la destruction de GestionMemoire, tout les blocs memoire contenu dans la table seront
	relacher (Free). 
	  
	\param item	[in] Pointeur vers un bloc d'adresse valide (allouer a l'aide de malloc)

	\return 0 si la fonction est executé avec succès
*/

void GestionMemoire::add(void* item)
{
	if ((++n)==MAXALLOCINGESTIONMEMOIRE)
	{
		tableaublock->next		= (GestionMemoireBlock*) malloc( sizeof(GestionMemoireBlock) );	
		tableaublock->next->tableau	= (void**) malloc(  MAXALLOCINGESTIONMEMOIRE * sizeof(void*) );
		tableaublock->next->next	= NULL;
		n				= 0;
		tableaublock			= tableaublock -> next;			
	}
	tableaublock->tableau[n]	=item;
}
/*! 
	\brief Destruction et "Garbage collector" de la classe GestionMemoire

	A la destruction de GestionMemoire, la toutes les adresses qui on été memorise
	sont relache (delete). 
*/

GestionMemoire::~GestionMemoire() 
{
	//Bon.. on efface tous jusqu'au dernier tableau
	GestionMemoireBlock *current= startblock;
	int MAX;
	while(current!=NULL)
	{
		//Est-ce le dernier tableau
		if (current->next!=NULL)
			MAX=MAXALLOCINGESTIONMEMOIRE;  //NON
		else
			MAX=n+1;  //OUI

		//Liberer la ram de chaque élément du tableau
		for(int i=0;i<MAX;i++)
		{
			if (current->tableau[i]!=NULL)
			{
				//if (UseMalloc)
				//{
					free(current->tableau[i]);				
				//}
				//else
					//SI SA BUG A CETTE LIGNE CI REGARDE LE CONTENU DE current->noligne[i]
					//Pour connaitre a quel ligne de code c'est fait l'affectation qui cause probleme
				//	memfreeIN(current->tableau[i]);
			}
		}
		
		//Libere les tableaux eux-meme
		free(current->tableau);

		//Avance au suivant
		current=current->next;
	}
	free(startblock);
}


// *************************************************************************** //
//		CREATION DES STRUCTURES NOEUDS
// *************************************************************************** //

/*! 
	\brief Initialise un vecteur de CIndSimul a partir d'une genealogie (ind,pere,mere)

		Cette fonctione initialise un vecteur de Noeud(CIndSimul) deja construit en incorporant les
		lien entre parent et enfant.
	  
	\param Noeud	[in] Adresse d'un vecteur de CIndSimul de iNind de int.

	\param noind	[in] Vecteur de taille iNind se composant des no des individu à étudier
	\param pere		[in] Vecteur de taille iNind se composant des no des peres
	\param mere		[in] Vecteur de taille iNind se composant des no des meres
	\param sex		[in] Vecteur de taille iNind se composant des no des sexes
	\param iNind	[in] Nombre d'individu dans la Genealogie. 

	\retval countchildren	[out] (peut etre NULL) Un pointeur vers un int qui s'il est different de NULL.
								  se vois assigne comme valeur, le nombre Total d'enfant
								  dans la genealogie.

	\retval Trie				[out] (peut etre NULL) Adresse d'un Vecteur de iNind CDuoPair.
									Au retour, si !=NULL, les CDuoPair contienne les noind 
									et leur indide dans ce meme vecteur. Mais le vecteur est 
									ordonne par ordre croissant de noind.
									
	\param ChildArrayStart		[in/out] Si !=NULL alors les parents/enfants seront chargé dans les noeud (CIndSimul).
								 Dans le cas contraire, seul les pointeurs vers les parents
								 seront construit.
								 ChildArray sera affecte la valeur vers le tableau d'enfants (n.b que l'on doit détruire avec détruit structure)
								 									
	\remark Cette fonction ne devrais pas etre utilise directement dans une fonction.
			Utiliser ReCreeStructure de preference car celle-ci n'est pas tres efficace.
			<br>Si LoadChildren!=0 alors DetruireStructure doit etre appeler pour evite les fuites de memoire.

	\sa DetruireStructure
*/
void CreeStructure(CIndSimul* Noeud,int* noind,int* pere, int* mere,int* sex, int iNind,int* countchildren,CDuoPair *Trie, Clist** ChildArrayStart)
{
	//CREER LA STRUCTURE
	int i;
	int nom;
	Clist *liste=NULL;
	Clist *current=NULL;
	int NombreEnfant=0;
	//int tst=0;

	for(i=0;i<iNind;i++)
	{
		Noeud[i].nom=noind[i];
		Noeud[i].mere=NULL;
		Noeud[i].pere=NULL;
		Noeud[i].sex=(sex?(sex_t)sex[i]:GEN_INCONNU);
		Noeud[i].fils=NULL;		
		Noeud[i].pGen=NULL;
		Noeud[i].etat=GENNONEXPLORER;
		Noeud[i].noind=i;
		Noeud[i].allele=0;									
	}

	//CREATION DE L'INDEX POUR LA RECHERCHE
	//SI AUCUN VECTEUR DE DUOPAIR EST FOURNI EN CONSTRUIRE 1
	INITGESTIONMEMOIRE;	
	if (Trie==NULL)
	 {
		Trie	=(CDuoPair*)memalloc(iNind,sizeof(CDuoPair));
	 }
	//Compte le nombre d'enfant
	for(i=0;i<iNind;i++)
	{
	    //VERIFICATION AU PASSAGE
	    Trie[i].nom=noind[i]; 
	    Trie[i].pos=i;
	    	
	    if (pere[i]!=0)	++NombreEnfant;
	    if (mere[i]!=0)	++NombreEnfant;
	}
	//TRIE DES DUOPAIR
	qsort(Trie,iNind,sizeof(CDuoPair),QsortDuoPaircompare);
	
	//ENFANT
	Clist* ChildArray=NULL;
	if (ChildArrayStart)
	{
		ChildArray= (Clist*) memallocIN(NombreEnfant,sizeof(Clist));
		*ChildArrayStart = ChildArray;
	}

	//Descendance
	int Code;
	for(i=0;i<iNind;i++)
	{
		if (mere[i]!=0)
		{
			Code=ReTrouverIndice(mere[i],Trie,iNind,&nom);
			if (Code==0)
			{
				Noeud[i].mere=&Noeud[nom];
				if (ChildArrayStart)
				{
					//Descendance mere
					liste=ChildArray++;
					liste->next=NULL;
					liste->noeud=&Noeud[i];

					if ( (*Noeud[i].mere).fils==NULL)
						(*Noeud[i].mere).fils=liste;				
					else
					{
						current=(*Noeud[i].mere).fils;
						while(current->next!=NULL)
							current=current->next;
						current->next=liste;
					}
				}
			}
		}
		if (pere[i]!=0)
		{
			Code=ReTrouverIndice(pere[i],Trie,iNind,&nom);
			if (Code==0)
			{
				Noeud[i].pere=&Noeud[nom];
				if (ChildArrayStart)
				{
					//Descendance pere
					liste=ChildArray++;					
					liste->next=NULL;
					liste->noeud=&Noeud[i];
	
					if ( (*Noeud[i].pere).fils==NULL)
						(*Noeud[i].pere).fils=liste;				
					else
					{
						current=(*Noeud[i].pere).fils;
						while(current->next!=NULL)
							current=current->next;
						current->next=liste;
					}
				}
			}

		}	

	}
	if (countchildren!=NULL)
		*countchildren=NombreEnfant;

	return;
}


/*! 
	\brief Recupere la memoire alloue par CreeStructure
		  
	\param ChildArray	[in, out] Clist à détruire

	\sa CreeStructure
*/
int DetruireStructure(Clist* ChildArray)
{       
	memfreeIN(ChildArray);
	return 0;
}

  
/*! 
	\brief INTERNE: Utiliser par CreeStructure

	Execute une recherche binaire p/r a nom dans un vecteur de CDuoPair
				
	\param nom [in] No d'individu a recherche

	\param Pair	 [in]  Adresse d'un vecteur de CIndSimul representant une genealogie

	\param iNind [in] Nombre d'individu dans la genealogie (Taille du vecteur Noeud)

	 \retval resultat [out] adresse d'un int.
							si succes, se nombre vaut l'indice correspondant au nom dans DuoPair
								Ce pointeur est fourni par la fonction ReCreeStructure
  	
	\return 0 si trouver, 1 si pas trouve
	
	 \sa CreeStructure
*/
static int ReTrouverIndice(int nom, CDuoPair* Pair, int iNind , int* resultat) 
{
	//Recherche binaire 	
	int uplimit=iNind; 
	int downlimit=-1;	 
	int pivot=iNind/2; 
	
	*resultat=-1;

	while(1) 
	{ 
		const int Value=Pair[pivot].nom;
		if (Value==nom) 
		{//Resultat trouver 
			*resultat = Pair[pivot].pos;
			return 0;
		} 
		else 
		{ 
			if (Value<nom) 
			{//La valeur est plus petite.. 
				downlimit=pivot; 
				pivot=(uplimit+pivot)/2; 
			
				if (pivot==downlimit) 
					return 1; //Valeur pas trouver 
			
			} 
			else 
			{//La valeur est plus grande 
				uplimit=pivot; 
				pivot=(downlimit+pivot)/2; 

				if (pivot==uplimit) 
					return 1; //Valeur pas trouver 

			} 
		
		} 
	} 
	
  
} 


// *************************************************************************** //
//		EBRANCHAGE ET ORDONNANCEMENT DES NOEUDS
// *************************************************************************** //

/*! 
	\brief Explore la genealogie et identifie l'utilite de chaque noeud.

	<br>Identifie tous les noeuds qui entre dans la 
		composition d'un ou plusieurs chemin qui vont d'un ancetre a un ou plusieur proposant.
		

		<br>Cette fonction s'utilise de la maniere suivante.
		initialisation :  Tous les Noeuds.etat = GENNONEXPLORER
		
		<br>1- Toute les proposants.etat = GENPROPOSANTINUTILE;
		<br>2- Toute les ancetre.etat = GENDEPART
		<br>3- Lancer ExploreArbre() sur chaque ancetre

  <br>Tous les noeuds qui font partie d'un ou plusieurs chemin d'un ancetre a un proposant verront
  leur propriete etat etre modifier de la facons suivante.

  <br>etat= 
  	<br>GENNONEXPLORER  : Ne fait partie de la descendance d'aucun des ancetres (donc inutile)
	<br>GENINUTILE	    : Fait partie de la descendance d'au moins un ancetre mais
						  ne fait pas partie d'un chemin entre un ancetre et un proposant		  
	<br>GENNOEUD		: Fait partie d'un chemin entre un ancetre et un proposant
	<br>GENDEPART		: Ce Noeud est un Ancetre 
	<br>GENPROPOSANTINUTILE : Ce noeud est un proposant mais il ne fait pas partie de la
							  descendance d'un des ancetres
	<br>GENPROPOSANT	: Ce Noeud est un proposant et fait parti de la descendance d'au moins un ancetre	
		  
	\param Noeud	[in] Adresse d'un vecteur de CIndSimul representa une genealogie



	\sa typenoeud_t
*/
int ExploreArbre(CIndSimul* Noeud)
{
	//Explore l'arbre et change l'état des noeud explorer selon qu'il sont utile ou non pour résoudre
	//le problème actuel
	//un noeud est utile s'il fait parti d'une chemin qui part d'un point de départ et vas à un proposant

	int isutile;
	Clist *current;

	//printf("\n%d == %s",Noeud->nom,stype[Noeud->etat]);

	//cas spéciaux
	switch (Noeud->etat)
	{
		case GENINUTILE:
				//Le noeud est déjà explorer inutile de le refaire
				return 0;
			break;

		case GENNOEUD:
		case GENPROPOSANT:
			//Le noeud est déjà explorer inutile de le refaire
				return 1;
			break;
		
		case GENDEPART:
			//Ce Noeud est un proposant ou un depart il est forcément utile mais il faut explorer les enfants
			current=Noeud->fils;
			if (current!=NULL)
			{
					do
					{
						ExploreArbre(current->noeud);				
						current=current->next;				
					}
					while(current!=NULL);
			}
			return 1;
			break;		

		case GENPROPOSANTINUTILE:
			//Ce Noeud est un proposant ou un depart il est forcément utile mais il faut explorer les enfants
			Noeud->etat=GENPROPOSANT;
			current=Noeud->fils;
			if (current!=NULL)
			{
					do
					{
						ExploreArbre(current->noeud);				
						current=current->next;				
					}
					while(current!=NULL);
			}
			
			return 1;
			break;
		
		case GENNONEXPLORER:
			//Noeud non-explorer. sont status dépend des ses enfants (s'Il y a un proposant dans sa descendance)
			isutile=0;
			current=Noeud->fils;
			if (current!=NULL)
			{
					do
					{
						isutile+=ExploreArbre(current->noeud);				
						current=current->next;				
					}
					while(current!=NULL);
			}

			if (isutile>0)
			{
				Noeud->etat=GENNOEUD;
				return 1;
			}
			else
			{
				Noeud->etat=GENINUTILE;
				return 0;
			}			
			break;
	}
			
	return 99;	
}



/*! 
	\brief Explore la genealogie et identifie tous les descendants d'un ancetre

	<br>Associe un nouvel etat a chaque noeud qui descend d'un ancetre
		
		<br>Cette fonction s'utilise de la maniere suivante.
		initialisation :  Tous les Noeuds.etat = GENNONEXPLORER
		
		<br>1- (facultatif) Toute les proposants.etat = GENPROPOSANTINUTILE; 
		<br>2- Toute les ancetre.etat = GENDEPART
		<br>3- Lancer ExploreArbreTousDescendant() sur chaque ancetre

  <br>Tous les noeuds qui font partie d'un ou plusieurs chemin d'un ancetre a un proposant verront
  leur propriete etat etre modifier de la facons suivante.

  <br>etat= 
  	<br>GENNONEXPLORER  : Ne fait partie de la descendance d'aucun des ancetres 
	<br>GENINUTILE	    : Ne fait partie de la descendance d'aucun des ancetres (ne devrais pas arrive)
						  	  
	<br>GENNOEUD		: Est un descendant de l'ancetre
	<br>GENDEPART		: Ce Noeud est un Ancetre 
	<br>GENPROPOSANTINUTILE : Ce noeud est un proposant et il ne fait pas partie de la
							  descendance de l'ancetres
	<br>GENPROPOSANT	: Ce Noeud est un proposant et fait parti de la descendance de l'ancetre
		  
	\param Noeud	[in] Adresse d'un vecteur de CIndSimul representa une genealogie
 */
void ExploreArbreTousDescendant(CIndSimul* Noeud)
{
	//Explore l'arbre et change l'état des noeud explorer selon qu'il sont utile ou non pour résoudre
	//le problème actuel
	//un noeud est utile s'il fait parti d'une chemin qui part d'un point de départ et vas à un proposant

	Clist *current;

	//cas spéciaux
	switch (Noeud->etat)
	{
		
		case GENNOEUD:
		case GENPROPOSANT:
			break;
		
		case GENDEPART:
			current=Noeud->fils;
			if (current!=NULL)
			{
					do
					{
						ExploreArbreTousDescendant(current->noeud);				
						current=current->next;				
					}
					while(current!=NULL);
			}
			return;
			break;		

		case GENPROPOSANTINUTILE:
			Noeud->etat=GENPROPOSANT;
			current=Noeud->fils;
			if (current!=NULL)
			{
					do
					{
						ExploreArbreTousDescendant(current->noeud);				
						current=current->next;				
					}
					while(current!=NULL);
			}
			break;

		case GENINUTILE:
		case GENNONEXPLORER:
			Noeud->etat=GENNOEUD;
			current=Noeud->fils;
			if (current!=NULL)
			{
					do
					{
						ExploreArbreTousDescendant(current->noeud);				
						current=current->next;				
					}
					while(current!=NULL);
			}
			break;
	}
			
	return;	
}


/*! 
	\brief Preparation d'un ordre de priorite

	Cette fonction prepare les Noeud pour demarrer la fonction PrepareSortPrioriteArbre<br>
	Utilise apres avoir lance ExploreArbre() sur chaque noeud<br>

	Initialise avec: Tous les Noeuds.bFlagSort = 0
		
	\param Noeud	[in] Adresse d'un vecteur de CIndSimul representa une genealogie

	\param iNind	[in] Nombre d'individu dans la Genealogie. 

	\sa StartSortPrioriteArbre 
	\sa ExploreArbre
*/
void PrepareSortPrioriteArbre(CIndSimul* Noeud,int iNind)
{
 	for(int i=0;i<iNind;i++)
	{
		
		if (Noeud[i].pere==NULL || (Noeud[i].pere->etat==GENNONEXPLORER || Noeud[i].pere->etat==GENINUTILE))
				Noeud[i].bFlagSort=-1;
		else
			if (Noeud[i].mere==NULL || (Noeud[i].mere->etat==GENNONEXPLORER || Noeud[i].mere->etat==GENINUTILE))
				Noeud[i].bFlagSort=-1;
			else
				if (Noeud[i].pere->etat==GENPROPOSANTINUTILE || Noeud[i].mere->etat==GENPROPOSANTINUTILE)
					Noeud[i].bFlagSort=-1;
				else
					Noeud[i].bFlagSort=0;
	}
	return;
}


/*! 
	\brief Ordonne une genealogie en respectant divers consigne.
	
	Genere un vecteur de pointeur de noeud de taille *index contenant seulement les noeud de type 
	GENNOEUD, GENDEPART, GENPROPOSANT.  

	<br>Les noeuds sont en ordre de dependance. Ce qui veut dire que tous les noeuds qui depende
	seulement d'un seul autre noeud vont etre un a la suite des autres.

    <br>Autrement dis, on peut dire que que l'ordre est une serie de petite liste de noeud qui sont 
	    interdependant.

	<br>La distance entre l'element et le dernier element de la chaine
	
	<br>Cette fonction s'utilise de la maniere suivante.
	initialisation :  Tous les Noeuds.etat = GENNONEXPLORER & Noeuds.bFlagSort=0
	
	<br>1- Modifier les .etat pour les proposant/ancetre
	<br>2- Lancer une exploration sur chaque ancetre : ExploreArbre()
	<br>3- Lancer PrepareSortPrioriteArbre() pour preparer	
	<br>4- Lancer StartSortPrioriteArbre sur chaque ancetre

	\param Noeud	 [in] Adresse d'un vecteur de CIndSimul representa une genealogie

  	\retval Ordre	 [out] ptr vers un vecteur de taille Nindividu, 
						   Les *index premiers element du vecteur pointeront vers les
						   *index Noeud qui compose l'ordre que l'on vient de calculer.						   

  	\retval index	 [out] Un pointeur vers un int
							En cas de succes, la valeur est le nombre d'element dans ordre 

	\retval TableSaut [out] Vecteur de int de taille *index qui indique
						   de combien il faut avancer dans le vecteur pour passer
						   a la liste de dependance suivante.

	\sa ExploreArbre PrepareSortPrioriteArbre 
*/
void StartSortPrioriteArbre(CIndSimul* Noeud,CIndSimul** Ordre,int* index,int* TableSaut)
{
	//Complete un ordre et un nombre de saut en tenant compte des dépendances		
	Clist *listedebut;	
	Clist *current;
	Clist *tmp;
	//int oldflag;

	//oldflag=Noeud->bFlagSort;
	Noeud->bFlagSort=5; //Ce noeud à déjà été traité


	// *******Première passe.. les fils à un parent "utile" (ce noeud-ci)*******
	current=Noeud->fils;
	if (current!=NULL)
	{
			listedebut=NULL;	
			SortPrioriteArbre(NULL,NULL,NULL,NULL,&listedebut);

			do
			{
				if (current->noeud->bFlagSort==-1)
					SortPrioriteArbre(current->noeud,Ordre,index,TableSaut);						
				current=current->next;				
			}
			while(current!=NULL);
			
			//DEPILAGE PILE
			tmp=NULL;
			while(listedebut!=NULL)
			{
				SortPrioriteArbre(listedebut->noeud,Ordre,index,TableSaut);						
				tmp=listedebut->next;
				memfreeIN(listedebut);
				listedebut=tmp;
			}			
	}

	// *******Deuxième passe.. les fils à deux parent "utile"***
	current=Noeud->fils;
	if (current!=NULL)
	{
			do
			{
				listedebut=NULL;	
				SortPrioriteArbre(NULL,NULL,NULL,NULL,&listedebut);

				switch(current->noeud->bFlagSort)
				{
				case -1:	//dejà fait
					break;
				case 0:
					//Dans ce cas on set le flag et on vas traiter ce noeud au passage suivant
					current->noeud->bFlagSort=1;
					break;
				case 1:						
					//On traite le noeud mais sont nombre de saut est ignorer
					//(*listenoeud)=(Clist*) memalloc(1,sizeof(Clist));
					//(*listenoeud)->noeud=current->noeud;
					//(*listenoeud)->next=NULL;
					//listenoeud=&((*listenoeud)->next);					
					SortPrioriteArbre(current->noeud,Ordre,index,TableSaut);		
					break;
				}
				current=current->next;				

				//Partie deux.. depilage de la liste..
				tmp=NULL;
				while(listedebut!=NULL)
				{
					SortPrioriteArbre(listedebut->noeud,Ordre,index,TableSaut);						
					tmp=listedebut->next;
					memfreeIN(listedebut);
					listedebut=tmp;
				}
			}
			while(current!=NULL);
	}

return;
}


/*! 
	\brief Fonction utiliser a l'interne par StartSortPrioriteArbre
			
		Fonction complementaire recursive qui sert a StartSortPrioriteArbre

	\param Noeud	 [in] Adresse d'un vecteur de CIndSimul representa une genealogie

  	\retval Ordre	 [out] ptr vers un vecteur de taille Nindividu, 
						   Les *index premiers element du vecteur pointeront vers les
						   *index Noeud qui compose l'ordre que l'on vient de calculer.						   

  	\retval index	 [out] Adresse d'un vecteur de CIndSimul representa une genealogie

	\retval TableSaut [out] vecteur de int de taille *index qui indiquombre d'individu dans la Genealogie. 
	
	\retval list	[out] List de noeud pour faire le traitement par apres. 

	\sa ExploreArbre PrepareSortPrioriteArbre 
*/
int SortPrioriteArbre(CIndSimul* Noeud,CIndSimul** Ordre,int *index,int *TableSaut,Clist **list)
{
	//Un essais pour optimiser les performances...
	static Clist **listenoeud=NULL;
	int Saut=0;
	Clist *current;
	int CurrentIndex;
	int oldflag;

	if (list!=NULL)
	{
		listenoeud=list;
		return 0;
	}


	//Les noeuds déjà traité sont ignorer ainsi que les autres types de noeuds
	if ((Noeud->etat!=GENPROPOSANT && Noeud->etat!=GENNOEUD) || Noeud->bFlagSort==5)
		return 0;
	
	//Ajouter dans le tableau d'ordre
	CurrentIndex=*index;
	Ordre[CurrentIndex]=Noeud;
	*index+=1;
	
	oldflag=Noeud->bFlagSort;
	Noeud->bFlagSort=5; //Ce noeud à déjà été traité

	//Première passe.. les fils à un parent "utile" (ce noeud-ci)
	current=Noeud->fils;
	if (current!=NULL)
	{
			do
			{
				if (current->noeud->bFlagSort==-1)
					Saut+=SortPrioriteArbre(current->noeud,Ordre,index,TableSaut);						
				current=current->next;				
			}
			while(current!=NULL);
	}
	if (TableSaut!=NULL)
		TableSaut[CurrentIndex]=Saut;	
	if (oldflag==-1)
		++Saut; 	

	//Deuxième passe.. les fils à deux parent "utile"
	current=Noeud->fils;
	if (current!=NULL)
	{
			do
			{
				switch(current->noeud->bFlagSort)
				{
				case -1:	//dejà fait
					break;
				case 0:
					//Dans ce cas on set le flag et on vas traiter ce noeud au passage suivant
					current->noeud->bFlagSort=1;
					break;
				case 1:						
					//On traite le noeud mais sont nombre de saut est ignorer
					(*listenoeud)=(Clist*) memallocIN(1,sizeof(Clist));
					(*listenoeud)->noeud=current->noeud;
					(*listenoeud)->next=NULL;
					listenoeud=&((*listenoeud)->next);
					//SortPrioriteArbre(current->noeud,Ordre,index,TableSaut);	
					break;
				}
				current=current->next;				
			}
			while(current!=NULL);
	}
	
return Saut;	
}


/*! 
	\brief Trie les noeuds par ordre Parent -> Enfant
			
	Trie les noeuds par ordre les parents toujours avant les enfants.
	Le resultat est retourne sous la forme d'un vecteur de pointeur de taille iNind.

	<br>Les noind de chaque noeud sont modifier pour indiqué leur nouvel ordre.
	Soit de 0 à (iNind-1) 

	\param Noeud	 [in]  Adresse d'un vecteur de CIndSimul representant une genealogie

  	\retval Ordre	 [out] Pointeur verr un vecteur de taille iNind
						   Si !=NULL et succes,un vecteur de pointeur
						   vers les noeud dans le bon ordre
						 
  	\param iNind	 [in]  Adresse d'un vecteur de CIndSimul representa une genealogie

	\param SensInverse [in] Si SensInverse!=0 alors l'ordre des noeuds est inversé.
							Auto pour les Noind que pour le vecteur ordre.
*/
int OrdonneStructure(CIndSimul* Noeud, CIndSimul** Ordre, int iNind, int SensInverse, int* profMax)
{
	try{
	//QUEL EST LE BON ORDRE...
	//noind 0.. (N-1) proposant dans l'ordre de priorite
	
	int i;
	//vecteur d'acceleration
	INITGESTIONMEMOIRE;
	int *next	=(int*)memalloc(iNind+1,sizeof(int));

	//initialisation
	for(i=0;i<iNind;i++)
	{
		Noeud[i].noind=-1; 
		Noeud[i].bFlagSort=-1;  //Flag qui indique a quel iteration le noeud a été place
		next[i]=i+1;
	}
	next[iNind-1]=-1;

	int n=0; //indice actuel a placer
	int itemchange;
	int debut=0;
	int indice;
	int oldindice;
	int indiceBoucle=-1; //Position du dernier implacable
	
	//No de l'iteration
	int iteration=0;
	while(n < iNind)
	{
	  itemchange =0;
	  iteration++;
	  indice     =debut;
	  oldindice  =-1;
	  
	  while(indice!=-1)
	  {
	    if ((Noeud[indice].pere==NULL || (Noeud[indice].pere->noind!=-1 && Noeud[indice].pere->bFlagSort!=iteration) )	&& \
		   (Noeud[indice].mere==NULL || (Noeud[indice].mere->noind!=-1 && Noeud[indice].mere->bFlagSort!=iteration) ) 	)
	    {
		 //Positionnnement valide
		Noeud[indice].noind=n++;
		Noeud[indice].bFlagSort=iteration;

		 //modification de la liste d'indice
		if (oldindice==-1)
			debut=next[indice];
		else
			next[oldindice]=next[indice];
		indice=next[indice];
		itemchange++;
	    }
	    else
	    {//Pas bon on continue
		 indiceBoucle=Noeud[indice].nom;
	      oldindice=indice;
		 indice=next[indice];
	    }
	  }//fin while chaque indice
	  
	  if (itemchange==0){ //La combinaison est invalide.. (existe des cycles)
//		GENError("The genealogy has at least one cycle (Number of individuals involved: %d    Number of an individual: %d )",iNind-n,indiceBoucle);
		char erreur[TAILLEDESCRIPTION];
		sprintf(erreur, "The genealogy has at least one cycle (Number of individuals involved: %d    Number of an individual: %d )",
				iNind-n,indiceBoucle);
		throw std::range_error(erreur);
		//GENError("LA GÉNÉALOGIE COMPORTE AU MOINS UN CYCLE (Nombre d'individus impliqués: %d    No d'un individu: %d )",iNind-n,indiceBoucle);
	  }
	}//jusqu'a la fin de la recherche

	//Inversion de sens 
	const int IndOne=iNind-1;
	if (SensInverse!=0) //mettre==0 si ancienne version
	{
		for(i=0;i<iNind;i++)
			Noeud[i].noind=IndOne-Noeud[i].noind; 
	}

	//Creation d'un vecteur dans le bon ordre
	if (Ordre!=NULL)
	{
		for(i=0;i<iNind;i++)
		{			
			Ordre[Noeud[i].noind] = &Noeud[i];		
		}
	}

	/*Profondeur maximale*/
	//if(*profMax)
		*profMax=iteration;

	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}

// *************************************************************************** //
//		DIVERS
// *************************************************************************** //


// *************************************************************************** //
//		utilitaire mathématique
// *************************************************************************** //

int interval(int x, int min, int max)
{
	if (x<min)
		return min;
	else
		if(x>max)
			return max;
		else
			return x;
}


/*! 
	\brief Arrondi un nombre(double) a l'entier le plus proche
				
	\param value [in]  Nombre double a arrondir
  	
	\return Le nombre intier le plus proche de value
*/
/*double round(double value)
{
	if (value-int(value)>.5) 
		return int(value)+1;
	
	return int(value);
}
*/

/*! 
	\brief Evalue/Calcule une puissance de 0.5

		Utilise une cache pour optimiser les performances
				
	\param y [in] Nombre double a arrondir
  	
	\return La puissance de (0.5)^y
*/
double pow2(int y)
{
   //Toute les puissance de 2 de 0.. pow2TailleCache precalculer
	
	const int pow2TailleCache = 150;

	static double arr[pow2TailleCache+1]={-99};
	if (arr[0]==-99)
	{		
		for(int a=0;a<=pow2TailleCache;a++)
			arr[a]=pow(0.5,a);
	}
	
	if (y>=0 && y<=pow2TailleCache)
	  return arr[y];
    
	double tmp=arr[pow2TailleCache];		
	for(int a=pow2TailleCache;a<y;a++)
		tmp*=0.5;		
	return tmp;
}



//TEST PRIMALITE MILLER RABIN 

/*! 
	\brief Algorithme miller Rabin évalue la primalite de n

	Evalue la primalite de n et ce pour t test.
	La probabilité d'erreur est de l'ordre de 2^t.
					
	\param n [in] Nombre entier a tester

	\param t [in] Nombre de test de primalite
  	
	\return 0 si le nombre est premier 1 dans le cas contraire..
*/
int millerRabin(unsigned int n, unsigned int t) 
{
	for (unsigned int i=0; i<t; i++) 
	{
		unsigned int a = irand(2, n-1);
		if (temoin(a,n)) 
		{
				return 1;          // n est composé
		}
	}
	return 0;             // n est tres probablement premier
} 

/*! 
	\brief Procedure Interne utiliser dans miller-rabin
				
	\param a [in] nombre tire aléatoirement entre 2..n-1
	\param n [in] Nombre donc on evalue la primalité
  	
	\return	1 s'il y a confirmation que n est composé. 0 dans le cas contraire
*/
static int temoin(unsigned int a, unsigned int n) 
{
    unsigned int m = n-1;
    unsigned int y = 1;

    while (m != 0) 
    {
		if (m%2 == 1) 
		{
			y =(unsigned int) (( ((XLONG) a)*( (XLONG) y)) % ((XLONG)n ));			//Possibilite overflow
			m = m-1;
		} 
		else 
		{
			unsigned int b = a;
			a = (unsigned int) (( ((XLONG) a)*( (XLONG) a)) % ((XLONG)n));			//Overflow			

			if (a==1 && b!=1 && b!=n-1) 
			{
				  // b est une racine carre non triviale de 1
				  return 1;        // n est composé
			}
			m = m/2;
		}
	}
    
    if (y != 1) 
		return 1;            // n est composé
    else 
		return 0;			// ?

}

/*! 
	\brief Tire un nombre aléatoire entre a et b [inclus]
				
	\param a [in] Limite inferieur
	\param b [in] Limite superieur
  	
	\return La puissance de (0.5)^y
*/
unsigned int irand(unsigned int a, unsigned int b) 
{
	#if defined _WIN32 || defined _WIN64
	unsigned seed2 = time(0);
	std::mt19937 gen(seed2);
	#else
//	boost::random_device gen;		// **chgt IGES**
	std::random_device gen;		// **chgt IGES**
	#endif
	//mt19937 gen(rd());
	//return a + (unsigned int)((double)gen()/(double)rd.max() * (b-a+1));
	//return a + (unsigned int)((double)rd()/(double)rd.max() * (b-a+1));
	return a + (unsigned int)((double)gen()/(double)gen.max() * (b-a+1));
    //return a + (unsigned int)(urand() * (b-a+1));
}


//***************************************************************************/
//
//     	OUTILS DE CLASSEMENT PAR GENERATION
//
//***************************************************************************/


//CLASSEMENT PAR GENERATION 

/*! 
	\brief Associe chaque individu a une generation

	Cette fonction attribue a chaque individu de la genealogie son numero de generation qu il enregistre dans l attribut bFlagSort
	Selon les parametres fournis, cette fonction peut en plus calculer le nombre d'individus par generation
	et/ou un tableau donnant pour chaque generation à partir de 0, la liste des individus qu elle contient
						
	\param Gen [in, out] Genealogie obtenue par LoadGenealogie ou chaque individu connait sa liste d enfants

	\param nbInd [in] Nombre de d individus de la genealogie (un des parametres retour de LoadGenealogie)
  						
	\param tab [out] optionnel, tableau de longeur le nombre de generations de la genealogie ou chaque int est initialise a 0

	\param tabind [out] optionnel, tableau de longeur le nombre de generations de la genealogie ou chaque CIndSimul est initialise a NULL
	<br> en sortie, la ieme case de ce tableau pointe sur un individu de generation i a partir duquel on peut avoir acces aux autres individus de la meme generation
	
	\remark Cette fonction modifie certains attributs :
	<br> bFlagSort	dans tous les cas	contient ne numero de generation (0 pour un proposant)
	<br> pGen		si tabind != NULL	pointe vers le "prochain" individu de la meme generation

	\return 0 
*/

int classeGen(CIndSimul* Gen, int nbInd, int* tab , CIndSimul** tabind)
{
	int i;

	for(i=nbInd-1; i>=0; i--)
	{
		// pour mettre bFlagSort au numero de generation
		if(Gen[i].fils==NULL)
			Gen[i].bFlagSort=0;
		else
		{
			int fmax = 0;
			Clist* tmp = Gen[i].fils;
			while (tmp != NULL)
			{
				int n = tmp->noeud->bFlagSort;
				if (fmax<n)
					fmax = n;
				tmp = tmp->next;
			}
			Gen[i].bFlagSort=fmax +1;	
		}

		// curpro : profondeur courante ie celle de l'individu traité
		const int curpro = Gen[i].bFlagSort;

		// incrémenter le nombre d'individus de la generation courante
		if(tab != NULL)
			tab[curpro]++;

		// ajouter l individu à la liste des individus de la generation courante
		if(tabind != NULL)
		{
			if(tabind[curpro]!=NULL)
				Gen[i].pGen = (double*)tabind[curpro];
			tabind[curpro]=Gen+i;
		}
	}

	return 0;
}


//CLASSEMENT PAR GENERATION MINIMALE

/*! 
	\brief Associe chaque individu a une generation minimale

	Cette fonction attribue a chaque individu de la genealogie son numero de generation minimale qu il enregistre dans l attribut bFlagSort
	Selon les parametres fournis, cette fonction peut en plus calculer le nombre d'individus par generation
	et/ou un tableau donnant pour chaque generation à partir de 0, la liste des individus qu elle contient
						
	\param Gen [in, out] Genealogie obtenue par LoadGenealogie ou chaque individu connait sa liste d enfants

	\param nbInd [in] Nombre de d individus de la genealogie (un des parametres retour de LoadGenealogie)
  						
	\param tab [out] optionnel, tableau de longueur le nombre de generations de la genealogie ou chaque int est initialise a 0

	\param tabind [out] optionnel, tableau de longeur le nombre de generations de la genealogie ou chaque CIndSimul est initialise a NULL
	<br> en sortie, la ieme case de ce tableau pointe sur un individu de generation i a partir duquel on peut avoir acces aux autres individus de la meme generation
	
	\remark Cette fonction modifie certains attributs :
	<br> bFlagSort	dans tous les cas	contient ne numero de generation (0 pour un proposant)
	<br> pGen		si tabind != NULL	pointe vers le "prochain" individu de la meme generation

	\return 0 
*/

int classeGenMin(CIndSimul* Gen, int nbInd, int* tab , CIndSimul** tabind)
{
	int i;

	for(i=nbInd-1; i>=0; i--)
	{
		// pour mettre bFlagSort au numero de generation
		if(Gen[i].fils==NULL)
			Gen[i].bFlagSort=0;
		else
		{
			int fmin = 0;
			Clist* tmp = Gen[i].fils;
			while (tmp != NULL)
			{
				int n = tmp->noeud->bFlagSort;
				if (fmin == 0)
					fmin = n;
				if (fmin>n)
					fmin = n;
				tmp = tmp->next;
			}
			Gen[i].bFlagSort=fmin +1;	
		}

		// curpro : profondeur courante ie celle de l'individu traité
		const int curpro = Gen[i].bFlagSort;

		// incrémenter le nombre d'individus de la generation courante
		if(tab != NULL)
			tab[curpro]++;

		// ajouter l individu à la liste des individus de la generation courante
		if(tabind != NULL)
		{
			if(tabind[curpro]!=NULL)
				Gen[i].pGen = (double*)tabind[curpro];
			tabind[curpro]=Gen+i;
		}
	}

	return 0;
}

//CLASSEMENT PAR GENERATION MOYENNE

/*! 
	\brief Associe chaque individu a une generation moyenne

	Cette fonction attribue a chaque individu de la genealogie son numero de generation moyenne qu il enregistre dans l attribut dFlagSort
						
	\param Gen [in, out] Genealogie obtenue par LoadGenealogie ou chaque individu connait sa liste d enfants

	\param nbInd [in] Nombre de d individus de la genealogie (un des parametres retour de LoadGenealogie)
  						
	\remark Cette fonction modifie certains attributs :
	<br> dFlagSort	dans tous les cas	contient ne numero de generation moyenne (0 pour un proposant)

	\return 0 
*/

int classeGenMoy(CIndSimul* Gen, int nbInd)
{
	int i;

	for(i=nbInd-1; i>=0; i--)
	{
		// pour mettre dFlagSort au numero de generation
		if(Gen[i].fils==NULL)
		{
			Gen[i].dFlagSort=(double)0;
			Gen[i].bFlagSort=1;
		}
		else
		{
			double fmoy = (double)0;
			Clist* tmp = Gen[i].fils;
			int recEnf = 0;
			while (tmp != NULL)
			{
				double n = tmp->noeud->dFlagSort;
				recEnf += tmp->noeud->bFlagSort;
				fmoy += n*tmp->noeud->bFlagSort;
				tmp = tmp->next;
			}
			Gen[i].dFlagSort=(fmoy/recEnf)+1;
			Gen[i].bFlagSort= recEnf;
		}

	}

	return 0;

}



//***************************************************************************/
//
//     	GENEALOGIE EN UN SEUL VECTEUR
//
//***************************************************************************/

///Fonction de comparaison pour QSORT: Compare des CDuoPair par ordre de nom
static int WINCDECL QsortDuoPaircompare(const void *p1, const void *p2)
{
	int i = ((CDuoPair*)p1)->nom;
        int j = ((CDuoPair*)p2)->nom;

        if (i > j)
		return 1;
        else 
		if (i < j)
			return -1;
		else
		        return (0);
}


static int WINCDECL QSORTintcmp(const void* t1,const void* t2)
{
    int i1 = *((int*)t1);
    int i2 = *((int*)t2);
    if (i1==i2)
      return 0;
    else
      if (i1>i2)
	 return 1;
      else
         return -1;
}

static int WINCDECL QsortCIndPereMereCmp(const void *p1, const void *p2)
{
	int i = ((CIndPereMere*)p1)->ind;
    int j = ((CIndPereMere*)p2)->ind;

    if (i > j)
		return 1;
    else 
		if (i < j)
			return -1;
		else
		    return (0);
}


int CompleteGenealogie(int* plIndividu,int* plPere,int* plMere,int* plSexe,int* fInd, int* fpere,int* fmere,int* fsexe,int* pNIndividu)
{
	//Copie du vecteur...petit dans gros
	for(int i=0;i<*pNIndividu;i++) 
	{
		fInd[i] =plIndividu[i];
		fpere[i]=plPere[i];
		fmere[i]=plMere[i];
		if (plSexe)
			fsexe[i]=plSexe[i];
	}

	//Trie des trois vecteurs (croissant) (dans les vecteurs d'origines)
	const int numelem=*pNIndividu;
	qsort(plIndividu,numelem,sizeof(int),QSORTintcmp);
	qsort(plPere,numelem,sizeof(int),QSORTintcmp);
	qsort(plMere,numelem,sizeof(int),QSORTintcmp);

	//Combine les 3 vecteurs et les completes
	int posind=0;
	int pospere=0;
	int posmere=0;
	int minimum=0;
	//Boucle a l'interieur de pere et mere
	while(posmere < numelem || pospere < numelem)
	{
		//Positionne au premiers element plus grand que l'element courant
		int SPot = 2; //Sexe du candidat
		while(posmere < numelem && plMere[posmere]<=minimum)
		      ++posmere;		
		while( pospere < numelem && plPere[pospere]<=minimum)
		      ++pospere;
		minimum=MIN(
			    (posmere<numelem?plMere[posmere]:INT_MAX),
			    (pospere<numelem?plPere[pospere]:INT_MAX) );
		if(posmere<numelem)SPot=2;
		if(pospere<numelem)SPot=1;
		//Avance l'element du vecteur individu
		while(posind < numelem && plIndividu[posind]<minimum)
		      ++posind;
		const int valind = (posind<numelem?plIndividu[posind]:INT_MAX);
		if (minimum<valind)  // Si egale, il est deja dans la liste
		{
			  //Pour chaque ajout;
			fInd[*pNIndividu]  =minimum;
			fpere[*pNIndividu] =0;
			fmere[*pNIndividu] =0;
			
			//Determine le sexe du nouveau candidat
			if (plSexe)
				fsexe[*pNIndividu]=SPot;
			++(*pNIndividu);
		}			       		
	}	
	return 0;
}

/*! 
	\brief Cette fonction creer un vecteur Genealogie (version 5) a partir d'une genealogie No Individu,No pere,No mere

	La fonction convertis une genealogie a 3 vecteur en une a un vecteur.
	le resultat est retourne dans le vecteur saveptr dont la taille doit avoir été calculée au prealabe
	a l'aide de TAILLEGENVERSION3
    	  
	<br><br>FORMAT DU FICHIER DE SAUVEGARDE
	<table>
	<tr><td>INDICE</TD><TD>VALEUR</TD></TR>
	<tr><td>0</TD><TD>'G'</TD></TR>
	<tr><td>1</TD><TD>'E'</TD></TR>
	<tr><td>2</TD><TD>'N'</TD></TR>
	<tr><td>3</TD><TD>5 //Version de la Genealogie</TD></TR>
	<tr><td>4-7</TD><TD>Signature MD5</TD></TR>
	<tr><td>8</TD><TD>Nombre d'individu</TD></TR>
	<tr><td>9</TD><TD>Nombre total d'enfant</TD></TR>
	<tr><td>10</TD><TD>Profondeur maximale</TD></TR>
	<tr><td>11 --> (loin)</TD><TD>Pour chaque individu dans l'ordre ORDONNE selon OrdonneStructure
	<br>On inscrit
	<br> No de Individu
	<br> Indice du pere (selon la genealogie ORDONNE)
	<br> Indice du mere (selon la genealogie ORDONNE)
	<br> Le nombre d'enfant direct de cet individu (0 si aucun)
	<br> L'indice de chaque enfant (selon la genealogie ORDONNE)</TD></TR>
	<tr><td></TD><TD>Vecteur de recherche binaire
	<br>Pour chaque individu par ordre croissant de No individu
	<br>Indice de l'individu</TD></TR>
	</table>
				
	\param plIndividu	[in] Vecteur taille lNIndividu, No(etiquette) de l'individu
	\param plPere		[in] Nombre entier a tester,	No du pere correspondant a l'individu
	\param plMere		[in] Nombre entier a tester,	No de la mere correspondant a l'individu
	\param plSexe		[in] Nombre entier a tester,	No correspondant au sexe de l'individu (1 ou 2)				

	\param lNIndividu [in] Nombre d'individu dans la genealogie

	\retval saveptr [out] Vecteur de int de taille calculer par la macro TAILLEGENVERSION3
						 Si Succes, contient la genealogie a un vecteur
  	
	\return 0 dans le cas contraire retourne un code d'erreur correspondant à sCodeErreurCreation

	\sa sCodeErreurCreation
*/
int CreerGenealogie(int* plIndividu, int* plPere, int* plMere, int* plSexe, int lNIndividu, int* saveptr)
{
	try{
	 //TAILLE GENEALOGIE    
		//ID(3 byte)+VERSION(1)+MD5(4)+NINDIVIDU+Nombred'enfant+
		// (colonne individu, pere , mere) + (nombre d'enfant) + listeindiceenfant+
		// Indexdespositionparordrecroissantdenoind+ (verification) + correction
	//const long TAILLESAUVEGARDE=3+1+1+1+lNIndividu*3+lNIndividu+NombreEnfant+lNIndividu+1; 
	
     //VARIABLE OPERATIONNEL	
	int NombreEnfant=0;
	int curseur;
	//NNiveau nombre de ligne dans le tableau
	INITGESTIONMEMOIRE;
	CIndSimul*  Noeud	=(CIndSimul*)  memalloc(lNIndividu,sizeof(CIndSimul)); 
	CIndSimul **Ordre	=(CIndSimul**)	memalloc(lNIndividu,sizeof(CIndSimul*));
	CDuoPair *Trie		=(CDuoPair*)	memalloc(lNIndividu,sizeof(CDuoPair));
	int i;
	
	//VALIDATION
	for(i=0;i<lNIndividu;i++) 
	{
	    Trie[i].nom=plIndividu[i]; 
	    Trie[i].pos=i; 	    

	    //VERIFICATION AU PASSAGE SI LES NO PERE ET NO MERE SONT VALIDE
	    if (plIndividu[i]<=0){
//			GENError("The index of an individual must be greater than zero.");
			throw std::range_error("The index of an individual must be greater than zero.");
			//GENError("L'indice d'un individu doit-être plus grand que zero");		
		}
		if (plPere[i]<0){
//			GENError("The father of individual %d must be greater than or equal to zero",plIndividu[i]);
			char erreur[TAILLEDESCRIPTION];
			sprintf(erreur, "The father of individual %d must be greater than or equal to zero",plIndividu[i]);
			throw std::range_error(erreur);
			//GENError("Le pere de l'individu %d doit être plus grand ou égal à zero",plIndividu[i]);			
		}
	    if (plMere[i]<0){
//			GENError("The mother of individual %d must be greater than or equal to zero",plIndividu[i]);
			char erreur[TAILLEDESCRIPTION];
			sprintf(erreur, "The mother of individual %d must be greater than or equal to zero",plIndividu[i]);
			throw std::range_error(erreur);
			//GENError("La mere de l'individu %d doit être plus grand ou égal à zero",plIndividu[i]);
		}
	    if (plPere[i]==plMere[i] && plPere[i]!=0){
//			GENError("Individual %d must have different mother and father",plIndividu[i]);
			char erreur[TAILLEDESCRIPTION];
			sprintf(erreur, "Individual %d must have different mother and father",plIndividu[i]);
			throw std::range_error(erreur);
			//GENError("L'individu %d doit avoir un pere et une mere différent",plIndividu[i]);
		}
	}
	
	//Partie 1: CREATION ET INITIALISATION DE LA STRUCTURE NOEUD ET L'ORDRE (SLOW)
	Clist* childrenArray=NULL;
	CreeStructure(Noeud,plIndividu,plPere,plMere,plSexe,lNIndividu,&NombreEnfant,Trie,&childrenArray); //Creer enfants
	//VALIDATION DU SEXE  0=INCONNU 1=M 2=F
	int oldnumber=-1;
	for(i=0;i<lNIndividu;i++) 
	{		
		if (Noeud[i].sex<GEN_INCONNU || Noeud[i].sex>GEN_FEM){
//				GENError("The sexe of individual %d is not valid (0=SEXE UNKNOWN, 1=MAN, 2=WOMAN)",plIndividu[i]);
				char erreur[TAILLEDESCRIPTION];
				sprintf(erreur, "The sexe of individual %d is not valid (0=SEXE UNKNOWN, 1=MAN, 2=WOMAN)",plIndividu[i]);
				throw std::range_error(erreur);
				//GENError("Le sexe de l'individu %d n'est pas une valeur valide (0=SEXE INCONNU, 1=HOMME, 2=FEMME)",plIndividu[i]);
		}
		if (Noeud[i].pere!=NULL) 
		{
			if (Noeud[Noeud[i].pere->noind].sex==GEN_FEM){
//				GENError("Individual %d is both mother and father to two different individuals\n\n",Noeud[i].pere->nom);
				char erreur[TAILLEDESCRIPTION];
				sprintf(erreur, "Individual %d is both mother and father to two different individuals\n\n",Noeud[i].pere->nom);
				throw std::range_error(erreur);
				//GENError("L'individu %d est a la fois un pere et une mere pour deux individu different\n\n",Noeud[i].pere->nom);
			}
			else
				Noeud[Noeud[i].pere->noind].sex=GEN_MASC; 
		}
		if (Noeud[i].mere!=NULL) 
		{
			if (Noeud[Noeud[i].mere->noind].sex==GEN_MASC){
//				GENError("Individual %d is both mother and father to two different individuals\n\n",Noeud[i].mere->nom);
				char erreur[TAILLEDESCRIPTION];
				sprintf(erreur, "Individual %d is both mother and father to two different individuals\n\n",Noeud[i].mere->nom);
				throw std::range_error(erreur);
				//GENError("L'individu %d est a la fois un pere et une mere pour deux individu different\n\n",Noeud[i].mere->nom);
			}
			else
				Noeud[Noeud[i].mere->noind].sex=GEN_FEM;
		}		
		//Test dupplicata individu
		if (Trie[i].nom==oldnumber){
//			GENError("Individual %d is duplicated in the genealogy",Trie[i].nom);
			char erreur[TAILLEDESCRIPTION];
			sprintf(erreur, "Individual %d is duplicated in the genealogy",Trie[i].nom);
			throw std::range_error(erreur);
			//GENError("L'individu %d est duppliqué dans la généalogie",Trie[i].nom);			
		}
		oldnumber=Trie[i].nom;
	}
	//VALIDATION DU SEXE (COMPLETE SI C'EST POSSIBLE)		
	int NIndMasc = (plSexe?0:-1);
	if(plSexe)
	{		
		for(i=0;i<lNIndividu;i++) 		
		{
			if (Noeud[i].sex==GEN_INCONNU){
//				GENError("The sexe of individual %d is unknown and must be set",plIndividu[i]);
				char erreur[TAILLEDESCRIPTION];
				sprintf(erreur, "The sexe of individual %d is unknown and must be set",plIndividu[i]);
				throw std::range_error(erreur);
				//GENError("Le sexe de l'individu %d est inconnu et doit-être fournie",plIndividu[i]);
			}
			//Compte le nombre d'individu masculin
			if (Noeud[i].sex==GEN_MASC)
				++NIndMasc;
		}	
	}

	//Trie de la genealogie selon parent toujours avant enfant
	int profondeurMax;	
	OrdonneStructure(Noeud,Ordre,lNIndividu,0/*SENSINVERSE*/,&profondeurMax); 

	//SAUVEGARDE ELLE MEME
	saveptr[0]='G'; 
	saveptr[1]='E'; 
	saveptr[2]='N'; 
	saveptr[3]=VERSIONGENEALOGIE;	//VERSION  
	//saveptr[4-7]= //CONSERVE POUR md5 plus tard
	saveptr[GENBIN_NIND_OFFSET]=lNIndividu;		//NOMBRE d'INDIVIDU 
	saveptr[GENBIN_NENFANT_OFFSET]=NombreEnfant;//NOMBRE d'ENFANT 
	saveptr[GENBIN_PROFMAX]=profondeurMax;		//PROFONDEUR MAXIMALE
	saveptr[GENBIN_NINDMASC]=NIndMasc;			//NOMBRE D'HOMME 

	//FAIRE LA SAUVEGARDE EN TENANT COMPTE DE LA PRIORITE
	curseur=GENBIN_NBIND_OFFSET;  //indice du curseur 
	for(i=0;i<lNIndividu;i++) 
	{        
	    //Sauvegarde donne principale
	    saveptr[curseur++]=Ordre[i]->nom; 
	    if (Ordre[i]->pere!=NULL) 
		saveptr[curseur++]=Ordre[i]->pere->noind; 
	    else 
		saveptr[curseur++]=-1; 
	    if (Ordre[i]->mere!=NULL) 
		saveptr[curseur++]=Ordre[i]->mere->noind; 
	    else 
		saveptr[curseur++]=-1; 
	    
		saveptr[curseur++]=Ordre[i]->sex; //Sexe de l'individu

	    //sauvegarde enfant
	    int indiceNEnfant=curseur++; 
	    int count=0;	 
	    Clist *current=Ordre[i]->fils; 
	    if (current!=NULL) 
	    { 
		do 
		{			 
			++count; 
			saveptr[curseur++]=current->noeud->noind; 
			current=current->next;				 
		} 
		while(current!=NULL); 
	    } 
	    saveptr[indiceNEnfant]=count; 
	} 

	//INDEX DE RECHERCHE
	for(i=0;i<lNIndividu;i++) 
		saveptr[curseur++]=Noeud[Trie[i].pos].noind; 	

	//FINAL
	saveptr[curseur++]=99999999; 	
	//NETTOYER LES ASSIGNEMENTS DE MÉMOIRES
	DetruireStructure(childrenArray);

	//GENERATION DE LA CHAINE POUR MD5
	const int TailleGen = TAILLEGENVERSION7(lNIndividu,NombreEnfant);	
	const int  ratio = sizeof(int)/sizeof(unsigned char);
	const int tailleGenInString = ((TailleGen-GENBIN_NIND_OFFSET)*ratio);
	const int taillestring = tailleGenInString+1024-(tailleGenInString%1024);
	unsigned char *md5string =(unsigned char*) memalloc(taillestring,sizeof(unsigned char));

	//Copie de la structure actuel dans la chaine MD5
	int md5cur=0;
	for(i=GENBIN_NIND_OFFSET;i<TailleGen;i++)
	{
		unsigned int val = (unsigned int) saveptr[i];
		md5string[md5cur++]=(unsigned char) (((0x000000FF&val)>>0 ) &0xFF);
		val = (unsigned int) saveptr[i];
		md5string[md5cur++]=(unsigned char) (((0x0000FF00&val)>>8 ) &0xFF);
		val = (unsigned int) saveptr[i];
		md5string[md5cur++]=(unsigned char) (((0x00FF0000&val)>>16) &0xFF);
		val = (unsigned int) saveptr[i];
		md5string[md5cur++]=(unsigned char) (((0xFF000000&val)>>24) &0xFF);
	}
	while(md5cur<taillestring)
		md5string[md5cur++]=0;
	
	//Hashage
	unsigned char digest[16];
	MD5_CTX context;
	md5_starts(&context);
	md5cur=0;
	while(md5cur<taillestring)		
	{
		md5_update(&context, md5string+md5cur, 64);
		md5cur+=64;
	}
	md5_finish(&context,digest);

	//return 10; //TailleGen;
	//Sauvegarde du résultat	
	saveptr[GENBIN_MD5_1]= (int) (
				(digest[ 0] << 0)  |
				(digest[ 1] << 8)  |
				(digest[ 2] << 16) |
				(digest[ 3] << 24));
	saveptr[GENBIN_MD5_2]= (int) (
				(digest[ 4] << 0)  |
				(digest[ 5] << 8)  |
				(digest[ 6] << 16) |
				(digest[ 7] << 24));
		
	saveptr[GENBIN_MD5_3]= (int) (
				(digest[ 8] << 0)  |
				(digest[ 9] << 8)  |
				(digest[10] << 16) |
				(digest[11] << 24));

	saveptr[GENBIN_MD5_4]= (int) (
				(digest[12] << 0)  |
				(digest[13] << 8)  |
				(digest[14] << 16) |
				(digest[15] << 24));					

	return 0; 
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


/*! 
	\brief Valide et retourne la taille d'une genealogie
	
	\param Genealogie [in] Une genealogie construite à l'aide de gen.genealogie

	\retval nenfant [out] Si !=NULL, alors retourne le nombre total d'enfant de la genealogie
  	
	\return Si succes, retourne le nombre d'individu de la genealogie
			dans le cas contraire retourne une valeur negative
*/
int LengthGenealogie(int* Genealogie, int* nenfant, int* nprofmax, int* nindmasc)
{
	try{
	//validation
	if (Genealogie[0]!='G' || Genealogie[1]!='E' || Genealogie[2]!='N' ){
//	  GENError("\nError: invalid genealogy given. Create one using gen.genealogie(ind,father,mother).");
	  throw std::range_error("Error: invalid genealogy given. Create one using gen.genealogie(ind,father,mother).");
	  //GENError("\nErreur: La généalogie fournie n'est pas valide. Créer en une à l'aide de gen.genealogie(ind,pere,mere).");
	}
	if (Genealogie[3]!=VERSIONGENEALOGIE){
//	  GENError("\nError: Given genealogy is not from current version.");
	  throw std::range_error("Error: Given genealogy is not from current version.");
	  //GENError("\nErreur: La généalogie fournie n'est pas de la version courante.");
	}
	//Initialisation
	const int iNind=Genealogie[GENBIN_NIND_OFFSET];

	if (nenfant!=NULL)
		*nenfant=Genealogie[GENBIN_NENFANT_OFFSET];
	if (nprofmax!=NULL)
		*nprofmax=Genealogie[GENBIN_PROFMAX];
	if (nindmasc!=NULL)
		*nindmasc=Genealogie[GENBIN_NINDMASC];

	return iNind;
 	} catch(std::exception &ex) {
 		forward_exception_to_r(ex);
 	} catch(...){
 		::Rf_error("c++ exception (unknown reason)"); 
 	} 
 	return 0;
}



int ValidateGenealogie(int* Genealogie)
{
	int nenfant;
	int nind = LengthGenealogie(Genealogie,&nenfant);  //Test pour certaine erreur

	//VALIDATION DE LA SIGNATURE MD5
	//GENERATION DE LA CHAINE POUR MD5
	INITGESTIONMEMOIRE;
	const int TailleGen = TAILLEGENVERSION7(nind,nenfant);	
	const int  ratio = sizeof(int)/sizeof(unsigned char);
	const int tailleGenInString = ((TailleGen-GENBIN_NIND_OFFSET)*ratio);
	const int taillestring = tailleGenInString+1024-(tailleGenInString%1024);
	unsigned char *md5string =(unsigned char*) memalloc(taillestring,sizeof(unsigned char));

	//Copie de la structure actuel dans la chaine MD5	
	int md5cur=0;
	for(int i=GENBIN_NIND_OFFSET;i<TailleGen;i++)
	{
		unsigned int val = (unsigned int) Genealogie[i];
		md5string[md5cur++]=(unsigned char) (((0x000000FF&val)>>0 ) &0xFF);
		val = (unsigned int) Genealogie[i];
		md5string[md5cur++]=(unsigned char) (((0x0000FF00&val)>>8 ) &0xFF);
		val = (unsigned int) Genealogie[i];
		md5string[md5cur++]=(unsigned char) (((0x00FF0000&val)>>16) &0xFF);
		val = (unsigned int) Genealogie[i];
		md5string[md5cur++]=(unsigned char) (((0xFF000000&val)>>24) &0xFF);		
	}
	while(md5cur<taillestring)
		md5string[md5cur++]=0;
	
	//Hashage
	unsigned char digest[16];
	MD5_CTX context;
	md5_starts(&context);
	md5cur=0;
	while(md5cur<taillestring)		
	{
		md5_update(&context, md5string+md5cur, 64);
		md5cur+=64;
	}
	md5_finish(&context,digest);

	//Comparaison du résultat
	int Md5_1= (int) (
				(digest[ 0] << 0)  |
				(digest[ 1] << 8)  |
				(digest[ 2] << 16) |
				(digest[ 3] << 24));
	int Md5_2= (int) (
				(digest[ 4] << 0)  |
				(digest[ 5] << 8)  |
				(digest[ 6] << 16) |
				(digest[ 7] << 24));
		
	int Md5_3= (int) (
				(digest[ 8] << 0)  |
				(digest[ 9] << 8)  |
				(digest[10] << 16) |
				(digest[11] << 24));

	int Md5_4= (int) (
				(digest[12] << 0)  |
				(digest[13] << 8)  |
				(digest[14] << 16) |
				(digest[15] << 24));
	if (		Genealogie[GENBIN_MD5_1]!=Md5_1 ||
			Genealogie[GENBIN_MD5_2]!=Md5_2 ||
			Genealogie[GENBIN_MD5_3]!=Md5_3 ||
			Genealogie[GENBIN_MD5_4]!=Md5_4 )
			return 0;   //Genealogie invalide
	else
			return 1;   //Genealogie valide
}


/*! 
	\brief 	Initialise un vecteur de CIndSimul pour qu'il represente une genealogie 
		
	  A partir d'une genealogie a un vecteur construite à l'aide de gen.genealogie.

	<br>Il faut fournir 2 vecteur a cette fonction.
	<br>Un vecteur de CIndSimul
	<br>Et un vecteur de Clist (optionnel)
	<br> La taille requise pour ces deux vecteur pour etre déterminé a l'aide de la fonction LengthGenealogie
	
	\param Genealogie [in] Une genealogie construite à l'aide de gen.genealogie

	\retval Noeud [out] Un vecteur CIndSimul de taille fournit par LengthGenealogie
					   Si succes, les noeuds sont initialiser de maniere a ce que le vecteur Noeud
					   represente la genealogie Genealogie
	
	\retval Children [out] Un vecteur de CList taille founit par  LengthGenealogie
							Si succes, le vecteur est utiliser pour genere les liste d'enfant pour chaque Noeud
							Si NULL, alors la liste fils de chaque NOEUD ne sera PAS initialisé.
	  	
	\retval IndexRecherche [out] pointeur vers un pointeur de int.
							      Si succes, le pointeur pointe vers le premier élément de la liste ordonné de
								  no etiquette de la genealogie
  	
	\return Si succes, retourne le nombre d'individu de la genealogie
			dans le cas contraire retourne une valeur negative

	\remark La raison de se faire fourni au depart des vecteurs de bonne taille permet d'utilisé le "gargage collecteur"
			facilement a partir de la procédure appelante.

	\sa LengthGenealogie
*/
int ReCreeStructure(int* Genealogie,CIndSimul* Noeud, Clist* Children, int** IndexRecherche)
{		
	int i;
	Clist **current = NULL;	
	//Clist *liste    = NULL;
	int nenfant=0;
	
	//Initialisation
	const int iNind=LengthGenealogie(Genealogie, &nenfant);

	for(i=0;i<iNind;i++)
	{		
		Noeud[i].pGen   = NULL;
		Noeud[i].etat   = GENNONEXPLORER;
		Noeud[i].noind  = i;
		Noeud[i].allele = 0;									
	}

	//Creation d'un tas de liste suffisant pour les enfants
	int curseur = GENBIN_NBIND_OFFSET;
	for(i=0; i<iNind; i++)
	{
		//Sauvegarde donne principale
		Noeud[i].nom=Genealogie[curseur++];
		if (Genealogie[curseur]!=-1)
			Noeud[i].pere=Noeud+Genealogie[curseur++];
		else
			{Noeud[i].pere=NULL;++curseur;}

		if (Genealogie[curseur]!=-1)
			Noeud[i].mere=Noeud+Genealogie[curseur++];
		else
			{Noeud[i].mere=NULL;++curseur;}

		Noeud[i].sex=(sex_t) Genealogie[curseur++];  //Sexe de l'individu
		//RECUPERATION ENFANT

		const int NEnfant=Genealogie[curseur++];
		
		if(Children!=NULL)
		{
			current = &(Noeud[i].fils);
			for(int a=0; a<NEnfant; a++)
			{
				*current = (Children++);		
				(*current)->noeud = Noeud + Genealogie[curseur++];
				current = &( (*current)->next );
			}
			*current=NULL;
		}
		else
		{
			//NE PAS CHARGE L'ENFANT
			Noeud[i].fils=NULL;
			curseur+=NEnfant;
		}
	}

	//Index de position pour la recherche
	if (IndexRecherche!=NULL)
	  *IndexRecherche=Genealogie+(curseur++);

	return 0;
}


/*! 
	\brief Retrouve rapidement l'indice d'un individu dans une genealogie

	Execute une recherche du no d'individu a partir de l'index de recherche de la genealogie.
				
	\param nom [in] No d'individu a recherche

	\param Noeud	 [in]  Adresse d'un vecteur de CIndSimul representant une genealogie

	\param IndexRecherche [in] Pointeur de int qui represente l'index de recherche.
								Ce pointeur est fourni par la fonction ReCreeStructure

	\param iNind [in] Nombre d'individu dans la genealogie (Taille du vecteur Noeud)
  	
	\return l'indice de nom dans la genealogie. -1 s'il n'est pas trouvé
	
	 \sa ReCreeStructure
*/
int ReTrouverIndiceStructure(int nom, CIndSimul* Noeud, int* IndexRecherche, int iNind)
{
	//Recherche binaire	
	int uplimit=iNind;
	int downlimit=-1;	
	int pivot=iNind/2;
       
	while(1)
	{
		const int Value=Noeud[IndexRecherche[pivot]].nom;
		
		if(Value==nom)
		{//Resultat trouver
			return IndexRecherche[pivot];
		}
		else
		{
			if (Value<nom)
			{//La valeur est plus petite..
				downlimit=pivot;
				pivot=(uplimit+pivot)/2;
			
				if (pivot==downlimit)
					return -1; //Valeur pas trouver
			
			}
			else
			{//La valeur est plus grande
				uplimit=pivot;
				pivot=(downlimit+pivot)/2;

				if (pivot==uplimit)
					return -1; //Valeur pas trouver
			}
		}
	}
}




void SortGenealogie3Vecteur(int* ind,int *pere, int* mere, int* sex, int nind)
{
	INITGESTIONMEMOIRE;
	CIndPereMere* vInd = (CIndPereMere*) memalloc(nind,sizeof(CIndPereMere));

	//Copie Vecteur -> Structure
	for(int a=0;a<nind;a++)
	{
		vInd[a].ind=ind[a];
		vInd[a].pere=pere[a];
		vInd[a].mere=mere[a];
		vInd[a].sex=(sex?sex[a]:0);
	}
	qsort(vInd,nind,sizeof(CIndPereMere),QsortCIndPereMereCmp);
	//Copie Structure -> Structure
	for(int a=0;a<nind;a++)
	{
		ind[a]=vInd[a].ind;
		pere[a]=vInd[a].pere;
		mere[a]=vInd[a].mere;
		if (sex)
			sex[a]=vInd[a].sex;
	}	
	return;
}

//***************************************************************************/
//
//     	FONCTION POUR LOADER ET CACHE GENEALOGIE,PROPOSANT & ANCETRE
//
//***************************************************************************/

//CACHE GENEALOGIE
int		g_CacheMD5Sign[4]	={0,0,0,0};
CIndSimul*  g_CacheGenArray	=NULL;		//Contient un tableau assigné
Clist*	g_CacheChildList	=NULL;		//Contient un tableau assigné
int*	g_CacheRecherche	=NULL;		//Contient un tableau assigné
int		g_CacheNInd			=0;
int		g_CacheProfMax		=0;
int		g_CacheNIndMasc		=0;		//Nombre d'individu Masculin

//CACHE PROPOSANT & ANCETRE
//enum ENUMBANQUE {PROPOSANT,ANCETRE};
CIndSimul** g_CacheVec[2]	={NULL,NULL};	//Contient deux tableau assigné
int	   g_CacheVecInd[2]	={-1,-1};

//CACHE GROUPE PROPOSANT & ANCETRE
int			g_CacheNbGroupe[2]	={-1,-1};		//Nombre total de groupe dans chaque banque
int*			g_CacheGrVecInd[2]	={NULL,NULL};	//Tableau nb individu dans groupe
CIndSimul***	g_CacheGroup[2]	={NULL,NULL};	//Tableau vers le début de chaque groupe

int LoadGenealogie(int* Genealogie,int loadChildren,int* NInd, CIndSimul **Noeudarr, int** IndexRecherche) 
{
	try{
/*	//Compare indice MD5 et la présence d'enfant
#ifndef USESPLUSALLOC
	//Detruit la genealogie précédente
	if(Genealogie[GENBIN_MD5_1]==g_CacheMD5Sign[0] && Genealogie[GENBIN_MD5_2]==g_CacheMD5Sign[1] &&
	   Genealogie[GENBIN_MD5_3]==g_CacheMD5Sign[2] && Genealogie[GENBIN_MD5_4]==g_CacheMD5Sign[3] &&
	   (loadChildren==GFALSE || g_CacheChildList) )				
	{	//C'est la meme genealogie pas besoin de reloader..	
#ifdef USESDEBUG		
		//printf("\nLa genealogie a été charge depuis la cache\n"); //#
#endif 			
		*NInd = g_CacheNInd;
		if (IndexRecherche)
			*IndexRecherche=g_CacheRecherche;
		*Noeudarr= g_CacheGenArray;
		return 1; 
	}	
#endif	
*/	//Dans le cas contraire, on flush la genealogie
	FlushGenealogie();

	//Validation genealogie et recupere nombre d'enfant et profondeurmax
	int nenfant;
	g_CacheNInd=LengthGenealogie(Genealogie,&nenfant,&g_CacheProfMax,&g_CacheNIndMasc);
	
	//Creation du tableau de Noeud
	g_CacheGenArray = (CIndSimul*)memallocIN(g_CacheNInd,sizeof(CIndSimul));	
	if (!g_CacheGenArray)
	{
		FlushGenealogie();
//		GENError("Not enough memory to load genealogy.");
		throw std::range_error("Not enough memory to load genealogy.");
		//GENError("Il n'y a pas assez de mémoire disponible pour charger la généalogie.");
	}	
	if (loadChildren)
	{
		g_CacheChildList=(Clist*) memallocIN(nenfant,sizeof(Clist));
		if (!g_CacheChildList)
		{
			FlushGenealogie();
//			GENError("Not enough memory to load genealogy.");
			throw std::range_error("Not enough memory to load genealogy.");
			//GENError("Il n'y a pas assez de mémoire disponible pour charger la généalogie.");
		}	
	}
	
	//Partie 1: CREATION ET INITIALISATION DE LA STRUCTURE NOEUD
	int* VecteurRecherche;
	if (loadChildren)
		ReCreeStructure(Genealogie, g_CacheGenArray, g_CacheChildList, &VecteurRecherche); //Avec enfant
	else
		ReCreeStructure(Genealogie, g_CacheGenArray, NULL		    , &VecteurRecherche);	//Sans enfant
	
	//Recupere le vecteur de recherche
	g_CacheRecherche=(int*) memallocIN(g_CacheNInd,sizeof(int));
	if (!g_CacheRecherche)
	{
		FlushGenealogie();
//		GENError("Not enough memory to load genealogy.");
		throw std::range_error("Not enough memory to load genealogy.");
		//GENError("Il n'y a pas assez de mémoire disponible pour charger la généalogie.");
	}	
	memcpy(g_CacheRecherche,VecteurRecherche,g_CacheNInd*sizeof(int));

	//Retourne la structure créer
	*NInd = g_CacheNInd;
	if (IndexRecherche)
		*IndexRecherche=g_CacheRecherche;
	*Noeudarr= g_CacheGenArray;
	
	//Sauvegarde la signature MD5
	g_CacheMD5Sign[0]=Genealogie[GENBIN_MD5_1];
	g_CacheMD5Sign[1]=Genealogie[GENBIN_MD5_2];
	g_CacheMD5Sign[2]=Genealogie[GENBIN_MD5_3];
	g_CacheMD5Sign[3]=Genealogie[GENBIN_MD5_4];
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}

int FlushGenealogie()
{	
	//Detruit la genealogie précédente
	g_CacheMD5Sign[0]=g_CacheMD5Sign[1]=g_CacheMD5Sign[2]=g_CacheMD5Sign[3]=0;	
	g_CacheNInd			=0;
	g_CacheProfMax		=0;	
	if (g_CacheGenArray)
	{
		memfreeIN( g_CacheGenArray );
		g_CacheGenArray=NULL;
	}

	if (g_CacheChildList)
	{
		memfreeIN(g_CacheChildList);
		g_CacheChildList=NULL;
	}
	if (g_CacheRecherche)
	{
		memfreeIN(g_CacheRecherche);
		g_CacheRecherche=NULL;
	}
	//Flush Cache proposant & ancetre....
	FlushProposantAncetre(PROPOSANT);
	FlushProposantAncetre(ANCETRE);

	//Detruit les groupes de proposants et ancetre

	return 0;
}

/**ATTENTION UTILISE LA GENEALOGIE PRESENTEMENT EN CACHE**/
int LoadProposant(int* Proposant, int nbProposant, CIndSimul*** pro) {return LoadVec(PROPOSANT, Proposant, nbProposant, pro);}

int LoadAncetre(int* Ancetre,int nbAncetre,CIndSimul*** anc) 	    {return LoadVec(ANCETRE,Ancetre,nbAncetre,anc);}

static int LoadVec(ENUMBANQUE banque,int* vec, int nb, CIndSimul*** NproAnc)
{
	try{
	const char* type[]={"proband","ancetre"};

	if (!g_CacheGenArray)
	{
		FlushProposantAncetre(banque);
//		GENError("Invalid use of LoadProposant or LoadAncetre function: start by LoadGenealogie");
		throw std::range_error("Invalid use of LoadProposant or LoadAncetre function: start by LoadGenealogie");
		//GENError("Utilisation invalide de la fonction LoadProposant ou LoadAncetre : il faut utiliser LoadGenealogie avant.");
	}

	//Bon.. creer le vecteur s'il est trop petit..
	CIndSimul** canRestaure = g_CacheVec[banque]; //Pointe sur le premier element (s'il existe)
	
	if(nb > g_CacheVecInd[banque])
	{
		FlushProposantAncetre(banque);
		g_CacheVec[banque] = (CIndSimul**)	memallocIN( nb, sizeof(CIndSimul*));
		if (!g_CacheVec[banque])
		{
			FlushGenealogie();
//			GENError("Not enough memory to load the %s.",type[banque]);
			char erreur[TAILLEDESCRIPTION];
			sprintf(erreur, "Not enough memory to load the %s.",type[banque]);
			throw std::range_error(erreur);
			//GENError("Il n'y a pas assez de mémoire disponible pour charger les %s.",type[banque]);
		}	
		g_CacheVecInd[banque]=nb;
		canRestaure = NULL;
		g_CacheVec[banque][0] = g_CacheGenArray;
		//return (long)(g_CacheVec[banque][0]);
	}
	//Boucle pour chaque proposant ou ancetre
	int count=0; //# prob ici ...
	for(int i=0; i<nb; i++)
	{	
		if (DescIEEEValue(vec+i)!=NULL)
		{
			FlushProposantAncetre(banque);
//			if(i>0) return 1;
//			GENError("Special IEEE caracter %s is not a valid %s", DescIEEEValue(vec+i),type[banque]);
			char erreur[TAILLEDESCRIPTION];
			sprintf(erreur, "Special IEEE caracter %s is not a valid %s", DescIEEEValue(vec+i),type[banque]);
			throw std::range_error(erreur);
			//GENError("Le caractère spécial IEEE %s n'est pas un %s valide",DescIEEEValue(vec+i),type[banque]);
		}

		//Tentative de recherche dans la cache		
		if (canRestaure && (*canRestaure)->nom==vec[i])
		{
			g_CacheVec[banque][i]=*(canRestaure++);
			count++; //#
		}
		else
		{
			int tmpindice = ReTrouverIndiceStructure(vec[i], g_CacheGenArray, g_CacheRecherche, g_CacheNInd);
			if (tmpindice==-1)
			{
				FlushProposantAncetre(banque);
				//printf("on retourne -1\n");
//				GENError("%s %d is not included in the genealogy ...",type[banque],vec[i]);
				char erreur[TAILLEDESCRIPTION];
				sprintf(erreur, "%s %d is not included in the genealogy ...",type[banque],vec[i]);
				throw std::range_error(erreur);
				//return -1;
				//printf("on devrait pas lire ça...\n");
				//GENError("Le %s %d n'est pas inclus dans la généalogie",type[banque],vec[i]);
			}
			//g_CacheVec[banque][i] = &(g_CacheGenArray[tmpindice]);
			//g_CacheVec[banque][i] = (CIndSimul*) memallocIN( 1, sizeof(CIndSimul));  (g_CacheGenArray + tmpindice) = 51426104
			g_CacheVec[banque][i] = g_CacheGenArray + tmpindice;
			//return (long)(g_CacheVec[banque][i]);
		}
	}
	
	//Retourne le vecteur de proposant
	*NproAnc = g_CacheVec[banque];
#ifdef USESDEBUG
	//printf("\nLoad %s from cache = %d\n",type[banque],count); //#
#endif
	return 0;
	 		} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


static int FlushProposantAncetre(ENUMBANQUE banque)
{
	if (g_CacheVec[banque])
		memfreeIN(g_CacheVec[banque]);
	g_CacheVec[banque]=NULL;
	g_CacheVecInd[banque]=-1;
	//Flush les groupe associés
	FlushGroupeProposantAncetre(banque);
	return 0;
}

//RECUPERE INFORMATION UTILITAIRE
///Recupere la profondeur maximale de la genealogie dans la cache
int LoadProfondeurMax() {return g_CacheProfMax;}
int LoadNIndMasc()		{return g_CacheNIndMasc;}


/**ATTENTION UTILISE LA GENEALOGIE PRESENTEMENT EN CACHE**/
int LoadGroupeProposant(int* Proposant, int* BorneGr, int nbGroupe, CIndSimul**** GRpro,int** nIndGr)
   {return LoadVecGroupe(PROPOSANT,BorneGr,nbGroupe,GRpro,nIndGr);}
int LoadGroupeAncetre  (int* Ancetre  , int* BorneGr, int nbGroupe, CIndSimul**** GRanc,int** nIndGr)
   {return LoadVecGroupe(ANCETRE,BorneGr,nbGroupe,GRanc,nIndGr);}

static int LoadVecGroupe(ENUMBANQUE banque,int* BorneGr,int nbGroupe, CIndSimul**** GrProAnc,int** nIndGr)
{	
	try{
	//Verifie si une généalogie a déjà été charge avec loadGenealogie
	if (!g_CacheGenArray)
	{
		FlushProposantAncetre(banque);
		FlushGroupeProposantAncetre(banque);
//		GENError("Invalid use of LoadGroupeProposant or LoadGroupeAncetre function: start with LoadGenealogie");
		throw std::range_error("Invalid use of LoadGroupeProposant or LoadGroupeAncetre function: start with LoadGenealogie");
		//GENError("Utilisation invalide de la fonction LoadGroupeProposant ou LoadGroupeAncetre : il faut utiliser LoadGenealogie avant.");
	}
	//Verifie si des proposants et ancetre ont déjà été chargé a l'aide de loadproposant ou loadancetre
	if (!g_CacheVec[banque])
	{
		FlushProposantAncetre(banque);
//		GENError("Invalid use of LoadGroupeProposant or LoadGroupeAncetre function: start with Loadproposant or loadancetre");
		throw std::range_error("Invalid use of LoadGroupeProposant or LoadGroupeAncetre function: start with Loadproposant or loadancetre");
		//GENError("Utilisation invalide de la fonction LoadGroupeProposant ou LoadGroupeAncetre: il faut utiliser Loadproposant ou loadancetre avant.");
	}

	//Creation des tableaux correspondants
	g_CacheNbGroupe[banque]=nbGroupe;
	//Tableau de pointeur vers les pointeurs de proposants
	g_CacheGroup[banque]	= (CIndSimul***) memallocIN(nbGroupe,sizeof(CIndSimul**));
	g_CacheGrVecInd[banque] = (int*)		 memallocIN(nbGroupe,sizeof(int));
	
	for(int i=0;i<nbGroupe;i++)
	{
		if (BorneGr[i]>g_CacheVecInd[banque])
		{
			FlushProposantAncetre(banque);
//			GENError("Invalid use of LoadGroupeProposant or LoadGroupeAncetre function: too many individuals in the group compared to those loaded by loadproposant");
			throw std::range_error("Invalid use of LoadGroupeProposant or LoadGroupeAncetre function: too many individuals in the group compared to those loaded by loadproposant");
			//GENError("Utilisation invalide de la fonction LoadGroupeProposant ou LoadGroupeAncetre : il y a trop d'individus dans le groupe p/r à ceux chargés par loadproposant");
		}
		g_CacheGroup[banque][i]=g_CacheVec[banque]+BorneGr[i];
		if (i<(nbGroupe-1))
			g_CacheGrVecInd[banque][i]=BorneGr[i+1]-BorneGr[i];
		else
			g_CacheGrVecInd[banque][i]=g_CacheVecInd[banque]-BorneGr[i];
	}


	//Attribution des valeurs de sortie
	*nIndGr=g_CacheGrVecInd[banque];
	*GrProAnc=g_CacheGroup[banque];
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}

static int FlushGroupeProposantAncetre(ENUMBANQUE banque)
{
	//Détruit les groupe précédemment chargé	
	if (g_CacheGroup[banque])
		memfreeIN(g_CacheGroup[banque]);
	if (g_CacheGrVecInd[banque])
		memfreeIN(g_CacheGrVecInd[banque]);

	g_CacheGroup[banque]=NULL;
	g_CacheGrVecInd[banque]=NULL;
	g_CacheNbGroupe[banque]=-1;
	return 0;
}


//***************************************************************************/
//
//     	FONCTION POUR LOADER UNE GENEALOGIE,PROPOSANT & ANCETRE MAIS SANS CACHE
//
//***************************************************************************/
// Utilise un objet GestionMemoire pour géné la memoire

int LoadGenealogieNC(int* Genealogie,int* NInd, CIndSimul **Noeudarr, int** IndexRecherche,	int* profMax,int* nenfant,int* nindmasc,
					 GestionMemoire& MemCheck) 
{
	//Validation genealogie et recupere nombre d'enfant et profondeurmax	
	*NInd=LengthGenealogie(Genealogie,nenfant,profMax,nindmasc);

	//Creation du tableau de Noeud en utilisant GestionMemoire
	Clist* CacheChildList=NULL;
	*Noeudarr = (CIndSimul*)	memalloc(*NInd,sizeof(CIndSimul));	
	if (nenfant)
		CacheChildList=(Clist*) memalloc(*nenfant,sizeof(Clist));
		
	//Partie 1: CREATION ET INITIALISATION DE LA STRUCTURE NOEUD
	if (nenfant)
		ReCreeStructure(Genealogie,*Noeudarr,CacheChildList ,IndexRecherche); //Avec enfant
	else
		ReCreeStructure(Genealogie,*Noeudarr,NULL		  ,IndexRecherche); //Sans enfant
	return 0;
}

/**ATTENTION UTILISE LA GENEALOGIE PRESENTEMENT EN CACHE**/
int LoadVectorNC(int* vec,int nb, CIndSimul*** NproAnc, CIndSimul* NoeudArr,int nbind, int* IndexRecherche, GestionMemoire& MemCheck)
{
	try{
	//Bon.. creer le vecteur s'il est trop petit..	
	*NproAnc = (CIndSimul**)	memalloc(nb,sizeof(CIndSimul*));		
	
	//Boucle pour chaque proposant ou ancetre du vecteur
	for(int i=0;i<nb;i++)
	{	
		if (DescIEEEValue(vec+i)!=NULL){
//			GENError("Special IEEE %s is not a valid proband",DescIEEEValue(vec+i));
			char erreur[TAILLEDESCRIPTION];
			sprintf(erreur, "Special IEEE %s is not a valid proband",DescIEEEValue(vec+i));
			throw std::range_error(erreur);
			//GENError("Le caractère spécial IEEE %s n'est pas un proposant valide",DescIEEEValue(vec+i));
		}
		int tmpindice = ReTrouverIndiceStructure(vec[i],NoeudArr,IndexRecherche, nbind);
		if (tmpindice==-1){
//			GENError("Proband %d is not in the genealogy",vec[i]);
			char erreur[TAILLEDESCRIPTION];
			sprintf(erreur, "Proband %d is not in the genealogy",vec[i]);
			throw std::range_error(erreur);
			//GENError("Le proposant %d n'est pas inclus dans la généalogie",vec[i]);
		}
		else
			(*NproAnc)[i]=NoeudArr+tmpindice;					
	}
	
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


