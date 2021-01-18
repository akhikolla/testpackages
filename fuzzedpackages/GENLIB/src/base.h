/*! \file base.h
\brief Library Genlib: Classe et Enumation de base

	Ce fichier comprend toutes les classes énumération qui sont d'interet a toute les fonctions
	Il automatise les modifications du code entre Unix et WIN32.
	De plus, ce fichier sert de fichier option et peut servir a activer les options suivante :
	1) L'utilisation de la gestion memoire integrée de SPLUS au lieu de celle standard (malloc)
	2) L'utilisation de la generation de nombre aleatoire à l'aide de SPLUS au lieu de rand,srand

\remark Sous windows, le detecteur de memory leak et automatiquement activer en mode debug

\author Sébastien Leclerc
\contributor Jean-François Lefebvre
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

*/


/*PARAMETRE DE COMPILATION*/

#define MODETEST 			//!< Si defini alors TOUS les blocs de code fesant appel a des fonctions specifique de splus seront déactiver
//#define USESPLUSALLOC		//!< Si defini alors tous les appels à memalloc seront explicitement converti en S_alloc
//#define USESDEBUG			//!< Si defini alors certain message d information seront affiché

#define USESPLUSRANDOM		//!< Si defini alors les fonctions utiliseront unif_rand au lieu de srand et rand
#define USESTIMER			//!< Si defini alors la variable .Last.timeSec equivaudra au temps d execution de la derniere fonction

#ifndef GENBASE
#define GENBASE

//#ifdef WIN32
	//#define MODETEST			//##A Effacer### defini automatiquement MODETEST (voir plus haut) sur plateforme WIN32
  	//#define _CRTDBG_MAP_ALLOC
//#endif

#ifndef MODETEST
//	#include <S.h>	
	#if defined _WIN32 || defined _WIN64
		#include "sconnect.h"
		#pragma warning(disable:4786)
	#endif
#else	
	#include<stdio.h>
	#include<stdlib.h>
#endif

#if defined _WIN32 || defined _WIN64
	#include "crtdbg.h"
	//#define WINDOWS_CONFLICT
	
	//Un XLONG est un entier unsigned a 64 peut importe la plateforme
	//#define XLONG		unsigned _int64  --JF
	#define XLONG		unsigned long long // --JF
	#define MAX_XLONG	_UI64_MAX

	//#define SXLONG		_int64	  --JF
	#define SXLONG		long long		//--JF
	#define MAX_SXLONG	_I64_MAX

	#define WINCDECL	__cdecl 
	#define FASTCALL    __fastcall
#else
//	#include <boost/random/random_device.hpp>
	#define XLONG unsigned long long int
	#define MAX_XLONG	ULLONG_MAX

	#define SXLONG long long int
	#define MAX_SXLONG	LLONG_MAX

	#define WINCDECL 
	#define FASTCALL
#endif

//#include<cstdlib>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/*COMPATIBILITE AVEC SPLUS ANSI QUE LES DEUX PLATEFORMES WIN32 ET LINUX*/
//#ifdef MODETEST
//	#undef USESPLUSALLOC
//	#undef USESPLUSRANDOM
//	#undef USEALTERNATEEXCEPTION
//#endif

//constante de base pour les algorithme
const int MAX_SUPPORTED_GENERATION=500;

/** @defgroup ALLOCMEM ALLOCATION MEMOIRE

	Les allocations memoire peuvent-être fait de deux manieres.
<OL>
  <li>Splus S_alloc : dans ce cas la memoire est automatiquement recuperer a la sortie de la fonction</li>
  <li>C malloc	   : Dans ce cas, un objet GestionMemoire est creer et il se charge d'effectuer
					 automatiquement la recuperation de la memoire a la sorti de la fonction</li>
</OL>
<br>Si USESPLUSALLOC est defini alors la methode d'allocation utiliser sera les S_alloc de splus.
		<br>Dans une fonction, avant de faire la premiere allocation memoire la macro
	INITGESTIONMEMOIRE doit etre faite.
<br>
	Cette macro defini (au besoin) un objet GestionMemoire qui est charge d'eliminer les fuites de memoire
<br><br>Pour les allocations memoires utiliser les fonctions suivantes
	<table>	
	<tr><td><b>INITGESTIONMEMOIRE</b></td><td>Precede la premiere utilisation de memalloc dans une fonction</td></tr>
	<tr><td><b>memalloc</b></td><td>Alloue dynamiquement de la memoire, la memoire sera automatiquement recupere au sortir de la fonction courante</td></tr>
	<tr><td><b>memallocIN</b></td><td>Alloue dynamiquement de la memoire mais n'utilise pas le gestionnaire memoire s'il y en a un
									 Permet d'allouer de la memoire qui "durera" plus longtemps que la fonction en cour</td></tr>
	<tr><td><b>memfreeIN</b></td><td>Recupere de la memoire alloue par memallocIN</td></tr>
	<tr><td><b>memfree</b></td><td>Obsolete comme la recuperation de la memoire est automatique</td></tr>
	</table>
		
*  @{
*/

/*#ifdef USESPLUSALLOC
	//ALLOCATION DANS SPLUS
	#ifdef WIN32
		#define memalloc(n,size)	S_alloc((n),(size))
		#define memallocIN(n,size)	S_alloc((n),(size))
	#else

		#define memalloc(n,size)	S_alloc((n),(size),S_evaluator)		
		#define memallocIN(n,size)	S_alloc((n),(size),S_evaluator)		
	#endif
	
	#define memfree(ad)			NULL
	#define memfreeIN(ad)		NULL
	#define INITGESTIONMEMOIRE	NULL;

#else*/
	///ALLOCATION CONVENTIONNEL
	#define INITGESTIONMEMOIRE	GestionMemoire MemCheck;
	#ifndef USESDEBUG
		#define memalloc(n,size)	MemCheck.alloc(n,size)	
	#else
		#define memalloc(n,size)	MemCheck.alloc(n,size)	
	#endif
	#define memfree(ad)			NULL
	#define memallocIN(n,size)	malloc((n)*(size))
	#define memfreeIN(ad)		free(ad)
//#endif
/** @} */ //fin du groupe





/** @defgroup GENNUM GENERATION DE NOMBRE ALEATOIRE  
	
	Il est possible de genere des nombres aleatoire avec SPLUS

	Je suggere pour la compatibilite que toutes les fonctions soit capable d'utilise l'un et l'autre
	Si possible utilise un 
	<br> #ifdef USESPLUSRANDOM  
	<br> #else
	<br> #endif
	<br> pour selectionne la methode a utiliser

	<br> Les fonctions utilisable pour la generation de nombre aleatoire sont les suivants
	<table>	
	<tr><td><b>initrand</b></td><td>Initialise le generateur de nombre aleatoire, devrais etre au debut de chaque fonction</td></tr>
	<tr><td><b>outrand</b></td><td>Termine le generateur, devrais etre utilise a la fin de la fonction qui utilise urand</td></tr>
	<tr><td><b>urand</b></td><td>Genere un nombre aleatoire a l'aide de splus (double)</td></tr>
	<tr><td><b>rand</b></td><td>Genere un nombre aleatoire en c classique (entier de 1.. RANDMAX)</td></tr>
	</table>

*  @{
*/
//#ifdef USESPLUSRANDOM
//	#define initrand()	seed_in((long*)NULL,S_evaluator)	///< Initialise le generateur de nombre aleatoire de SPLUS	
//	#define outrand()	seed_out((long*)NULL,S_evaluator)	///< Termine le generateur de nombre aleatoire de SPLUS	
//	#define urand()	(rand()/RAND_MAX1) //unif_rand(S_evaluator)	***JFL***	///< Genere un nombre aleatoire de 0-1 (double) a l'aide de SPLUS
//#else

	#define initrand()		srand(time(NULL))	///< Initialise le generateur de nombre aleatoire rand a l'aide de srand
	#define outrand()		NULL				///< Termine le generateur de nombre aleatoire
	#define urand()		(rand()/RAND_MAX1)	///< Genere un nombre aleatoire de 0-1 (double) a l'aide de SPLUS			
//#endif

/** @} */ //fin du groupe


// ************************************************************/
//      DECLARATION DES STRUCTURES POUR LES INDIVIDUS
// ************************************************************//


//Definit tous les types d'individu possible
/*
	Une fois que tout les individu on été transforme en CIndSimul (Appele Noeud)
	il ne sont pas tous utile, souvent seulement ceux entre un certain ancetre et
	un proposant sont utile. (utile au sens, ce noeud peut influence le resultat du calcul courant)

	typenoeud sert a identifier quel et l'utilite d'un noeud dans le probleme courant

  \remark la valeur par defaut est GENNONEXPLORER
  
  \sa ExploreArbre
 */

enum typenoeud_t {
	GENNONEXPLORER=0, ///< Utilite n'as pas encore ete determiner (DEFAULT)
	GENINUTILE,		  ///< Noeud Inutile, on peut l'ignorer dans le calculs
	GENNOEUD,		  ///< Noeud Utile, doit etre considere dans les calculs
	GENDEPART,		  ///< Ce Noeud est un Ancetre 
	GENPROPOSANT,     ///< Proposant, ce noeud est un proposant et il est confirme utile
	GENPROPOSANTINUTILE	 ///< Le noeud est un proposant mais il n'est pas encore confirme utile
};

enum sex_t {GEN_INCONNU=0,GEN_MASC=1,GEN_FEM=2};

struct Clist;
struct CIndSimul;

///Structure qui represente un individu (NOEUD) 
/**
	Soit
	<br>allele,
	<br>prob[]
	<br>pGen
	<br>Noind
	<br>bFlagSort
	<br>Ces parametres n'ont pas de signification permanente, ils changent de signification dependamment
	 de la fonction dans laquelle ils sont utilises.

	<br>Par consequent, ils ne sont pas initialises par les fonctions de creation de genealogie.
		 
	<br>Ils peuvent etre utilises a volonte pour differentes utilisations. les utilisations qui suivent ne sont qu'a titre d'exemple

		\remark si Splus sous unix avait supporte les templates, ce type de variable n'aurait pas ete necessaire
*/
struct CIndSimul 
{
	//Essentiel
	int nom;			///< No de l individu
	sex_t sex;		///< Sexe de l'individu
	CIndSimul *pere;    ///< Ptr vers le Noeud qui represente son Pere. NULL si omis
	CIndSimul *mere;	///< Ptr vers le Noeud qui represente sa Mere. NULL si omis	 
	Clist* fils;		///< Liste des enfants de ce Noeud
	int noind;	///< EX : No indice du Noeud dans un vecteur de priorite (d ordre) 

	//Accessoire (peuvent-être modifier)
	int allele;			///<EX : Nombre d allele que l individu possede 	
	int allele2Pos[2];		///<EX : Allele 1 (chiffre) et allele 2 (chiffre) que l individu possede 
	int alleleAttendu;		///<EX : Nombre d allele que l individu devrait posseder
	
	union 
	{
		double prob[3];	///<EX : Probabilite d avoir 0,1,2 allele 
		void* ptr[6];
	};
	union 
	{
		double *pGen;		///<EX : Pointeur vers un tableau de valeur double
		int iind;			///<EX : Indicatif du no d'individu
	};
	
	int bFlagSort;	        ///< EX : Indique si ce noeud a 0-1-2-++ enfant qui sont classe utile 	
	double dFlagSort;
	typenoeud_t etat;	///< Etat/Utilite de ce noeud			 
}; 

///Liste de Noeud. Utilise pour les enfants d un individu (CIndSimul)
struct Clist
{
	Clist*	   next;	
	CIndSimul* noeud;	
};

// ************************************************************/
// MESSAGE & CODE  D ERREUR
// ************************************************************// 
const int GFALSE = 0;
const int GTRUE  = 1;
#include "exception.h"

// ************************************************************/
// GESTIONNAIRE DE MEMOIRE
// ************************************************************//

struct GestionMemoireBlock;
///Gestionnaire de memoire 
/**
	Gestionnaire de memoire qui permet de faire des allocations style malloc
	Au moment ou l'objet GestionMemoire est detruit il appelle delete automatiquement
	sur toute les allocations memoire qu'il a fait.
*/
class GestionMemoire
{
private:
	GestionMemoireBlock *tableaublock; ///< Liste de tableau de void utilise a l'interne pour memorise toute les allocations
	GestionMemoireBlock *startblock;   ///< Ptr vers le premier tableau
	int n;							   ///< Nombre d'allocation sur le tableau courant
	char UseMalloc;					   ///< voir constructeur
public:
	GestionMemoire(char UseStdMalloc=0);
	~GestionMemoire(); 
	void* alloc(int n,size_t size);
	void    add(void* item);
};

// ************************************************************/
// DECLARATION 
// ************************************************************//

/// Probabilite avec un pere qui a INDICE1 allele et une mere INDICE2 allele d obtenir INDICE3 allele
extern double TransGenCum[3][3][3]; 

/// Probabilite avec un pere qui a INDICE1 allele et une mere INDICE2 allele d obtenir INDICE3 allele
extern double TransGen[3][3][3]; 

extern const char *stype[]; 
///Constante utilise pour genere des nombres aleatoires entre 0 et <1
const double RAND_MAX1 = (double)RAND_MAX+1;

#define MIN(x,y) ( (x)<(y)?(x):(y) )
#define MAX(x,y) ( (x)>(y)?(x):(y) )

#endif



