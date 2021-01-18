/*! \file statanal.cc
\brief Implementation des outils pour l analyse statistique


\author Claire Gires

*/ 

/// Autorise l'affichage du niveau de progression sur la sortie standard
/** 
	Si ALLOWPRINTPROGRESS est defini, les fonctions suivantes peuvent
	indiquer leur niveau de progression sur la sortie standard stdout
*/
#define ALLOWPRINTPROGRESS


#include"base.h" 
#include"outils.h" 
#include"statanal.h" 



// ********************************************************************
//
//			MACRO
//
// ********************************************************************
/**	\brief	macro pour d�terminer quel est le plus grand des deux

	\param	x [in]

	\param	y [in]

	\return	x ou y 

	
*/
#define GLMAX(x,y) ( (x)>(y) ? x : y )


// ********************************************************************
//
//			VARIABLE GLOBALE
//
// ********************************************************************

/**	\brief
	de type int*, g_Complet est une variable globale qui sert au d�compte des individus par g�n�ration
	dans les fonctions addCompl et initComplet
*/
static int* g_Complet;

/**	\brief
	de type int*, g_Fond est une variable globale qui sert au d�compte des fondateurs par g�n�ration
	dans les fonctions addFond et FondParGen
*/
static int* g_Fond;

/**	

	\brief	outil pour initComplet
	\internal
	\param	ind		[in] individu de la g�n�alogie
	\param	prof	[in] profondeur courante ie num�ro de g�n�ration

	\return	static void

	
*/
static void addCompl(CIndSimul* ind, int prof)
{

	g_Complet[prof]++;
	
	if (ind->pere != NULL)
		addCompl(ind->pere, prof+1);
	if (ind->mere != NULL)
		addCompl(ind->mere, prof+1);
}

//------------------------------------------------------------------

/**	

	\brief	outils pour FondParGen
	\internal
	\param	ind		[in] individu de la g�n�alogie
	\param	prof	[in] profondeur courante ie num�ro de g�n�ration

	\return	static void

	
*/
static void addFond(CIndSimul* ind, int prof)
{

	if(ind->pere==NULL && ind->mere==NULL)
		g_Fond[prof]++;
	
	if (ind->pere != NULL)
		addCompl(ind->pere, prof+1);
	if (ind->mere != NULL)
		addCompl(ind->mere, prof+1);
}


/*! 
	\brief Calcule le nombre de fondateurs par generation

	Consid�rant les proposants donn�s de g�n�ration 0, on effectue une remont�, g�n�ration par
	g�n�ration afin de d�terminer le nombre de fondateurs � chaque g�n�ration. Les r�sultats
	sont pr�sent�s dans le tableau de retour.

	\param Genealogie	[in] Une genealogie binaire (gen.genealogie)

	\param prop			[in] Vecteur des individus de depart � consid�rer

	\param nbProp		[in] Le nombre de ces individus
	
	\retval retour		[out] Un pointeur vers un double

	\return 0 si la fonction est execut� avec succ�s

	\remark lors de l'appel, retour doit �tre imp�rativement de taille nbgen(Genealogie)
			et doit �tre initialise avec des zeros dans chaque case
*/
int FondParGen(int* Genealogie,int* prop,int nbProp, int* retour)
{
	//CREATION DE TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GFALSE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(prop,nbProp,&NoeudPro);
	
	g_Fond= retour;

	for(int i=0; i<nbProp; i ++)
	{
		addFond(NoeudPro[i], 0);
	}

	return 0;

}


