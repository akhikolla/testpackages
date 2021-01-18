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
/**	\brief	macro pour déterminer quel est le plus grand des deux

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
	de type int*, g_Complet est une variable globale qui sert au décompte des individus par génération
	dans les fonctions addCompl et initComplet
*/
static int* g_Complet;

/**	\brief
	de type int*, g_Fond est une variable globale qui sert au décompte des fondateurs par génération
	dans les fonctions addFond et FondParGen
*/
static int* g_Fond;

/**	

	\brief	outil pour initComplet
	\internal
	\param	ind		[in] individu de la généalogie
	\param	prof	[in] profondeur courante ie numéro de génération

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
	\param	ind		[in] individu de la généalogie
	\param	prof	[in] profondeur courante ie numéro de génération

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

	Considérant les proposants donnés de génération 0, on effectue une remonté, génération par
	génération afin de déterminer le nombre de fondateurs à chaque génération. Les résultats
	sont présentés dans le tableau de retour.

	\param Genealogie	[in] Une genealogie binaire (gen.genealogie)

	\param prop			[in] Vecteur des individus de depart à considérer

	\param nbProp		[in] Le nombre de ces individus
	
	\retval retour		[out] Un pointeur vers un double

	\return 0 si la fonction est executé avec succès

	\remark lors de l'appel, retour doit être impérativement de taille nbgen(Genealogie)
			et doit être initialise avec des zeros dans chaque case
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


