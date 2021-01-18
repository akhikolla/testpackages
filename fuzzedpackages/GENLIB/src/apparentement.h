/*! \file apparentement.h
\brief Interface des fonctions de calcul de Phi


\author Sébastien Leclerc
\contributor Jean-François Lefebvre

*/

#ifndef GENAPPARENTEMENT
#define GENAPPARENTEMENT

#include "outilsalloc.h"
#include "hashtable.h"
#include <Rcpp.h>
#define R_NO_REMAP

/// Structure de gestion memoire mémoire utilise par kinship2
/**
	Utiliser pour effectuer une serie d'assignement memoire et ensuite laisser la voie libre au
	gestionnaire de memoire pour qu'il effectue lui-meme la "gargage collection"
	  
	Normalement: la structure est initialiser comme suit
		<br>&nbsp;&nbsp;&nbsp;&nbsp; MemGest : On lui assigne un gestionnaire memoire de porte locale
		<br>&nbsp;&nbsp;&nbsp;&nbsp; Tableau : NULL
		<br>&nbsp;&nbsp;&nbsp;&nbsp; ctableau: Un gros nombre, habituellement LONG_MAX
	
*/
struct Kinship4Struct
{	
private:	
	short NiveauMax; //Niveau maximal voulu 0..Niveaumax = NiveauMax+1 Case...
	BlockAlloc<double> memblock; 
	SMPile<double*,MAX_SUPPORTED_GENERATION> PileCosan; //L'adresse de la cosanguinité a écrire

	//SUPPORT MT		
	static CSema m_acces; //Semaphore pour l'acces a allele qui sert de variable partage	
public:
	//CREATION ET INITIALISATION
	Kinship4Struct();	
	Kinship4Struct(short NiveauMax, double* resultat); //Resultat tableau de taille NiveauMax+1		
	void Initialise(short NiveauMax, double* Resultat);
	static void InitMT();
	static void ReleaseMT();
	inline void remplace(double* resultat) {PileCosan.pop();PileCosan.push(resultat);} 
	
    //FONCTION AMIS
	friend void FASTCALL Kinship4(CIndSimul* Ind1, CIndSimul* Ind2, short ttl1, short ttl2,Kinship4Struct &T);
	friend void FASTCALL Kinship4MT(CIndSimul* Ind1, CIndSimul* Ind2, short ttl1, short ttl2,Kinship4Struct &T);	
};

int Phis(int* Genealogie, int* proposant, int NProposant,int NiveauMin,int NiveauMax, double* pdMoyenne,double *MatrixArray, int printprogress);
int PhisMT(int* Genealogie, int* proposant, int NProposant,int NiveauMin,int NiveauMax, double* pdMoyenne,double *MatrixArray, int printprogress);


double FASTCALL Kinship(CIndSimul* Ind1,CIndSimul* Ind2,short ttl1,short ttl2);
void FASTCALL Kinship4(CIndSimul* Ind1, CIndSimul* Ind2, short ttl1, short ttl2,Kinship4Struct &T);
void FASTCALL Kinship4MT(CIndSimul* Ind1, CIndSimul* Ind2, short ttl1, short ttl2,Kinship4Struct &T);
		
#endif




