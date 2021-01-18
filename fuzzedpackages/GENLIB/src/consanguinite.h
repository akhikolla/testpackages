/*! \file consanguinite.h
\brief Interface des fonctions de calcul de F

Interface de toutes les fonctions en rapport avec la consanguinit�

\author S�bastien Leclerc
\contributor Jean-Fran�ois Lefebvre

*/

#ifndef GENCONSANG
#define GENCONSANG

int consan(int* Genealogie,
	int* proposant, int NProposant,int Niveau,
	double* pdConsan, int printprogress);


int consanFs(int* Genealogie,
	int* proposant, int NProposant,int NiveauMin,int NiveauMax,
	double* pdDeepConsan, int printprogress);

#endif


