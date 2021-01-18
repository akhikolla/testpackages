/*! \file probsimul.h
\brief Interface library Genlib: fonction diverse

non déterminé

\author Sébastien Leclerc, Claire Gires
\contributor Jean-Francois Lefebvre

*/

#ifndef GENPROBSIMUL
#define GENPROBSIMUL

int CountChild(int* Genealogie, int* plProposant,int NProposant, int* retour);
int TakeInd(int* Genealogie,int FonOuPro, int* retour, int* taille);

//int testEbranche(int * Genealogie, int * plProposant, int lNProposant, int * plAncetre, int lNAncetre, int * NouvelGenealogie, int * tailleNouvelGenealogie);
int ebranche(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre, int* NouvelGenealogie,int* tailleNouvelGenealogie);

int compareGen(int* Gen1, int* Gen2, int* retour, int mustprint);
int numeroGen(int* Genealogie, int* plProposant, int NProposant, int* retour);
int numeroGenMin(int* Genealogie, int* plProposant, int NProposant, int* retour);
int numeroGenMoy(int* Genealogie, int* plProposant, int NProposant, double* retour);
#endif



