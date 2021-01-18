/*! \file outils.h
\brief Interface library Genlib: outils généraux

Interface de toutes les fonctions qui sont d'un intéret général pour chaque menu

\author Sébastien Leclerc
\contributor Jean-Francois Lefebvre
*/

#ifndef GENOUTILS
#define GENOUTILS

struct CDuoPair;
typedef unsigned int HugeInt;

//Creation classique de la structure (surtout pour gen.genealgie)
void CreeStructure(CIndSimul* Noeud,int* indice,int* pere, int* mere,int* sex,int iNind,int* countchildren,CDuoPair *Trie,Clist** ChildArray);
int DetruireStructure(Clist* ChildArray);
void SortGenealogie3Vecteur(int* ind,int *pere, int* mere,int* sex, int nind);

//Construction et utilisation de genealogie "binaire"
///Version courante du format de genealogie a un vecteur
const int VERSIONGENEALOGIE=1; 
///Calcule la taille du vecteur requis pour inserer une genalogie a 1 vecteur pour NIndividu et Nenfant
#define TAILLEGENVERSION7(NIndividu,Nenfant) ((NIndividu)*6+9+(Nenfant)+4)

int LengthGenealogie(int* Genealogie,int* nenfant=NULL,int* nprofmax=NULL,int* nindmasc=NULL);
int ValidateGenealogie(int* Genealogie);

int CompleteGenealogie(int* ind,int* pere,int*mere,int*sex,int* newind, int* newpere,int* newmere,int* newsex,int *pNIndividu);

int CreerGenealogie(int* plIndividu,int* plPere, int* plMere,int* plSexe,int lNIndividu,int* saveptr);

int ReCreeStructure(int* Genealogie,CIndSimul* Noeud, Clist* Children, int** IndexRecherche);


//Recherche l'indice d'un individu dans le tableau de Noeud à l'aide d'un index de recherche
int ReTrouverIndiceStructure(int nom, CIndSimul* Noeud, int* IndexRecherche, int iNind);

//Ebranchage et etiquetage des noeuds
int ExploreArbre(CIndSimul* Noeud);
void ExploreArbreTousDescendant(CIndSimul* Noeud);

//Ordonnancement de la structure
void PrepareSortPrioriteArbre(CIndSimul* Noeud,int iNind);
void StartSortPrioriteArbre(CIndSimul* Noeud,CIndSimul** Ordre,int *index,int *TableSaut);
int SortPrioriteArbre(CIndSimul* Noeud,CIndSimul** Ordre,int *index,int *TableSaut,Clist **list=NULL);
int OrdonneStructure(CIndSimul* Noeud,CIndSimul** Ordre,int iNind,int SensInverse=0,int* profMax=NULL);
int classeGen(CIndSimul* Gen, int nbInd, int* tab = NULL, CIndSimul** tabind=NULL);
int classeGenMin(CIndSimul* Gen, int nbInd, int* tab = NULL, CIndSimul** tabind=NULL);
int classeGenMoy(CIndSimul* Gen, int nbInd);

//Accessoire
int interval(int x, int min, int max);
//long round(double value); // enleve par JFL : conflit avec double round(double value) de mathcalls.h
double pow2(int y);

//Utilitaire cryptographique
int millerRabin(unsigned int n, unsigned int t);
unsigned int irand(unsigned int a, unsigned int b);


//FONCTION UTILISANT LA CACHE
//Fonction de haut-niveau de creation de genealogie
int LoadGenealogie(int* Genealogie,int loadChildren,int* NInd, CIndSimul **Noeudarr, int** IndexRecherche=NULL);
int LoadProposant(int* Proposant,int nbProposant,CIndSimul*** pro);
int LoadAncetre(int* Ancetre,int nbAncetre,CIndSimul*** anc);
int LoadProfondeurMax();
int LoadNIndMasc();
int FlushGenealogie();
//Fonction de haut-niveau pour l'assignation de groupe a une série de proposant
int LoadGroupeProposant(int* Proposant, int* BorneGr, int nbGroupe, CIndSimul**** GRpro,int** nIndGr);
int LoadGroupeAncetre  (int* Ancetre  , int* BorneGr, int nbGroupe, CIndSimul**** GRanc,int** nIndGr);


//Fonction de haut-niveau de creation de genealogie mais sans la cache
//Utilise l'objet MemCheck Creer par INITGESTIONMEMOIRE
int LoadGenealogieNC(int* Genealogie,int* NInd, CIndSimul **Noeudarr, 
					 int** IndexRecherche,	int* profMax,int* nenfant,int* nindmasc,
					 GestionMemoire& MemCheck);
int LoadVectorNC(int* vec,int nb, CIndSimul*** NproAnc,
				 CIndSimul* NoeudArr,int nbind, int* IndexRecherche,
				 GestionMemoire& MemCheck);

//Fonction spécifique pour la gestion des valeurs spéciale IEEE
char* DescIEEEValue(int* val);

#endif





