
/*! \file apparentement.cc
\brief Implementation des fonctions de calcul de Phi

Calcul et Analyse de diverse valeur dérivé de Phi et Phi moyen

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
#include "basemt.h"
#include "userInterface.h"
 
//Temporaire
#include <fstream>
#include <assert.h>
#include <ctime>
#include <cstdlib>
#include <Rcpp.h>
using namespace std;
//#include "windows.h"

#define R_NO_REMAP

#define DIV 1024

#define WIDTH 7



//Fichier de sortie pour le debuggage
//#define TESTRACEFILE
#include <string.h>
#include <math.h>
#include <limits.h>


// ********************************************************************
//
//			CONSTANTE
//
// ********************************************************************


// ********************************************************************
//
//			PUBLIC
//
// ********************************************************************

/*! 
	\brief Calcule les matrices Phis pour chaque profondeur de NiveauMin a NiveauMax (OPTIMISER)

		La fonction calcule l'apparentement entre chaque proposant qui lui est fourni.
		Et ce pour chaque profondeur de NiveauMin à NiveauMax soit ((NiveauMax - NiveauMin)+1) Matrices.
		
	\param Genealogie		[in] Une genealogie construite à l'aide de gen.genealogie

	\param proposant		[in] Vecteur avec les numeros des proposants à étudier
	\param NProposant		[in] Nombre d'éléments du vecteur proposant
  
	\param NiveauMin		[in] Profondeur minimale
	\param NiveauMax		[in] Profondeur maximale

	\retval pdMoyenne		[out] Adresse d'un tableau de ((NiveauMax - NiveauMin)+1) element
							En cas de succes, contient la valeur du Phi moyen pour chaque niveau
							<br>Ex: NiveauMin=5  NiveauMax:10 <br>
							alors <br>
								pdMoyenne[0]= Phis Moyen profondeur 5<br>
								pdMoyenne[1]= Phis Moyen profondeur 6<br>
								..<br>
								pdMoyenne[5]= Phis Moyen profondeur 10


	\retval MatrixArray		[out] Adresse d'un tableau de NProposant x NProposant x ((NiveauMax - NiveauMin)+1) element
							En cas de Succes, pdMatricePhi contient la matrices Phi de chaque niveau.
		  					<br>Ex: <br>
								soit:									<br>
								    NiveauMin: 5<br>
									NiveauMax: 10<br>
									Ind1 : Indice du proposant 1 dans la liste proposant<br>
									Ind2 : Indice du proposant 2 dans la liste proposant<br>
									ProfondeurVoulu: 6										<br>							
							  Apparentement entre proposant 1 et proposant 2 <br>
								= pdMatricePhi[((ProfondeurVoulu-NiveauMin)*NProposant*NProposant) +Ind1 x NProposant + Ind2 ]									
						
	\param printprogress	[in] Imprime une serie de message indiquant la progression du calcul

	\return 0 si la fonction est executé avec succès
	
	\remark Version optimise de Phis(), au lieu de lancer une fois kinship pour chaque niveau possible.
			elle utilise kinship2 qui ne fait qu'une seule remonter et note la valeur du kinship pour chaque profondeur

*/
int Phis(int* Genealogie, int* proposant, int NProposant,int NiveauMin,int NiveauMax, double* pdMoyenne, double *MatrixArray,int printprogress)
{
 	try{
	//TEST D'ERREUR DE BASE
	if (NProposant<2){
//		GENError("At least two probands are required for this function");
		throw std::range_error("At least two probands are required for this function");
	}
	if (NiveauMin<0){
//		GENError("depthmax and depthmin must be greater than zero.");
		throw std::range_error("depthmax and depthmin must be greater than zero.");
	}
	if (NiveauMax<NiveauMin){
//		GENError("depthmax must be greater or equal to depthmin");
		throw std::range_error("depthmax must be greater or equal to depthmin");
	}
	if (NiveauMax>SHRT_MAX){
//		GENError("depthmax must be smaller than %d",SHRT_MAX);
		char erreur[TAILLEDESCRIPTION];
		sprintf(erreur, "depthmax must be smaller than %d\n", SHRT_MAX);
		throw std::range_error(erreur);
	}
	//CREATION DU TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GFALSE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(proposant,NProposant,&NoeudPro);

	//Creation du tableau de Noeud
	const int LongEcart		= NiveauMax-NiveauMin+1;
	INITGESTIONMEMOIRE;
	double*    Phideep		= (double*)	memalloc(NiveauMax+1,sizeof(double)); 	
	
	//Initialisation de l'algorithme (valeur par niveau)
	for(int Niveau=0;Niveau<LongEcart;Niveau++)
		pdMoyenne[Niveau]=0.0;		
	
	//Initialisation des noeuds
	for(int i=0;i<lNIndividu;i++)
		Noeud[i].pGen=NULL;
	
	//Iterateur pour toute la matrice
	const int NItem	 = ((NProposant*NProposant)-NProposant)/2;
	const int taillematrix	 = NProposant*NProposant;

	//CALCUL LUI MEME	
	const short niveauMax = short(NiveauMax);
	Kinship4Struct Element(niveauMax,Phideep);	

	//Barre de progression
//	CREATE_PROGRESS_BAR_MATRIX(NProposant,printprogress)
	for(int cPro1=0;cPro1<NProposant;++cPro1)
	{
		for(int cPro2=cPro1;cPro2<NProposant;++cPro2)
		{
			/*
			if (cPro2==cPro1)
			{
				for(int a=0;a<LongEcart;a++)
					MatrixArray[(a*taillematrix)+cPro1*NProposant+cPro2]=0.5;					
			}
			else
			{*/
				//Partie 2: EBRANCHER L'ARBRE, ORGANISER LE NOEUD EN VECTEUR EN RESPECTANT PRECEDENCE
  				//Calcul de la valeur de retours
				 
				//remise à zero de phi deep
				for(int a=0;a<=niveauMax;a++)	Phideep[a]=0.0;

				//calcul de phi				
				Kinship4(NoeudPro[cPro1],NoeudPro[cPro2],niveauMax,niveauMax,Element);				

				//Recupere le résultat de phideep				
				for(int a=0;a<LongEcart;a++)
				{
					if (Phideep[a]<0.5)
						pdMoyenne[a]+=Phideep[a+NiveauMin];
					MatrixArray[(a*taillematrix)+cPro1*NProposant+cPro2]=Phideep[a+NiveauMin];
					MatrixArray[(a*taillematrix)+cPro2*NProposant+cPro1]=Phideep[a+NiveauMin];					
				}	

				//Affichage des progress
//				INCREMENT_PROGRESS_BAR()
			//}			

		}//Fin itérateur proposant 2
		
	}//Fin itérateur proposant 1
		
	//Completer les matrices

	//Calcul de la moyenne........
	for(int Niveau=0;Niveau<LongEcart;Niveau++)
		pdMoyenne[Niveau]/=NItem;
		
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			} catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}


// ********************************************************************
//
//			PUBLIC :VERSION MT 
//
// ********************************************************************

struct CBASEMTPhisMT : public CMtGlobalMessage
{
	CIndSimul* ind1;
	CIndSimul* ind2;
	short niveauMax;		//Nombre de remonte dans la généalogie

	//Indice ou placer la reponse dans le tableau de reponse
	int indice1;			
	int indice2;			
	
	// C'est lui qui contient les valeurs de retours
	Kinship4Struct elem; //Doit-être initialisé
	double* ret;
};

BASEMT_CREATE_GLOBALMESSAGE(CBASEMTPhisMT,1)

BASEMT_DEBUT_HELPERFCT(CBASEMTPhisMT,1)
		//LANCEMENT DU CALCUL ET RECUPERE LE RESULTAT
		Kinship4MT(BASEMT_HLPMES.ind1,BASEMT_HLPMES.ind2,
				 BASEMT_HLPMES.niveauMax,BASEMT_HLPMES.niveauMax,BASEMT_HLPMES.elem); 
BASEMT_FIN_HELPERFCT() 

int PhisMT(int* Genealogie, int* proposant, int NProposant,int NiveauMin,int NiveauMax,
	double* pdMoyenne, double *MatrixArray,int printprogress)
{
	try{
	//TEST D'ERREUR DE BASE
	if (NProposant<2){
//		GENError("At least two probands are required for this function");
		throw std::range_error("At least two probands are required for this function");
		//GENError("Il faut au minimum 2 proposant pour utilise cette fonction");
	}
	if (NiveauMin<0){
//		GENError("depthmax and depthmin must be greater than zero.");
		throw std::range_error("depthmin and depthmin must be greater than zero.");
		//GENError("Le niveau minimum et le niveau maximum doivent-être supérieur à zéro");
	}
	if (NiveauMax<NiveauMin){
//		GENError("depthmax must be greater or equal to depthmin");
		throw std::range_error("depthmax must be greater or equal to depthmin");
		//GENError("Le niveau maximum doit-être supérieur ou égal au niveau minimum");
	}
	if (NiveauMax>SHRT_MAX){
//		GENError("depthmax must be smaller than %d",SHRT_MAX);
		char erreur[TAILLEDESCRIPTION];
		sprintf(erreur, "depthmin must be smaller than %d",SHRT_MAX);
		throw std::range_error(erreur);
		//GENError("Le niveau maximum doit-être inférieur à %d",SHRT_MAX);
	}
	//CREATION DU TABLEAU D'INDIVIDU
	int lNIndividu;
	CIndSimul *Noeud=NULL;
	LoadGenealogie(Genealogie,GFALSE,&lNIndividu,&Noeud);

	//CREATION D'UN VECTEUR DE PROPOSANT
	CIndSimul **NoeudPro=NULL;
	LoadProposant(proposant,NProposant,&NoeudPro);

	//Creation du tableau de Noeud
	const short niveauMax = short(NiveauMax);
	const int LongEcart		= NiveauMax-NiveauMin+1;
	INITGESTIONMEMOIRE;
		
	//Initialisation de l'algorithme (valeur par niveau)
	for(int Niveau=0;Niveau<LongEcart;Niveau++)
		pdMoyenne[Niveau]=0.0;		
	
	//Initialisation des noeuds
	for(int i=0;i<lNIndividu;i++)
		Noeud[i].pGen=NULL;					
	
	
	//Initialisation d'une structure multithread
	BASEMT_DEBUTBOUCLE_INITMESSAGE(1) 		
		//initialisation des données Phis de la structure
		BASEMT_MESSAGE(1).ind1=NULL;
		BASEMT_MESSAGE(1).ind2=NULL;
		BASEMT_MESSAGE(1).niveauMax=niveauMax;

		//Indice ou placer la reponse dans le tableau de reponse
		BASEMT_MESSAGE(1).indice1	= -1;
		BASEMT_MESSAGE(1).indice2	= -1;	 			
		
		//Initialise la structure associe pour le retour de résultat
		BASEMT_MESSAGE(1).ret = (double*) memalloc(NiveauMax+1,sizeof(double)); 
		BASEMT_MESSAGE(1).elem.Initialise(niveauMax, BASEMT_MESSAGE(1).ret);
		//Initialisation
		for(int a=0;a<=niveauMax;a++)
			BASEMT_MESSAGE(1).ret[a]=0.0;

	BASEMT_FINBOUCLE_INITMESSAGE(1)

	//Iterateur pour toute la matrice
	const int NItem	 = ((NProposant*NProposant)-NProposant)/2;
	const int taillematrix	 = NProposant*NProposant;
	
	//Barre de progression
	//Tres important.. initialise la sémaphore
	Kinship4Struct::InitMT();
	CREATE_PROGRESS_BAR_MATRIX(NProposant,printprogress)
	for(int cPro1=0;cPro1<NProposant;++cPro1)
	{
		for(int cPro2=cPro1;cPro2<NProposant;++cPro2)
		{
			/*
			if (cPro2==cPro1)
			{
				for(int a=0;a<LongEcart;a++)
					MatrixArray[(a*taillematrix)+cPro1*NProposant+cPro2]=0.5;					
			}
			else
			{*/
				//Partie 2: EBRANCHER L'ARBRE, ORGANISER LE NOEUD EN VECTEUR EN RESPECTANT PRECEDENCE
  				//Calcul de la valeur de retours
				 
				BASEMT_DEBUT_REQUETE(1) 
					//Résupere le résultat de phi
					if (BASEMT_MESSAGE(1).indice1!=-1)
					{
						//Recupere le résultat de phideep				
						double* retour	= BASEMT_MESSAGE(1).ret;
						const int ind1	= BASEMT_MESSAGE(1).indice1;
						const int ind2	= BASEMT_MESSAGE(1).indice2;
						
						for(int a=0;a<LongEcart;a++)
						{
							if (retour[a]<0.5)
								pdMoyenne[a]+=retour[a+NiveauMin];
							
							MatrixArray[(a*taillematrix)+ind1*NProposant+ind2]=retour[a+NiveauMin];
							MatrixArray[(a*taillematrix)+ind2*NProposant+ind1]=retour[a+NiveauMin];					
						}	
						
						//remise à zero de phi deep
						for(int a=0;a<=niveauMax;a++)
							retour[a]=0.0;
					}
					//Parametre pour un nouveau calcul de phi
					BASEMT_MESSAGE(1).indice1 = cPro1;
					BASEMT_MESSAGE(1).indice2 = cPro2;
					BASEMT_MESSAGE(1).ind1=NoeudPro[cPro1];
					BASEMT_MESSAGE(1).ind2=NoeudPro[cPro2];						
				BASEMT_FIN_REQUETE(1)

				//Affichage des progress
				INCREMENT_PROGRESS_BAR()
			//}			

		}//Fin itérateur proposant 2
		
	}//Fin itérateur proposant 1	

	//Fermeture des threads
	BASEMT_DEBUT_FERMETURE(1)
		//RECUPERE LES DERNIERS RESULTATS DE PHI S'IL SONT VALIDE
		if (BASEMT_MESSAGE(1).indice1!=-1)
		{
			//Recupere le résultat de phideep				
			double* retour	= BASEMT_MESSAGE(1).ret;
			const int ind1	= BASEMT_MESSAGE(1).indice1;
			const int ind2	= BASEMT_MESSAGE(1).indice2;
			
			for(int a=0;a<LongEcart;a++)
			{
				if (retour[a]<0.5)
					pdMoyenne[a]+=retour[a+NiveauMin];
				
				MatrixArray[(a*taillematrix)+ind1*NProposant+ind2]=retour[a+NiveauMin];
				MatrixArray[(a*taillematrix)+ind2*NProposant+ind1]=retour[a+NiveauMin];					
			}	
			
			//remise à zero de phi deep
			for(int a=0;a<=niveauMax;a++)
				retour[a]=0.0;
		}
	BASEMT_FIN_FERMETURE(1)

	//Calcul de la moyenne........
	for(int Niveau=0;Niveau<LongEcart;Niveau++)
		pdMoyenne[Niveau]/=NItem;
	
	//Retourne la sémaphore
	Kinship4Struct::ReleaseMT();
	return 0;
 			} catch(std::exception &ex) {
 				forward_exception_to_r(ex);
 			}
 			  catch(...){
 				::Rf_error("c++ exception (unknown reason)"); 
 			} 
 			return 0;
}

//**********************************************************************************/
//
//					FONCTION RECURSIVE POUR LE KINSHIP
//
//**********************************************************************************/

/*! 
	\brief Calcule l'apparentement entre deux noeud a l'interieur d'une genealogie

		La fonction calcule l'apparentement entre chaque proposant qui lui est fourni.
		Et ce pour chaque profondeur de NiveauMin à NiveauMax soit ((NiveauMax - NiveauMin)+1) Matrices.

	\attention prob[1] doit-être initialise a -1 avant de lancer cette fonction pour la 1e fois.

	\param Ind1				[in] Ptr vers le Noeud de l'individu 1 et fesant partie d'une genealogie valide
	\param Ind2				[in] Ptr vers le Noeud de l'individu 2 et fesant partie d'une genealogie valide

	\param ttl1				[in] Nombre de remonte maximal a partie de l'individu 1
	\param ttl2				[in] Nombre de remonte maximal a partie de l'individu 2
 
	\return La valeur de l'apparentement calculer
	
	\remark Normalement, ttl1 et ttl2 sont égale lorsque l'on lance un calcul
*/
double FASTCALL Kinship(CIndSimul* Ind1,CIndSimul* Ind2,short ttl1,short ttl2)
{
//#define testtrace  
 #ifdef testtrace
    static int deep=0;
    static int profondeur=0;
    static int tableau[22]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  #endif

	if (Ind1==Ind2)
	{				
		double val;
		{
		  
			if (Ind1->mere && Ind1->pere )
			{
				const int max=(ttl1>ttl2?ttl1:ttl2);
				if (max>0)
				{
					#ifdef testtrace
						++deep;
					#endif
					val=Kinship(Ind1->pere,Ind1->mere,max-1,max-1);
					#ifdef testtrace
						--deep;
					#endif
				}
				else
					val=0;					
			}
			else
				val=0;
#ifdef testtrace
    if (deep==0)
    {
        //printf("\nPISTE= ");
        //for(int a=0;a<profondeur;a++)
	       //printf("%d, ",tableau[a]);	
	//printf(" === prof=%d  ->  COSAN=%f  --> VALEUR= %f \n",profondeur+1, val, pow2(profondeur+1)*(1+val) );
    }
#endif
		}			
		return (1+val)*0.5;
	}
		
#ifdef testtrace
    if (deep==0) tableau[profondeur]=Ind1->nom;
#endif

	//Trie improvisé: Toujours le plus bas de la hérarchie en premiers
	if (Ind2->noind > Ind1->noind)
	{
		CIndSimul* tmp=Ind1;
		Ind1=Ind2;
		Ind2=tmp;
		short tmp2=ttl1;
		ttl1=ttl2;
		ttl2=tmp2;
	}      

	if(Ind1->mere==NULL && Ind1->pere==NULL)
	{	
	   return 0;
	}
	
	double valpere=0;
	double valmere=0;
	if (ttl1>0)
	{
		//MERE
		if (Ind1->mere!=NULL )
		{
		  #ifdef testtrace
		       if (deep==0) ++profondeur;
		  #endif
		  valmere=Kinship(Ind1->mere,Ind2,ttl1-1,ttl2);
		  #ifdef testtrace
		       if (deep==0) --profondeur;
		  #endif
		}

		//PERE
		if (Ind1->pere!=NULL)
		{
		  #ifdef testtrace
		       if (deep==0) ++profondeur;
		  #endif		  
		  valpere=Kinship(Ind1->pere,Ind2,ttl1-1,ttl2);
		  #ifdef testtrace
		       if (deep==0) --profondeur;
		  #endif
		}

	}
	return (valpere+valmere)/2.;
}


/**************************************************************************************************************************************
//
//  TEST NOUVELLE VERSION KINSHIP
****************************************************************************************************************************************/

//TABLEAU POUR LES UTILISATION MULTIPROCESSOR


//DEFINITION DES FONCTIONS DE LA CLASSE KINSHIP4STRUCT UTILISÉ PAR KINSHIP4
CSema  Kinship4Struct::m_acces;

Kinship4Struct::Kinship4Struct() {m_acces=NULL;}
Kinship4Struct::Kinship4Struct(short NiveauMax, double* Resultat) {m_acces=NULL;Initialise(NiveauMax,Resultat);}

void Kinship4Struct::InitMT() {CSema_init(m_acces,1);}
void Kinship4Struct::ReleaseMT() {CSema_destroy(m_acces);}

void Kinship4Struct::Initialise(short NiveauMax, double* Resultat) 
{
	this->NiveauMax=NiveauMax;	
	memblock.setTaille(NiveauMax+1);	
	PileCosan.push(Resultat);
}

/*! 
	\brief Calcule l'apparentement entre deux noeud a l'interieur d'une genealogie pour toute une serie de profondeur

		La fonction calcule l'apparentement entre chaque proposant qui lui est fourni.
		Et ce pour chaque profondeur de NiveauMin à NiveauMax soit ((NiveauMax - NiveauMin)+1) Matrices.

	\attention pGen doit etre initialiser a Null avant de lancer cette fonction la premiere fois
			   allele doit etre initialiser a 0 avant de lancer cette fonction la premiere fois
			   bFlagsort doit etre initialiser a 0 avant de lancer cette fonction la premiere fois
	
	\param Ind1				[in] Ptr vers le Noeud de l'individu 1 et fesant partie d'une genealogie valide
	\param Ind2				[in] Ptr vers le Noeud de l'individu 2 et fesant partie d'une genealogie valide

	\param ttl1				[in] Nombre de remonte maximal a partie de l'individu 1
	\param ttl2				[in] Nombre de remonte maximal a partie de l'individu 2
 
	\retval TabResultat		[out] Ptr vers un tableau de (NiveauMax-NiveauMin+1)
									En cas de succes, contient l'apparentement entre ind1 et in2 pour chaque profondeur

	\param NiveauMin		[in] Doit-être egal a NiveauMax de kinship4struct lors du lancement initial
	\param NiveauMax		[in] Doit-être egal a NiveauMax de kinship4struct  lors du lancement initial

	\param InputStruct		[in] Ptr vers une structure Kinship2Struct remplie correctement.
									La structure Kinship2Struct contient un gestionnaire de memoire qui permet d'assigne
									Dynamiquement de la memoire a chaque noeud si c'est necessaire.
									<br>Ex:									<br>
										GestionMemoire ParentFMem(1);       <br> 
										Kinship2Struct iKinship2;<br>
										iKinship2.MemGest	= &ParentFMem;<br>
										iKinship2.Tableau	= NULL;<br>
										iKinship2.ctableau	= LONG_MAX;

	\return (NE PAS UTILISER) Devrais etre la valeur de l'apparentement pour un niveau quelconque
	
	\remark Normalement, ttl1 et ttl2 sont égale a NiveauMax lorsque l'on lance le calcul
			
*/
void FASTCALL Kinship4(CIndSimul* Ind1, CIndSimul* Ind2, short ttl1, short ttl2,Kinship4Struct &T)
{	

	//REQUIS ET PREALABLEE
	//pGen = NULL
	
	if (Ind1==Ind2)	 
	{//Ind1 et Ind2 croisement.. reste a calculer l'apparentement	       				
		const short max=(ttl1>ttl2?ttl1:ttl2);
		const short min=(ttl1<ttl2?ttl1:ttl2);
		const short niveauMax = T.NiveauMax;	
		
		//Si requis calcule le facteur de cosanguinite
		if (Ind1->mere && Ind1->pere && Ind1->pGen==NULL)  
		{					
				//alloue un nouveau block de mémoire
				Ind1->pGen=T.memblock.Alloc();
			
				//Calcule la consanguinité
				T.PileCosan.push(Ind1->pGen); //Fera un autre essais....							
					Kinship4(Ind1->pere,Ind1->mere,niveauMax,niveauMax,T);	//C'est pas niveauMax-1 surtout pas
				T.PileCosan.pop(); //Retour en arriere					
													
		}//if (Ind1->mere && Ind1->pere && pgen==null)

		//Il y a de la consanguinité et on la traite...							
		double *retour = T.PileCosan.top(); //Valeur de retour

		//Ce calcul est correct
		const short DeepMin=niveauMax-max;
		const short DeepMax=niveauMax-min;
		const double tmpvaleur=pow2(DeepMin+DeepMax+1); //calcul de la valeur en fct de la distance
				
		//Pour chaque niveau possible
		//Ajoute la valeur de phi au calcul courant...
		const double Niv_min_touchable = niveauMax-min;		
		if (Ind1->pGen)			
			for(short niv=niveauMax,sim=0;niv>=Niv_min_touchable;niv--,sim--)
			{								
				const short post= (max-1)+sim; //manque peut-être une legere correction...
				if (post>=0)
				{
					assert(Ind1->pGen[post]!=-1);
					retour[niv]+=tmpvaleur*(1+Ind1->pGen[post]); 
				}
				else
					retour[niv]+=tmpvaleur;
			}
		else //ignore la cosanguinite
			for(short niv=niveauMax;niv>=Niv_min_touchable;niv--) 
				retour[niv]+=tmpvaleur;		
		return; 
	}//fin if end1 == end2
	
	//Inversion pour la priorite	
	//Trie improvisé: Toujours le plus bas de la hérarchie en premiers

	if (Ind2->noind > Ind1->noind) //Pour vérifier....
	{
		//ON SWAP
		//si on remonte pas plus loin
		if(Ind2->mere==NULL && Ind2->pere==NULL)	
			return;
			
		if (ttl2>0)
		{
			//MERE
			if (Ind2->mere!=NULL )
					Kinship4(Ind2->mere,Ind1,ttl2-1,ttl1,T);		
						
			//PERE
			if (Ind2->pere!=NULL)
					Kinship4(Ind2->pere,Ind1,ttl2-1,ttl1,T);							
		}
	}
	else
	{	
		//si on remonte pas plus loin
		if(Ind1->mere==NULL && Ind1->pere==NULL)	
			return;
			
		if (ttl1>0)
		{
			//MERE
			if (Ind1->mere!=NULL )
					Kinship4(Ind1->mere,Ind2,ttl1-1,ttl2,T);		
						
			//PERE
			if (Ind1->pere!=NULL)
					Kinship4(Ind1->pere,Ind2,ttl1-1,ttl2,T);							
		}
	}
	return;
}

void FASTCALL Kinship4MT(CIndSimul* Ind1, CIndSimul* Ind2, short ttl1, short ttl2,Kinship4Struct &T)
{	
	//REQUIS ET PREALABLEE
	//pGen = NULL
	
	if (Ind1==Ind2)	 
	{//Ind1 et Ind2 croisement.. reste a calculer l'apparentement	       				
		const short max=(ttl1>ttl2?ttl1:ttl2);
		const short min=(ttl1<ttl2?ttl1:ttl2);
		const short niveauMax = T.NiveauMax;	
		
		//Synchronisation MT
		CSema_wait(T.m_acces);

		//Si requis calcule le facteur de cosanguinite
		if (Ind1->mere && Ind1->pere && Ind1->pGen==NULL)  
		{					
				//alloue un nouveau block de mémoire
				Ind1->pGen=T.memblock.Alloc();
			
				//Calcule la consanguinité
				T.PileCosan.push(Ind1->pGen); //Fera un autre essais....
					//On appelle la version non-mt pour évité deadlock...
					Kinship4(Ind1->pere,Ind1->mere,niveauMax,niveauMax,T);	//C'est pas niveauMax-1 surtout pas
				T.PileCosan.pop(); //Retour en arriere					
													
		}//if (Ind1->mere && Ind1->pere && pgen==null)
		
		//Synchronisation MT
		CSema_post(T.m_acces);

		//Il y a de la consanguinité et on la traite...							
		double *retour = T.PileCosan.top(); //Valeur de retour

		//Ce calcul est correct
		const short DeepMin=niveauMax-max;
		const short DeepMax=niveauMax-min;
		const double tmpvaleur=pow2(DeepMin+DeepMax+1); //calcul de la valeur en fct de la distance
				
		//Pour chaque niveau possible
		//Ajoute la valeur de phi au calcul courant...
		const double Niv_min_touchable = niveauMax-min;		
		if (Ind1->pGen)			
			for(short niv=niveauMax,sim=0;niv>=Niv_min_touchable;niv--,sim--)
			{								
				const short post= (max-1)+sim; //manque peut-être une legere correction...
				if (post>=0)
				{
					assert(Ind1->pGen[post]!=-1);
					retour[niv]+=tmpvaleur*(1+Ind1->pGen[post]); 
				}
				else
					retour[niv]+=tmpvaleur;
			}
		else //ignore la cosanguinite
			for(short niv=niveauMax;niv>=Niv_min_touchable;niv--) 
				retour[niv]+=tmpvaleur;		
		return; 
	}//fin if end1 == end2
	
	//Inversion pour la priorite	
	//Trie improvisé: Toujours le plus bas de la hérarchie en premiers

	if (Ind2->noind > Ind1->noind) //Pour vérifier....
	{
		//ON SWAP
		//si on remonte pas plus loin
		if(Ind2->mere==NULL && Ind2->pere==NULL)	
			return;
			
		if (ttl2>0)
		{
			//MERE
			if (Ind2->mere!=NULL )
					Kinship4MT(Ind2->mere,Ind1,ttl2-1,ttl1,T);		
						
			//PERE
			if (Ind2->pere!=NULL)
					Kinship4MT(Ind2->pere,Ind1,ttl2-1,ttl1,T);							
		}
	}
	else
	{	
		//si on remonte pas plus loin
		if(Ind1->mere==NULL && Ind1->pere==NULL)	
			return;
			
		if (ttl1>0)
		{
			//MERE
			if (Ind1->mere!=NULL )
					Kinship4MT(Ind1->mere,Ind2,ttl1-1,ttl2,T);		
						
			//PERE
			if (Ind1->pere!=NULL)
					Kinship4MT(Ind1->pere,Ind2,ttl1-1,ttl2,T);							
		}
	}
	return;
}
