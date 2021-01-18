
#ifndef GENBASEMT
#define GENBASEMT

#include "hal.h"

/*
	EXEMPLE D'UTILISATION voir phiMatrixMT

struct maclasse : public CMtGlobalMessage
{

};

//Classe contenant les éléments nécessaire a l'appel d'une fonction de calcul
ainsi que le retour des résultats

BASEMT_CREATE_GLOBALMESSAGE(maclasse,1)

//DANS LA FONCTION DE CALCUL 
	//Initialisation du multithread
	BASEMT_DEBUTBOUCLE_INITMESSAGE(1) 
		//utiliser BASEMT_MESSAGE(1) pour accédé a l'intance maclasse courante
	BASEMT_FINBOUCLE_INITMESSAGE(1)

	//Dans la boucle de service
		BASEMT_DEBUT_REQUETE(1) 
			//utiliser BASEMT_MESSAGE(1) pour accédé a l'intance maclasse courante
			//sauvegarde la solution courante si valide
			//Attribut un nouveau calcul
		BASEMT_FIN_REQUETE(1)

	BASEMT_DEBUT_FERMETURE(1)
	  //utiliser BASEMT_MESSAGE(1) pour accédé a l'intance maclasse courante
		//sauvegarde le résultat du dernier calcul
	BASEMT_FIN_FERMETURE(1)

//POUR CREER LA FONCTION HELPER
	BASEMT_DEBUT_HELPERFCT(maclasse,1)
		//utiliser BASEMT_HLPMES pour accede a l'intance ma classe actuel
		//Lance le calcul et recupere le résultat
	BASEMT_FIN_HELPERFCT() \
*/




//Tas de macro utile
//DEVRAIS CORRESPONDRE AU NOMBRE DE CHAMP ptr de CIndSimul
const int MAXPROCESSOR = 6; 


//THREADRETURN PhiMatrixThreadHelper(void* ptr);
enum PM_STATUS { PM_GO, PM_WAIT,PM_END };
struct CMtGlobalMessage
{
  	
  //TRUC STANDARS POUR UN ALGORITHME MULTITHREAD
  PM_STATUS control;
  CSema mtArret;	//Pour indique la position d'arret du thread
  CSema mtData;		//Est-ce que touche est disponible
  CSema *sm;		//Semaphone pour la demande de resource au parent
};


#define BASEMT_CREATE_GLOBALMESSAGE(NomClasseDerive,IdInstance)  \
	static NomClasseDerive g_Message##IdInstance[MAXPROCESSOR];\
	static CSema g_smMustGo##IdInstance; \
	static THREADRETURN ThreadHelper##IdInstance(void* ptr);



#define BASEMT_DEBUTBOUCLE_INITMESSAGE(IdInstance) \
	/*Calcul du nombre de processus à utilisé*/ \
	int BASEMTThreadCount = processorCount(); /*4;*/\
	if (BASEMTThreadCount>MAXPROCESSOR)\
		BASEMTThreadCount = MAXPROCESSOR;\
	Cthread *BASEMTThreadId = (Cthread*) memalloc(BASEMTThreadCount,sizeof(Cthread));\
	/*INITIATION DU SYSTEME MULTITHREAD	*/ \
	{/*int goNextBoucle=-1;*/\
	CSema_init(g_smMustGo##IdInstance,BASEMTThreadCount);\
	for(int i=0;i<BASEMTThreadCount;i++)\
	{ \
	  /*Initialisation de la partie controle de la structure*/ \
	  CSema_init(g_Message##IdInstance[i].mtData,1);\
	  CSema_init(g_Message##IdInstance[i].mtArret,0);\
	  g_Message##IdInstance[i].sm=&g_smMustGo##IdInstance;\
	  g_Message##IdInstance[i].control=PM_WAIT;\
	  Cthread_create(BASEMTThreadId[i],ThreadHelper##IdInstance,&(g_Message##IdInstance[i]));

#define BASEMT_MESSAGE(IdInstance) g_Message##IdInstance[i]	  

#define BASEMT_FINBOUCLE_INITMESSAGE(IdInstance) }/*finfor*/}/*finblock*/

	  
#define BASEMT_DEBUT_REQUETE(IdInstance) \
		{CSema_wait(g_smMustGo##IdInstance); /*count = nombre de thread requerant service*/ \
		int goNextBoucle;\
		do /*En cas de desync*/\
		{\
			goNextBoucle=-1;\
			for(int i=0;i<BASEMTThreadCount && goNextBoucle==-1;i++)\
			{\
			  CSema_wait(g_Message##IdInstance[i].mtData);	/*protection control*/ \
			  if (g_Message##IdInstance[i].control==PM_WAIT) \
			  {	

#define BASEMT_FIN_REQUETE(IdInstance) \
				/*redémarrage du thread*/\
				g_Message##IdInstance[i].control=PM_GO;\
				goNextBoucle = i;\
				CSema_post(g_Message##IdInstance[i].mtArret);\
				}/*fin if wait*/ \
			  CSema_post(g_Message##IdInstance[i].mtData); /*relache control*/ \
			} /*Fin du pour chaque thread*/ \
		  }\
		  while(goNextBoucle==-1);}





#define BASEMT_DEBUT_FERMETURE(IdInstance) \
	/*RECUPERATION DES DERNIERES DONNES ET FERMETURE DES THREAD*/\
	{int goNextBoucle=-1;\
	for(int eachi=0;eachi<BASEMTThreadCount;eachi++)\
	{\
		CSema_wait(g_smMustGo##IdInstance); /*count = nombre de thread requerant service*/ \
		goNextBoucle=-1;\
		for(int i=0;i<BASEMTThreadCount && goNextBoucle==-1;i++)\
		{\
			CSema_wait(g_Message##IdInstance[i].mtData);	/*Bloque control*/\
			if (g_Message##IdInstance[i].control==PM_WAIT)\
			{


#define BASEMT_FIN_FERMETURE(IdInstance) \
				/*TERMINE LA THREAD EN QUESTION & JOIN LE PROCESSUS*/ \
				g_Message##IdInstance[i].control=PM_END;\
				goNextBoucle = i;\
				CSema_post(g_Message##IdInstance[i].mtArret);  /*REDEMARRE thread*/\
				Cthread_join(BASEMTThreadId[i]);	\
			}\
			CSema_post(g_Message##IdInstance[i].mtData);  /*relache control */ \
		} /*Fin du pour chaque thread*/ \
	} /* fin du pour chaque thread*/ \
	CSema_destroy(g_smMustGo##IdInstance);\
	for(int i=0;i<BASEMTThreadCount;i++) \
	{\
		Cthread_destroy(BASEMTThreadId[i]);\
		CSema_destroy(g_Message##IdInstance[i].mtArret);\
		CSema_destroy(g_Message##IdInstance[i].mtData);\
	}\
	}/*fin du block*/

#define BASEMT_DEBUT_HELPERFCT(NomClasseDerive,IdInstance) \
	THREADRETURN ThreadHelper##IdInstance(void* psMes)\
	{\
	  NomClasseDerive &t = *( (NomClasseDerive*)psMes) ;\
	  do\
	  {\
		CSema_wait(t.mtArret);\
		switch(t.control)\
		{\
		case PM_GO:

#define BASEMT_HLPMES t

#define BASEMT_FIN_HELPERFCT() \
		/*Mise en attente*/\
		CSema_wait(t.mtData);\
		t.control=PM_WAIT; /*mise en attente*/\
		CSema_post(t.mtData);\
		CSema_post(*t.sm); /*Incremente la semaphore comme un thread requere du service*/ \
		break;\
		case PM_END:\
			Cthread_exit();\
		break;\
		case PM_WAIT:\
			;\
		}\
	  }\
	  while (1); \
	  THREADONEXIT;\
	} /*fin fct*/

#endif

