//contributor Jean-Francois Lefebvre

#ifndef GENUSERINTERFACE
#define GENUSERINTERFACE

#ifdef USESTIMER
	#define STARTTIMER TimerOnStart();
	#define STOPTIMER  TimerOnStop();
#else
	#define STARTTIMER 
	#define STOPTIMER  
#endif



//CHRONOMETRE POUR LA GESTION DU TEMPS
void TimerOnStart(); 
void TimerOnStop();	
int getLastTimer();

//PROGRESSBAR
//Maintenant la progress bar inclus un chronometre qui arrete 
//la procédure si remarque qu'un certain temps est sur le points d'etre dépassé
class CTextProgressBar
{
private:
	//Curseur
	int m_afficheBar;			//Si !=0 alors affiche la barre de progression
	int m_Last;	
	SXLONG m_max;
	enum {FIRST,ADD,FINISH} m_status;
	SXLONG m_cur;
	

	//Pour l'estimation du temps
	int m_tempsdebut;
		
public:
	//Constructeur
	CTextProgressBar(SXLONG max,int affiche);

	//Utilisation
	void operator++();	
};

class CTextProgressBarFloat
{
private:
	//Curseur
	int m_afficheBar;			//Si !=0 alors affiche la barre de progression
	int m_Last;	
	enum {FIRST,ADD,FINISH} m_status;

	//Pour l'intervalle
	double m_max;	
	double *m_pcur; //Pointeur vers la valeur courante
	
	//Pour l'estimation du temps
	int m_tempsdebut;
	SXLONG smallit;
		
public:
	//Constructeur
	CTextProgressBarFloat(double max,double* curseur,int affiche);
	void End();
	//Utilisation
	void operator++();	
};


//Raccourci pour dans les fonctions seulement si ALLOWPRINTPROGRESS EST DEFINI
const double PR_MAXINTERVAL = 50000;
const double PR_MAXNOINTERVAL = 200000;
#ifdef ALLOWPRINTPROGRESS
	#define CREATE_FLOAT_PROGRESS_BAR(maxValue,Valueptr,mustprint) \
			CTextProgressBarFloat Prbar((maxValue),(Valueptr),(mustprint));

	#define CREATE_PROGRESS_BAR(nbelem,mustprint) \
			const SXLONG PrMaxInterval= (SXLONG) MIN( ceil(double(nbelem)/PR_MAXINTERVAL) , PR_MAXNOINTERVAL ); \
			SXLONG PrInterval = 0;\
			CTextProgressBar Prbar(nbelem/PrMaxInterval,mustprint);		
#else
	#define CREATE_FLOAT_PROGRESS_BAR(maxValue,Valueptr,mustprint) \
			CTextProgressBarFloat Prbar((maxValue),(Valueptr),0);
	
	#define CREATE_PROGRESS_BAR(nbelem,mustprint) \
			const SXLONG PrMaxInterval= (SXLONG) ceil(double(nbelem)/PR_MAXINTERVAL); \
			SXLONG PrInterval = 0;	\
			CTextProgressBar Prbar(nbelem/PrMaxInterval,0);
#endif

#define INCREMENT_FLOAT_PROGRESS_BAR() ++Prbar;
#define END_FLOAT_PROGRESS_BAR() Prbar.End();

#define INCREMENT_PROGRESS_BAR()	\
			if ((++PrInterval)==PrMaxInterval) {++Prbar;PrInterval=0;}

#define CREATE_PROGRESS_BAR_MATRIX(nbelem,mustprint) CREATE_PROGRESS_BAR( (((nbelem*nbelem)-nbelem)/2) ,mustprint)



//GESTION DU TEMPS ET DE L'ARRET AUTOMATIQUE DE LA PROCEDURE
void setCurrentMaxTime(double min);
void getCurrentMaxTime(double* min);


#endif

