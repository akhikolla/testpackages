//contributor Jean-Francois Lefebvre

#include "base.h"
#include "userInterface.h"
#include "hal.h"
#include <float.h>
#include <string.h>
#include <stdio.h> //jfl
#include <cstdlib>
#include <Rcpp.h>
using namespace std;
#define R_NO_REMAP

//Implementation d'un chronometre
static int g_TimerStart=-1;
static int g_TimerLast=-1;
void TimerOnStart() {g_TimerStart=thetime();g_TimerLast=-1;}
void TimerOnStop()	{g_TimerLast =thetime()-g_TimerStart;}
int getLastTimer() {return g_TimerLast;}

//Implémentation d'un système d'arret automatique si une procédure s'annonce trop longue
const int g_tempsEchantillon = 30;		   //Temps d'échantillonnage pour le temps d'execution (seconde)
const int g_sautEchantillon = 3;			   //Nombre d'incrément entre chaque test
static double g_maxEstimatedProcessingTime=20*60;  //En seconde: Valeur par défaut: 20 min

void setCurrentMaxTime(double minute) 
{
	if (minute==0)
		g_maxEstimatedProcessingTime=DBL_MAX;
	else
		g_maxEstimatedProcessingTime=minute*60;
}
void getCurrentMaxTime(double* minute) 
{
	if (g_maxEstimatedProcessingTime==DBL_MAX)
		*minute=0;
	else
		*minute=g_maxEstimatedProcessingTime/60;
}



//Implémentation d'une progress bar en mode texte...
const char  ProgresBARLength=62;
//const char  sProgressMarge[]="   ";
//const char  sProgressBAR[]  ="0..........20..........40...........60..........80..........100";
//const char  sProgresSymbol[]="XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
const int  CTextProgressBar_Tampon = 2;

CTextProgressBar::CTextProgressBar(SXLONG max,int affiche) 
{
	m_max=max;
	m_status=FIRST;
	m_Last=0;
	m_cur=0;
	m_afficheBar=affiche;
}

//L'opérateur magique
void CTextProgressBar::operator++() 
{
	try{
	switch(m_status)
	{
	case FIRST:
		m_status=ADD;		
		++m_cur;
		//DEBUT
		/*if (m_afficheBar)
		{
			//printf("\n%s%s\n%s",sProgressMarge,sProgressBAR,sProgressMarge); 
			fflush(stdout);		
		}*/
		//DÉMARRÉ CHRONOMETRE
		m_tempsdebut=thetime()+g_tempsEchantillon;

	case ADD:
		{
			++m_cur;			
			//TEST SI LE TEMPS EST TROP LONG
			int tst;
			if (m_tempsdebut!=-1 && m_cur%g_sautEchantillon==0 && (tst=thetime())>m_tempsdebut)
			{				
				//Temps estimer				
				const double testime=(double(m_max)/double(m_cur))*((tst-m_tempsdebut)+g_tempsEchantillon+1); //1 marge
				m_tempsdebut=-1;			
				if (testime>g_maxEstimatedProcessingTime){
//					GENError("Execution time exceeded maximum allowed: ESTIMATED: %d min MAXIMUM: %d min\nSee gen.maxexetime() definition",
//							int(testime/60),int(g_maxEstimatedProcessingTime/60));
					char erreur[TAILLEDESCRIPTION];
					sprintf(erreur, "Execution time exceeded maximum allowed: ESTIMATED: %d min MAXIMUM: %d min\nSee gen.maxexetime() definition",
							int(testime/60),int(g_maxEstimatedProcessingTime/60));
					throw std::range_error(erreur);
					//GENError("Le temps d'execution estimé est trop long: ESTIME: %d min   MAXIMUM: %d min"\n Regarde la definition de la fonction gen.maxexetime()",
				}
			}
			if (m_afficheBar)
			{
				//Affichage des progress			
				const int CurrentPos  = int((m_cur*ProgresBARLength)/m_max);				
				if ((CurrentPos-m_Last)>=CTextProgressBar_Tampon)
				{
					  //printf("%s",(sProgresSymbol+ProgresBARLength-CurrentPos+m_Last+1));
					  //fflush(stdout);
					  m_Last=CurrentPos;
				}
				if (m_cur>=m_max)
				{
					m_status=FINISH;
					//printf("\n");
				}
			}
		}			
	case FINISH: // Rien faire
	{}
	}//fin du switch
	} catch(...){
		::Rf_error("c++ exception (unknown reason)"); 
	}
}

//Implémentation d'une progress bar en mode texte...
//Variante pour des constantes floatante
CTextProgressBarFloat::CTextProgressBarFloat(double max,double* curseur,int affiche)
{	
	m_status=FIRST;
	m_Last=0;
	smallit=0;
	m_afficheBar=affiche;
	
	m_max=max;
	m_pcur=curseur;	
}

void CTextProgressBarFloat::End()
{
	//Finis la barre de tache
	m_pcur=&m_max;
	++(*this);
}
//L'opérateur magique
void CTextProgressBarFloat::operator++() 
{
	try{
	switch(m_status)
	{
	case FIRST:
		m_status=ADD;		
		//DEBUT
		/*if (m_afficheBar)
		{
			//printf("\n%s%s\n%s",sProgressMarge,sProgressBAR,sProgressMarge); 
			fflush(stdout);		
		}*/
		//DÉMARRÉ CHRONOMETRE
		m_tempsdebut=thetime()+g_tempsEchantillon;

	case ADD:
		{					
			//TEST SI LE TEMPS EST TROP LONG
			int tst;
			++smallit;
			if (m_tempsdebut!=-1 && smallit%g_sautEchantillon==0 && (tst=thetime())>m_tempsdebut)
			{				
				//Temps estimer				
				const double testime=m_max/ (*m_pcur) * ((tst-m_tempsdebut)+g_tempsEchantillon+1); //1 marge
				m_tempsdebut=-1;			
				if (testime>g_maxEstimatedProcessingTime){			
//					GENError("Execution time exceeded maximum allowed: ESTIMATED: %.10G min MAXIMUM: %.10G min\nSee gen.maxexetime() definition",
//						testime/60.,g_maxEstimatedProcessingTime/60.);
					//GENError("Le temps d'execution estimé est trop long: ESTIME: %.10G min   MAXIMUM: %.10G min\n Regarde la definition de la fonction gen.maxexetime()",
					char erreur[TAILLEDESCRIPTION];
					sprintf(erreur, "Execution time exceeded maximum allowed: ESTIMATED: %.10G min MAXIMUM: %.10G min\nSee gen.maxexetime() definition",
								  testime/60.,g_maxEstimatedProcessingTime/60.);
					throw std::range_error(erreur);
				}
			}
			if (m_afficheBar)
			{
				//Affichage des progress						
				const int CurrentPos  = int( *m_pcur / m_max*ProgresBARLength );									
				if ((CurrentPos-m_Last)>=CTextProgressBar_Tampon)
				{
					  //printf("%s",(sProgresSymbol+ProgresBARLength-CurrentPos+m_Last+1));
					  //fflush(stdout);
					  m_Last=CurrentPos;
				}
				if (*m_pcur>=m_max)
				{
					m_status=FINISH;
					//printf("\n");
				}
			}
		}			
	case FINISH: // Rien faire
	{}
	}//fin du switch
	} catch(...){
 		::Rf_error("c++ exception (unknown reason)"); 
 	}
}
