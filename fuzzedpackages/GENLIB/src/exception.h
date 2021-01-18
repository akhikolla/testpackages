
#ifndef GENEXCEPTION
#define GENEXCEPTION

///Fonction pour récupéré le dernier message d'erreur
const char* getLastMessage();


//Utilisation interne
const int TAILLEDESCRIPTION=1024;
extern char g_LastMessage[];
void ErrorHandler();


//DECLARATION DES HANDLERS D'EXCEPTION
template<class T1,class T2,class T3,class T4> 
void GENError(const char* format,T1 item1,T2 item2,T3 item3,T4 item4)
{
	char buffer[TAILLEDESCRIPTION];
	sprintf(buffer,"%s\n",format);
	sprintf(buffer,format,item1,item2,item3,item4);
	sprintf(g_LastMessage,"\nError: %s \n",buffer);
	ErrorHandler();
}
template<class T1,class T2,class T3> 
void GENError(const char* format,T1 item1,T2 item2,T3 item3)
{
	char buffer[TAILLEDESCRIPTION];
	sprintf(buffer,"%s\n",format);
	sprintf(buffer,format,item1,item2,item3);
	sprintf(g_LastMessage,"\nError: %s \n",buffer);
	ErrorHandler();
}
template<class T1,class T2> 
void GENError(const char* format,T1 item1,T2 item2)
{
	char buffer[TAILLEDESCRIPTION];
	sprintf(buffer,"%s\n",format);
	sprintf(buffer,format,item1,item2);
	sprintf(g_LastMessage,"\nError: %s \n",buffer);
	ErrorHandler();
}
template<class T1> 
void GENError(const char* format,T1 item1)
{
	char buffer[TAILLEDESCRIPTION];
	sprintf(buffer,"%s\n",format);
	sprintf(buffer,format,item1);
	sprintf(g_LastMessage,"\nError: %s \n",buffer);
	ErrorHandler();
}
void GENError(const char* format);


#endif


