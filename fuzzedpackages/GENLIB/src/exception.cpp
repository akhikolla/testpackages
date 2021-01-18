
#include "base.h"
#include "exception.h"
#include <stdlib.h>
#include <cstdlib>
#include <string.h>
#include <Rcpp.h>
using namespace Rcpp;
#define R_NO_REMAP
/*

N.b La motivation d'utilise un systeme d'erreur par ErrorHandler au lieu d'un systeme par try..catch

c'est que les exception C++ sont tres mal gere par Splus version Unix
Pour etre precis elle sont completement ignorer

*/


//Variable globale pour regle le probleme d'interception d'exception sous Solaris


char g_LastMessage[TAILLEDESCRIPTION+13];
const char* getLastMessage() {return g_LastMessage;}

void ErrorHandler()
{
#ifdef MODETEST
	//printf("\n%s",getLastMessage() );
	//abort();
#else
	//PROBLEM "%s",getLastMessage();
     //RECOVER(NULL_ENTRY); 
#endif 
}

void GENError(const char* format)
{
	char buffer[TAILLEDESCRIPTION];
	sprintf(buffer, "%s\n", format);
	sprintf(g_LastMessage,"\nError: %s \n",buffer);
	ErrorHandler();
}
