/***************************
	Fichier qui contient les fonctions très spécifique au hardware
\contributor Jean-François Lefebvre
  ***********************************/
#ifndef GENLIBHAL
#define GENLIBHAL
#define WIN32_LEAN_AND_MEAN 

#if defined _WIN32 || defined _WIN64
 #include <Rcpp.h>
 #include <Windows.h>
#else
 #include <R_ext/RS.h>
#endif
#define R_NO_REMAP
/*FONCTION CLASSIQUE*/
int processorCount();
int thetime();

//FONCTION POUR LES SEMAPHORES
struct Opa_Cema;
typedef Opa_Cema *CSema;
int CSema_init(CSema& Semaphore,int compteinitial);
int CSema_wait(CSema& Semaphore);
int CSema_post(CSema& Semaphore);
int CSema_destroy(CSema& Semaphore);

//FONCTION POUR LE MULTITHREAD

#if defined _WIN32 || defined _WIN64
	#define THREADRETURN DWORD WINAPI  	//unsigned int
	#define THREADONEXIT return 0;
#else	
	#define THREADRETURN void*
	#define THREADONEXIT return NULL;
#endif
struct Opa_Thread;
typedef Opa_Thread *Cthread;
int  Cthread_create(Cthread& thread, THREADRETURN (*start_routine)(void*),void* arg);
void Cthread_destroy(Cthread& thread);
void Cthread_exit();
void Cthread_join(Cthread& thread);



#endif


