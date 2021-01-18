/***************************
	Fichier qui contient les fonctions très spécifique au hardware


\contributor Jean-François Lefebvre

  ***********************************/
#include "hal.h"
#include <time.h>	

/* VERSION WINDOWS  */
#if defined _WIN32 || defined _WIN64
	/*#define _CRTDBG_MAP_ALLOC
	#include<stdio.h>
	#include<stdlib.h>
	#include<crtdbg.h> */
	#include <cstdlib>
	#include <unistd.h>
	
	int thetime()
	{
		return time(NULL);
	}

	int processorCount()
	{
	   SYSTEM_INFO siSysInfo;
	   GetSystemInfo(&siSysInfo); 
	   return siSysInfo.dwNumberOfProcessors;
	}

	
	//MULTITHREAD
	struct Opa_Thread
	{
		HANDLE tid;
	};

	int Cthread_create(Cthread& thread, THREADRETURN (*start_routine)(void*),void* arg)
	{
		DWORD dwThreadId = 1;
		thread = (Opa_Thread*) malloc(sizeof(Opa_Thread));
		thread->tid = CreateThread( 
			NULL,                        // default security attributes 	
			0,                           // use default stack size  		
			(LPTHREAD_START_ROUTINE) start_routine,               // thread function 			
			arg,		                 // argument to thread function 	
			0,                           // use default creation flags 	
			&dwThreadId);                // returns the thread identifier 	
			//NULL);                // returns the thread identifier 
		return 0;
			
	}

	void Cthread_destroy(Cthread& thread)
	{	
		free(thread);
		thread=NULL;
	}

	void Cthread_exit()
	{	
		ExitThread(0);	
	}

	void Cthread_join(Cthread& thread)
	{
		WaitForSingleObject(thread->tid, INFINITE);
	}


	//ENCAPSULATION D'UNE SEMAPHORE

	struct Opa_Cema
	{
		HANDLE sema;
	};
	int CSema_init(CSema& Semaphore,int compteinitial)
	{
		Semaphore= (Opa_Cema*) malloc(sizeof(Opa_Cema));
		Semaphore->sema = CreateSemaphore( 
		NULL,   // no security attributes
		compteinitial,   // initial count
		99999,   // maximum count
		NULL);  // unnamed semaphore
		return 0;
	}

	int CSema_destroy(CSema& Semaphore)
	{
		CloseHandle(Semaphore->sema);
		free(Semaphore);
		Semaphore=NULL;
		return 0;
	}

	int CSema_wait(CSema& Semaphore)
	{
		
		WaitForSingleObject(
			Semaphore->sema,
			INFINITE
			);
		return 0;
	}

	int CSema_post(CSema& Semaphore)
	{
		ReleaseSemaphore( 
			Semaphore->sema,  // handle to semaphore
			1,           // increase count by one
			NULL);       // not interested in previous count
			
		return 0;
	}

	//CAS SPECIAL WINDOWS (EQUIVALENT DU MAIN)

	int FlushGenealogie();
	BOOL APIENTRY DllMain( HANDLE hModule, 
						   DWORD  ul_reason_for_call, 
						   LPVOID lpReserved
						 )
	{
		switch (ul_reason_for_call)
		{
			
			case DLL_THREAD_ATTACH:				  
			case DLL_THREAD_DETACH:
					break;
			case DLL_PROCESS_ATTACH:							
					FlushGenealogie();
					break;
			case DLL_PROCESS_DETACH:
					//Evite un gros memory leak sur windows 95-98
					FlushGenealogie();
					//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_CHECK_CRT_DF | _CRTDBG_DELAY_FREE_MEM_DF | _CRTDBG_LEAK_CHECK_DF);					
					//_CrtDumpMemoryLeaks();
				break;
		}
		return TRUE;
	}


/*VERSION UNIX */
#else

	#include <stdlib.h>
//	#include <iostream>
	#include <unistd.h>
	#include <sys/unistd.h> 
	//#include <sys/processor.h> 
	#include <sys/types.h> 
	#include <stdio.h> 
	#include <cstdio> 
	#include <errno.h> 

	#include <pthread.h>
	#include <semaphore.h>

	int thetime()
	{
		return time(NULL);
	}

	int processorCount()
	{
	  int n = sysconf(_SC_NPROCESSORS_ONLN);  //Seulement ca devrais suffire normalement
	  return (n<=0?1:n);
	}

	//MULTITHREAD
	/*struct Opa_Thread
	{
		thread_t tid;
	};
	int Cthread_create(Cthread& thread, THREADRETURN (*start_routine)(void*),void* arg)
	{
		thread= (Opa_Thread*) malloc(sizeof(Opa_Thread));
		return thr_create(NULL,NULL,start_routine,arg,NULL,&(thread->tid));
	}

	void Cthread_destroy(Cthread& thread)
	{	
		free(thread);
		thread=NULL;
	}

	void Cthread_exit()
	{	
		thr_exit(NULL);
	}

	void Cthread_join(Cthread& thread)
	{		
		thr_join(thread->tid,NULL,NULL); 
	}

	//ENCAPSULATION D'UNE SEMAPHORE//
	struct Opa_Cema
	{
		sema_t sema;
	};
	int CSema_init(CSema& Semaphore,int compteinitial)
	{
		Semaphore= (Opa_Cema*) malloc(sizeof(Opa_Cema));
		return sema_init(&Semaphore->sema,compteinitial,USYNC_PROCESS,NULL);		
	}

	int CSema_destroy(CSema& Semaphore)
	{
		int tmp = sema_destroy(&Semaphore->sema);
		free(Semaphore);
		Semaphore=NULL;
		return tmp;
	}

	int CSema_wait(CSema& Semaphore)
	{
		return sema_wait(&Semaphore->sema);
	}

	int CSema_post(CSema& Semaphore)
	{
		return sema_post(&Semaphore->sema);
	}
*/



	//MULTITHREAD
	struct Opa_Thread
	{
		pthread_t tid;
	};
	int Cthread_create(Cthread& thread, THREADRETURN (*start_routine)(void*),void* arg)
	{
		thread= (Opa_Thread*) malloc(sizeof(Opa_Thread));
		return pthread_create(&(thread->tid), NULL, start_routine, arg);
	}

	void Cthread_destroy(Cthread& thread)
	{	
		free(thread);
		thread=NULL;
	}

	void Cthread_exit()
	{	
		pthread_exit(NULL);
	}

	void Cthread_join(Cthread& thread)
	{		
		pthread_join(thread->tid, NULL);		
	}

	//ENCAPSULATION D'UNE SEMAPHORE
	struct Opa_Cema
	{
		sem_t sema;
	};
	int CSema_init(CSema& Semaphore,int compteinitial)
	{
		Semaphore= (Opa_Cema*) malloc(sizeof(Opa_Cema));
		return sem_init(&Semaphore->sema,1,compteinitial);		
	}

	int CSema_destroy(CSema& Semaphore)
	{
		int tmp = sem_destroy(&Semaphore->sema);
		free(Semaphore);
		Semaphore=NULL;
		return tmp;
	}

	int CSema_wait(CSema& Semaphore)
	{
		return sem_wait(&Semaphore->sema);
	}

	int CSema_post(CSema& Semaphore)
	{
		return sem_post(&Semaphore->sema);
	}


#endif

/*TOUTE VERSION*/


