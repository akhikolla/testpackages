#ifndef TIMER_H
#define TIMER_H

#include <ctime>

#ifdef _WIN64
#undef VOID
// AMR: We have to make some adjustments for a Windows compilation
// http://www.mzan.com/article/11588765-using-rcpp-with-windows-specific-includes.shtml
#ifdef Realloc
#undef Realloc
#endif
#ifdef Free
#undef Free
#endif
#include <windows.h>
#elif _WIN32
#undef VOID
// AMR: We have to make some adjustments for a Windows compilation
// http://www.mzan.com/article/11588765-using-rcpp-with-windows-specific-includes.shtml
#ifdef Realloc
#undef Realloc
#endif
#ifdef Free
#undef Free
#endif
#include <windows.h>
#elif __APPLE__
    #include "TargetConditionals.h"
    #if TARGET_OS_MAC
    #include <unistd.h>
    #include <sys/time.h>
    #include <netinet/in.h>
    #endif
#elif __linux
#include <unistd.h>
#include <sys/time.h>
#include <netinet/in.h>
#endif // end of checking platforms

#define CLK_PS 1000
//#define CLK_PS CLOCKS_PER_SEC

/*Define the namespace*/
namespace ROPTLIB{

	unsigned long getTickCount();
} /*end of ROPTLIB namespace*/
#endif
