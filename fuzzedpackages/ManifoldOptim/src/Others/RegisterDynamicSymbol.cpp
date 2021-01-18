/* 
 * Date Created: 3/31/17
 * Description:
 * This is an internal R routine to register functions internally.*/

// RegisteringDynamic Symbols
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
void R_init_ManifoldOptim(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}
