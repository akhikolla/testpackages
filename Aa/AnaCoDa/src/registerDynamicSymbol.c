// RegisteringDynamic Symbols
//http://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_AnaCoDa(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

