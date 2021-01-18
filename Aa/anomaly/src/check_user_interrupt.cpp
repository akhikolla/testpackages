#include <R.h>
#include <Rinternals.h>

static void check_user_interrupt_handler(void* dummy) {R_CheckUserInterrupt();}

bool check_user_interrupt() {return R_ToplevelExec(check_user_interrupt_handler,NULL) == FALSE;}
