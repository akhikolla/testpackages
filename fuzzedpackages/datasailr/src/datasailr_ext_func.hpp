#ifndef DATASAILR_EXT_FUNC_H
#define DATASAILR_EXT_FUNC_H

extern "C" {
#include "sailr.h"
#include "sailr_ext.h"
}

int sailr_external_push_row( arg_list* arglist , unsigned int num_args, vm_stack* vmstack );
int sailr_external_discard_row( arg_list* arglist , unsigned int num_args, vm_stack* vmstack );

// Test functions
int sailr_external_println( arg_list* arglist , unsigned int num_args, vm_stack* vmstack);
int sailr_external_add2( arg_list* arglist , unsigned int num_args, vm_stack* vmstack);

#endif // DATASAILR_EXT_FUNC_H
