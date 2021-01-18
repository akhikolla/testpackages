#ifndef VM_CALL_FUNC_H
#define VM_CALL_FUNC_H

#include "vm_stack.h"

int call_func(vm_stack* vmstack, char* fname, int num_args, ptr_table** table);

#endif /* CALL_FUNC_H */
