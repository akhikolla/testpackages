#ifndef VM_ERROR_H
#define VM_ERROR_H

#include "vm_stack.h"
#include <stdbool.h>

void vm_error_raise(vm_stack* vmstack);
unsigned int vm_error_num(vm_stack* vmstack);
bool vm_error_exist(vm_stack* vmstack);

#endif /* VM_ERROR_H */
