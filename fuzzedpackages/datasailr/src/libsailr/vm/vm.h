#ifndef VM_H
#define VM_H

#include "vm_code.h"
#include "ptr_table.h"
#include "vm_stack.h"
#include "func/ext_func/ext_func_hash.h"

int vm_exec_code( vm_inst* code , int num_insts ,int start_inst_idx, ptr_table* table , vm_stack* vmstack, ext_func_hash* extfunc_hash);

#endif /* VM_H */
