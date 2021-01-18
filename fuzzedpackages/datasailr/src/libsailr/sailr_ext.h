#ifndef SAILR_EXT_H
#define SAILR_EXT_H

#include "vm/vm.h"
#include "vm/vm_stack.h"
#include "string/common_string.h"
#include "simple_re/simple_re.h"
#include "simple_date/simple_date.h"
#include "vm/func/ext_func/ext_func_hash.h"
#include "sailr.h"

// External Functions

ext_func_hash_object* sailr_ext_func_hash_init();
void sailr_ext_func_hash_add(ext_func_hash_object** hash, const char* fname, unsigned int num_args, int (* func)(arg_list*, unsigned int, vm_stack*));
void sailr_ext_func_hash_free(ext_func_hash_object** hash);

const char* sailr_ext_func_hash_get_last_executed(ext_func_hash_object** hash);
void sailr_ext_func_hash_reset_last_executed(ext_func_hash_object** hash);

void sailr_ext_vm_stack_end(vm_stack_object* vmstack);

#endif

