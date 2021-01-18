#ifndef VM_CALC_H
#define VM_CALC_H

#include "vm_stack.h"
#include "ptr_table.h"

int vm_calc_addx(vm_stack*, ptr_table** table);
int vm_calc_mulx(vm_stack*);
int vm_calc_subx(vm_stack*);
int vm_calc_powx(vm_stack*);
int vm_calc_modx(vm_stack*);
int vm_calc_divx(vm_stack*);
int vm_calc_factorial(vm_stack*);
int vm_calc_uminus(vm_stack*);

int vm_calc_and(vm_stack*);
int vm_calc_or(vm_stack*);
int vm_calc_eq(vm_stack*);
int vm_calc_neq(vm_stack*);
int vm_calc_gt(vm_stack*);
int vm_calc_lt(vm_stack*);
int vm_calc_ge(vm_stack*);
int vm_calc_le(vm_stack*);
int vm_calc_neg(vm_stack*);
#endif

