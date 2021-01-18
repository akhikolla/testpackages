#ifndef GEN_CODE_UTIL_H
#define GEN_CODE_UTIL_H

#include "vm/vm_code.h"

typedef vm_inst vm_inst_list;

vm_inst* new_vm_inst_command( VM_CMD );
vm_inst* new_vm_inst_push_ival( int );
vm_inst* new_vm_inst_push_dval( double );
vm_inst* new_vm_inst_push_pp_ival( char* ptr_key );
vm_inst* new_vm_inst_push_pp_dval( char* ptr_key );
vm_inst* new_vm_inst_push_pp_str( char* ptr_key );
vm_inst* new_vm_inst_push_pp_rexp( char* ptr_key );
vm_inst* new_vm_inst_push_null( char* ptr_key );
vm_inst* new_vm_inst_label( char* label );
vm_inst* new_vm_inst_fjmp( char* label );
vm_inst* new_vm_inst_jmp( char* label );

vm_inst* vm_inst_list_last( vm_inst_list* );
vm_inst* vm_inst_list_reach_last( vm_inst_list* );
vm_inst* vm_inst_list_get(vm_inst_list* , int );
int vm_inst_list_size(vm_inst_list* );
vm_inst_list* vm_inst_list_cat( vm_inst_list* , vm_inst_list*  );
void vm_inst_list_show_all( vm_inst_list* );
// vm_inst_list* vm_inst_list_insert( vm_inst_list*, int, vm_inst_list* );
int vm_inst_list_free( vm_inst_list* );

vm_inst* vm_inst_list_to_code( vm_inst_list* );

void vm_inst_set_loc_to_last( struct script_loc, vm_inst* );
// Used to set script location to vm instruction

#endif /* GEN_CODE_UTIL */
