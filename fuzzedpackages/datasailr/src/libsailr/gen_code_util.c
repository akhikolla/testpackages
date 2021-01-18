#include <R_ext/Print.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "gen_code_util.h"
#include "vm/vm_code.h"
#include "helper.h"

vm_inst*
new_vm_inst_command( VM_CMD cmd )
{
	vm_inst* temp_inst = (vm_inst*)malloc(sizeof(vm_inst));
	temp_inst->cmd = cmd;
	temp_inst->prev = NULL;
	temp_inst->next = NULL;
	temp_inst->last = temp_inst;
	temp_inst->loc = loc_init();
	return temp_inst;
}

vm_inst*
new_vm_inst_push_ival( int ival )
{
	vm_inst* temp_inst = new_vm_inst_command( VM_PUSH_IVAL );
	temp_inst->ival = ival;
	return temp_inst;
}

vm_inst*
new_vm_inst_push_dval( double dval )
{
	vm_inst* temp_inst = new_vm_inst_command( VM_PUSH_DVAL );
	temp_inst->dval = dval;
	return temp_inst;
}

vm_inst*
new_vm_inst_push_pp_ival( char* ptr_key )
{
    vm_inst* temp_inst = new_vm_inst_command( VM_PUSH_PP_NUM );

	int length = strlen(ptr_key);
	char* key_copied = (char*) malloc(length*sizeof(char)+1);
	strcpy(key_copied, ptr_key);
	temp_inst->ptr_key = key_copied;

	return temp_inst;
}

vm_inst*
new_vm_inst_push_pp_dval( char* ptr_key )
{
    vm_inst* temp_inst = new_vm_inst_command( VM_PUSH_PP_NUM );

	int length = strlen(ptr_key);
	char* key_copied = (char*) malloc(length*sizeof(char)+1);
	strcpy(key_copied, ptr_key);
	temp_inst->ptr_key = key_copied;

	return temp_inst;
}

vm_inst*
new_vm_inst_push_pp_str( char* ptr_key )
{
	vm_inst* temp_inst = new_vm_inst_command( VM_PUSH_PP_STR );

	int length = strlen(ptr_key);
	char* key_copied = (char*) malloc(length*sizeof(char)+1);
	strcpy(key_copied, ptr_key);
	temp_inst->ptr_key = key_copied;

	return temp_inst;
}

vm_inst*
new_vm_inst_push_pp_rexp( char* ptr_key )
{
	vm_inst* temp_inst = new_vm_inst_command( VM_PUSH_PP_REXP );

	int length = strlen(ptr_key);
	char* key_copied = (char*) malloc(length*sizeof(char)+1);
	strcpy(key_copied, ptr_key);
	temp_inst->ptr_key = key_copied;

	return temp_inst;
}

vm_inst*
new_vm_inst_push_null( char* ptr_key )
{
	vm_inst* temp_inst = new_vm_inst_command( VM_PUSH_NULL );

	int length = strlen(ptr_key);
	char* key_copied = (char*) malloc(length*sizeof(char)+1);
	strcpy(key_copied, ptr_key);
	temp_inst->ptr_key = key_copied;

	return temp_inst;
}

vm_inst*
new_vm_inst_label( char* label )
{
	vm_inst* temp_inst = new_vm_inst_command( VM_LABEL );

	int length = strlen(label);
	char* label_copied = (char*) malloc(length*sizeof(char)+1);
	strcpy(label_copied, label);
	temp_inst->label = label_copied;

	return temp_inst;
}


vm_inst* 
new_vm_inst_fjmp( char* label )
{
	vm_inst* temp_inst = new_vm_inst_command( VM_FJMP );

	int length = strlen(label);
	char* label_copied = (char*) malloc(length*sizeof(char)+1);
	strcpy(label_copied, label);
	temp_inst->label = label_copied;

	return temp_inst;
}


vm_inst*
new_vm_inst_jmp( char* label )
{
	vm_inst* temp_inst = new_vm_inst_command( VM_JMP );

	int length = strlen(label);
	char* label_copied = (char*) malloc(length*sizeof(char)+1);
	strcpy(label_copied, label);
	temp_inst->label = label_copied;

	return temp_inst;
}


vm_inst*
vm_inst_list_last( vm_inst_list* list )
{
	
	vm_inst* last_inst = list->last ;
	return last_inst;
}


vm_inst*
vm_inst_list_reach_last( vm_inst_list* list )
{
	vm_inst* curr_inst = list;
	while(curr_inst->next != NULL)
	{
		curr_inst = curr_inst->next;
	}
	return curr_inst;
}

vm_inst*
vm_inst_list_get(vm_inst_list* list , int zero_indexed)
{
	int idx = zero_indexed;
	vm_inst* curr_inst = list;
	do{
		if(idx != 0){
			idx = idx - 1;
			curr_inst = curr_inst->next;
		} else {
			return curr_inst;
		}
	}while(curr_inst->next != NULL);

	Rprintf("index is out of bound. The returned inst is the last one.\n");
	return curr_inst;
}

int
vm_inst_list_size(vm_inst_list* list )
{
	int size = 1;
	vm_inst* curr_inst = list;
	while(curr_inst->next != NULL)
	{
		curr_inst = curr_inst->next;
		size = size + 1;
	}
	return size;
}

vm_inst_list*
vm_inst_list_cat( vm_inst_list* list1, vm_inst_list* list2 )
{
	vm_inst* list1_last = vm_inst_list_last(list1);

	list1_last->next = list2;
	list2->prev = list1_last;

	list1->last = vm_inst_list_last(list2);

	return list1;
}

void
vm_inst_list_show_all( vm_inst_list* list )
{
	vm_inst* curr_inst;
	curr_inst = list;
    DEBUG_PRINT("For ptr_table record, VM instructions just holds 'key name'.\n");
    DEBUG_PRINT("For values, VM instructions holds values themselves.\n");

	do{
		vm_inst_show( curr_inst );
		curr_inst = curr_inst->next;
	}while(curr_inst != NULL);
}

/*
vm_inst_list*
vm_inst_list_insert( vm_inst_list* list1, int idx, vm_inst_list* list2)
{
	if(idx <= 0 )
		Rprintf("idx should be larger than zero.\n");
	if(idx > (vm_inst_list_size(list1) - 1) )
		Rprintf("idx is out of bound.\n");

	vm_inst* list1_inst_before_idx = vm_inst_list_get(list1, idx - 1);
	vm_inst* list1_inst_at_idx = vm_inst_list_get(list1, idx );
	vm_inst* list2_head = list2;
	vm_inst* list2_tail = vm_inst_list_tail( list2 );

	list1_inst_before_idx->next = list2_head;
	list2_head->prev = list1_inst_before_idx;

	list1_inst_at_idx->prev = list2_tail;
	list2_tail->next = list1_inst_at_idx;

	return list1;
}
*/

int
vm_inst_list_free( vm_inst_list* inst_list )
{
	vm_inst* curr_inst = inst_list;
	vm_inst* next_inst = curr_inst->next;
	vm_inst_free( curr_inst );
	if(next_inst != NULL){
		vm_inst_list_free(next_inst);
	}
	return 1;
}

vm_inst*
vm_inst_list_to_code( vm_inst_list* list)
{
	int size = vm_inst_list_size(list);
	vm_inst* vm_code_start = (vm_inst*)malloc(sizeof(vm_inst) * size);
	vm_inst* vm_code_ptr = vm_code_start;
	vm_inst* vm_inst_list_ptr = list;
	int idx;
	for( idx = 0 ; idx < size ; idx = idx + 1){
		memcpy( vm_code_ptr, vm_inst_list_ptr, sizeof(vm_inst) );
		vm_inst_list_ptr = vm_inst_list_ptr->next;
		vm_code_ptr = vm_code_ptr + 1;
	}
	return vm_code_start;
}


// Used to set script location to vm instruction
void
vm_inst_set_loc_to_last( struct script_loc loc, vm_inst* insts)
{
	insts->last->loc = loc;
}



// main function to check behavior 
// gcc gen_code_util.c vm/vm_code.o vm/vm_cmd.o -Ivm

/*
#define cat_insts(a, b) vm_inst_list_cat(a, b) 

int
main(int argc, char** argv)
{
	vm_inst* list = new_vm_inst_push_pp_ival("x");
	cat_insts( list, new_vm_inst_push_ival(11));
	cat_insts( list, new_vm_inst_push_ival(22));
	cat_insts( list, new_vm_inst_command(VM_ADDX));
	cat_insts( list, new_vm_inst_command(VM_DISP));
	cat_insts( list, new_vm_inst_command(VM_STO));
	cat_insts( list, new_vm_inst_command(VM_END));

	Rprintf("---Show all---\n");
	vm_inst_list_show_all( list );
	Rprintf("---Free all---\n");
	vm_inst_list_free( list );
}
*/

