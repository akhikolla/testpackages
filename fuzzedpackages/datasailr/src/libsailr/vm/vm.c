#include <R_ext/Print.h>
#include <stdio.h>
#include <stdbool.h>
#include "ptr_table.h"
#include "helper.h"

// Using information of vm_code
// VM does
// 1. Stack manipulation
// 2. Calculation

#include "vm.h"
#include "vm_cmd.h"
#include "vm_code.h"
#include "vm_stack.h"
#include "vm_calc.h"
#include "vm_assign.h"
#include "vm_rexp.h"
#include "vm_error.h"

#include "script_loc.h"

int vm_run_inst (vm_inst* , ptr_table* , vm_stack* , ext_func_hash* );

int
vm_exec_code( vm_inst* code , int num_insts , int start_inst_idx, ptr_table* table , vm_stack* vmstack, ext_func_hash* extfunc_hash )
{
	int exec_result = 1; // 0:fail, 1:success, 2:suspend
	int inst_idx = start_inst_idx;
	int move_forward = 0;	
	stack_item* top_item;
	vm_inst* inst;

	for(  ; inst_idx < num_insts ; ++inst_idx ){
		DEBUG_PRINT("\n===== Current Stack (sp=%d) =====\n", vmstack->sp );
		vm_stack_display_all(vmstack);

		inst = &(code[inst_idx]);
		if(inst->cmd == VM_JMP){
			move_forward = vm_code_jmp(code, inst_idx, inst->label , num_insts);
			DEBUG_PRINT("----- VM INSTRUCTION ----- \n");
			DEBUG_PRINT("VM_JMP: move forward by %d .\n", move_forward);
			inst_idx = inst_idx + move_forward ;
		}else if(inst->cmd == VM_FJMP){
			top_item = &( vmstack->stack[vmstack->sp] );
			if( (top_item->type == BOOLEAN) && (top_item->boolean == false) ){
				move_forward = vm_code_jmp(code, inst_idx, inst->label , num_insts);
				DEBUG_PRINT("----- VM INSTRUCTION ----- \n");
				DEBUG_PRINT("VM_JMP: move forward by %d .\n", move_forward);
				inst_idx = inst_idx + move_forward ;

				//The top item on stack should be removed.
				vm_stack_clean_and_pop( vmstack , 1 ) ;
			} else if( (top_item->type == BOOLEAN) && (top_item->boolean == true) ){
				DEBUG_PRINT("Top item of the current stack is boolean. and it says 'true' \n");
				//The top item on stack should be removed.
				vm_stack_clean_and_pop( vmstack , 1 ) ;

			} else {
				Rprintf("ERROR: Top item of the current stack is not boolean... \n");
				exec_result = 0; // fail
				break;
			}
		}else {
			DEBUG_PRINT("----- VM INSTRUCTION ----- \n");
			DEBUG_PRINT("%s(%d/%d)\n", vm_cmd_to_string(inst->cmd), inst_idx, num_insts);
			exec_result = vm_run_inst(inst, table , vmstack, extfunc_hash);
			if(exec_result == 0){ // fail
				Rprintf("ERROR: current vm instruction causing some problem.\n");
				break;
			}else if(exec_result == 2){ // suspend
				vm_stack_set_code_position(vmstack, inst_idx + 1 ); // Set next code postion, which is important to resume
				break;
			}
		}

		// During execution of virtual machine, two mechanisms report errors.
		// One is using function return value (0: success, 1: fail), and the other is using vm_error_*** functions defined in vm_error.c. 
		// Which one is used depends on each function, meaning both are used.
		// The following code checks the errors raised by the latter mechanism.
		if( (inst->cmd != VM_END) && vm_error_exist( vmstack ) ){
			Rprintf("ERROR REPORT: %d runtime error(s) raised", vm_error_num( vmstack ));
			exec_result = 0; // fail
			break;
		}
	}

	if(exec_result == 0){ // fail
		loc_show( inst->loc );
	}

	return exec_result;
}

int
vm_run_inst (vm_inst* inst, ptr_table* table, vm_stack* vmstack , ext_func_hash* extfunc_hash )
{
	int result = 1;
	ext_func_elem* ext_func = NULL;
	switch (inst->cmd){
	case VM_PUSH_IVAL:
		result = vm_stack_push_ival(vmstack, inst->ival);
		break;
	case VM_PUSH_DVAL:
		result = vm_stack_push_dval(vmstack, inst->dval);
		break;
	case VM_PUSH_PP_IVAL:
//		result = vm_stack_push_pp_ival(vmstack, &table, inst->ptr_key);
        Rprintf("ERROR: This instruction is not used. Use VM_PUSH_PP_NUM.");
		result = 0;
		break;
	case VM_PUSH_PP_DVAL:
//		result = vm_stack_push_pp_dval(vmstack, &table, inst->ptr_key);
        Rprintf("ERROR: This instruction is not used. Use VM_PUSH_PP_NUM.");
		result = 0;
		break;
	case VM_PUSH_PP_NUM:
		result = vm_stack_push_pp_num(vmstack, &table, inst->ptr_key);
		break;
	case VM_PUSH_PP_STR:
		result = vm_stack_push_pp_str(vmstack, &table, inst->ptr_key);
		break;
	case VM_PUSH_PP_REXP:
		result = vm_stack_push_pp_rexp(vmstack, &table, inst->ptr_key);
		break;
    case VM_PUSH_NULL:
    	result = vm_stack_push_corresp_item(vmstack, &table, inst->ptr_key);
	    break;
	case VM_POP:
		if( vm_stack_pop(vmstack) == NULL){
			result = 0;
		}
		break;
	case VM_FJMP:
	case VM_JMP:
		Rprintf("ERROR: This code should never be run. ");
		result = 0;
		break;
	case VM_END:
		result = vm_stack_end(vmstack);
		break;
	case VM_DISP:
		result = vm_stack_display_item(vmstack, vmstack->sp);
		break;
	case VM_STO:
		result = vm_stack_store_val(vmstack);
		break;
	case VM_FCALL:
		if( extfunc_hash != NULL ){
			if((ext_func = ext_func_hash_find(&extfunc_hash, inst->fname))){
				DEBUG_PRINT("External function is to be executed: %s\n", inst->fname);
				DEBUG_PRINT("External function pointer: %p \n", ext_func);
				result = ext_func_elem_apply(&extfunc_hash, ext_func, vmstack);
			}else{
				DEBUG_PRINT("Function is not found in external function list: %s\n", inst->fname);
				DEBUG_PRINT("Function is to be executed as an internal function \n");
				result = vm_stack_fcall(vmstack, inst->fname, inst->num_arg, &table );
			}
		}else{
			DEBUG_PRINT("Function is to be executed as an internal function: %s \n", inst->fname);
			result = vm_stack_fcall(vmstack, inst->fname, inst->num_arg, &table );
		}
		break;
	case VM_ADDX:
		result = vm_calc_addx(vmstack, &table);
		break;
	case VM_MULX:
		result = vm_calc_mulx(vmstack);
		break;
	case VM_SUBX:
		result = vm_calc_subx(vmstack);
		break;
	case VM_DIVX:
		result = vm_calc_divx(vmstack);
		break;
	case VM_POWX:
		result = vm_calc_powx(vmstack);
		break;
	case VM_FAC:
		result = vm_calc_factorial(vmstack);
		break;
	case VM_UMINUS:
		result = vm_calc_uminus(vmstack);
		break;
	case VM_AND:
		result = vm_calc_and(vmstack);
		break;
	case VM_OR:
		result = vm_calc_or(vmstack);
		break;
	case VM_EQ:
		result = vm_calc_eq(vmstack);
		break;
	case VM_NEQ:
		result = vm_calc_neq(vmstack);
		break;
	case VM_LT:
		result = vm_calc_lt(vmstack);
		break;
	case VM_LE:
		result = vm_calc_le(vmstack);
		break;
	case VM_GT:
		result = vm_calc_gt(vmstack);
		break;
	case VM_GE:
		result = vm_calc_ge(vmstack);
		break;
	case VM_NEG:
		result = vm_calc_neg(vmstack);
		break;
	case VM_REXP_MATCH:
		result = vm_rexp_match(vmstack);
		break;
	case VM_LABEL:
		// Do nothing.
		break;
	default:
		Rprintf("ERROR: undefined VM command specified. \n");
		result = 0;
		break;
	}

//	ptr_table_show_all(&table);
	return result;
}

// Sample Code
/*
vm_code vm_code1 = {
    { .cmd = VM_PUSH_PP_IVAL, .ptr_key = "x" },
    { .cmd = VM_PUSH_IVAL, .ival = 11 },
    { .cmd = VM_PUSH_IVAL, .ival = 22 },
    { .cmd = VM_ADDX },
	{ .cmd = VM_DISP },
    { .cmd = VM_STO },
    { .cmd = VM_END }
};


int main(int args, char** argv){
	vm_stack* vmstack = vm_stack_init();
	ptr_table* table = ptr_table_init();
	int* address_for_x = (int *)malloc(sizeof(int));
	*address_for_x = 0;
	ptr_table_add( &table, "x", (void**) &address_for_x, PTR_INT, GC_NO);
	Rprintf("---- Show ptr_table ----\n");
	ptr_table_show_all(&table);
	Rprintf("\n");

	vm_exec_code(vm_code1, sizeof(vm_code1)/sizeof(vm_code1[0]), table , vmstack);
	Rprintf("\n");
	Rprintf("---- Show ptr_table ----\n");
	ptr_table_show_all(&table);
	Rprintf("\n");

	free(address_for_x);
	return 0;
}
*/

