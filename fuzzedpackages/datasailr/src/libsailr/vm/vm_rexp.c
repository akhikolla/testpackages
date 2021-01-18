#include <R_ext/Print.h>
#include "vm_rexp.h"
#include "simple_re/simple_re.h"
#include "string/common_string.h"
#include "vm_item_pp2val.h"
#include "vm_error.h"
#include <stdio.h>
#include <stdbool.h>

int
vm_rexp_match(vm_stack* vmstack)
{
	stack_item* top_item = vm_stack_top(vmstack);
	stack_item* sec_item = vm_stack_second(vmstack);
	stack_item_pp2value( top_item );
	stack_item_pp2value( sec_item );

	stack_item* str_item = NULL;
	stack_item* rexp_item = NULL;

	string_object* str_obj;
	simple_re* rexp_obj;
	int matched_pos;
	bool result_bool;

	if( (top_item->type == PP_STR ) && (sec_item->type == PP_REXP)){
		str_item = top_item;
		rexp_item = sec_item;
	} else if ( (top_item->type == PP_REXP ) && (sec_item->type == PP_STR)){
		str_item = sec_item;
		rexp_item = top_item;
	} else {
		Rprintf("ERROR: REXP_MATCH should have REXP and STR on each side respectively.\n");
		vm_error_raise(vmstack);
		return 0;
	}

	str_obj = *(str_item->pp_str);
	rexp_obj = *(rexp_item->pp_rexp);
	simple_re** ptr_last_rexp_field = vm_stack_get_ptr_last_rexp_field(vmstack); // The last regular expression executed is tracked as a vm_stack information.

#ifdef DEBUG
	simple_re*  ptr_last_rexp_obj = *ptr_last_rexp_field;
	Rprintf("Just before simple_re_match on vm_stack: Address to last_rexp field of vm_stack info %p\n", &(vmstack->stack[0].p_vm_stack_info->last_rexp));
	Rprintf("Just before simple_re_match on vm_stack: Address to last_rexp field obtained through vm_stack_get_ptr_last_rexp_field(vmstack) %p\n",  ptr_last_rexp_field);
	Rprintf("Just before simple_re_match on vm_stack: Address to last_rexp object %p\n", vmstack->stack[0].p_vm_stack_info->last_rexp);
	Rprintf("Just before simple_re_match on vm_stack: Address to last_rexp object obtained through vm_stack_get_ptr_last_rexp_field(vmstack) %p\n",  ptr_last_rexp_obj);
#endif

	matched_pos = simple_re_match( rexp_obj , string_read(str_obj), ptr_last_rexp_field);

#ifdef DEBUG
	ptr_last_rexp_field = vm_stack_get_ptr_last_rexp_field(vmstack); 
	ptr_last_rexp_obj = *ptr_last_rexp_field;
	Rprintf("Just before simple_re_match on vm_stack: Address to last_rexp field of vm_stack info %p\n", &(vmstack->stack[0].p_vm_stack_info->last_rexp));
	Rprintf("Just before simple_re_match on vm_stack: Address to last_rexp field obtained through vm_stack_get_ptr_last_rexp_field(vmstack) %p\n",  ptr_last_rexp_field);
	Rprintf("Just before simple_re_match on vm_stack: Address to last_rexp object %p\n", vmstack->stack[0].p_vm_stack_info->last_rexp);
	Rprintf("Just before simple_re_match on vm_stack: Address to last_rexp object obtained through vm_stack_get_ptr_last_rexp_field(vmstack) %p\n",  ptr_last_rexp_obj);
#endif

	if(matched_pos > 0 ){
		result_bool = true;
	}else{
		result_bool = false;
	}

	vm_stack_clean_and_pop(vmstack, 2); 
	vm_stack_push_boolean(vmstack, result_bool);

	return 1;
}

	




