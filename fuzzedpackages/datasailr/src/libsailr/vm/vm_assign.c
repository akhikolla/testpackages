#include <R_ext/Print.h>
#include <string.h>
#include <stdio.h>
#include "string/common_string.h"
#include "vm_stack.h"
#include "vm_assign.h"
#include "vm_item_pp2val.h"
#include "ptr_table.h"
#include "simple_re/simple_re.h"
#include "helper.h"
#include "vm_error.h"

int vm_stack_assign_numval_to_ptr_dbl_record(ptr_record* left_record, stack_item* rvalue);
int vm_stack_assign_numval_to_ptr_int_record(ptr_record* left_record, stack_item* rvalue);
int vm_stack_assign_temp_str_to_record(ptr_record* left_record, stack_item* rvalue);
int vm_stack_assign_copy_str_to_record(ptr_record* left_record, stack_item* rvalue);
int vm_stack_assign_temp_rexp_to_record(ptr_record* left_record, stack_item* rvalue);
int vm_stack_assign_copy_rexp_to_record(ptr_record* left_record, stack_item* rvalue);

// Needed?
int vm_stack_convert_str_item_into_void(ptr_record* left_record, stack_item* rvalue);

#define PTR_TABLE_NULL_UPDATED( ptr_to_ptr_record , ptr_type ) do{ ptr_table* table = ptr_record_obtain_table( ptr_to_ptr_record ); ptr_table_info_change_null_updated_by_type(&table, ptr_type ); } while(0)

int
vm_stack_store_val(vm_stack* vmstack)
{
	stack_item* lvalue = vm_stack_second(vmstack);
	stack_item* rvalue = vm_stack_top(vmstack);
	stack_item_pp2value(rvalue);
	
	// lvalue should point to heap memory to store objects or values, b/c they should be returned to library user.
	// If they are not stored on heap, there are possibilities that they are automatically destroyed.
	ptr_record* left_record = NULL;
	if( (lvalue->type != NULL_ITEM) && (lvalue->type != PP_IVAL ) && (lvalue->type != PP_DVAL ) && (lvalue->type != PP_STR ) && (lvalue->type != PP_REXP )){
		Rprintf("ERROR: lvalue should be pointer to pointer, such as PP_IVAL, PP_DVAL, PP_STR or PP_REXP, or NULL_ITEM.\n");
		vm_error_raise(vmstack);
		return 0;
	} else {
		left_record = (ptr_record*)lvalue->p_record;
	}

	if(lvalue->type == NULL_ITEM){ // Unknown variables
		// ------------------------------
		// Unknown and undefined lvalue
		// Library users do not know the types of these variables beforehand, and they are not defined yet.
		// All the memory that is pointed by this address shoud be GC_YES. 
		// ------------------------------
		if(left_record->type == PTR_NULL){  // Undefined variables
			if(rvalue->type == IVAL){ // Undefined varialbes become defined ones.
				DEBUG_PRINT("Thin lvalue is an unknown and undefined variable, which becomes PTR_INT.\n");
				// change type from PTR_NULL to PTR_INT on ptr_table.
				PTR_TABLE_NULL_UPDATED(left_record, PTR_INT);
				left_record->type = PTR_INT;
				left_record->address = malloc(sizeof(int));
				left_record->gc = GC_YES;
				left_record->ex_type = PTR_DBL;
				left_record->ex_addr = malloc(sizeof(double));
				left_record->ex_gc = GC_YES;
				// assign value to the newly allocated memory.
				// memcpy( left_record->address, &(rvalue->ival), sizeof(int));
				*((int*)left_record->address) = rvalue->ival ;
				*((double*)left_record->ex_addr) = 0.0;
			}else if( rvalue->type == DVAL){
				DEBUG_PRINT("Thin lvalue is an unknown and undefined variable, which becomes PTR_DBL.\n");
				// change type from PTR_NULL to PTR_DBL on ptr_table.
				PTR_TABLE_NULL_UPDATED(left_record, PTR_DBL);
				left_record->type = PTR_DBL;
				left_record->address = malloc(sizeof(double));
				left_record->gc = GC_YES;
				left_record->ex_type = PTR_INT;
				left_record->ex_addr = malloc(sizeof(int));
				left_record->ex_gc = GC_YES;
				// assign value to the newly allocated memory.
				// memcpy( left_record->address, &(rvalue->dval), sizeof(double));
				*((double*)left_record->address) = rvalue->dval ;
				*((int*)left_record->ex_addr) = 0;
			}else if( rvalue->type == PP_STR){
				DEBUG_PRINT("Thin lvalue is an unknown and undefined variable, which becomes PTR_STR.\n");
				// change type from PTR_NULL to PTR_STR on ptr_table.
				PTR_TABLE_NULL_UPDATED(left_record, PTR_STR);
				left_record->type = PTR_STR;
				if( vm_stack_item_is_temp(rvalue) ){ // If rvalue is temporary, use the object.
					left_record->address = (string_object*) *(rvalue->pp_str);
					free(rvalue->pp_str);
					rvalue->pp_str = NULL;
					rvalue->type = VOID_ITEM;
				}else{  // If rvalue is not tempoary, create a new string and manage it 
					left_record->address = (string_object*) string_new(string_read((string_object*) *(rvalue->pp_str)));
				}
				left_record->gc = GC_YES;

			}else if( rvalue->type == PP_REXP){
				DEBUG_PRINT("Thin lvalue is an unknown and undefined variable, which becomes PTR_REXP.\n");
				// change type from PTR_NULL to PTR_STR on ptr_table.
				PTR_TABLE_NULL_UPDATED(left_record, PTR_REXP);
				left_record->type = PTR_REXP;
				if( vm_stack_item_is_temp(rvalue) ){ // If rvalue is temporary, use the object.
					Rprintf("LVALUE is unknown and rvalue is temp rexp.\n");
					left_record->address = (simple_re*) *(rvalue->pp_rexp);
					free(rvalue->pp_rexp);
					rvalue->pp_rexp = NULL;
					rvalue->type = VOID_ITEM;
				}else{  // If rvalue is not tempoary, create a new regular expression object and manage it 
					// printf("LVALUE is unknown and rvalue is ptr_table rexp.\n");
					left_record->address = (simple_re*) simple_re_compile( (*(rvalue->pp_rexp))->pattern , (*(rvalue->pp_rexp))->encoding );
				}
				left_record->gc = GC_YES;
			}else{
				Rprintf("ERROR: Only IVAL, DVAL, PP_STR or PP_REXP can be rvalue for assignment operator.\n");
				vm_error_raise(vmstack);
            }

		// ------------------------------
		// Unknown but defined lvalue
		// Library users do not know the types of these variables beforehand, but they are already defined during execution.
		// ------------------------------
		// Unknown but defined variable as PTR_INT
		}else if(left_record->type == PTR_INT){ 
			DEBUG_PRINT("This lvalue is an originally unknown but is now defined variable, PTR_INT.\n");
			if( vm_stack_assign_numval_to_ptr_int_record(left_record, rvalue) != 1){
				vm_error_raise(vmstack);
			}
		// Unknown but defined variable as PTR_DBL
		}else if(left_record->type == PTR_DBL){
			DEBUG_PRINT("This lvalue is an originally unknown but is now defined variable, PTR_DBL.\n");
			if( vm_stack_assign_numval_to_ptr_dbl_record(left_record, rvalue) != 1){
				vm_error_raise(vmstack);
			}
		// Unknown but defined variable as PTR_STR
		}else if(left_record->type == PTR_STR){
			if(rvalue->type == PP_STR){
				if( vm_stack_item_is_temp(rvalue) ){ // If rvalue is temporary, use the object.
					ptr_record_free_gc_required_memory( left_record );
					if( vm_stack_assign_temp_str_to_record(left_record, rvalue) != 1){
						vm_error_raise(vmstack);
					}
				}else{  // If rvalue is not tempoary, create a new string and manage it
					if(rvalue->p_record == lvalue->p_record){
					// Assign to oneself (e.g.) cyl = cyl
					}else{
						ptr_record_free_gc_required_memory( left_record );
						if( vm_stack_assign_copy_str_to_record(left_record, rvalue) != 1 ){
							vm_error_raise(vmstack);
						}
					}
				}
			}else {
				Rprintf("ERROR: Object other than PP_STR is trying to be assigned to PTR_STR.\n");
				vm_error_raise(vmstack);
			}
		// Unknown but defined variable as PTR_REXP
		}else if(left_record->type == PTR_REXP){
			if(rvalue->type == PP_REXP){
				ptr_record_free_gc_required_memory( left_record );
				if( vm_stack_item_is_temp(rvalue) ){ // If rvalue is temporary, use the object.
					if( vm_stack_assign_temp_rexp_to_record(left_record, rvalue) != 1){
						vm_error_raise(vmstack);
					}
				}else{  // If rvalue is not tempoary, create a new regular expression and manage it
					if( vm_stack_assign_copy_rexp_to_record(left_record, rvalue) != 1 ){
						vm_error_raise(vmstack);
					} 
				}
			}else {
				Rprintf("ERROR: Object other than PP_REXP is trying to be assigned to PTR_REXP.\n");
				vm_error_raise(vmstack);
			}
		}

	// ------------------------------
	// Known and defined lvalue
	// Library users know the types of these variables beforehand.
	// ------------------------------
	// as PP_IVAL. 
	}else if(lvalue->type == PP_IVAL){
		if(left_record->type != PTR_INT){
			Rprintf("ERROR: ptr record should be PTR_INT. This branch should never be executed. \n");
			vm_error_raise(vmstack);
		}else{
			if( vm_stack_assign_numval_to_ptr_int_record(left_record, rvalue) != 1){
				vm_error_raise(vmstack);
			}
		}

	// as PP_DVAL
	}else if(lvalue->type == PP_DVAL){
		if(left_record->type != PTR_DBL){
			Rprintf("ERROR: ptr record should be PTR_DBL. This branch should never be executed. \n");
			vm_error_raise(vmstack);
		}else{
			if( vm_stack_assign_numval_to_ptr_dbl_record(left_record, rvalue) != 1){
				vm_error_raise(vmstack);
			}
		}

	// as PP_STR
	}else if(lvalue->type == PP_STR){
		if(left_record->type != PTR_STR){
			Rprintf("ERROR: ptr record should be PTR_STR. This branch should never be executed. \n");
			vm_error_raise(vmstack);
		}else{
			if(rvalue->type == PP_STR){
				if( vm_stack_item_is_temp(rvalue) ){ // If rvalue is temporary, use the object.
					ptr_record_free_gc_required_memory( left_record );
					if( vm_stack_assign_temp_str_to_record(left_record, rvalue) != 1){
						vm_error_raise(vmstack);
					}
				}else{  // If rvalue is not tempoary, create a new string and manage it
					if(rvalue->p_record == lvalue->p_record){
					// Assign to oneself (e.g.) cyl = cyl
					}else{
						ptr_record_free_gc_required_memory( left_record );
						if( vm_stack_assign_copy_str_to_record(left_record, rvalue) != 1){
							vm_error_raise(vmstack);
						}
					}
				}
			}else {
				Rprintf("ERROR: Object other than PP_STR is trying to be assigned to PTR_STR.\n");
				vm_error_raise(vmstack);
			}
		}

	// as PP_REXP
	}else if(lvalue->type == PP_REXP){
		if(left_record->type != PTR_REXP){
			Rprintf("ERROR: ptr record should be PTR_REXP. This branch should never be executed. \n");
			vm_error_raise(vmstack);
		}else{
			if(rvalue->type == PP_REXP){
				ptr_record_free_gc_required_memory( left_record );
				if( vm_stack_item_is_temp(rvalue) ){ // If rvalue is temporary, use the object.
					if( vm_stack_assign_temp_rexp_to_record(left_record, rvalue) != 1){
						vm_error_raise(vmstack);
					}
				}else{  // If rvalue is not tempoary, create a new regular expression and manage it
					if( vm_stack_assign_copy_rexp_to_record(left_record, rvalue) != 1){
						vm_error_raise(vmstack);
					}
				}
			}else {
				Rprintf("ERROR: Object other than PP_REXP is trying to be assigned to PTR_REXP.\n");
				vm_error_raise(vmstack);
			}
		}
	}

	vm_stack_clean_and_pop( vmstack , 2 );
	return 1;
}

int
vm_stack_assign_numval_to_ptr_dbl_record(ptr_record* left_record, stack_item* rvalue)
{
			if(rvalue->type == IVAL){// Type mismatch
				DEBUG_PRINT("lvalue is PTR_DBL and rvalue is IVAL, and assign the value after converting the lvalue into PTR_INT.");
				// Main address type is now PTR_INT
				ptr_record_swap_addresses(left_record); 
				*((int*)left_record->address) = rvalue->ival;
				left_record->gc = left_record->gc; // This heap area is usually prepared by library user.
			}else if( rvalue->type == DVAL){ // Type compatible
				// continue to be left_record->type == PTR_DBL
				*((double*)left_record->address) = rvalue->dval;
				left_record->gc = left_record->gc; // This heap area is usually prepared by library user.
			}else {
				Rprintf("ERROR: Object other than IVAL an DVAL is trying to be assigned to PTR_DBL.\n");
				return 0;
			}
			return 1;
}

int
vm_stack_assign_numval_to_ptr_int_record(ptr_record* left_record, stack_item* rvalue)
{
			if(rvalue->type == IVAL){ // Type compatible
				DEBUG_PRINT("lvalue is PTR_INT and rvalue is IVAL, and just assign the value.");
				// continue to be left_record->type == PTR_INT
				*((int*)left_record->address) = rvalue->ival;
				left_record->gc = left_record->gc; // This heap area is usually prepared by library user.
			}else if( rvalue->type == DVAL){ // Type mismatch
				DEBUG_PRINT("lvalue is PTR_INT and rvalue is DVAL, and assign the value after converting the lvalue into PTR_DBL.");
				// Main address type is now PTR_DBL
				ptr_record_swap_addresses(left_record); 
				*((double*)left_record->address) = rvalue->dval;
				left_record->gc = left_record->gc; // This heap area is usually prepared by library user.
			}else {
				Rprintf("ERROR: Object other than IVAL an DVAL is trying to be assigned to PTR_INT.\n");
				return 0;
			}
			return 1;
}

int
vm_stack_assign_temp_str_to_record(ptr_record* left_record, stack_item* rvalue)
{
	left_record->address = (void*) *(rvalue->pp_str);
	left_record->gc = GC_YES;
	free(rvalue->pp_str);
	rvalue->pp_str = NULL;
	rvalue->type = VOID_ITEM;
	return 1;
}

int
vm_stack_assign_copy_str_to_record(ptr_record* left_record, stack_item* rvalue)
{
	left_record->address = (void*) string_new(string_read((string_object*) *(rvalue->pp_str)));
	left_record->gc = GC_YES;
	return 1;
}

int
vm_stack_assign_temp_rexp_to_record(ptr_record* left_record, stack_item* rvalue)
{
	left_record->address = (void*) *(rvalue->pp_rexp);
	left_record->gc = GC_YES;
	free(rvalue->pp_rexp);
	rvalue->pp_rexp = NULL;
	rvalue->type = VOID_ITEM;
	return 1;
}

int
vm_stack_assign_copy_rexp_to_record(ptr_record* left_record, stack_item* rvalue)
{
	left_record->address = (simple_re*) simple_re_compile( (*(rvalue->pp_rexp))->pattern , (*(rvalue->pp_rexp))->encoding );
	left_record->gc = GC_YES;
	return 1;
}


