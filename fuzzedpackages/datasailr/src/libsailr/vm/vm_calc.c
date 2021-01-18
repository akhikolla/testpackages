#include <R_ext/Print.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>

#include "vm_calc.h"
#include "string/common_string.h"
#include "vm_item_pp2val.h"
#include "ptr_table.h"
#include "helper.h"
#include "vm_error.h"

#define INT_OVERFLOW


// Helper functions --------------


bool
within_int_limits( double num )
{
	if( (INT_MIN < num) && (num < INT_MAX) ){
		return true;
	}else{
		return false;
	}
}

INT_OVERFLOW int
int_add( int x, int y)
{
	return (x + y);
}

double
dbl_add( double x, double y)
{
	return (x + y);
}

INT_OVERFLOW int
int_mul( int x, int y)
{
	return (x * y);
}

double
dbl_mul( double x, double y)
{
	return (x * y) ;
}

INT_OVERFLOW int
int_sub( int x, int y)
{
	return (x - y);
}

double
dbl_sub( double x, double y)
{
	return (x - y) ;
}

INT_OVERFLOW int
int_pow( int x, int y)
{
	double result = pow((double)x, (double)y) + 0.5;
	int int_result = (int) result;
	return int_result;
}

double
dbl_pow( double x, double y)
{
	double result = pow(x, y);
	// TODO Check range of double ? -> I think when double no buffer overflow happens doesn't it.
	return result;
}

int
int_mod( int x, int y)
{
	return ( x % y );
}

double
dbl_mod( double x, double y)
{
	double result = modf( x , &y );
	return  result;
}

// div operator is going to produce double type.
/*
int
divisible_int_div( int x, int y )
{
	int result = x / y;
	return result;
}
*/

double
int_div( int x, int y )
{
	double result;
	result = ((double) x) / ((double) y);
	return result;
}

double
dbl_div( double x, double y )
{
	double result = x / y;
	return result;
}

double
int_factorial( int x )
{
	if( x < 0 )
		return -1;

	double result;
	if( x == 1){
		result = 1.0;
	}else{
		result = x * (int) (int_factorial( x - 1 ));
	}
	return result ;
}


bool
item_is_num (stack_item* item)
{
	if((item->type == IVAL) || (item->type == DVAL) || (item->type == PP_IVAL) || (item->type == PP_DVAL) )
		return true;
	else
		return false;
}

bool
item_is_str (stack_item* item)
{
	if( item->type == PP_STR) 
		return true;
	else
		return false;
}

bool
item_is_bool (stack_item* item)
{
	if( item->type == BOOLEAN) 
		return true;
	else
		return false;
}

bool
item_is_nan (stack_item* item)
{
	if( item->type == DVAL ){
		if( isnan( item->dval )){
			return true;
		}
	}
	return false;
}

// -------------------------------


int
vm_calc_addx(vm_stack* vmstack, ptr_table** table)
{
	stack_item* top_item = vm_stack_top(vmstack);
	stack_item* sec_item = vm_stack_second(vmstack);
	stack_item_pp2value( top_item );
	stack_item_pp2value( sec_item );

	string_object* str1;
	string_object* str2;
	string_object* new_p_str;
	string_object** new_pp_str;

	int result_ival;
	double result_dval;
	double temp_dval;
	// Add calcuation
	// For ivals
	if((top_item->type == IVAL) && (sec_item->type == IVAL)){
		temp_dval = dbl_add ( (double) sec_item->ival , (double) top_item->ival ) ; 
		if( within_int_limits( temp_dval )){ 
			result_ival = int_add ( sec_item->ival , top_item->ival ) ; 
			vm_stack_clean_and_pop( vmstack, 2);
			vm_stack_push_ival( vmstack, result_ival );
		}else{ 
			result_dval = temp_dval ;  
			vm_stack_clean_and_pop( vmstack, 2);
			vm_stack_push_dval( vmstack, result_dval );
		} 
	} else {
		// For dvals
		if((top_item->type == IVAL) && (sec_item->type == DVAL)){
			result_dval = sec_item->dval + (double)top_item->ival ;
		}else if((top_item->type == DVAL) && (sec_item->type == IVAL)){
			result_dval = (double)sec_item->ival + top_item->dval ;
		}else if((top_item->type == DVAL) && (sec_item->type == DVAL)){
			result_dval = sec_item->dval + top_item->dval ;
		// For strings
		}else if((top_item->type == PP_STR ) && (sec_item->type == PP_STR )){
			str1 = (* (sec_item->pp_str)); /* Left on ADDX operator */
			str2 = (* (top_item->pp_str)); /* Right on ADDX operator */

			new_p_str = (string_object*) string_ptr_concat( str1, str2 );
			new_pp_str = (string_object**) malloc(sizeof(string_object*));
			*new_pp_str = new_p_str;
			vm_stack_clean_and_pop( vmstack, 2);
			vm_stack_push_temp_pp_str( vmstack, new_pp_str);
			return 1;
		}else{
			Rprintf("ERROR: ADDX should be applied to 'num and num' or 'str and str' on stack.\n");
			vm_error_raise(vmstack);
			return 0;
		}
		vm_stack_clean_and_pop( vmstack, 2);
		vm_stack_push_dval( vmstack, result_dval );
	}
	return 1;
}

#define DEFINE_BINARY_OPERATOR( INTFUNC , DBLFUNC ) \
do { \
	stack_item* top_item = vm_stack_top(vmstack); \
	stack_item* sec_item = vm_stack_second(vmstack); \
	stack_item_pp2value( top_item ); \
	stack_item_pp2value( sec_item ); \
 \
	int result_ival; \
	double result_dval; \
	double temp_dval; \
	/* OP calcuation*/ \
	/* For ivals*/ \
	if((top_item->type == IVAL) && (sec_item->type == IVAL)){ \
		temp_dval = DBLFUNC ( (double) sec_item->ival , (double) top_item->ival ) ; \
		if( within_int_limits( temp_dval )){ \
			result_ival = INTFUNC ( sec_item->ival , top_item->ival ) ; \
			vm_stack_clean_and_pop( vmstack, 2); \
			vm_stack_push_ival( vmstack, result_ival ); \
		}else{ \
			result_dval = temp_dval ; \
			vm_stack_clean_and_pop( vmstack, 2); \
			vm_stack_push_dval( vmstack, result_dval ); \
		} \
	} else { \
		/* For dvals */ \
		if((top_item->type == IVAL) && (sec_item->type == DVAL)){ \
			result_dval = DBLFUNC ( sec_item->dval , (double)top_item->ival ) ; \
		}else if((top_item->type == DVAL) && (sec_item->type == IVAL)){ \
			result_dval = DBLFUNC ( (double)sec_item->ival , top_item->dval ) ; \
		}else if((top_item->type == DVAL) && (sec_item->type == DVAL)){ \
			result_dval = DBLFUNC ( sec_item->dval , top_item->dval ) ; \
		}else{ \
			Rprintf("ERROR: This VM_CMD should be applied to num and num on stack.\n"); \
			vm_error_raise(vmstack); \
			return 0; \
		} \
		vm_stack_clean_and_pop( vmstack, 2); \
		vm_stack_push_dval( vmstack, result_dval ); \
	} \
	return 1; \
} while(0)

int
vm_calc_mulx(vm_stack* vmstack)
{
  DEFINE_BINARY_OPERATOR( int_mul, dbl_mul );
}

int
vm_calc_subx(vm_stack* vmstack)
{
  DEFINE_BINARY_OPERATOR( int_sub, dbl_sub );
}

int
vm_calc_powx(vm_stack* vmstack)
{
  DEFINE_BINARY_OPERATOR( int_pow, dbl_pow );
}

int
vm_calc_modx(vm_stack* vmstack)
{
  DEFINE_BINARY_OPERATOR( int_mod, dbl_mod );
}


int
vm_calc_divx(vm_stack* vmstack)
{
	stack_item* top_item = vm_stack_top(vmstack);
	stack_item* sec_item = vm_stack_second(vmstack);
	stack_item_pp2value( top_item );
	stack_item_pp2value( sec_item );

//	int result_ival;  // divx always return double
	double result_dval;
	// For ivals
	if((top_item->type == IVAL) && (sec_item->type == IVAL)){

// Previously, division returned int as much as possible. This caused unintentional errors, and now division always returns double.
//		if( sec_item->ival % top_item->ival == 0 ){
//			result_ival = divisible_int_div( sec_item->ival , top_item->ival) ;
//			vmstack->sp = vmstack->sp - 1; 
//			sec_item->type = IVAL;
//			sec_item->ival = result_ival;
//		}else{
//			DEBUG_PRINT("calculation dividing int by int is not divisable. %d / %d \n", sec_item->ival, top_item->ival );

			result_dval = int_div( sec_item->ival, top_item->ival ) ;
			vmstack->sp = vmstack->sp - 1; 
			sec_item->type = DVAL;
			sec_item->dval = result_dval;
//		}
	} else {
		// For dvals
		if((top_item->type == IVAL) && (sec_item->type == DVAL)){
			result_dval = dbl_div( sec_item->dval , (double)top_item->ival ) ;
		}else if((top_item->type == DVAL) && (sec_item->type == IVAL)){
			result_dval = dbl_div( (double)sec_item->ival , top_item->dval ) ;
		}else if((top_item->type == DVAL) && (sec_item->type == DVAL)){
			result_dval = dbl_div( sec_item->dval , top_item->dval ) ;
		}else{
			Rprintf("ERROR: DIVX should be applied to num and num on stack.\n");
			vm_error_raise(vmstack);
			return 0;
		}
		vmstack->sp = vmstack->sp - 1; 
		sec_item->type = DVAL;
		sec_item->dval = result_dval;
	}
	return 1;
}

// factorial (single operator)
int
vm_calc_factorial(vm_stack* vmstack)
{
	stack_item* top_item = vm_stack_top(vmstack);
	stack_item_pp2value( top_item );

	int result_ival;
	double result_dval;
	double temp_dval;

	// Factorial calcuation
	// For ivals
	if(top_item->type == IVAL){
		temp_dval = int_factorial( top_item->ival) ;
		if( within_int_limits(temp_dval) ){
			result_ival = (int) temp_dval;
			vmstack->sp = vmstack->sp - 0; 
			top_item->type = IVAL;
			top_item->ival = result_ival;
		}else{
			result_dval = temp_dval;
			vmstack->sp = vmstack->sp - 0; 
			top_item->type = DVAL;
			top_item->dval = result_dval;
		}	
	} else if(top_item->type == DVAL) {
		// For dval
			result_ival = int_factorial( (int) top_item->dval );
			vmstack->sp = vmstack->sp - 0; 
			top_item->type = IVAL;
			top_item->ival = result_ival;		
	} else {
		Rprintf("ERROR: FACT should be applied to num and num on stack.\n");
		vm_error_raise(vmstack);
		return 0;
	}
	return 1;
}

// unitary minus (single operator)
int vm_calc_uminus(vm_stack* vmstack)
{
	stack_item* top_item = vm_stack_top(vmstack);
	stack_item_pp2value( top_item );

	int result_ival;
	double result_dval;

	// For ivals
	if(top_item->type == IVAL){
			result_ival = (-1) * ( top_item->ival) ;
			vmstack->sp = vmstack->sp - 0; 
			top_item->type = IVAL;
			top_item->ival = result_ival;		
	} else if(top_item->type == DVAL) {
		// For dval
			result_dval = (-1) * ( top_item->dval );
			vmstack->sp = vmstack->sp - 0; 
			top_item->type = DVAL;
			top_item->dval = result_dval;		
	} else {
		Rprintf("ERROR: uminus should be applied to num and num on stack.\n");
		vm_error_raise(vmstack);
		return 0;
	}
	return 1;
}

// Logical Calculations

int
vm_calc_and(vm_stack* vmstack)
{
	stack_item* top_item = vm_stack_top(vmstack);
	stack_item* sec_item = vm_stack_second(vmstack);
	bool result_bool;

	if( (top_item->type == BOOLEAN) && (sec_item->type == BOOLEAN)){
		result_bool = ( ( sec_item->boolean ) && ( top_item->boolean )) ;
		vmstack->sp = vmstack->sp - 1; 
		sec_item->type = BOOLEAN ;
		sec_item->boolean = result_bool;
	}else{
		Rprintf("ERROR: AND should be applied to boolean and boolean.\n");
		vm_error_raise(vmstack);
		return 0;
	}
	return 1;
}

int
vm_calc_or(vm_stack* vmstack)
{
	stack_item* top_item = vm_stack_top(vmstack);
	stack_item* sec_item = vm_stack_second(vmstack);
	bool result_bool;

	if( (top_item->type == BOOLEAN) && (sec_item->type == BOOLEAN)){
		result_bool = ( ( sec_item->boolean ) || ( top_item->boolean )) ;
		vmstack->sp = vmstack->sp - 1; 
		sec_item->type = BOOLEAN ;
		sec_item->boolean = result_bool;
	}else{
		Rprintf("ERROR: AND should be applied to boolean and boolean.\n");
		vm_error_raise(vmstack);
		return 0;
	}
	return 1;
}

int
vm_calc_eq(vm_stack* vmstack)
{
	stack_item* top_item = vm_stack_top(vmstack);
	stack_item* sec_item = vm_stack_second(vmstack);
	stack_item_pp2value( top_item ); 
	stack_item_pp2value( sec_item ); 
	bool result_bool;

	if( item_is_nan(sec_item) && item_is_nan(top_item)){
		result_bool = true;
		goto eq_compared;
	}else if(item_is_nan(sec_item) || item_is_nan(top_item)){
		result_bool = false;
		goto eq_compared;
	}else{
		// None of the items are nan.
	}

	if( item_is_num(sec_item) && item_is_num(top_item)){
		stack_item_pp2value( top_item );
		stack_item_pp2value( sec_item );
		if((top_item->type == IVAL) && (sec_item->type == IVAL)){
			result_bool = ( sec_item->ival == top_item->ival ) ;
		}else if((top_item->type == IVAL) && (sec_item->type == DVAL)){
			result_bool = ( sec_item->dval == top_item->ival ) ;
		}else if((top_item->type == DVAL) && (sec_item->type == IVAL)){
			result_bool = ( sec_item->ival == top_item->dval ) ;
		}else if((top_item->type == DVAL) && (sec_item->type == DVAL)){
			result_bool = ( sec_item->dval == top_item->dval ) ;
		}else{
			result_bool = false;
			Rprintf("ERROR: This branch should not be executed.\n");
			vm_error_raise(vmstack);
			return 0;
		}
	}else if(item_is_str(sec_item) && item_is_str(top_item)){
		result_bool = string_compare( *(sec_item->pp_str), *(top_item->pp_str) ) ;
	}else{
		result_bool = false;
		Rprintf("ERROR: Types are invalied for VM_EQ command.\n");
		vm_error_raise(vmstack);
		return 0;
	}

	eq_compared:

	vmstack->sp = vmstack->sp - 1; 
	sec_item->type = BOOLEAN ;
	sec_item->boolean = result_bool;
	return 1;
}


int
vm_calc_neq(vm_stack* vmstack)
{
	stack_item* top_item = vm_stack_top(vmstack);
	stack_item* sec_item = vm_stack_second(vmstack);
	stack_item_pp2value( top_item ); 
	stack_item_pp2value( sec_item ); 
	bool result_bool;

	if( item_is_nan(sec_item) && item_is_nan(top_item)){
		result_bool = false;
		goto neq_compared;
	}else if(item_is_nan(sec_item) && item_is_nan(top_item)){
		result_bool = true;
		goto neq_compared;
	}else{
		// None of the items are nan.
	}

	if( item_is_num(sec_item) && item_is_num(top_item)){
		stack_item_pp2value( top_item );
		stack_item_pp2value( sec_item );
		if((top_item->type == IVAL) && (sec_item->type == IVAL)){
			result_bool = ( sec_item->ival != top_item->ival ) ;
		}else if((top_item->type == IVAL) && (sec_item->type == DVAL)){
			result_bool = ( sec_item->dval != top_item->ival ) ;
		}else if((top_item->type == DVAL) && (sec_item->type == IVAL)){
			result_bool = ( sec_item->ival != top_item->dval ) ;
		}else if((top_item->type == DVAL) && (sec_item->type == DVAL)){
			result_bool = ( sec_item->dval != top_item->dval ) ;
		}else{
			result_bool = false;
			Rprintf("ERROR: This branch should not be executed.\n");
			vm_error_raise(vmstack);
			return 0;
		}
	}else if(item_is_str(sec_item) && item_is_str(top_item)){
		result_bool = !(string_compare( *(sec_item->pp_str), *(top_item->pp_str) )) ;
	}else{
		result_bool = true;  // Currently comparing incompatible types returns true, though this should be an error case. 
		Rprintf("ERROR: Types are invalied for VM_EQ command.\n");
		vm_error_raise(vmstack);
		return 0;
	}

	neq_compared:

	vmstack->sp = vmstack->sp - 1; 
	sec_item->type = BOOLEAN ;
	sec_item->boolean = result_bool;
	return 1;
}

#define DEFINE_LOGICAL_OPERATOR( OP ) \
do { \
	stack_item* top_item = vm_stack_top(vmstack); \
	stack_item* sec_item = vm_stack_second(vmstack); \
	stack_item_pp2value( top_item ); \
	stack_item_pp2value( sec_item ); \
	bool result_bool; \
 \
	if( item_is_num(sec_item) && item_is_num(top_item)){ \
		if((top_item->type == IVAL) && (sec_item->type == IVAL)){ \
			result_bool = ( sec_item->ival OP top_item->ival ) ; \
		}else if((top_item->type == IVAL) && (sec_item->type == DVAL)){ \
			result_bool = ( sec_item->dval OP top_item->ival ) ; \
		}else if((top_item->type == DVAL) && (sec_item->type == IVAL)){ \
			result_bool = ( sec_item->ival OP top_item->dval ) ; \
		}else if((top_item->type == DVAL) && (sec_item->type == DVAL)){ \
			result_bool = ( sec_item->dval OP top_item->dval ) ; \
		}else{ \
			result_bool = false; \
			Rprintf("ERROR: This branch should not be executed.\n"); \
			vm_error_raise(vmstack); \
			return 0; \
		} \
	}else if(item_is_str(sec_item) && item_is_str(top_item)){ \
		Rprintf("ERROR: String is not supported for OP calculation.\n"); \
		vm_error_raise(vmstack); \
		return 0; \
	}else{ \
		result_bool = false; \
		Rprintf("ERROR: Types are invalied for OP calculation.\n"); \
		vm_error_raise(vmstack); \
		return 0; \
	} \
	vmstack->sp = vmstack->sp - 1;  \
	sec_item->type = BOOLEAN ; \
	sec_item->boolean = result_bool; \
	return 1; \
} while(0)

int
vm_calc_gt(vm_stack* vmstack)
{
	DEFINE_LOGICAL_OPERATOR( > ) ;
}

int
vm_calc_lt(vm_stack* vmstack)
{
	DEFINE_LOGICAL_OPERATOR( < ) ;
}

int
vm_calc_ge(vm_stack* vmstack)
{
	DEFINE_LOGICAL_OPERATOR( >= ) ;
}

int
vm_calc_le(vm_stack* vmstack)
{
	DEFINE_LOGICAL_OPERATOR( <= ) ;
}


int
vm_calc_neg(vm_stack* vmstack)
{
	stack_item* top_item = vm_stack_top(vmstack);
	bool result_bool;

	if( item_is_bool(top_item)){
		result_bool = ( ! top_item->boolean );
		vmstack->sp = vmstack->sp - 0; 
		top_item->type = BOOLEAN ;
		top_item->boolean = result_bool;
	}else{
		Rprintf("ERROR: Type is invalied for VM_NEG command.\n");
		vm_error_raise(vmstack);
		return 0;
	}
	return 1;
}


