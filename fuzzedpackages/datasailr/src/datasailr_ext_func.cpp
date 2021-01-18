#include <iostream>
#include <Rcpp.h>

#include "datasailr_ext_func.hpp"

int
sailr_external_push_row( arg_list* arglist , unsigned int num_args, vm_stack* vmstack )
{
	return 2; // suspend
}

int
sailr_external_discard_row( arg_list* arglist , unsigned int num_args, vm_stack* vmstack )
{
	return 2; // suspend is used to exit vm loop.
}

int
sailr_external_println( arg_list* arglist , unsigned int num_args, vm_stack* vmstack)
{
	arg_item* arg = arglist;
	if(arg_item_confirm_string(arg)){
		string_object* obj = arg_item_string_obj( arg );
		Rcpp::Rcout << string_read(obj) << std::endl;
	}
	arg_list_finalize( vmstack, num_args , arglist);
	return 1;
}

int
sailr_external_add2( arg_list* arglist , unsigned int num_args, vm_stack* vmstack)
{
	arg_item* arg = arglist;

	double x1, x2, total ;
	if(arg_item_confirm_int(arg)){
		x1 = (double) arg_item_int_value( arg );
	}else if(arg_item_confirm_double(arg)){
		x1 = arg_item_double_value( arg );
	}

	arg_item_next(&arg);
	if(arg_item_confirm_int(arg)){
		x2 = (double) arg_item_int_value( arg );
	}else if(arg_item_confirm_double(arg)){
		x2 = arg_item_double_value( arg );
	}

	total = x1 + x2;
	arg_list_finalize( vmstack, num_args , arglist);
	vm_stack_push_dval( vmstack , total );
	return 1;
}
