#include <R_ext/Print.h>
#include "vm_call_func.h"
#include "func/c_func/c_func.h"
#include "vm_stack.h"
#include <string.h>
#include <stdio.h>

#define FUNC_NAME_IS( a , b ) ( strcmp( ( a ) , ( b ) ) ==  0 )

int
call_func( vm_stack* vmstack, char* fname, int num_args, ptr_table** table )
{
	int result = 1;
    if(FUNC_NAME_IS(fname, "print")){
        result = sailr_func_print(vmstack, num_args);
    }else if(FUNC_NAME_IS(fname, "num_to_str")){
        result = sailr_func_num_to_str(vmstack, num_args, table);
    }else if(FUNC_NAME_IS(fname, "str_strip")){
        result = sailr_func_str_strip(vmstack, num_args, table);
    }else if(FUNC_NAME_IS(fname, "str_lstrip")){
        result = sailr_func_str_lstrip(vmstack, num_args, table);
    }else if(FUNC_NAME_IS(fname, "str_rstrip")){
        result = sailr_func_str_rstrip(vmstack, num_args, table);
    }else if(FUNC_NAME_IS(fname, "str_concat")){
        result = sailr_func_str_concat(vmstack, num_args, table);
    }else if(FUNC_NAME_IS(fname, "str_repeat")){
        result = sailr_func_str_repeat(vmstack, num_args, table);
    }else if(FUNC_NAME_IS(fname, "str_subset")){
        result = sailr_func_str_subset(vmstack, num_args, table);
    }else if(FUNC_NAME_IS(fname, "str_to_num")){
        result = sailr_func_str_to_num(vmstack, num_args);
    }else if(FUNC_NAME_IS(fname, "rexp_matched")){
        result = sailr_func_rexp_matched(vmstack, num_args, table);
    }else if(FUNC_NAME_IS(fname, "date_ymd")){
        result = sailr_func_date_ymd(vmstack, num_args);
    }else if(FUNC_NAME_IS(fname, "date_ym_weekday_nth")){
        result = sailr_func_date_ym_weekday_nth(vmstack, num_args);
    }else if(FUNC_NAME_IS(fname, "date_add_n_years")){
        result = sailr_func_date_add_n_years(vmstack, num_args);
    }else if(FUNC_NAME_IS(fname, "date_add_n_months")){
        result = sailr_func_date_add_n_months(vmstack, num_args);
    }else if(FUNC_NAME_IS(fname, "date_add_n_days")){
        result = sailr_func_date_add_n_days(vmstack, num_args);
    }else if(FUNC_NAME_IS(fname, "date_format")){
        result = sailr_func_date_format(vmstack, num_args, table);
	}else{
        Rprintf("ERROR: Function, %s , cannot be found. \n", fname );
		result = 0;
    }
//	vm_stack_display_all(vmstack);
	return result;
}
