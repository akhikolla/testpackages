#ifndef C_FUNCS_H
#define C_FUNCS_H

#include "vm/vm_stack.h"

int sailr_func_print(vm_stack* vmstack, int num_args);
int sailr_func_num_to_str(vm_stack* vmstack, int num_args, ptr_table** table );
int sailr_func_str_strip( vm_stack* vmstack, int num_args, ptr_table** table ); // (remove white space from the head/tail)
int sailr_func_str_lstrip( vm_stack* vmstack, int num_args, ptr_table** table ); 
int sailr_func_str_rstrip( vm_stack* vmstack, int num_args, ptr_table** table );
int sailr_func_str_concat( vm_stack* vmstack, int num_args, ptr_table** table );  // str1, str2 ...

int sailr_func_str_repeat( vm_stack* vmstack, int num_args, ptr_table** table );
int sailr_func_str_subset( vm_stack* vmstack, int num_args, ptr_table** table );  // index starts from zero.

int sailr_func_str_to_num( vm_stack* vmstack, int num_args );

int sailr_func_rexp_matched( vm_stack* vmstack, int num_args , ptr_table** table );

int sailr_func_date_ymd(vm_stack* vmstack, int num_args ); 
int sailr_func_date_ym_weekday_nth(vm_stack* vmstack, int num_args ); 
int sailr_func_date_add_n_years(vm_stack* vmstack, int num_args ); 
int sailr_func_date_add_n_months(vm_stack* vmstack, int num_args ); 
int sailr_func_date_add_n_days(vm_stack* vmstack, int num_args );  
int sailr_func_date_format( vm_stack* vmstack, int num_args , ptr_table** table); 

#endif /* C_FUNCS_H */
