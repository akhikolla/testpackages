#include <R_ext/Print.h>

#include "c_func.h"
#include "c_func_helper.h"
#include "vm/vm_stack.h"
#include "ptr_table.h"
#include "vm/vm_error.h"

#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>


// Tips to implemnt new functions
// 
// Use the following functions for stack manipulation
// int vm_stack_push_ival( vm_stack* , int );
// int vm_stack_push_dval( vm_stack* , double );
// int vm_stack_push_pp_str( vm_stack* , ptr_table**, char* );
// Before pushing string, register it as an annonymous string into ptr_table.
// ptr_record* ptr_table_create_anonym_string(ptr_table** table, string_object** strptr);
// 
// When error occurs, call vm_error_raise( vmstack ). 
// Be sure that macro functions hard code the variable name, vmstack, so function parameters should include `vm_stack* vmstack` parameter.
//

#define ASSIGN_INT_VALUE( VAR , ARG_PTR , ERROR_MSG ) \
do { \
	if(arg_item_confirm_int( ARG_PTR )){ \
		VAR = arg_item_int_value( ARG_PTR ); \
	}else{ \
		VAR = 0; /* This branch should not be executed */ \
		Rprintf( ERROR_MSG ); \
		vm_error_raise( vmstack ); \
	} \
} while (0)

#define ASSIGN_DOUBLE_VALUE( VAR , ARG_PTR , ERROR_MSG ) \
do { \
	if(arg_item_confirm_double( ARG_PTR )){ \
		VAR = arg_item_double_value( ARG_PTR ); \
	}else{ \
		VAR = 0.0; /* This branch should not be executed */ \
		Rprintf( ERROR_MSG ); \
		vm_error_raise( vmstack ); \
	} \
} while (0)

#define ASSIGN_STRING_PTR( VAR , ARG_PTR , ERROR_MSG ) \
do { \
	if(arg_item_confirm_string( ARG_PTR )){ \
		VAR = arg_item_string_obj( ARG_PTR ); \
	}else{ \
		VAR = NULL; /* This branch should not be executed */ \
		Rprintf( ERROR_MSG ); \
		vm_error_raise( vmstack ); \
	} \
} while (0)


// This fucntion is used by print() and str_concat(), to concatenate multiple objects as string.
void
append_arg_list_as_string( string_object* new_str, arg_list* arglist, vm_stack* vmstack )
{
  arg_item* argitem = arglist; 

  string_object* tmp_str;
  int tmp_int;
  double tmp_double;
  simple_re* tmp_re;
  bool tmp_bool;

  for( argitem = arglist ; argitem != NULL; arg_item_next(&argitem) ){
//	printf("ARG ITEM TYPE: %c \n", arg_item_interpret_type(argitem));

	switch(arg_item_interpret_type(argitem)){ 
    case 's': // string
      ASSIGN_STRING_PTR( tmp_str , argitem , "ERROR: This should be string.\n" );
//	  printf("%s\n", string_read(tmp_str));
      string_append_string(new_str, tmp_str);
    break;
    case 'i': // int
      ASSIGN_INT_VALUE( tmp_int , argitem , "ERROR: This should be int.\n" );
//	  printf("%d\n", tmp_int);
	  tmp_str = string_new_int2str(tmp_int);
      string_append_string(new_str, tmp_str);
	  string_free(tmp_str);
    break;
    case 'd': // double
      ASSIGN_DOUBLE_VALUE( tmp_double , argitem , "ERROR: This should be int.\n" );
//	  printf("%f\n", tmp_double);
	  tmp_str = string_new_double2str(tmp_double);
      string_append_string(new_str, tmp_str);
	  string_free(tmp_str);
    break;
    case 'r': // regexp
      tmp_re = arg_item_rexp_obj( argitem );
//	  printf("%s\n", simple_re_read_pattern(tmp_re));
      string_append_cstring(new_str, simple_re_read_pattern(tmp_re));      
    break;
    case 'b': // boolan
      tmp_bool = arg_item_bool_value( argitem );
      tmp_int = tmp_bool ? 1 : 0 ;
//	  printf("%d\n", tmp_int);
	  tmp_str = string_new_int2str(tmp_int);
      string_append_string(new_str, tmp_str);
	  string_free(tmp_str);
    break;
    case 'n': // null
      // Nothing to be added.
//	  printf("null\n");
    break;
    default:  // shoud never enter this branch
      Rprintf("ERROR: This should never be executed (sailr_func_print()).\n");      
    break;
    }
  }
}

int
sailr_func_print(vm_stack* vmstack, int num_args)
{
  if( ! arg_num_should_be_larger_than( num_args, 0 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  string_object* new_str = string_new("");

  append_arg_list_as_string( new_str, arglist, vmstack);

  Rprintf("%s", string_read(new_str));

  string_free(new_str);

  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  return 1;
}

int
sailr_func_num_to_str(vm_stack* vmstack, int num_args , ptr_table** table)
{
  if( ! arg_num_should_be( num_args, 1 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 

  string_object* p_str ;
  string_object** pp_str = (string_object**) malloc(sizeof(string_object*));

  if(arg_item_confirm_int(argitem)){
    p_str = (string_object*) string_new_int2str( arg_item_int_value(argitem) );
    *pp_str = p_str;  // Never do "pp_str = &p_str;"
  }else if(arg_item_confirm_double(argitem)){
    p_str = string_new_double2str( arg_item_double_value(argitem) );
    *pp_str = p_str;  // Never do "pp_str = &p_str;"
  }else{
    Rprintf("ERROR: For argument, number shouble be specified.\n");
  }

  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_temp_pp_str( vmstack, pp_str);

//  The following is a code previously used. This results in creating new anonymous object for each row.
//  ptr_record* pr = ptr_table_create_anonym_string(table, pp_str);
//  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
//  vm_stack_push_pp_str( vmstack , table , pr->key );

  return 1;
}

int
sailr_func_str_func_ptr (vm_stack* vmstack, int num_args , ptr_table** table , string_object* (*str_func_ptr)(string_object*) )
{
  if( ! arg_num_should_be( num_args, 1 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 

  string_object* p_str ; 
  string_object** pp_str = (string_object**) malloc(sizeof(string_object*));
  string_object* tmp_str ;

  ASSIGN_STRING_PTR( tmp_str, argitem,  "ERROR: For argument, string shouble be specified.\n");
  p_str = (*str_func_ptr)( tmp_str );


  *pp_str = p_str;  // Never do "pp_str = &p_str;"
  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_temp_pp_str( vmstack, pp_str);

//  The following is a code previously used. This results in creating new anonymous object for each row.
//  ptr_record* pr = ptr_table_create_anonym_string(table, pp_str);
//  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
//  vm_stack_push_pp_str( vmstack , table , pr->key );
  return 1;
}

int
sailr_func_str_strip( vm_stack* vmstack, int num_args, ptr_table** table )
{
  return sailr_func_str_func_ptr(vmstack, num_args, table, &string_strip);
}

int
sailr_func_str_lstrip( vm_stack* vmstack, int num_args, ptr_table** table )
{
  return sailr_func_str_func_ptr(vmstack, num_args, table, &string_lstrip); 
}
 
int
sailr_func_str_rstrip( vm_stack* vmstack, int num_args, ptr_table** table )
{
  return sailr_func_str_func_ptr(vmstack, num_args, table, &string_rstrip);
}

int
sailr_func_str_concat( vm_stack* vmstack, int num_args, ptr_table** table )
{
  if( ! arg_num_should_be_larger_than( num_args, 0 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );

  string_object* p_str = string_new(""); 
  string_object** pp_str = (string_object**) malloc(sizeof(string_object*));

  append_arg_list_as_string( p_str, arglist, vmstack );

  *pp_str = p_str;   // Never do "pp_str = &p_str;"
  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_temp_pp_str( vmstack, pp_str);

  return 1;
}


int
sailr_func_str_repeat( vm_stack* vmstack, int num_args, ptr_table** table )
{
  if( ! arg_num_should_be( num_args, 2 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 

  int repeat_num ;

  string_object* p_str ;
  string_object** pp_str = (string_object*) malloc(sizeof(string_object*));
  string_object* ori_str;

  ASSIGN_STRING_PTR( ori_str , argitem , "ERROR: For 1st argument, string shouble be specified.\n" );

  arg_item_next(&argitem);

  ASSIGN_INT_VALUE( repeat_num , argitem , "ERROR: For 2nd argument, int value shouble be specified.\n" ); 

  p_str = string_repeat( ori_str , repeat_num );
  *pp_str = p_str; // Never do "pp_str = &p_str;"
  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_temp_pp_str( vmstack, pp_str);

  return 1;
}


// Arguments for str_subset() function are one-indexed.
int sailr_func_str_subset( vm_stack* vmstack, int num_args, ptr_table** table )
{
  if( ! arg_num_should_be( num_args, 3 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 

  int from_idx, to_idx ;

  string_object* p_str ;
  string_object** pp_str = (string_object*) malloc(sizeof(string_object*));
  string_object* ori_str;

  ASSIGN_STRING_PTR( ori_str , argitem , "ERROR: For 1st argument, string shouble be specified.\n" );

  arg_item_next(&argitem);

  ASSIGN_INT_VALUE( from_idx , argitem , "ERROR: For 2nd argument, int value shouble be specified.\n" ); 

  arg_item_next(&argitem);

  ASSIGN_INT_VALUE( to_idx , argitem , "ERROR: For 3rd argument, int value shouble be specified.\n" ); 

  // Arguments for Sailr str_subset() are one-indexed.
  // Arguments for C/C++ string_subset() are zero-indexed.
  if(from_idx >= 1){
    from_idx = from_idx - 1;
  }else{
    from_idx = 0;
  }

  if(to_idx >= 1){
    to_idx = to_idx - 1;
  }else{
    to_idx = 0;
  }

  p_str = string_subset( ori_str , from_idx, to_idx , vm_stack_get_encoding(vmstack));
  *pp_str = p_str; // Never do "pp_str = &p_str;"
  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_temp_pp_str( vmstack, pp_str);

  return 1;
}


int
sailr_func_str_to_num ( vm_stack* vmstack, int num_args )
{
  if( ! arg_num_should_be( num_args, 1 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 

  if(arg_item_confirm_string(argitem)){
    // The argument should be string.
  }else{
    Rprintf("ERROR: For argument, string shouble be specified.\n");
  }

  int ival;
  double dval;
  char ival_or_dval;

  string_object* str = arg_item_string_obj(argitem);
  if( string_has_char( str, '.')  == 1 ){  // This suggests double.
	dval = string_str2double( str );
    ival_or_dval = 'd';
  } else {                                 // This suggests integer.
	ival = string_str2int( str );
    ival_or_dval = 'i';
  }

  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 

  if(ival_or_dval == 'd'){
    vm_stack_push_dval( vmstack , dval);
  } else if(ival_or_dval == 'i'){
    vm_stack_push_ival( vmstack , ival);
  }
  return 1;
}

int
sailr_func_rexp_matched( vm_stack* vmstack, int num_args , ptr_table** table )
{
  if( ! arg_num_should_be( num_args, 1 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 

  simple_re** ptr_rexp_last_matched = vm_stack_get_ptr_last_rexp_field( vmstack );

#ifdef DEBUG
  simple_re* ptr_last_rexp_obj = *ptr_rexp_last_matched;
  Rprintf("Just before back referencing the last regexp match information: Address to last_rexp object %p\n",  ptr_last_rexp_obj);
#endif

  int matched_group_index;

  ASSIGN_INT_VALUE( matched_group_index , argitem , "ERROR: For the argument, int value shouble be specified.\n" ); 

  char* matched_str = simple_re_matched_str( *ptr_rexp_last_matched , matched_group_index );
  string_object* str_obj = string_new(matched_str);
  free(matched_str);

  string_object** pp_str = (string_object**) malloc(sizeof(string_object*));
  *pp_str = str_obj;
  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_temp_pp_str( vmstack, pp_str);

  return 1;
}

int
sailr_func_date_ymd(vm_stack* vmstack, int num_args ) 
{
  if( ! arg_num_should_be( num_args, 3 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 
  
  int y;
  ASSIGN_INT_VALUE( y , argitem , "ERROR: For 1st argument, int value shouble be specified.\n" ); 

  arg_item_next(&argitem);
  int m;
  ASSIGN_INT_VALUE( m , argitem , "ERROR: For 2nd argument, int value shouble be specified.\n" ); 

  arg_item_next(&argitem);
  int d;
  ASSIGN_INT_VALUE( d , argitem , "ERROR: For 3rd argument, int value shouble be specified.\n" ); 

  int idate = simple_date_ymd( y, m , d);

  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_ival( vmstack , idate );

  return 1;
}

int
sailr_func_date_ym_weekday_nth(vm_stack* vmstack, int num_args ) 
{
  if( ! arg_num_should_be( num_args, 4 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 
  
  int y;
  ASSIGN_INT_VALUE( y , argitem , "ERROR: For 1st argument, int value shouble be specified.\n" ); 

  arg_item_next(&argitem);
  int m;
  ASSIGN_INT_VALUE( m , argitem , "ERROR: For 2nd argument, int value shouble be specified.\n" ); 

  arg_item_next(&argitem);
  string_object* str_wd;
  ASSIGN_STRING_PTR( str_wd , argitem , "ERROR: For 3rd argument, string object shouble be specified.\n" );

  arg_item_next(&argitem);
  int n_th;
  ASSIGN_INT_VALUE( n_th , argitem , "ERROR: For 4th argument, int value shouble be specified.\n" ); 

  int idate = simple_date_ym_weekday_nth( y, m , string_read( str_wd ), n_th);

  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_ival( vmstack , idate );

  return 1;
}


int
sailr_func_date_add_n_unit(vm_stack* vmstack, int num_args, int (*add_n_unit)(int , int) ) 
{
  if( ! arg_num_should_be( num_args, 2 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 
  
  int unix_date;
  ASSIGN_INT_VALUE( unix_date , argitem , "ERROR: For 1st argument, int value shouble be specified.\n" ); 

  arg_item_next(&argitem);
  int add_period;
  ASSIGN_INT_VALUE( add_period , argitem , "ERROR: For 2nd argument, int value shouble be specified.\n" ); 

  int idate = (*add_n_unit)( unix_date, add_period);

  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_ival( vmstack , idate );

  return 1;
}


int
sailr_func_date_add_n_years(vm_stack* vmstack, int num_args ) 
{
  return sailr_func_date_add_n_unit(vmstack, num_args, &simple_date_add_n_years);
}

int
sailr_func_date_add_n_months(vm_stack* vmstack, int num_args ) 
{
  return sailr_func_date_add_n_unit(vmstack, num_args, &simple_date_add_n_months);
}

int
sailr_func_date_add_n_days(vm_stack* vmstack, int num_args ) 
{
  return sailr_func_date_add_n_unit(vmstack, num_args, &simple_date_add_n_days);
}



int
sailr_func_date_format(vm_stack* vmstack, int num_args , ptr_table** table )
{
  if( ! arg_num_should_be( num_args, 2 ) ){
    return 0;
  }
  arg_list* arglist = arg_list_initialize( vmstack, num_args );
  arg_item* argitem = arglist; 

  int unix_date ;

  string_object* p_str ;
  string_object** pp_str = (string_object*) malloc(sizeof(string_object*));
  string_object* fmt_str;

  ASSIGN_INT_VALUE( unix_date , argitem , "ERROR: For 1st argument, int value shouble be specified.\n" ); 

  arg_item_next(&argitem);

  ASSIGN_STRING_PTR( fmt_str , argitem , "ERROR: For 2nd argument, string shouble be specified.\n" );

  char* formatted_cstr = simple_date_new_cstr_format ( unix_date, string_read(fmt_str));
  p_str = string_new( formatted_cstr );
  free(formatted_cstr);

  *pp_str = p_str; // Never do "pp_str = &p_str;"
  arg_list_finalize( vmstack, num_args, arglist ); // Move sp pointer & deallocate arg_list* 
  vm_stack_push_temp_pp_str( vmstack, pp_str);

  return 1;
}




