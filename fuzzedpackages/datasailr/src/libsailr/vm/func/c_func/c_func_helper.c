#include <R_ext/Print.h>
#include "c_func_helper.h"

// Definitions
bool
arg_num_should_be(int num_args1, int num_args2)
{
  if( num_args1 == num_args2){
    return true;
  } else {
    Rprintf("ERROR: number of args is not specified correctly. Specified: %d , Intended: %d\n", num_args1, num_args2);
    return false;
  }
}

bool
arg_num_should_be_larger_than(int num_args1, int num_args2)
{
  if( num_args1 > num_args2){
    return true;
  } else {
    Rprintf("ERROR: number of args is not specified correctly. Specified: %d , Intended: %d\n", num_args1, num_args2);
    return false;
  }
}

bool
arg_num_should_be_smaller_than(int num_args1, int num_args2)
{
  if( num_args1 < num_args2){
    return true;
  } else {
    Rprintf("ERROR: number of args is not specified correctly. Specified: %d , Intended: %d\n", num_args1, num_args2);
    return false;
  }
}

arg_item*
arg_item_new( stack_item* item )
{
  arg_item* arg = (arg_item*) malloc(sizeof( arg_item ));
  arg->head = arg;
  arg->tail = arg;
  arg->next = NULL;
  arg->item = item;
  return arg;
}

void
arg_item_push_back(arg_item** pp_items, arg_item** pp_new_item)
{
  arg_item* p_items = *pp_items;
  arg_item* p_new_item = *pp_new_item;

  (p_items->tail)->head = p_items;
  (p_items->tail)->next = p_new_item;
  (p_items->tail)->tail = p_new_item;

  p_new_item->head = p_items;
  p_new_item->next = NULL;
  p_new_item->tail = p_new_item;

  p_items->tail = p_new_item;
}

arg_list*
arg_list_initialize(vm_stack* vmstack , int num_args )
{
  int idx ;
  stack_item* curr_item ;
  arg_item* new_arg = NULL ;
  arg_item* args = NULL ;
  for( idx = num_args ; idx >= 1 ; idx = idx - 1 ){
    curr_item = vm_stack_nth( vmstack , idx );
    new_arg = arg_item_new( curr_item );
    if( idx == num_args ){
      args = new_arg;
    } else {
      arg_item_push_back( &args, &new_arg );
    }
  }
  return args;
}

int
arg_list_finalize(vm_stack* vmstack, int num_args, arg_list* alist)
{
  arg_list_free(alist);
  return vm_stack_clean_and_pop(vmstack, num_args);
}

int
arg_item_next( arg_item** item )
{
  *item = (*item)->next;
  if(*item != NULL){
//    printf("The next arg_item is not null.\n");
    return 1;
  }else{
//    printf("The next arg_item is null.\n");
    return 0;
  }
}

bool
arg_num_confirm_size( int num_args_from_code , int num_args_at_c_level )
{
  bool result;
  if(num_args_from_code == num_args_at_c_level){
    result = true;
  }else if(num_args_from_code < num_args_at_c_level){
    result = false;
  }else{ // num_args_from_code > num_args_at_c_level
    result = false;
  }
  return result;
}

bool
arg_item_confirm_type( arg_item* arg, ItemType type )
{
  /* Possible types: IVAL, DVAL, BOOLEAN, PP_IVAL, PP_DVAL, PP_STR, PP_REXP, NULL_ITEM */
  if ( arg->item->type == type )
    return true;
  else
    return false;
}

// confirm int
bool
arg_item_confirm_int( arg_item* arg)
{
  if( arg_item_confirm_type( arg, IVAL) || arg_item_confirm_type( arg, PP_IVAL))
    return true;
  else
    return false;
}

// confirm double
bool
arg_item_confirm_double( arg_item* arg)
{
  if( arg_item_confirm_type( arg, DVAL) || arg_item_confirm_type( arg, PP_DVAL))
    return true;
  else
    return false;
}

// confirm string
bool
arg_item_confirm_string( arg_item* arg)
{
  return arg_item_confirm_type( arg, PP_STR);
}

// confirm regular expression
bool
arg_item_confirm_rexp( arg_item* arg)
{
  return arg_item_confirm_type( arg, PP_REXP);
}

// return type using char
char
arg_item_interpret_type( arg_item* arg)
{
  if(arg_item_confirm_int(arg)){
    return 'i';
  }else if(arg_item_confirm_double(arg)){
    return 'd';
  }else if(arg_item_confirm_string(arg)){
    return 's';
  }else if(arg_item_confirm_rexp(arg)){
    return 'r';
  }else if(arg_item_confirm_type(arg, BOOLEAN)){
    return 'b';
  }else if(arg_item_confirm_type(arg, NULL_ITEM)){
    return 'n';
  }else{
    return 'x';
  }
}

// Obtain int
int
arg_item_int_value( arg_item* arg )
{
  int value;
  if(arg_item_confirm_type(arg, IVAL) ){
    value = arg->item->ival;
  }else if (arg_item_confirm_type(arg, PP_IVAL) ){
    value = **(arg->item->pp_ival);
  }else{
    value = 0; // This branch should never be executed.
    Rprintf("ERROR: the stack item does not hold int value. \n");
  }
  return value;
}

// Obtain double
double
arg_item_double_value( arg_item* arg )
{
  double value;
  if(arg_item_confirm_type(arg, DVAL) ){
    value = arg->item->dval;
  }else if (arg_item_confirm_type(arg, PP_DVAL) ){
    value = **(arg->item->pp_dval);
  }else{
    value = 0.0; // This branch should never be executed.
    Rprintf("ERROR: the stack item does not hold double value. \n");

  }
  return value;
}

// Obtain string object
string_object*
arg_item_string_obj( arg_item* arg)
{
  string_object* temp_obj = NULL;
  string_object** pp_obj;
  if(arg_item_confirm_type(arg, PP_STR) ){
    pp_obj = (arg->item)->pp_str;
    temp_obj = *pp_obj;
  }else{
    Rprintf("ERROR: the stack item does not hold string value. \n");
  }
  return (string_object*) temp_obj;
}


// Obtain const char*
const char* /* const char * is a pointer to a const char */
arg_item_string_char( arg_item* arg )
{
  const char* str = "" ;
  if(arg_item_confirm_type(arg, PP_STR) ){
    str = string_read( (string_object*) *( arg->item->pp_str)) ;
  }else{
    Rprintf("ERROR: the stack item does not hold string value. \n");
  }
  return str;
}

// Obtain boolean
bool
arg_item_bool_value( arg_item* arg )
{
  bool result = false;
  if(arg_item_confirm_type(arg, BOOLEAN) ){
    result = arg->item->boolean;
  }else{
    Rprintf("ERROR: the stack item does not hold boolean value. \n");
  }
  return result;
}

// Obtain regular expression
simple_re*
arg_item_rexp_obj( arg_item* arg)
{
  simple_re* temp_obj = NULL;
  simple_re** pp_obj;
  if(arg_item_confirm_type(arg, PP_REXP) ){
    pp_obj = (arg->item)->pp_rexp;
    temp_obj = *pp_obj;
  }else{
    Rprintf("ERROR: the stack item does not hold rexp value. \n");
  }
  return (simple_re*) temp_obj;
}


// Free
int
arg_list_free( arg_list* head)
{
  arg_item* temp;
  arg_item* curr = head;
  do{
    temp = curr->next;
    free(curr);
    curr = temp;
  }while(temp != NULL);
  return 1;
}

// Show
/*
int
arg_list_show( arg_list* head )
{

}
*/

