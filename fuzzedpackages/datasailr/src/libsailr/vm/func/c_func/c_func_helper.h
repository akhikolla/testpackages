/* This helper file is recommended to be included when defining new functions. */
#ifndef FUNC_HELPER_H
#define FUNC_HELPER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "vm/vm_stack.h"
#include "string/common_string.h"
#include "simple_re/simple_re.h"
#include "simple_date/simple_date.h"

// Define iterator for function arguments.

struct _arg_item{
  struct _arg_item* head;
  struct _arg_item* next;
  struct _arg_item* tail;
  stack_item* item;
};
typedef struct _arg_item  arg_item;
typedef arg_item arg_list;

// Prototypes
arg_item* arg_item_new( stack_item* item );
bool arg_num_should_be( int num_args1, int num_args2 );
bool arg_num_should_be_larger_than( int num_args1, int num_args2 );
bool arg_num_should_be_smaller_than( int num_args1, int num_args2 );
void arg_item_push_back(arg_item** items, arg_item** new_item);
arg_list* arg_list_initialize(vm_stack* vmstack , int num_args );
int arg_list_finalize(vm_stack* vmstack, int num_args, arg_list* alist);
int arg_item_next( arg_item** item );
bool arg_num_confirm_size( int num_args_from_code , int num_args_at_c_level );
bool arg_item_confirm_type( arg_item* arg, ItemType type );
bool arg_item_confirm_int( arg_item* arg);
bool arg_item_confirm_double( arg_item* arg);
bool arg_item_confirm_string( arg_item* arg);
char arg_item_interpret_type( arg_item* arg);
int arg_item_int_value( arg_item* arg );
double arg_item_double_value( arg_item* arg );
string_object* arg_item_string_obj( arg_item* arg);
simple_re* arg_item_rexp_obj(arg_item* arg);
const char* arg_item_string_char( arg_item* arg );
bool arg_item_bool_value( arg_item* arg );
int arg_list_free( arg_list* head);


#endif /* FUNC_HELPER_H */

