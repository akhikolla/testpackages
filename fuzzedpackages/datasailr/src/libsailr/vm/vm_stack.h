#ifndef VM_STACK_H
#define VM_STACK_H

#include <stdbool.h>

#include "vm_code.h"
#include "ptr_table.h"
#include "helper.h"
#include "simple_re/simple_re.h"
#include "string/common_string.h"

#define DEFAULT_VM_CHARACTER_ENCODING  SAILR_DEFAULT_VM_CHARACTER_ENCODING
#define MAXSTACKSIZE 1000

#define JUST_A_VALUE  NULL
#define TEMP_OBJECT  NULL
#define NOT_ON_PTR_TABLE  NULL

// Stack type and name
// 
// NULL_ITEM : The corresponding variable is not defined yet (type of the variable is still unknown.)
// VOID_ITEM : When the stack item is no longer used, it is marked as VOID_ITEM.
// INFO_ITEM : The zero-index item holds information of the vm_stack.

#define VM_STACK_ITEM_NAME_TABLE \
	Y(IVAL, "IVAL") /* The stack is holding an int. */ \
	Y(DVAL, "DVAL") /* The stack is holding a double. */ \
	Y(BOOLEAN, "BOOLENAN") \
	Y(PP_IVAL, "PP_IVAL") \
	Y(PP_DVAL, "PP_DVAL") \
	Y(PP_STR, "PP_STR") \
	Y(PP_REXP, "PP_REXP") \
	Y(NULL_ITEM, "NULL_ITEM")\
	Y(VOID_ITEM, "VOID_ITEM")\
	Y(INFO_ITEM, "INFO_ITEM")

#define Y(a, b) a,
enum _ItemType {
  VM_STACK_ITEM_NAME_TABLE
};
#undef Y

typedef enum _ItemType ItemType;

char* display_item_type(ItemType type);

// Structure to store VM stack information
struct _vm_stack_info{
  const char* characterEncoding;
  int max_size;
  simple_re* last_rexp; // Pointer to hold the last executed regular expression.
  int vm_code_pos;
};
typedef struct _vm_stack_info vm_stack_info;

// Stack structure for VM

struct _stack_item {
	ItemType type ; 
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)  /* C11 */
	union {
#endif
		int ival;
		double dval;
		bool boolean;
		int** pp_ival;
		double** pp_dval;
		string_object** pp_str; 
		simple_re** pp_rexp;
		vm_stack_info* p_vm_stack_info;
		void* ptr; // Used by NULL_ITEM
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
	};
#endif
    ptr_record* p_record;  /* Pointer to ptr_record. If the object is just temporary, NOT_ON_PTR_TABLE(=NULL) is assigned. */
};
typedef struct _stack_item stack_item ;

struct _vm_stack{
	int sp ; 
	stack_item stack[ MAXSTACKSIZE ];
	unsigned int error;
};
typedef struct _vm_stack vm_stack;


// Function Prototypes 

// public
// int run_vm ( vm_stack* , vm_code* , ptr_table* ); Deprecated

// private: manipulate stack information
vm_stack* vm_stack_init();
int vm_stack_set_encoding(vm_stack* , const char*);
const char* vm_stack_get_encoding(vm_stack* );
simple_re** vm_stack_get_ptr_last_rexp_field(vm_stack* vmstack);
void vm_stack_clear_last_rexp_hisotry(vm_stack* vmstack);

// public
void vm_stack_set_code_position(vm_stack* vmstack, int pos);
int vm_stack_get_code_position(vm_stack* vmstack);

// private: manipulate stack
int vm_stack_push_item( vm_stack* , stack_item* );
int vm_stack_push_ival( vm_stack* , int );
int vm_stack_push_dval( vm_stack* , double );
int vm_stack_push_pp_ival( vm_stack* , ptr_table**, char* );
int vm_stack_push_pp_dval( vm_stack* , ptr_table**, char* );
int vm_stack_push_pp_num( vm_stack* , ptr_table**, char* );
int vm_stack_push_pp_str( vm_stack* , ptr_table**, char* );
int vm_stack_push_temp_pp_str( vm_stack*, string_object** );
int vm_stack_push_pp_rexp( vm_stack* , ptr_table**, char* );
int vm_stack_push_null( vm_stack* , ptr_table**, char* );
int vm_stack_push_corresp_item( vm_stack* , ptr_table** , char* );
int vm_stack_push_boolean( vm_stack* , bool );

int vm_stack_fcall( vm_stack* , char* , int , ptr_table** );

stack_item* vm_stack_pop( vm_stack* );
int vm_stack_clean_top(vm_stack* vmstack);
int vm_stack_clean_and_pop( vm_stack* , int n);
int vm_stack_clean_items_from_zero_to_top( vm_stack* );

bool vm_stack_item_is_temp( stack_item* item );
int vm_stack_display_item( vm_stack*, int );
void vm_stack_display_all( vm_stack* );
int vm_stack_end( vm_stack* );
int vm_stack_is_full( vm_stack* );
int vm_stack_is_empty( vm_stack* );
int vm_stack_size( vm_stack* );
int vm_stack_free( vm_stack* );


stack_item* vm_stack_top( vm_stack*);
stack_item* vm_stack_second( vm_stack* );
stack_item* vm_stack_third( vm_stack* );
stack_item* vm_stack_nth( vm_stack*, int );

#endif

