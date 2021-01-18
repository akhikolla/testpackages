// ptr_table stores address information for variables and strings. 
// Virtual stack machine does not directly manipulate these information on its stack.

#ifndef PTR_TABLE_H
#define PTR_TABLE_H

#include "uthash.h"
#include "string/common_string.h"
#include <stdbool.h>

#define MAX_KEY_LEN 511
#define ANONYM_KEY_WIDTH 16

enum _PtrType{
	PTR_INT = 0,
	PTR_DBL = 1,
	PTR_STR = 2,
	PTR_REXP = 3,
	PTR_NULL = 4,
	PTR_INFO = 5
};
typedef enum _PtrType PtrType ;

enum _GCReq {
	GC_NO,
	GC_YES
};
typedef enum _GCReq GCReq ;
// GC_NO: 0, GC_YES: 1

struct _ptr_record {
	char key[ MAX_KEY_LEN ];  /* Usually this should be a variable name. For string, temporary name is generated. */
	void* address;  /* value1 */
	PtrType type; /* value2 */
	GCReq gc; /* value3: whether garbage collected. 0: GC_NO ,  1: GC_YES */
	/* Mainly used for numbers which requires memory for double and int. */
	void* ex_addr; 
	PtrType ex_type;
	GCReq ex_gc; 
	int anonym; /* When the record is created to be anonymous, this is set to be 1. Otherwise 0. */
	UT_hash_handle hh; /* This macro makes this structure hashable */
};
typedef struct _ptr_record ptr_record ;
typedef ptr_record ptr_table;

struct _ptr_table_info {
	int str_counter;
	int rexp_counter;
	unsigned int null_updated;
};
typedef struct _ptr_table_info ptr_table_info ;

ptr_table*  ptr_table_init();
ptr_record* ptr_table_add(ptr_table** table, const char* key, void** address, PtrType type, GCReq gc);

ptr_record* ptr_table_create_int(ptr_table** table, const char* key, int ival);
ptr_record* ptr_table_create_int_from_ptr(ptr_table** table, const char* key, int** iptr, double** dptr);
int ptr_table_update_int(ptr_table** table, const char* key, int ival);

ptr_record* ptr_table_create_double(ptr_table** table, const char* key, double dval);
ptr_record* ptr_table_create_double_from_ptr(ptr_table** table, const char* key, double** dptr, int** iptr);
int ptr_table_update_double(ptr_table** table, const char* key, double dval);

int ptr_record_update_extra_address(ptr_record* pr, void** ptr_ex_addr, PtrType ex_type, GCReq ex_gc );
int ptr_record_swap_addresses(ptr_record* pr);


void ptr_record_set_anonym( ptr_record* pr, int val);
int ptr_record_get_anonym( ptr_record* pr);

ptr_record* ptr_table_create_anonym_string(ptr_table** table, string_object** strptr);
ptr_record* ptr_table_create_string_from_cstring(ptr_table** table, const char* key, const char* str);
// (DEPRECATED) ptr_record* ptr_table_create_string_from_ptr(ptr_table** table, const char* key, string_object** strptr);

string_object* ptr_table_get_ptr_string(ptr_table** table, const char* key);
const char* ptr_table_read_string(ptr_table** table, const char* key);
int ptr_table_update_string(ptr_table** , const char* key, string_object** );
int ptr_record_update_string(ptr_record* pr , string_object** pp_str, GCReq gc);

ptr_record* ptr_table_create_anonym_rexp(ptr_table** table, const char* pattern, const char* enc);
int ptr_record_reset_rexp(ptr_record* pr);

ptr_record* ptr_table_create_null(ptr_table** table, const char* key);
int ptr_table_del_record(ptr_table** table, const char* key);
void ptr_table_del_records_except(ptr_table** table, const char** keys, int key_num );
void ptr_table_del_all(ptr_table** table);

void ptr_table_show_all(ptr_table** table);
void ptr_record_show(ptr_record* pr);
ptr_table* ptr_record_obtain_table(ptr_record* pr);
// int ptr_table_info_set_null_updated(ptr_table** table, int updated_value);
int ptr_table_info_change_null_updated_by_type(ptr_table** table, PtrType type);
int ptr_table_info_get_null_updated(ptr_table** table);
int ptr_table_info_reset_null_updated(ptr_table** table);

PtrType ptr_table_get_type(ptr_table** table, const char* key);
int ptr_record_is_ptr_null(ptr_table** table, const char* key);
void** ptr_table_get_pptr(ptr_table** table, const char* key);

ptr_record* ptr_table_first_record(ptr_table** table);
ptr_record* ptr_record_next(ptr_record* pr);

void ptr_record_free_gc_required_memory(ptr_record*);

ptr_record* ptr_table_find(ptr_table** table, const char* key);

// private
ptr_record* ptr_table_insert(ptr_table** table, ptr_record* pr);
int ptr_record_update(ptr_record* pr, void* address, PtrType type, GCReq gc);
void ptr_record_free(ptr_record* );
int ptr_table_free(ptr_table**);

bool ptr_table_points_to_header(ptr_table** table);

#endif /* PTR_TABLE_H */
