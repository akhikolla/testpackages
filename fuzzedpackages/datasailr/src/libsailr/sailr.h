#ifndef SAILR_H
#define SAILR_H

typedef void ptr_table_object;
typedef void ptr_record_object;
typedef void parser_state_object;
typedef void vm_inst_object;
typedef void vm_stack_object;
typedef void string_type_object;
typedef void ext_func_hash_object;

ptr_table_object*     sailr_ptr_table_init();
parser_state_object*  sailr_new_parser_state(const char* , ptr_table_object* );
int                   sailr_parser_state_set_source_encoding(parser_state_object* ps , const char* source_encoding);
const char*           sailr_parser_state_get_source_encoding(parser_state_object* ps );
int                   sailr_run_parser (const char* code , parser_state_object* ps); 
int                   sailr_parser_state_free(parser_state_object* ps);
void                  sailr_tree_dump( parser_state_object* ps ); 
void                  sailr_tree_free( parser_state_object* ps ); 
vm_inst_object*       sailr_gen_code( parser_state_object* ps , ptr_table_object*); // VM Code is generated.
vm_inst_object*       sailr_vm_inst_list_to_code( vm_inst_object* );
int                   sailr_vm_inst_list_size( vm_inst_object* );
void                  sailr_vm_inst_list_show_all( vm_inst_object* );
void                  sailr_vm_inst_list_free( vm_inst_object* inst_list );
void                  sailr_vm_inst_code_free( vm_inst_object* vmcode );
vm_stack_object*      sailr_vm_stack_init();
int                   sailr_vm_stack_set_encoding(vm_stack_object* , const char*);
const char*           sailr_vm_stack_get_encoding(vm_stack_object* );
int                   sailr_vm_exec_code( vm_inst_object* code , int num_insts , ptr_table_object* table , vm_stack_object* vmstack , ext_func_hash_object* extfunc_hash);
int                   sailr_vm_resume_code( vm_inst_object* code , int num_insts , int start_inst_idx, ptr_table_object* table , vm_stack_object* vmstack, ext_func_hash_object* extfunc_hash);
int                   sailr_vm_stack_get_code_position( vm_stack_object*);

// Create
ptr_record_object* sailr_ptr_table_create_int_from_ptr(ptr_table_object** table, const char* key, int** i_pp, double** d_pp);
ptr_record_object* sailr_ptr_table_create_double_from_ptr(ptr_table_object** table, const char* key, double** d_pp, int** i_pp);
ptr_record_object* sailr_ptr_table_create_anonym_string(ptr_table_object** table, const char* str);
ptr_record_object* sailr_ptr_table_create_string_from_cstring(ptr_table_object** table, const char* key, const char* str);
ptr_record_object* sailr_ptr_table_create_null(ptr_table_object** table, const char* key);

// (Deprecated) ptr_record_object* sailr_ptr_table_create_string_from_ptr(ptr_table_object** table, const char* key, string_type_object** pp); 

// Read
char sailr_ptr_table_get_type(ptr_table_object** table, const char* key);
char sailr_ptr_record_get_type(ptr_record_object* pr);
int sailr_ptr_record_is_anonym(ptr_record_object* pr);
int sailr_ptr_record_is_ptr_null(ptr_table_object** table, const char* key);
void** sailr_ptr_table_get_pptr(ptr_table_object** table, const char* key);
const char* sailr_ptr_table_read_string(ptr_table_object** table, const char* key);

ptr_record_object* sailr_ptr_table_find( ptr_table_object** table, const char* key );
ptr_record_object* sailr_ptr_table_first_record(ptr_table_object** table);
ptr_record_object* sailr_ptr_record_next(ptr_record_object* pr );

// (Deprecated) string_type_object* sailr_ptr_table_get_ptr_string(ptr_table_object** table, const char* key);

// Update
int sailr_ptr_table_update_int(ptr_table_object** table, const char* key, int ival);
int sailr_ptr_table_update_double(ptr_table_object** table, const char* key, double dval);
int sailr_ptr_table_update_string(ptr_table_object** table, const char* key, string_type_object** str);
int sailr_ptr_record_reset_rexp(ptr_record_object* pr);
void sailr_ptr_table_free_objects(ptr_table_object** table, const char* key);
void sailr_ptr_record_free_objects(ptr_record_object* pr);

// Delete
void sailr_ptr_table_del_records_except(ptr_table_object** table, const char** keys, int key_num );
void sailr_ptr_table_del_all(ptr_table_object** table);

// Utility
void sailr_ptr_table_show_all(ptr_table_object** table);
int sailr_ptr_table_info_get_null_updated( ptr_table_object** table);
int sailr_ptr_table_info_reset_null_updated( ptr_table_object** table);

string_type_object* sailr_new_string(const char* str);

char** sailr_varnames(parser_state_object* );
char** sailr_rhs_varnames(parser_state_object* );
char** sailr_lhs_varnames(parser_state_object* );
int sailr_varnames_num(parser_state_object* psobj);
int sailr_rhs_varnames_num(parser_state_object* psobj);
int sailr_lhs_varnames_num(parser_state_object* psobj);
void sailr_varnames_free( char** varnames , int num);

#endif



