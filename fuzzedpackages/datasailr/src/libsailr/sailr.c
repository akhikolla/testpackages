#include <R_ext/Print.h>
/*
    libsailr - virtual machine library for arithmetic calculation and string manipulation 
    Copyright (C) 2018-2019 Toshi Umehara <toshi@niceume.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include <stdio.h>

#include "sailr.h"
#include "ptr_table.h"
#include "var_hash.h"
#include "parser_state.h"
#include "string/common_string.h"
#include "node.h"
#include "tree_dump.h"
#include "tree_free.h"
#include "gen_code.h"
#include "gen_code_util.h"
#include "vm/vm_code.h"
#include "vm/vm_stack.h"
#include "vm/vm.h"
#include "vm/func/ext_func/ext_func_hash.h"
#include "y.tab.h"

// To prevent implicit declaration warnings
int yylex_init(void*);
void* yy_scan_string(const char*, void*);
int yylex_destroy(void*);


ptr_table_object*
sailr_ptr_table_init()
{
	return (ptr_table_object*) ptr_table_init();
}

parser_state_object*
sailr_new_parser_state(const char* fname, ptr_table_object* table)
{
	return (parser_state_object*) new_parser_state(fname, (ptr_table*) table);
}

int
sailr_parser_state_set_source_encoding(parser_state_object* ps , const char* source_encoding)
{
	((parser_state*) ps)->rexp_encoding = source_encoding;
	return 1;
}

const char*
sailr_parser_state_get_source_encoding(parser_state_object* ps )
{
	return (((parser_state*) ps)->rexp_encoding);
}

int
sailr_run_parser (const char* code , parser_state_object* ps)
{
	/* %option reentrant in lex.l enables reentrant scanner (=tokenizer) */
	/* We need to track the tokenizer object, and pass it to parser appropriately.*/
	/* The tokenizer is initialized by yylex_init() */
	/* The parser function brings code into tokenizer's buffer by yy_scan_string() */
	/* The parser function works with tokenizer in yyparse*() */
	/* The tokenizer is destroyed by yylex_destroy() */
	void* scanner;
	yylex_init(&scanner);
	yy_scan_string(code, scanner);
	int result = yyparse( (parser_state*) ps, scanner );
    yylex_destroy(scanner);
	return result;
}

int
sailr_parser_state_free(parser_state_object* ps)
{
	return parser_state_free((parser_state*) ps);
}

void
sailr_tree_dump( parser_state_object* ps )
{
	tree_dump( ((parser_state*) ps)->tree, 0);
}

void
sailr_tree_free(parser_state_object* ps)
{
	tree_free( ((parser_state*) ps)->tree, 0);
}

vm_inst_object*
sailr_gen_code( parser_state_object* ps , ptr_table_object* table)
{
	vm_inst* result = gen_code( ((parser_state*) ps)->tree, (ptr_table*) table );
	return (vm_inst_object*) result; 
}

vm_stack_object*
sailr_vm_stack_init()
{
	return (vm_stack_object*) vm_stack_init();
}

int
sailr_vm_stack_set_encoding(vm_stack_object* vmstack , const char* encoding)
{
	return vm_stack_set_encoding( (vm_stack*) vmstack , encoding);
}

const char*
sailr_vm_stack_get_encoding(vm_stack_object* vmstack)
{
	return vm_stack_get_encoding( (vm_stack*) vmstack);
}


vm_inst_object*
sailr_vm_inst_list_to_code( vm_inst_object* insts)
{
	return (vm_inst_object*) vm_inst_list_to_code((vm_inst_list*) insts);
}


int
sailr_vm_inst_list_size( vm_inst_object* insts)
{
	return vm_inst_list_size((vm_inst_list*) insts);
}

void
sailr_vm_inst_list_show_all( vm_inst_object* insts )
{
	vm_inst* inst = (vm_inst*) insts;
    vm_inst_list_show_all( inst );
    
//	do{
//		vm_inst_show(inst);		
//		inst = inst->next;
//	}while(inst != NULL);
}

void
sailr_vm_inst_list_free( vm_inst_object* inst_list )
{
  vm_inst_list_free( (vm_inst_list*) inst_list ); 
}

void
sailr_vm_inst_code_free( vm_inst_object* vmcode )
{
  free( (vm_inst*) vmcode );
}

// Set NULL for ext_func_hash_object* extfunc_hash argument, when no external functions are added. 
int
sailr_vm_exec_code( vm_inst_object* code , int num_insts , ptr_table_object* table , vm_stack_object* vmstack, ext_func_hash_object* extfunc_hash )
{
	return vm_exec_code((vm_inst*)code, num_insts, 0, (ptr_table*)table, (vm_stack*)vmstack, (ext_func_hash*) extfunc_hash );
}

int
sailr_vm_resume_code( vm_inst_object* code , int num_insts , int start_inst_idx, ptr_table_object* table , vm_stack_object* vmstack, ext_func_hash_object* extfunc_hash)
{
	return vm_exec_code((vm_inst*)code, num_insts, start_inst_idx , (ptr_table*)table, (vm_stack*)vmstack, (ext_func_hash*) extfunc_hash );
}

int
sailr_vm_stack_get_code_position( vm_stack_object* vmstack)
{
	return vm_stack_get_code_position((vm_stack*) vmstack); 
}

ptr_record_object*
sailr_ptr_table_create_int_from_ptr(ptr_table_object** table, const char* key, int** i_pp, double** d_pp)
{
	return (ptr_record_object*) ptr_table_create_int_from_ptr((ptr_table**)table, key, i_pp, d_pp);
}

ptr_record_object*
sailr_ptr_table_create_double_from_ptr(ptr_table_object** table, const char* key, double** d_pp, int** i_pp)
{
	return (ptr_record_object*) ptr_table_create_double_from_ptr((ptr_table**)table, key, d_pp, i_pp);
}

ptr_record_object*
sailr_ptr_table_create_null(ptr_table_object** table, const char* key)
{
	return (ptr_record_object*) ptr_table_create_null((ptr_table**)table, key);
}

ptr_record_object*
sailr_ptr_table_create_anonym_string(ptr_table_object** table, const char* str)
{
	string_object* new_str = string_new(str);
	return (ptr_record_object*) ptr_table_create_anonym_string((ptr_table**)table, &new_str);
}


ptr_record_object*
sailr_ptr_table_create_string_from_cstring(ptr_table_object** table, const char* key, const char* str)
{
	return (ptr_record_object*) ptr_table_create_string_from_cstring((ptr_table**)table, key, str);
}


// Deprecated: hard to mange string object.
// ptr_record_object*
// sailr_ptr_table_create_string_from_ptr(ptr_table_object** table, const char* key, string_type_object** pp)
// {
// 	return (ptr_record_object*) ptr_table_create_string_from_ptr((ptr_table**)table, key, (string_object**)pp);
// 
// }
//

char
sailr_ptr_table_get_type(ptr_table_object** table, const char* key)
{
	ptr_record* record = ptr_table_find((ptr_table**) table, key);
	return sailr_ptr_record_get_type(record);
}


char
sailr_ptr_record_get_type(ptr_record_object* pr)
{
	ptr_record* record = (ptr_record*) pr;
	switch(record->type){
		case PTR_INT:
			return 'i';
		break;
		case PTR_DBL:
			return 'd';
		break;
		case PTR_STR:
			return 's';
		break;
		case PTR_REXP:
			return 'r';
		break;
		case PTR_NULL:
			return 'n';
		break;
		case PTR_INFO:
			return 'f';
		break;
		default:
			Rprintf("ERROR: Unintended type of ptr_record. %c \n", record->type);
			return 'x';
		break;
	}
}

int
sailr_ptr_record_is_anonym(ptr_record_object* pr)
{
	if( ptr_record_get_anonym((ptr_record*) pr) == 1 ){
		return 1;
	}else{
		return 0;
	}
}

int sailr_ptr_record_is_ptr_null(ptr_table_object** table, const char* key)
{
	return ptr_record_is_ptr_null((ptr_table**) table, key);
}

void**
sailr_ptr_table_get_pptr(ptr_table_object** table, const char* key)
{
	return ptr_table_get_pptr((ptr_table**) table, key);
}

// (Deprecated) This string object may be deleted without notifying user. Not recommended to use this function.
//string_type_object*
//sailr_ptr_table_get_ptr_string(ptr_table_object** table, const char* key){
//  return (string_type_object*) ptr_table_get_ptr_string((ptr_table**)table, key);
//}

const char*
sailr_ptr_table_read_string(ptr_table_object** table, const char* key){
  return  ptr_table_read_string((ptr_table**)table, key);
}

ptr_record_object*
sailr_ptr_table_find( ptr_table_object** table, const char* key )
{
	return (ptr_record_object*) ptr_table_find((ptr_table**) table, key);
}

ptr_record_object*
sailr_ptr_table_first_record(ptr_table_object** table)
{
	return (ptr_record_object*) ptr_table_first_record((ptr_table**)table);
}

// If this is the last ptr_record, this function returns NULL.
ptr_record_object*
sailr_ptr_record_next(ptr_record_object* pr )
{
	ptr_record* pr_next = ptr_record_next((ptr_record*) pr);
	return (ptr_record_object*) pr_next ; 
}

int
sailr_ptr_table_update_int(ptr_table_object** table, const char* key, int ival)
{
	return ptr_table_update_int((ptr_table**) table, key, ival);
}

int
sailr_ptr_table_update_double(ptr_table_object** table, const char* key, double dval)
{
	return ptr_table_update_double((ptr_table**) table, key, dval);
}

int
sailr_ptr_table_update_string(ptr_table_object** table, const char* key, string_type_object** str)
{
	return ptr_table_update_string((ptr_table**) table, key, (string_object**) str);
}

int
sailr_ptr_record_reset_rexp(ptr_record_object* pr)
{
	return ptr_record_reset_rexp((ptr_record*) pr);
}

void
sailr_ptr_table_free_objects(ptr_table_object** table, const char* key)
{
	ptr_record* record = ptr_table_find((ptr_table**) table, key);
	ptr_record_free_gc_required_memory(record);
}

void
sailr_ptr_record_free_objects(ptr_record_object* pr)
{
	ptr_record_free_gc_required_memory(pr);
}

void
sailr_ptr_table_del_records_except(ptr_table_object** table, const char** keys, int key_num )
{
	ptr_table_del_records_except((ptr_table**)table, keys, key_num );
}

void
sailr_ptr_table_del_all(ptr_table_object** table)
{
	ptr_table_del_all((ptr_table**) table);
}


// Utility

void
sailr_ptr_table_show_all(ptr_table_object** table)
{
	ptr_table_show_all((ptr_table**) table);
}

int // 0 : non-null_updated.  1: null_updated.
sailr_ptr_table_info_get_null_updated( ptr_table_object** table)
{
	return ptr_table_info_get_null_updated((ptr_table**) table);
}

int
sailr_ptr_table_info_reset_null_updated( ptr_table_object** table)
{
	return ptr_table_info_reset_null_updated((ptr_table**) table);
}


string_type_object*
sailr_new_string(const char* str)
{
	string_object* new_str = string_new(str);
	return (string_type_object*) new_str;
}


// Functions fro var_hash, varnames

char**
sailr_varnames(parser_state_object* psobj)
{
	parser_state* ps = (parser_state*) psobj;
	return var_hash_names(&(ps->vars));
}

char**
sailr_rhs_varnames(parser_state_object* psobj)
{
	parser_state* ps = (parser_state*) psobj;
	return var_hash_names(&(ps->rhsvars));
}

char**
sailr_lhs_varnames(parser_state_object* psobj)
{
	parser_state* ps = (parser_state*) psobj;
	return var_hash_names(&(ps->lhsvars));
}

int
sailr_varnames_num(parser_state_object* psobj)
{
	parser_state* ps = (parser_state*) psobj;
	return var_hash_size(&(ps->vars));
}

int
sailr_rhs_varnames_num(parser_state_object* psobj)
{
	parser_state* ps = (parser_state*) psobj;
	return var_hash_size(&(ps->rhsvars));
}

int
sailr_lhs_varnames_num(parser_state_object* psobj)
{
	parser_state* ps = (parser_state*) psobj;
	return var_hash_size(&(ps->lhsvars));
}

void
sailr_varnames_free( char** varnames , int size)
{
	var_hash_names_free( varnames, size );
}

