#include <R_ext/Print.h>
#include <stdio.h>
#include <string.h>
#include "node.h"
#include "vm/vm_cmd.h"
#include "vm_label.h"
#include "gen_code_util.h"

#include "helper.h"

/* Wrapper macro to add script location information */

#define ADD_LOC( VM_CODE , LOC ) vm_inst_set_loc_to_last( LOC, VM_CODE )


/* Helper Functions ---------------- */

vm_inst*
gen_code_stmt(vm_inst* code1, vm_inst* code2)  // OK
{
	vm_inst* code;
	code = vm_inst_list_cat( code1, code2 );
	return code;
}

vm_inst*
gen_code_stmt_nosib(vm_inst* code1)  // OK. Add VM_END instruction.
{
	vm_inst* code;
	code = code1;
	return code;
}

vm_inst* 
gen_code_int(TreeNode* nd)  // Terminal OK
{
	int ival = nd->e1.ival;
	vm_inst* code = new_vm_inst_push_ival( ival );
	return code; 
}

vm_inst*
gen_code_double(TreeNode* nd)  // Terminal OK
{
	double dval = nd->e1.dval;
	vm_inst* code = new_vm_inst_push_dval( dval );
	return code; 
}

vm_inst*
gen_code_str(TreeNode* nd )  // Terminal OK
{
//	printf("instruction for VM to push pp to STR: %s. \n", nd->e1.str_key);
	char* key = nd->e1.str_key;
	vm_inst* code = new_vm_inst_push_pp_str(key);
	return code; 
}

vm_inst*
gen_code_rexp(TreeNode* nd )  // Terminal OK
{
	char* key = nd->e1.rexp_key;
	vm_inst* code = new_vm_inst_push_pp_rexp(key);
	return code; 
}

vm_inst* 
gen_code_ident(TreeNode* nd, ptr_table* table)  // Terminal OK
{
	vm_inst* code = NULL;
	char* id_name = nd->e1.id;
	ptr_record* record = ptr_table_find( &table, id_name );
	if(record->type == PTR_INT){
		code = new_vm_inst_push_pp_ival(id_name);
	} else if (record->type == PTR_DBL){
		code = new_vm_inst_push_pp_dval(id_name);
	} else if (record->type == PTR_STR){
		code = new_vm_inst_push_pp_str(id_name);
	} else if (record->type == PTR_NULL){
		code = new_vm_inst_push_null(id_name);
	} else {
		Rprintf("ERROR: Inappropriate type is specified for varialbe. \n");
	}
	return code; 
}


VM_CMD
convert_op(char* op_name)
{
	if ( strcmp( op_name, "PLUS") == 0 ) {
		return VM_ADDX;
	} else if ( strcmp( op_name, "SUB") == 0 ) {
		return VM_SUBX;
	} else if ( strcmp( op_name, "MULT") == 0 ) {
		return VM_MULX;
	} else if ( strcmp( op_name, "DIV") == 0 ) {
		return VM_DIVX;
	} else if ( strcmp( op_name, "MOD") == 0 ) {
		return VM_MODX;
	} else if ( strcmp( op_name, "POWER") == 0 ) {
		return VM_POWX;
	} else if ( strcmp( op_name, "FACTOR") == 0 ) {
		return VM_FAC;
	} else if ( strcmp( op_name, "UMINUS") == 0 ) {
		return VM_UMINUS;
	} else if ( strcmp( op_name, "AND") == 0 ) {
		return VM_AND;
	} else if ( strcmp( op_name, "OR") == 0 ) {
		return VM_OR;
	} else if ( strcmp( op_name, "EQ") == 0 ) {
		return VM_EQ;
	} else if ( strcmp( op_name, "NEQ") == 0 ) {
		return VM_NEQ;
	} else if ( strcmp( op_name, "GT") == 0 ) {
		return VM_GT;
	} else if ( strcmp( op_name, "LT") == 0 ) {
		return VM_LT;
	} else if ( strcmp( op_name, "GE") == 0 ) {
		return VM_GE;
	} else if ( strcmp( op_name, "LE") == 0 ) {
		return VM_LE;
	} else if ( strcmp( op_name, "NEG") == 0 ) {
		return VM_NEG;
	} else if ( strcmp( op_name, "REXP_MATCH") == 0 ) {
		return VM_REXP_MATCH;
	} else {
		Rprintf("ERROR: node op has undefined oprator!!\n");
		return VM_NOP;
	}
}

vm_inst*
gen_code_op(VM_CMD cmd, vm_inst* code1, vm_inst* code2) // OK operator
{
	vm_inst* code;
	code = vm_inst_list_cat(code1, code2);

	vm_inst* op_code;
	op_code = new_vm_inst_command(cmd);

	code = vm_inst_list_cat(code, op_code);
	return code;
}

vm_inst*
gen_code_unitary_op(VM_CMD cmd, vm_inst* code1) // OK operator
{
	vm_inst* code;
	vm_inst* op_code;

	op_code = new_vm_inst_command(cmd);

	code = vm_inst_list_cat(code1, op_code);
	return code;
}

vm_inst*
gen_code_let(vm_inst* code1, vm_inst* code2)
{
	vm_inst* code;
	vm_inst* let_code;
	let_code = new_vm_inst_command(VM_STO);
	code = vm_inst_list_cat(code1, code2);
	code = vm_inst_list_cat(code, let_code);
	return code;
}

vm_inst*
gen_code_fcall( char* fname, int num_arg, vm_inst* farg_code)
{
	vm_inst* code;
    vm_inst* fcall_inst;

    fcall_inst = new_vm_inst_command(VM_FCALL);
    int fname_len = strlen(fname);
    if( fname_len < MAX_FUNC_NAME_LEN){
        memcpy( fcall_inst->fname, fname, fname_len + 1 );
    } else { 
        Rprintf("ERROR: function name is too long. over %d.", MAX_FUNC_NAME_LEN);
    }  
    fcall_inst->num_arg = num_arg;

    if(farg_code != NULL){
      code = vm_inst_list_cat( farg_code, fcall_inst);
    } else {
      code = fcall_inst;
    }
    return code;
}

/* END Helper Functions ---------------- */



/* Main Logic */

vm_inst* gen_code(TreeNode* nd, ptr_table* table){
  vm_inst* tmp_code1 = NULL;
  vm_inst* tmp_code2 = NULL;
  vm_inst* tmp_code3 = NULL;
  vm_inst* nd_code = NULL;
  VM_CMD cmd = 0;

  switch (nd->type){
  case NODE_PRGM:
    DEBUG_PRINT("NODE_PRGM\n");
    tmp_code1 = gen_code(nd->e1.nd, table);
    vm_inst* termin_code = new_vm_inst_command(VM_END);
    nd_code = vm_inst_list_cat( tmp_code1 , termin_code );
    return nd_code;
    break;

  case NODE_STMT:
    DEBUG_PRINT("NODE_STMT\n");
    tmp_code1 = gen_code(nd->e1.nd, table);
    if(nd->e3.sibling != NULL){
		tmp_code3 = gen_code(nd->e3.sibling, table);
		nd_code = gen_code_stmt(tmp_code1, tmp_code3);
	}else{
		nd_code = gen_code_stmt_nosib(tmp_code1);
	}
    return nd_code;
    break;

  case NODE_INT:  // terminal node
    DEBUG_PRINT("NODE_INT\n");
    nd_code = gen_code_int(nd);
    ADD_LOC(nd_code, nd->loc);
    return nd_code;
    break;

  case NODE_DBL:  // terminal node
    DEBUG_PRINT("NODE_DBL\n");
    nd_code = gen_code_double(nd);
    ADD_LOC(nd_code, nd->loc);
    return nd_code;
    break;

  case NODE_STR:  // terminal node
    DEBUG_PRINT("NODE_STR\n");
    nd_code = gen_code_str(nd);
    ADD_LOC(nd_code, nd->loc);
    return nd_code;
    break;

  case NODE_REXP:  // terminal node
    DEBUG_PRINT("NODE_REXP\n");
    nd_code = gen_code_rexp(nd);
    ADD_LOC(nd_code, nd->loc);
    return nd_code;
    break;

  case NODE_IDENT:  // terminal node
    DEBUG_PRINT("NODE_IDENT\n");
    nd_code = gen_code_ident(nd, table);
    ADD_LOC(nd_code, nd->loc);
    return nd_code;
    break;

  case NODE_OP:
    DEBUG_PRINT("NODE_OP\n");
    cmd = convert_op(nd->e1.op);
    tmp_code2 = gen_code(nd->e2.nd, table);
    tmp_code3 = gen_code(nd->e3.nd, table);
    nd_code = gen_code_op( cmd , tmp_code2, tmp_code3 );
    ADD_LOC(nd_code, nd->loc);
    return nd_code;
    break;

  case NODE_UNIOP:
    DEBUG_PRINT("NODE_UNIOP\n");
    cmd = convert_op(nd->e1.op);
    tmp_code2 = gen_code(nd->e2.nd, table);
    nd_code = gen_code_unitary_op( cmd , tmp_code2 );
    ADD_LOC(nd_code, nd->loc);
    return nd_code;
    break;

  case NODE_LET:
    DEBUG_PRINT("NODE_LET\n");
    tmp_code1 = gen_code(nd->e1.nd, table);
    tmp_code2 = gen_code(nd->e2.nd, table);
    nd_code = gen_code_let( tmp_code1, tmp_code2 );
    ADD_LOC(nd_code, nd->loc);
    return nd_code;
    break;

  case NODE_FCALL:
    DEBUG_PRINT("NODE_FCALL\n");
    char* fname = nd->e1.nd->e1.id;
    int num_arg;
    if( nd->e3.nd->type == NODE_FARG ){
      DEBUG_PRINT("function, %s \n", fname);
      tmp_code3 = gen_code(nd->e3.nd, table); /* Code for FARG */
      num_arg = count_num_farg( nd );
      DEBUG_PRINT("function, %s , with %d args \n", fname, num_arg);
      nd_code = gen_code_fcall(fname, num_arg, tmp_code3 );
      ADD_LOC(nd_code, nd->loc);
    }else if( nd->e3.nd->type == NODE_NULL ){
      nd_code = gen_code_fcall(fname, 0, NULL);
      ADD_LOC(nd_code, nd->loc);
    }else{
      DEBUG_PRINT("ERROR: Unintended node under NODE_FCALL\n");
    }
    return nd_code;
    break;

  case NODE_FARG:
    DEBUG_PRINT("NODE_FARG\n");
    tmp_code1 = gen_code(nd->e1.nd, table);
    if( nd->e3.sibling != NULL ){
      DEBUG_PRINT("Small sibling exists \n");
	  tmp_code3 = gen_code(nd->e3.sibling, table);
      nd_code = vm_inst_list_cat( tmp_code1, tmp_code3);
    }else{
      DEBUG_PRINT("No siblings exist \n");
      nd_code = tmp_code1;
    }
    return nd_code;
    break;

  case NODE_IF:
    DEBUG_PRINT("NODE_IF\n");
    tmp_code1 = gen_code(nd->e1.nd, table);
    nd_code = tmp_code1;

	char* label_L1 = new_vm_label();
	char* label_L2 = new_vm_label();

	// fjmp L1
	// Always generate code
	nd_code = vm_inst_list_cat( nd_code, new_vm_inst_fjmp(label_L1));

	// Code for "then node"
	if(nd->e2.nd != NULL){
	    tmp_code2 = gen_code (nd->e2.nd, table);
		nd_code = vm_inst_list_cat( nd_code, tmp_code2 );
	}

	// If e3 exists. jmp L2
	if(nd->e3.nd != NULL){ 
	    nd_code = vm_inst_list_cat( nd_code, new_vm_inst_jmp(label_L2) );
	}

	// label L1
	// Always generate 
	nd_code = vm_inst_list_cat( nd_code, new_vm_inst_label(label_L1));

	// If e3 exists
	// Code for "else node"
	if(nd->e3.nd != NULL){
	    tmp_code3 = gen_code (nd->e3.nd, table);
		nd_code = vm_inst_list_cat( nd_code, tmp_code3 );
	}

	// If e3 exists
	if(nd->e3.nd != NULL){ // If e3 exists.
	    nd_code = vm_inst_list_cat( nd_code, new_vm_inst_label(label_L2));
	}

	free_vm_label(label_L1);
	free_vm_label(label_L2);

	return nd_code;
    break;

  case NODE_NULL:
  default:
    DEBUG_PRINT("This part should not be executed.\n");

	return NULL;
    break;
  }
  DEBUG_FLUSH();
}

