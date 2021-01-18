#include <R_ext/Print.h>
#include "node.h"
#include "ptr_table.h"
#include "helper.h"

int node_cnt = 0;

#define NEW_NODE_HEADER \
do { \
  node_cnt++; \
} \
while (0);


/* ----- Helper functions to convert char* to int/double. ----- */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <math.h>

int
char_to_int(char* num)
{
	char *ep;
	int ival; 
	long longval; 

	errno = 0; 
	longval = strtol(num, &ep, 10); 
	if (num[0] == '\0' || *ep != '\0') 
		Rprintf("ERROR: Not a valid int number. \n");
	if ((errno == ERANGE && (longval == LONG_MAX || longval == LONG_MIN)) || (longval > INT_MAX || longval < INT_MIN)) 
		Rprintf("ERROR: Invalid number out of int range. \n");
	ival = longval;
	return ival;
}

double
char_to_double(char* num)
{
	char* ep;
	double dval; 

	errno = 0; 
	dval = strtod(num, &ep); 
	if (num[0] == '\0' || *ep != '\0') 
		Rprintf("ERROR: Not a valid decimal number. \n");
	if (errno == ERANGE && (dval == DBL_MAX || dval == DBL_MIN)) 
		Rprintf("ERROR: Invalid number out of double range. \n");

	return dval;
}
/* ---------- */


TreeNode* new_node_prgm( TreeNode* tree_node )
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_PRGM;
  nd->e1.nd = tree_node;
  return nd;
}

TreeNode* new_node_stmt( TreeNode* tree_node )
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_STMT;
  nd->e1.nd = tree_node;
  nd->e2.last_sibling = nd; // point to self
  nd->e3.sibling  = (TreeNode*)NULL;
  return nd;
}

TreeNode* pushback_node_stmt( TreeNode* node_stmts , TreeNode* node_stmt2 )
{
  node_stmts->e2.last_sibling->e3.sibling = node_stmt2;
  node_stmts->e2.last_sibling = node_stmt2;
  return node_stmts;
}

TreeNode*
new_node_int( char* num )
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_INT;
  nd->e1.ival = char_to_int( num );
  return nd;
}

TreeNode*
new_node_double( char* num )
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_DBL;
  nd->e1.dval = char_to_double( num );
  return nd;
}

TreeNode*
new_node_nan_double ( )
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_DBL;
  nd->e1.dval = sqrt(-1) ;  // creating nan
  return nd;
}

TreeNode*
new_node_str( string_object* str , ptr_table* table )
{
  NEW_NODE_HEADER
  ptr_record* record = ptr_table_create_anonym_string(&table, &str);
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_STR;
  nd->e1.str_key = record->key ;
  // printf("table's pointer is %p \n", table);
  // ptr_table_show_all(&table);
  return nd;
}

TreeNode*
new_node_rexp( string_object* pattern , ptr_table* table , const char* rexp_encoding)
{
  NEW_NODE_HEADER
  ptr_record* record = ptr_table_create_anonym_rexp(&table, string_read(pattern), rexp_encoding );
  string_free(pattern);
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_REXP;
  nd->e1.rexp_key = record->key ;
  // printf("table's pointer is %p \n", table);
  // ptr_table_show_all(&table);
  return nd;
}

TreeNode*
new_node_ident( char* ident )
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_IDENT;
  nd->e1.id = ident;
  return nd;
}

TreeNode* new_node_fcall( TreeNode* ident_node , TreeNode* arg )
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_FCALL;
  nd->e1.nd = ident_node;
  nd->e3.nd = arg;
  return nd;
}

TreeNode* new_node_farg( TreeNode* arg)
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_FARG;
  nd->e1.nd = arg;
  nd->e3.sibling = (TreeNode*)NULL;
  return nd;
}

TreeNode* pushback_node_farg( TreeNode* first_arg, TreeNode* new_arg)
{
  TreeNode* node_ptr = first_arg->e3.sibling;
  if(node_ptr == NULL){
    first_arg->e3.sibling = new_arg;
  }else{
    while( node_ptr->e3.sibling != NULL ){
       node_ptr = node_ptr->e3.sibling;
    }
    node_ptr->e3.sibling = new_arg;
  }
  return first_arg;
}

TreeNode* new_node_op( char* op_type, TreeNode* t1, TreeNode* t2)
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_OP;
  nd->e1.op = op_type;
  nd->e2.nd = t1;
  nd->e3.nd = t2;
  return nd;
}

TreeNode* new_node_uniop( char* op_type, TreeNode* t1)
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_UNIOP;
  nd->e1.op = op_type;
  nd->e2.nd = t1;
  return nd;
}

TreeNode* new_node_let(TreeNode* lval, TreeNode* rval)
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_LET;
  nd->e1.nd = lval;
  nd->e2.nd = rval;
  return nd;
}

TreeNode* new_node_if(TreeNode* cond, TreeNode* then_node, TreeNode* else_node)
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_IF;
  nd->e1.nd = cond;
  nd->e2.nd = then_node;
  nd->e3.nd = else_node;
  return nd;
}

TreeNode* new_node_null()
{
  NEW_NODE_HEADER
  TreeNode* nd = (TreeNode*)malloc(sizeof(TreeNode));
  nd->type = NODE_NULL;
  nd->e1.nd = NULL;
  nd->e2.nd = NULL;
  nd->e3.nd = NULL;
  return nd;
}


int
count_num_farg(TreeNode* fcall_node)
{
  int count;
  if(fcall_node->e3.nd->type == NODE_NULL){
    // printf("NODE_NULL");
    return 0;
  }else if(fcall_node->e3.nd->type == NODE_FARG){
    // printf("NODE_FARG");
    count = 0;
  }else{
    return -1;
  }

  TreeNode* farg_node;
  farg_node = fcall_node->e3.nd;
  while( farg_node != NULL ){
    count = count + 1;
    farg_node = farg_node->e3.sibling;
  } 
  return count;
}



