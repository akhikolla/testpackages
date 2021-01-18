#ifndef NODE_H
#define NODE_H

#include <stdlib.h>
#include <string.h>
#include "ptr_table.h"
#include "string/common_string.h"
#include "script_loc.h"

extern int node_cnt;

/* node types */

typedef enum {
  NODE_PRGM,
  NODE_STMT,   // e1: <nd>, e3: <sibling>
  NODE_INT,    // e1: <ival>
  NODE_DBL,  // e1: <dval>
  NODE_STR,    // e1: <str_key> string
  NODE_REXP,   // e1: <rexp_key> pattern string
  NODE_IDENT,  // e1: <id> variable name
  NODE_FCALL,  // e1: <id> function name, e3: <nd> argument
  NODE_FARG,   // e1: <nd> node, e3: <sibling> next argument
  NODE_OP,     // e1: <id> op e2:<nd>left of operator, e3:<nd> right of operator
  NODE_UNIOP,  // e1: <id> op e2:<nd> argument
  NODE_LET,    // e1: lvalue, e2 rvalue
  NODE_IF,     // e1: condtion expression, e2: then stmt, e3: else stme
  NODE_NULL
} NodeType;


/* Basic node type */

typedef struct TreeNode_ {
  NodeType type;
  char* print_value;
  union {
    struct TreeNode_ *nd;
    char *id;
    char *op;
    char *str_key ;
    char *rexp_key ;
    int ival ; 
    double dval;
  } e1;
  union {
    struct TreeNode_ *nd;
	struct TreeNode_ *last_sibling;
  } e2;
  union {
    struct TreeNode_ *nd;
    struct TreeNode_ *sibling;
  } e3;

  struct script_loc loc; // hold the corresponding Sailr script location information 
} TreeNode ;

/* ********************************************* */
/* function prototypes for each NodeType        */
/* ********************************************* */

TreeNode* new_tree();
void init_tree(TreeNode*);

/* Create each type of node */

extern TreeNode* new_node_prgm( TreeNode* );
extern TreeNode* new_node_stmt( TreeNode* );
extern TreeNode* pushback_node_stmt( TreeNode*, TreeNode* );
extern TreeNode* new_node_int(char*);
extern TreeNode* new_node_double(char*);
extern TreeNode* new_node_nan_double();
extern TreeNode* new_node_str( string_object*, ptr_table* );
extern TreeNode* new_node_rexp( string_object*, ptr_table* , const char* rexp_encoding);
extern TreeNode* new_node_ident( char* );
extern TreeNode* new_node_fcall( TreeNode* , TreeNode*  );
extern TreeNode* new_node_farg( TreeNode* );
extern TreeNode* pushback_node_farg( TreeNode* , TreeNode* );
extern TreeNode* new_node_op( char* , TreeNode* , TreeNode* );
extern TreeNode* new_node_uniop( char* , TreeNode*);
extern TreeNode* new_node_let(TreeNode*, TreeNode*);
extern TreeNode* new_node_if(TreeNode* , TreeNode* , TreeNode* );
extern TreeNode* new_node_null();

extern int count_num_farg(TreeNode* fcall_node);

#endif /* NODE_H */


