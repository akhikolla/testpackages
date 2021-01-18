#include <stdio.h>
#include <stdlib.h>
#include "tree_free.h"
#include "tree_dump.h"
#include "node.h"
#include "helper.h"

void
tree_free( TreeNode* nd, int level)
{
  if(nd == NULL){
    DEBUG_PRINT( "Pointer to TreeNode is NULL (nothing to free)\n" );
    return;
  }

  switch(nd->type){
  case NODE_PRGM:
    tree_free(nd->e1.nd, level + 1);
    DEBUG_PRINT( "%s%s\n", rep_spaces(level), "Free NODE_PRGM" );
	free(nd);
    break;
  case NODE_STMT:
    if( nd->e3.sibling != NULL)
      tree_free(nd->e3.sibling, level + 1);
    if( nd->e1.nd != NULL )
      tree_free(nd->e1.nd, level + 1);
    DEBUG_PRINT( "%s%s\n", rep_spaces(level), "Free NODE_STMT" );
    free(nd);
    break;
  case NODE_INT:
    DEBUG_PRINT( "%s%s(%d)\n", rep_spaces(level), "Free NODE_INT", nd->e1.ival );
    free(nd);
    break;
  case NODE_DBL:
    DEBUG_PRINT( "%s%s(%f)\n", rep_spaces(level), "Free NODE_FLOAT", nd->e1.dval );
    free(nd);
    break;
  case NODE_STR:
    DEBUG_PRINT( "%s%s(%s)\n", rep_spaces(level), "Free NODE_STR", nd->e1.str_key );
    free(nd);
    break;
  case NODE_REXP:
    DEBUG_PRINT( "%s%s(%s)\n", rep_spaces(level), "Free NODE_REXP", nd->e1.rexp_key );
	free(nd);
    break;
  case NODE_IDENT:
    DEBUG_PRINT("%s%s(%s)\n", rep_spaces(level),"Free NODE_IDENT", nd->e1.id);
	free(nd->e1.id);
    free(nd);
    break;
  case NODE_FCALL:
    if( nd->e3.nd != NULL ){
      tree_free(nd->e3.nd, level + 1);
	}
    tree_free(nd->e1.nd, level);
    DEBUG_PRINT("%s%s(%s)\n", rep_spaces(level), "Free NODE_FCALL",nd->e1.nd->e1.id);
    free(nd);
    break;
  case NODE_FARG:
    if( nd->e3.sibling != NULL )
      tree_free(nd->e3.sibling, level );
    if( nd->e1.nd != NULL )
      tree_free(nd->e1.nd, level + 1 );
    DEBUG_PRINT( "%s%s\n", rep_spaces(level), "Free NODE_FARG" );
    free(nd);
    break;
  case NODE_OP:
    if( nd->e3.nd != NULL )
      tree_free(nd->e3.nd, level + 1);
    if( nd->e2.nd != NULL )
      tree_free(nd->e2.nd, level + 1);
    DEBUG_PRINT( "%s%s(%s)\n", rep_spaces(level), "Free NODE_OP", nd->e1.op );
    free(nd);
    break;
  case NODE_UNIOP:
    if( nd->e2.nd != NULL )
      tree_free(nd->e2.nd, level + 1);
    DEBUG_PRINT( "%s%s(%s)\n", rep_spaces(level), "Free NODE_UNIOP", nd->e1.op );
    free(nd);
    break;
  case NODE_LET:
    if( nd->e2.nd != NULL)
      tree_free(nd->e2.nd, level + 1);
    if( nd->e1.nd != NULL)
      tree_free(nd->e1.nd, level + 1);
    DEBUG_PRINT( "%s%s\n", rep_spaces(level), "Free NODE_LET" );
    free(nd);
    break;
  case NODE_IF:
    if( nd->e3.nd != NULL)
      tree_free(nd->e3.nd, level + 1);
    if( nd->e2.nd != NULL)
      tree_free(nd->e2.nd, level + 1);
    if( nd->e1.nd != NULL)
      tree_free(nd->e1.nd, level + 1);
    DEBUG_PRINT( "%s%s\n", rep_spaces(level), "Free NODE_IF" );
    free(nd);
    break;
  case NODE_NULL:
    DEBUG_PRINT( "%s%s\n", rep_spaces(level), "Free NODE_NULL" );
    free(nd);
    break;
  }
}

