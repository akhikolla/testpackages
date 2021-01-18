#include <stdlib.h>
#include <stdio.h>
#include "parser_state.h"
#include "var_hash.h"
#include "helper.h"

parser_state*
new_parser_state(const char* fname, ptr_table* table)
{
	parser_state* ps = (parser_state*)malloc(sizeof(parser_state));
	ps->fname = fname;
	ps->tree = NULL;
	ps->lineno = 0;
	ps->tline = 0;
	ps->yynerrs = 0;
	ps->ptrtable = table;
	ps->vars = var_hash_init();
	ps->lhsvars = var_hash_init();
	ps->rhsvars = var_hash_init();
	ps->rexp_encoding = SAILR_DEFAULT_REXP_ENCODING;
	return ps;
}

int
parser_state_free(parser_state* ps)
{
	// fname is constant. Should not be freed.
	// Tree is already freed.
	// Ptrtable is already freed.
	DEBUG_PRINT("Going to free vars\n");
	var_hash_free(&(ps->vars)); // Free vars
	DEBUG_PRINT("Going to free lhsvars\n");
	var_hash_free(&(ps->lhsvars)); // Free lhsvars
	DEBUG_PRINT("Going to free rhsvars\n");
	var_hash_free(&(ps->rhsvars)); // Free rhsvars
	DEBUG_PRINT("Going to free parser state object\n");
	free(ps); // Free parser_state
	return 1;
}
