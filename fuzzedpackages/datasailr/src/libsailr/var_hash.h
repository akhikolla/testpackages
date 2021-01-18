#ifndef VAR_HASH_H
#define VAR_HASH_H

#include "uthash.h"
#include "ptr_table.h"

struct _var_hash {
	char name[MAX_KEY_LEN];
	int dummy;
	UT_hash_handle hh; /* This macro makes this structure hashable */
};
typedef struct _var_hash var_hash ;

typedef var_hash var_elem;

var_hash* var_hash_init();
var_elem* var_hash_add_name ( var_hash**, char* );
var_elem* var_hash_insert (var_hash** , var_elem* );
var_elem* var_hash_find (var_hash** , char* );
int var_hash_size( var_hash** );
char** var_hash_names( var_hash** );
void var_hash_names_free( char** hash_names, int size );
void var_hash_free( var_hash** );

#endif // VAR_HASH_H
