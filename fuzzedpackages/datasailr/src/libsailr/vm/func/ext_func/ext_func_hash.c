#include <R_ext/Print.h>
#include <stdio.h>
#include "ext_func_hash.h"

ext_func_hash*
ext_func_hash_init()
{
	ext_func_hash* hash = NULL;
	ext_func_elem* new_elem ;
	new_elem = (ext_func_elem*) malloc(sizeof(ext_func_elem));
	strncpy( new_elem->fname , "_HEAD_OF_EXT_FUNC_UTHASH_", MAX_FUNC_NAME_LEN - 1); 
	new_elem->num_args = 0;
	new_elem->last_executed = NULL;
	new_elem->func = NULL;

	ext_func_hash_insert(&hash, new_elem);
	return hash;
}

ext_func_elem*
ext_func_hash_find (ext_func_hash** hash, const char* fname)
{
	ext_func_elem* temp;
	HASH_FIND_STR(*hash, fname, temp);
	return temp;
}


ext_func_elem*
ext_func_hash_add (ext_func_hash** hash, const char* fname, unsigned int num_args, int (* func)(arg_list*, unsigned int, vm_stack*) )
{
	if (strlen(fname) >= MAX_FUNC_NAME_LEN - 1){
		Rprintf("ERROR: Function name is too long: %s \n", fname);
	}

	ext_func_elem* result = NULL;
	result = ext_func_hash_find(hash, fname);
	if(result == NULL){
		ext_func_elem* new_elem = (ext_func_elem*) malloc(sizeof(ext_func_elem));
		strncpy( new_elem->fname, fname, MAX_FUNC_NAME_LEN - 1) ;
		new_elem->num_args = num_args;
		new_elem->func = func;
		new_elem->last_executed = NULL; // Not used. Only used by the first elemnt.
		ext_func_hash_insert(hash, new_elem);
//		printf("%s added", new_elem->fname);
		return new_elem;
	}else{
		return result;
	}
}

const char*
ext_func_hash_get_last_executed(ext_func_hash** hash)
{
  ext_func_elem* first_elem = (ext_func_elem*)(*hash);
  return (first_elem->last_executed);
}

void
ext_func_hash_reset_last_executed(ext_func_hash** hash)
{
  ext_func_elem* first_elem = (ext_func_elem*)(*hash);
  first_elem->last_executed = NULL;
}

void
ext_func_hash_free( ext_func_hash** hash )
{
	ext_func_elem* elem;
	ext_func_elem* temp;

	/* When deleting elements of UTHASH, use both HASH_DEL() and free(). */
	/* Without HASH_DEL(), some memory leak happens*/
	HASH_ITER(hh, *hash, elem, temp) {
//		printf("Delete %s\n", elem->fname );
		HASH_DEL(*hash, elem);
		free(elem);  /* Free structure & memory */
	}
}

int
ext_func_elem_apply(ext_func_hash** hash, ext_func_elem* elem, vm_stack* vmstack)
{
	arg_list* arglis;
	if(elem->num_args > 0){
		arglis = arg_list_initialize( vmstack , elem->num_args );
	}else{
		arglis = NULL;
	}
	int result = (*(elem->func))( arglis , elem->num_args, vmstack);

	ext_func_hash_set_last_executed(hash, elem->fname);
	return result;
}

/* PRIVATE */

ext_func_elem*
ext_func_hash_insert (ext_func_hash** hash, ext_func_elem* new_elem )
{
	HASH_ADD_STR(*hash, fname, new_elem );
	return new_elem ;
}

void
ext_func_hash_set_last_executed(ext_func_hash** hash, const char* fname)
{
  ext_func_elem* first_elem = (ext_func_elem*)(*hash);
  first_elem->last_executed = fname;
}

