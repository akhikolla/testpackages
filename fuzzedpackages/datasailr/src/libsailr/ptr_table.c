#include <R_ext/Print.h>
#include "ptr_table.h"
#include <stdio.h>
#include <string.h>
#include "string/common_string.h"
#include "simple_re/simple_re.h"
#include "helper.h"


ptr_table*
ptr_table_init (){
	ptr_table *table = NULL;
	ptr_record* new_ptr_record;
	ptr_table_info* new_ptr_table_info;
	new_ptr_record = (ptr_record *)malloc(sizeof(ptr_record));
	strncpy( new_ptr_record->key , "_HEAD_OF_UTHASH_", MAX_KEY_LEN) ;
	new_ptr_table_info = (ptr_table_info *) malloc(sizeof(ptr_table_info));
	new_ptr_table_info->str_counter = 0;
	new_ptr_table_info->rexp_counter = 0;
//	new_ptr_table_info->null_updated = 0b0000 ; // GNU extension
	new_ptr_table_info->null_updated = 0x0 ;
	new_ptr_record->address = new_ptr_table_info;
	new_ptr_record->type = PTR_INFO;
	new_ptr_record->gc = GC_YES;
	new_ptr_record->ex_addr = (void*) NULL;
	new_ptr_record->ex_type = PTR_NULL;
	new_ptr_record->ex_gc = GC_NO;
	new_ptr_record->anonym = 0;
	ptr_table_insert(&table, new_ptr_record );
	return table;
}

// type specified should be PTR_INT, PTR_DBL, PTR_STR or PTR_NULL
ptr_record*
ptr_table_add (ptr_table** table, const char* key, void** address, PtrType type, GCReq gc )
{
    ptr_record * result = NULL;
    result = ptr_table_find (table, key);
    if(result == NULL){
        // Create new key/value.
        ptr_record* new_ptr_record;
        new_ptr_record = (ptr_record *)malloc(sizeof(ptr_record));
        if( (strlen(key) + 1 ) <= MAX_KEY_LEN ){
            strcpy( new_ptr_record->key , key);
        }else{
            memcpy( new_ptr_record->key , key, MAX_KEY_LEN );
            new_ptr_record->key[ MAX_KEY_LEN - 1] = '\0';
        }
        if(type != PTR_NULL){
            new_ptr_record->address = *address;
        }else{
            new_ptr_record->address = NULL;
        }
        new_ptr_record->type = type;
        new_ptr_record->gc = gc ;
        // If necessary, extra address and its related information should be updated.
        new_ptr_record->ex_addr = NULL;
        new_ptr_record->ex_type = PTR_NULL;
        new_ptr_record->ex_gc = GC_NO ; 
        new_ptr_record->anonym = 0;
        // Insert new ptr_record on ptr_table
        ptr_table_insert( table, new_ptr_record );
        result = new_ptr_record;
    } else {
        if(type != PTR_NULL){
            ptr_record_update( result, *address, type, gc );
        }else{
            ptr_record_update( result, NULL, type, gc );
        }
	}
	return result;
}

ptr_record*
ptr_table_create_int(ptr_table** table, const char* key, int ival)
{
	ptr_record* result = NULL;
	int* new = (int*)malloc(sizeof(int));
	*new = ival;
	result = ptr_table_add(table, key, (void**) &new, PTR_INT, GC_YES);
	return result;
}

ptr_record*
ptr_table_create_int_from_ptr(ptr_table** table, const char* key, int** iptr, double** dptr)
{
	ptr_record* result = NULL;
	result = ptr_table_add(table, key, (void**) iptr, PTR_INT, GC_NO);
	ptr_record_update_extra_address(result, (void**) dptr, PTR_DBL, GC_NO);
	return result;
}

int
ptr_table_update_int(ptr_table** table, const char* key, int ival)
{
	ptr_record* result = ptr_table_find(table, key);
	int* p_int = (int*) result->address;
	*p_int = ival;
	return 1;
}

string_object*
ptr_table_get_ptr_string(ptr_table** table, const char* key)
{
  ptr_record* result = ptr_table_find(table, key);
  void* ptr_address = result->address;
  return ptr_address ;
}

const char*
ptr_table_read_string(ptr_table** table, const char* key)
{
  ptr_record* result = ptr_table_find(table, key);
  void* ptr_address = result->address;
//  printf("The pointer referred is : %p", ptr_address);
  return string_read ((string_object*) ptr_address);
}

ptr_record*
ptr_table_create_double(ptr_table** table, const char* key, double dval)
{
	ptr_record* result = NULL;
	double* new = (double*)malloc(sizeof(double));
	*new = dval;
	result = ptr_table_add(table, key, (void**) &new, PTR_DBL, GC_YES);
	return result;
}

ptr_record*
ptr_table_create_double_from_ptr(ptr_table** table, const char* key, double** dptr, int** iptr)
{
	ptr_record* result = NULL;
	result = ptr_table_add(table, key, (void**) dptr, PTR_DBL, GC_NO);
	ptr_record_update_extra_address(result, (void**) iptr, PTR_INT, GC_NO);
	return result;
}

int
ptr_table_update_double(ptr_table** table, const char* key, double dval)
{
	ptr_record* result = ptr_table_find(table, key);
	double* p_dbl = (double*) result->address;
	*p_dbl = dval;
	return 1;
}

int
ptr_record_update_extra_address(ptr_record* pr, void** ptr_ex_addr, PtrType ex_type, GCReq ex_gc )
{
	pr->ex_addr = *ptr_ex_addr;
	pr->ex_type = ex_type;
	pr->ex_gc = ex_gc;
	return 1;
}

int
ptr_record_swap_addresses(ptr_record* pr)
{
	void* temp_address;
	PtrType temp_type;
	GCReq temp_gc;

	temp_address = pr->address;
	temp_type = pr->type;
	temp_gc = pr->gc;

	pr->address = pr->ex_addr;
	pr->type = pr->ex_type;
	pr->gc = pr->ex_gc;

	pr->ex_addr = temp_address;
	pr->ex_type = temp_type;
	pr->ex_gc = temp_gc;
	return 1;
}


char*
create_new_str_key(ptr_table** table){
	char* new_str = (char *)malloc(sizeof(char)* ANONYM_KEY_WIDTH );
	ptr_table_info* info = (ptr_table_info*) ((*table)->address) ;
	info->str_counter = (info->str_counter) + 1;
    const char* prefix = "STR%0*d";
	sprintf(new_str, prefix , ANONYM_KEY_WIDTH -3 -1 ,info->str_counter);
	return new_str ; 
}

char*
create_new_rexp_key(ptr_table** table){
	char* new_str = (char *)malloc(sizeof(char)* ANONYM_KEY_WIDTH );
	ptr_table_info* info = (ptr_table_info*) ((*table)->address) ;
	info->rexp_counter = (info->rexp_counter) + 1;
    const char* prefix = "REXP%0*d";
	sprintf(new_str, prefix , ANONYM_KEY_WIDTH -4 -1 , info->rexp_counter);
	return new_str ; 
}

void
ptr_record_set_anonym( ptr_record* pr, int val)
{
	pr->anonym = val;
}

int
ptr_record_get_anonym( ptr_record* pr)
{
	return pr->anonym;
}

ptr_record*
ptr_table_create_anonym_string(ptr_table** table, string_object** strptr)
{
	char* new_key;
	new_key = create_new_str_key(table);
	ptr_record* new_ptr_record;
	new_ptr_record = ptr_table_add(table, new_key, (void**)strptr, PTR_STR, GC_YES);
	ptr_record_set_anonym( new_ptr_record, 1);
	free(new_key);
	return new_ptr_record ;
}

int
ptr_record_reset_rexp(ptr_record* pr)
{
	if(pr->type == PTR_REXP){
		simple_re* rexp = (simple_re*) pr->address; 
		simple_re_reset( rexp );
		return 0;
	}else{
		DEBUG_PRINT( "This record is not regular expression ptr_record.\n");
		return -1;
	}
}

ptr_record*
ptr_table_create_string_from_cstring(ptr_table** table, const char* key, const char* str)
{
	ptr_record* new_ptr_record;
	string_object* p_str = string_new(str);
	new_ptr_record = ptr_table_add(table, key, (void**) &p_str, PTR_STR, GC_YES);
	return new_ptr_record ;

//	ptr_record* new_ptr_record;
//	string_object** pp_str = malloc(sizeof(string_object*));
//	pp_str* = string_new(str);
//	new_ptr_record = ptr_table_add(table, key, (void**) pp_str, PTR_STR, GC_YES);
//   free(pp_str);
//	return new_ptr_record ;
}

// (Deprecated)
//ptr_record*
//ptr_table_create_string_from_ptr(ptr_table** table, const char* key, string_object** strptr)
//{
//	ptr_record* result = ptr_table_create_string(table, key, strptr);
//	return result;
//}
//

int
ptr_table_update_string(ptr_table** table, const char* key, string_object** strptr)
{
    ptr_record* to_be_updated = ptr_table_find(table, key);
    if(to_be_updated->type != PTR_STR){
        Rprintf("ERROR: Record with non-string is trying to be updated with string.");
        return -1;
    }

    if(to_be_updated->gc == GC_YES){
        free(to_be_updated->address);
    }
    
    to_be_updated->address = *strptr;
    return 1;
}

int
ptr_record_update_string(ptr_record* pr , string_object** pp_str, GCReq gc)
{
	pr->address = *pp_str;
	pr->gc = gc;
    return 1;
}

ptr_record*
ptr_table_create_anonym_rexp(ptr_table** table, const char* pattern, const char* enc)
{
	char* new_key;
	new_key = create_new_rexp_key(table);
	simple_re* new_re;
	new_re = simple_re_compile( pattern, enc );
	ptr_record* new_ptr_record;
	new_ptr_record = ptr_table_add(table, new_key, (void**) &new_re, PTR_REXP, GC_YES);
	ptr_record_set_anonym( new_ptr_record, 1);
	free(new_key);
	return new_ptr_record ;
}

ptr_record*
ptr_table_create_null(ptr_table** table, const char* key)
{
	ptr_record* result = NULL;
    void** ppv = NULL;
	result = ptr_table_add(table, key, ppv, PTR_NULL, GC_NO);
	return result;
}

int
ptr_table_del_record(ptr_table** table, const char* key)
{
	ptr_record* to_be_deleted = ptr_table_find(table, key);
	if (to_be_deleted != NULL){
		HASH_DEL( *table, to_be_deleted);  /* Just removed from hash. (Structure exists) */
		ptr_record_free(to_be_deleted);  /* Free structure and momory pointed by address */
		return 1;
	} else {
		Rprintf("Cannot find record to be deleted.\n");
		return -1;
	}
}

void
ptr_record_free_gc_required_memory(ptr_record* pr)
{
//	ptr_record_show(pr);
	if(pr->gc == GC_YES){
		switch( pr->type ){	
			case PTR_INT:
			case PTR_DBL:
				free(pr->address);
				break;
			case PTR_STR:
				string_free((string_object*)pr->address);
				break;
			case PTR_REXP:
				simple_re_free((simple_re*)pr->address);
				break;
			case PTR_NULL:
				// Nothing to be freed
				break;
			case PTR_INFO:
				free((ptr_table_info*)pr->address);
				break;
			default:
				free(pr->address);
				break;
		}
		pr->address = NULL;
		pr->gc = GC_NO;
	}	
	if(pr->ex_gc == GC_YES){
		switch( pr->ex_type ){	
			case PTR_INT:
			case PTR_DBL:
				free(pr->ex_addr);
				break;
			case PTR_STR:
				string_free((string_object*)pr->ex_addr);
				break;
			case PTR_REXP:
				simple_re_free((simple_re*)pr->ex_addr);
				break;
			case PTR_NULL:
				break;
			case PTR_INFO:
				free((ptr_table_info*)pr->ex_addr);
				break;
			default:
				free(pr->ex_addr);
				break;
		}
		pr->ex_addr = NULL;
		pr->ex_gc = GC_NO;
	}
}

void
ptr_record_free_struct(ptr_record* pr)
{
	free(pr);
}

void
ptr_record_free(ptr_record* pr)
{
	ptr_record_free_gc_required_memory(pr);  /* Free memory pointed by address. */
	ptr_record_free_struct(pr); /* Free this record itself. */
}

PtrType
ptr_table_get_type(ptr_table** table, const char* key)
{
	ptr_record* temp = ptr_table_find(table, key);
	return temp->type;
}

int
ptr_record_is_ptr_null(ptr_table** table, const char* key)
{
	PtrType ptrtype = ptr_table_get_type(table, key);
	if(ptrtype == PTR_NULL){
		return 1;
	} else {
		return 0;
	}
}

void**
ptr_table_get_pptr(ptr_table** table, const char* key)
{
	ptr_record* temp = ptr_table_find(table, key);
	void** temp_pptr = (void**) &(temp->address);
	return temp_pptr;
}

ptr_record*
ptr_table_first_record(ptr_table** table)
{
    ptr_record *pr;
    pr = *table;
	return pr;
}

ptr_record*
ptr_record_next(ptr_record* pr){
	return pr->hh.next;
}

void
ptr_table_del_records_except(ptr_table** table, const char** keys, int key_num )
{
	/* keys is array of pointers to chars. */
	ptr_record *current_record;
	ptr_record *temp_record;
	char* current_record_key;
	const char* key_name;
	int idx;
	int matched;

	for( idx = 0; idx < key_num ; ++idx ){
		key_name = keys[idx];
		Rprintf("* %s\n", key_name);
	}

	for(current_record = *table; current_record != NULL; current_record=temp_record) {
		current_record_key = current_record->key;
		temp_record = current_record->hh.next;
		matched = 0;
		for( idx = 0; idx < key_num ; ++idx ){
			key_name = keys[idx];
			if(strcmp( current_record_key , key_name ) == 0){ /*matched*/
				matched = 1;
			}
		}
		if(matched == 1){
			/* Don't delete the curent record. */
		}else{
			/* Delete the current record. */
    		HASH_DEL( *table , current_record);  /* Just removed from hash. (Structure exists) */
			ptr_record_free( current_record );  /* Free structure and momory pointed by address */
		}
	}
}

void
ptr_table_del_all(ptr_table** table)
{
  ptr_record *current_record;
  ptr_record *temp_record;

  for(current_record = *table; current_record != NULL; current_record=temp_record) {
    temp_record = current_record->hh.next;
    HASH_DEL( *table , current_record);  /* Just removed from hash. (Structure exists) */
	ptr_record_free( current_record );  /* Free structure and momory pointed by address */
  }
}

void
ptr_table_show_all(ptr_table** table)
{
    ptr_record *pr;
    for( pr = *table ; pr != NULL; pr=(ptr_record *)pr->hh.next) {
		ptr_record_show(pr);
    }
}

void
ptr_record_show(ptr_record* pr)
{
		if(pr->type == PTR_INT){
			if( pr->address != NULL ){
	        	Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t VAL:%d\t (EXTR_ADR:%p\t TYPE:%d\t GC:%d\t VAL:%lf) [Anonym:%d]\n", 
				pr->key, pr->address, pr->type, pr->gc, *((int*)(pr->address)),
				pr->ex_addr, pr->ex_type, pr->ex_gc, *((double*)(pr->ex_addr)), pr->anonym );
			}else{
	        	Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t VAL:(NULL)\t (EXTR_ADR:%p\t TYPE:%d\t GC:%d\t VAL:(NULL)) [Anonym:%d]\n", 
				pr->key, pr->address, pr->type, pr->gc,
				pr->ex_addr, pr->ex_type, pr->ex_gc, pr->anonym );
			}
		}else if(pr->type == PTR_DBL){
			if( pr->address != NULL ){
	        	Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t VAL:%lf\t (EXTR_ADR:%p\t TYPE:%d\t GC:%d\t VAL:%d) [Anonym:%d]\n", 
				pr->key, pr->address, pr->type, pr->gc, *((double*)(pr->address)),
				pr->ex_addr, pr->ex_type, pr->ex_gc, *((int*)(pr->ex_addr)), pr->anonym);
			}else{
	        	Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t VAL:(NULL)\t (EXTR_ADR:%p\t TYPE:%d\t GC:%d\t VAL:(NULL)) [Anonym:%d]\n", 
				pr->key, pr->address, pr->type, pr->gc, 
				pr->ex_addr, pr->ex_type, pr->ex_gc, pr->anonym );
			}
		}else if(pr->type == PTR_STR){
			if( pr->address != NULL ){
				Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t VAL:%s\t (EXTR_ADR:%p (Not used for string)) [Anonym:%d]\n", 
				pr->key, pr->address, pr->type, pr->gc, string_read((string_object*)(pr->address)),
				pr->ex_addr, pr->anonym );
			}else{
				Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t VAL:(NULL)\t (EXTR_ADR:%p (Not used for string)) [Anonym:%d]\n", 
				pr->key, pr->address, pr->type, pr->gc, pr->ex_addr, pr->anonym );
			}
		}else if(pr->type == PTR_REXP){
			if( pr->address != NULL ){
        		Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t VAL:%s\t (EXTR_ADR:%p) [Anonym:%d]\n", 
				pr->key, pr->address, pr->type, pr->gc, simple_re_read_pattern((simple_re*)(pr->address)),
				pr->ex_addr, pr->anonym );
			}else{
        		Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t VAL:(NULL)\t (EXTR_ADR:%p) [Anonym:%d]\n", 
				pr->key, pr->address, pr->type, pr->gc, 
				pr->ex_addr, pr->anonym );
			}
		}else if(pr->type == PTR_NULL){
        	Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t (EXTR_ADR:%p) [Anonym:%d]\n", 
			pr->key, pr->address, pr->type, pr->gc, pr->ex_addr, pr->anonym);
		}else if(pr->type == PTR_INFO){
        	Rprintf("KEY:%s\t ADR:%p\t TYPE:%s\t GC:%d\t (EXTR_ADR:%p)", 
			pr->key, pr->address, "INFO" , pr->gc, pr->ex_addr);
			ptr_table_info* pti = (ptr_table_info*) (pr->address);
			Rprintf("\t str_counter %d, rexp_counter %d, null_updated %d \n", pti->str_counter, pti->rexp_counter, pti->null_updated);
		}else{
        	Rprintf("KEY:%s\t ADR:%p\t TYPE:%d\t GC:%d\t (EXTR_ADR:%p) [Anonym:%d]\n", 
			pr->key, pr->address, pr->type, pr->gc,	pr->ex_addr, pr->anonym);
		}
}


ptr_table*
ptr_record_obtain_table(ptr_record* pr)
{
	ptr_record* temp ;
	ptr_record* pre_temp = NULL;
	ptr_table** table;

    for( temp = pr ; temp != NULL; temp = (ptr_record*) (temp->hh.prev)) {
		pre_temp = temp;
    }
	table = &pre_temp;
	if( ptr_table_points_to_header(table) ){
		return *table;
	}else{
		Rprintf("ERROR: The function cannot find header of UTHASH.\n");
		return NULL;
	}
}

/*
int
ptr_table_info_set_null_updated(ptr_table** table, int updated_value)
{
	ptr_record* pr;
	if(ptr_table_points_to_header(table)){
		pr = (ptr_record*) *table;
		((ptr_table_info*) (pr->address))->null_updated = updated_value;
		return 1;
	}else{
		Rprintf("ERROR: The pointer passed is not pointing to valid ptr_table.");
		return 0;
	}
}
*/

int
ptr_table_info_change_null_updated_by_type(ptr_table** table, PtrType type)
{
	ptr_record* pr;
	unsigned int current_null_updated;
	unsigned int bitmask;
	if(ptr_table_points_to_header(table)){
		pr = (ptr_record*) *table;
		current_null_updated = ((ptr_table_info*) (pr->address))->null_updated ;
		
		if( PTR_INT <= type && type <= PTR_REXP ){
//			bitmask = ( 0b0001 << ((int)type));  // GNU extension
			bitmask = ( 0x1 << ((int)type));
			DEBUG_PRINT( "new type is %d, new bitmask is %d \n", type, bitmask);
			DEBUG_PRINT( "current null_update value is %d, new value is %d \n", current_null_updated, (current_null_updated | bitmask));
		}else{
			bitmask = 0 ; // This branch should never be executed.
			Rprintf("ERROR: Null may be converted to unintentional type on ptr_table." );
		}
		((ptr_table_info*) (pr->address))->null_updated = current_null_updated | bitmask ;
		return 1;
	}else{
		Rprintf("ERROR: The pointer passed is not pointing to valid ptr_table.");
		return 0;
	}
}

int
ptr_table_info_get_null_updated(ptr_table** table)
{
	ptr_record* pr;
	if(ptr_table_points_to_header(table)){
		pr = (ptr_record*) *table;
		ptr_table_info* pti = (ptr_table_info*) (pr->address);
		return pti->null_updated;
	}else{
		Rprintf("ERROR: The pointer passed is not pointing to valid ptr_table. This branch works, but should never executed. ");
		return 0; 
	}
}

int
ptr_table_info_reset_null_updated(ptr_table** table)
{
	ptr_record* pr;
	if(ptr_table_points_to_header(table)){
		pr = (ptr_record*) *table;
//		((ptr_table_info*) (pr->address))->null_updated = 0b0000 ;  // GNU extension
		((ptr_table_info*) (pr->address))->null_updated = 0x0 ;
		return 1;
	}else{
		Rprintf("ERROR: The pointer passed is not pointing to valid ptr_table.");
		return 0;
	}
}


// private
ptr_record*
ptr_table_find(ptr_table** table, const char* key)
{
	ptr_record* temp;
	HASH_FIND_STR(*table, key, temp);
	return temp;
}

ptr_record*
ptr_table_insert(ptr_table** table, ptr_record* new_ptr_record)
{
	HASH_ADD_STR(*table, key , new_ptr_record );
	return new_ptr_record;
}

int
ptr_record_update(ptr_record* pr, void* address, PtrType type, GCReq gc )
{
	pr->address = address;
	pr->type = type;
	pr->gc = gc;
	return 1;
}

bool
ptr_table_points_to_header(ptr_table** table)
{
	ptr_record* pr = (ptr_record*) (*table);
	if( (pr->type == PTR_INFO) &&  (strcmp(pr->key , "_HEAD_OF_UTHASH_") == 0) ){
		return true;
	}else{
		return false;
	}
}

/*
int main(int argc, char** argv){
	struct _person {
		char* name;
		int age;
	};
	typedef struct _person person;

	person* ind1 = (person*)malloc(sizeof(person));
	ind1->name = "toshi";
	ind1->age = 35;
	struct_string* hello = new_struct_string( "Hello", strlen("Hello"));

	ptr_table* table;
	table = ptr_table_init();
	// printf("table\t %p \n");
	
	ptr_table_add(&table, ind1->name, (void**)&ind1, PTR_OBJ, GC_YES);
	Rprintf("toshi\t %p \n", ind1);
	// printf("table\t %p \n");

	ptr_table_add_anonym_str(&table, &hello);
	Rprintf("hello\t %p \n", hello);
	// printf("table\t %p \n");

	Rprintf("\n");
	ptr_table_show_all(&table);

	ptr_table_del_record(&table, "toshi");
	Rprintf("\n");
	ptr_table_show_all(&table);

	ptr_table_del_record(&table, "STR000000000001");
	Rprintf("\n");
	ptr_table_show_all(&table);	
}
*/



