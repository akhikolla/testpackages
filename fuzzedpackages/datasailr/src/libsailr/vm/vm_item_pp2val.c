#include <R_ext/Print.h>
#include "simple_re/simple_re.h"
#include "string/common_string.h"
#include "vm_item_pp2val.h"
#include <stdio.h>

#include "ptr_table.h"
#include "helper.h"

void
stack_item_pp2value(stack_item* item)
{
//	int result = 1;

	ptr_record* null_record ;
	switch(item->type){
		case PP_IVAL:
			item->type = IVAL;
			int tmp_int = **(item->pp_ival) ;
			item->ival = tmp_int;
			item->p_record = TEMP_OBJECT;
			DEBUG_PRINT("PP_IVAL on stack is converted to IVAL.\n");
			break;
		case PP_DVAL:
			item->type = DVAL;
			double tmp_double = **(item->pp_dval) ;
			item->dval = tmp_double;
			item->p_record = NOT_ON_PTR_TABLE;
			DEBUG_PRINT("PP_DVAL on stack is converted to DVAL.\n");
			break;
		case IVAL:
			DEBUG_PRINT("IVAL on stack is already a value. ( The item itself holds a value, %d )\n", item->ival );
			break;
		case DVAL:
			DEBUG_PRINT("DVAL on stack is already a value. ( The item itself holds a value, %f )\n", item->dval );
			break;
		case NULL_ITEM:
			null_record = (ptr_record*) item->p_record;
			if(null_record->type == PTR_NULL){
				Rprintf("ERROR: The variable, %s, should not be null. ", null_record->key );
				Rprintf("Variable of null value cannot be used for calculation. \n");
//				result = 0;
			}else{
				if(null_record->type == PTR_INT){
					item->type = IVAL;
					int tmp_int = *((int*)(null_record->address)) ;
					item->ival = tmp_int;
					item->p_record = null_record;
					DEBUG_PRINT("NULL_ITEM is converted to IVAL.\n");
				}else if(null_record->type == PTR_DBL){
					item->type = DVAL;
					double tmp_double = *((double*)(null_record->address)) ;
					item->dval = tmp_double;
					item->p_record = null_record;
					DEBUG_PRINT("NULL_ITEM is converted to DVAL.\n");
				}else if(null_record->type == PTR_STR){
					item->type = PP_STR;
					string_object** tmp_pp_str = (string_object**) &( null_record->address) ;
					item->pp_str = tmp_pp_str;
					item->p_record = null_record;
					DEBUG_PRINT("NULL_ITEM is converted to PP_STR.\n");
				}else if(null_record->type == PTR_REXP){
					item->type = PP_REXP;
					simple_re** tmp_pp_rexp = (simple_re**) &(null_record->address);
					item->pp_rexp = tmp_pp_rexp;
					item->p_record = null_record;
					DEBUG_PRINT("NULL_ITEM is converted to PP_REXP.\n");
				}else{
					Rprintf("ERROR: NULL_ITEM points to a ptr_record with unintended type: %s", null_record->key );
//					result = 0;
				}
			}
			break;
		default:
			DEBUG_PRINT("This item is not PP_IVAL, PP_DVAL, IVAL, DVAL or NULL_ITEM. Need not to be converted. \n");
//			result = 1;
	}
//	return result;
}

