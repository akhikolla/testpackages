#ifndef VM_CODE_H
#define VM_CODE_H

#include "vm_cmd.h"
#include "script_loc.h"

// VM codes consist of a command and an optional argument.

struct _vm_inst{
	VM_CMD cmd;
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)  /* C11 */
	union{ 
#endif
		int ival;
		double dval;
		char* ptr_key;
		char* label;
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L) 
	};
#endif
	char fname[MAX_FUNC_NAME_LEN + 1];
	int num_arg;
	struct _vm_inst* prev; // Only used for linked list structure
	struct _vm_inst* next; // Only used for linked list structure
	struct _vm_inst* last; // Only used for linked list structure. The first element should have a pointer to the last element for performance purpose.

	struct script_loc loc; // Location information of raw Sailr script. (Not essential for VM itself)
};
typedef struct _vm_inst vm_inst;

void vm_inst_show( vm_inst* );
int vm_inst_free( vm_inst* );

// vm_code
typedef vm_inst vm_code[];

int vm_code_jmp(vm_code code, int idx, char* label, int max_line);
// fjmp code is implemented at vm.c level.

#endif
