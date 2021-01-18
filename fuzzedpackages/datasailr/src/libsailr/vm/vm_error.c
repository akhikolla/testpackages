#include "vm_error.h"

void
vm_error_raise(vm_stack* vmstack)
{
	vmstack->error = vmstack->error + 1;
}

unsigned int
vm_error_num(vm_stack* vmstack)
{
	return vmstack->error;
}

bool
vm_error_exist(vm_stack* vmstack)
{
	if( vm_error_num(vmstack) > 0 ){
		return true;
	}else{
		return false;
	}
}
