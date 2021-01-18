#include "vm_cmd.h"

#define X(a, b) b,
char *vm_cmd_name[] = {
  VM_CMD_TABLE
};
#undef X

char* vm_cmd_to_string(VM_CMD cmd){
	return vm_cmd_name[cmd];
}

