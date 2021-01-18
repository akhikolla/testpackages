#ifndef VM_CMD_H
#define VM_CMD_H

// VM commands.

#define VM_CMD_TABLE \
/* Stack manipulation */ \
	X(VM_PUSH_IVAL, "VM_PUSH_IVAL") /* For rvalue */ \
	X(VM_PUSH_DVAL, "VM_PUSH_DVAL")  /* For rvalue */ \
	X(VM_PUSH_PP_IVAL, "VM_PUSH_PP_IVAL") \
	X(VM_PUSH_PP_DVAL, "VM_PUSH_PP_DVAL") \
	X(VM_PUSH_PP_NUM, "VM_PUSH_PP_NUM") \
	X(VM_PUSH_PP_STR, "VM_PUSH_PP_STR") \
	X(VM_PUSH_PP_REXP, "VM_PUSH_PP_REXP") \
	X(VM_PUSH_NULL, "VM_PUSH_NULL") \
	X(VM_POP, "VM_POP") \
	X(VM_END, "VM_END") \
	X(VM_DISP, "VM_DISP") \
\
/* Code jump */ \
	X(VM_FJMP, "VM_FJMP") \
	X(VM_JMP, "VM_JMP") \
	X(VM_LABEL, "VM_LABEL") \
\
/* Value assign */ \
	X(VM_STO, "VM_STO") \
\
/* Calculations on stack */ \
/* Function Call */ \
	X(VM_FCALL, "VM_FCALL") \
\
/* Arithmetic Calculation */ \
	X(VM_ADDX, "VM_ADDX") \
	X(VM_SUBX, "VM_SUBX") \
	X(VM_MULX, "VM_MULX") \
	X(VM_DIVX, "VM_DIVX") \
	X(VM_MODX, "VM_MODX") \
	X(VM_POWX, "VM_POWX") \
	X(VM_FAC, "VM_FAC") \
	X(VM_UMINUS, "VM_UMINUS") \
\
/* Regular Expression Matching */ \
	X(VM_REXP_MATCH, "VM_REXP_MATCH") \
\
/* Logical Calculation */ \
	X(VM_AND, "VM_AND") \
	X(VM_OR, "VM_OR") \
	X(VM_EQ, "VM_EQ") \
	X(VM_NEQ, "VM_NEQ") \
	X(VM_GT, "VM_GT") \
	X(VM_LT, "VM_LT") \
	X(VM_GE, "VM_GE") \
	X(VM_LE, "VM_LE") \
	X(VM_NEG, "VM_NEG") \
\
/* No operation */ \
	X(VM_NOP, "VM_NOP") 

#define X(a, b) a,
enum _VM_CMD {
  VM_CMD_TABLE
};
#undef X
typedef enum _VM_CMD VM_CMD;

#define MAX_FUNC_NAME_LEN 511

char* vm_cmd_to_string(VM_CMD cmd);

#endif

