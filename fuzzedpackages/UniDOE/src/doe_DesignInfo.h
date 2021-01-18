
/* variable information*/
#ifndef DESIGNINFO_H
#define DESIGNINFO_H

typedef struct
{
	char *isbalanced;
	int *nxvl;  // number of value levels
	int **xid; //sorted index for x (column-wise)
	int **xvl; //different value level(not the same as group level)
	int **nle;  //number of entries for each value level
	//xid<->xvl
}XINFO;


typedef struct
{
	int *levels; //number of group level
	int *entries; //number of entries in each group level
	char *isequal; //is group level equal to value level
}XLEVEL;


extern XINFO xinfo; // changed from XINFO xinfo
extern XLEVEL xlevel; // changed from XLEVEL xlevel

//int TEST_INTEGER_VAR;

void create_xinfo(double **x,int nnew1,int np1,int nv1);
void create_xlevel(double **x,int *levels);
char checkSymm();
void free_xinfo();
void free_xlevel();
void xinfo_snap(int ncol);
void xinfo_reset(int ncol);
void xinfo_full_snap();
void xinfo_full_reset();

void xinfo_cp(int ncol,int idx1,int idx2);
void xinfo_pm(int ncol,int npm,int *idx1,int *idx2);
void xinfo_symm_pm(int ncol,int npm,int *idx1,int *idx2);
void xinfo_symm_cp(int ncol,int idx1,int idx2);

#endif
