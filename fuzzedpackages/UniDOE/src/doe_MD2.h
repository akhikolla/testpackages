void create_mxcl2(double **x,int nnew1, int np1,int nv1,CRITOPT *critopt);
void free_mxcl2();
double mxcl2_set(double **x0);
double mxcl2_cp_set(int ncol,int ncp,int *idx1,int *idx2);
double mxcl2_cp(int ncol,int ncp,int *idx1,int *idx2);
double mxcl2_cp1(int ncol,int i1,int i2);
double mxcl2_pm(int ncol, int npm,int *idx,int *idxp);
double mxcl2_pm_set(int ncol, int npm,int *idx,int *idxp);
double mxcl2();
double **mxcl2_x(double **xnew);
void mxcl2_snap(int ncol);
void mxcl2_reset(int ncol);
double mxcl2_eval(double **x);
void mxcl2_full_snap();
void mxcl2_full_reset();
void mxcl2_global_x(double **xnew);




