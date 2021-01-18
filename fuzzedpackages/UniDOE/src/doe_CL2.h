
void create_discrcl2(double **x,int nnew1, int np1,int nv1,CRITOPT *critopt);
void free_discrcl2();
double discrcl2_set(double **x0);
double discrcl2_cp_set(int ncol,int ncp,int *idx1,int *idx2);
double discrcl2_cp(int ncol,int ncp,int *idx1,int *idx2);
double discrcl2_cp1(int ncol,int i1,int i2);
double discrcl2_pm(int ncol, int npm,int *idx,int *idxp);
double discrcl2_pm_set(int ncol, int npm,int *idx,int *idxp);
double discrcl2();
double **discrcl2_x(double **xnew);
void discrcl2_snap(int ncol);
void discrcl2_reset(int ncol);
double discrcl2_eval(double **x);
void discrcl2_full_snap();
void discrcl2_full_reset();
void discrcl2_global_x(double **xnew);




