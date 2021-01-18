void create_wdl2(double **x,int nnew1, int np1,int nv1,CRITOPT *critopt);
void free_wdl2();
double wdl2_set(double **x0);
double wdl2_cp_set(int ncol,int ncp,int *idx1,int *idx2);
double wdl2_cp(int ncol,int ncp,int *idx1,int *idx2);
double wdl2_cp1(int ncol,int i1,int i2);
double wdl2_pm(int ncol, int npm,int *idx,int *idxp);
double wdl2_pm_set(int ncol, int npm,int *idx,int *idxp);
double wdl2();
double **wdl2_x(double **xnew);
void wdl2_snap(int ncol);
void wdl2_reset(int ncol);
double wdl2_eval(double **x);
void wdl2_full_snap();
void wdl2_full_reset();
void wdl2_global_x(double **xnew);




