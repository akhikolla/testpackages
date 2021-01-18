typedef struct
{
	char type;
	char ismax;
	double *scale;
	double goal;
	int npars[3];
	double *pars[3];
	char func[50];
}CRITOPT;
void create_criteria(double **x,int nnew,int np,int nv,CRITOPT *critopt);
void free_criteria(void);
double criteria();
double **criteria_x(double **x);
double criteria_set(double **x);
double criteria_cp_set(int ncol,int ncp,int *idx1,int *idx2);
double criteria_cp(int ncol,int ncp,int *idx1,int *idx2);
double criteria_cp1(int ncol,int idx1,int idx2);
double criteria_pm(int ncol, int npm,int *idx,int *idxp);
double criteria_pm_set(int ncol, int npm,int *idx,int *idxp);
void criteria_snap(int ncol);
void criteria_reset(int ncol);
double criteria_eval(double **x);
double criteria_min();
void criteria_full_snap();
void criteria_full_reset();
void criteria_global_x(double **x);
