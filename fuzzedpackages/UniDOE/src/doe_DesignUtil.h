#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void nrerror(const char* error_text);
float *vector();
float **matrix();
float **submatrix();
float **convert_matrix();
float ***f3tensor();
double *dvector(long,long);
double **dmatrix(long,long,long,long);
int *ivector(long,long);
int **imatrix(long,long,long,long);
unsigned char *cvector(long,long);
unsigned long *lvector(long,long);
void free_vector();
void free_dvector(double *,long,long);
void free_ivector(int*,long,long);
void free_cvector(unsigned char *,long,long);
void free_lvector(unsigned long *,long,long);
void free_matrix();
void free_submatrix();
void free_convert_matrix();
void free_dmatrix(double **,long,long,long,long);
void free_imatrix(int**,long,long,long,long);
void free_f3tensor();
