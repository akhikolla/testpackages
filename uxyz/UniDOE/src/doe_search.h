#include <vector>
typedef struct
{
char isnorm;
//char issymm;
char israndcol;
char israndpairs;
char isperm;
int maxiter; // max iteration without improvement
int maxcol; //max column exhanges per iteration
int maxpairs;
double maxtime;
double tol;
double th0;
double factor; //change factor of threshold
double hits_ratio;
int *levels;
double *colweight;
}SEARCHOPT;

void create_search(double **x,int nnew1,int np1,int nv1,SEARCHOPT *options);
void free_search();
void get_new_options(SEARCHOPT *options);
long get_ntotal();
std::vector<double> search(double **x);
