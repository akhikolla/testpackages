#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define SWAP(a,b,itemp) (itemp)=(a);(a)=(b);(b)=(itemp)
#define ABS(x) ((x)>0?(x):(-(x)))
#define MINIDOUBLE 1.0e-200
#define MAXDOUBLE 1.0e200
#define EPS 1.0e-10
#define EPS2 1.0e-15

double **normalize(int type,double **x,int nsamp,int nv);//,int *level);
void unnormalize(double **x,double **rang,int nsamp,int nv);
int mymin(int *v,int n);
int iCheckValue(int min,int max,int def,int value);
unsigned int uCheckValue(unsigned int min,unsigned int max,unsigned int def,unsigned int value);
char cCheckValue(char min,char max,char def,char value);
double dCheckValue(double min,double max,double def,double value);
void permute(int *x,int n);
double mult(double x,int p);


