typedef char ** CMATRIX;
typedef int ** IMATRIX;
typedef double ** DMATRIX;
typedef float ** FMATRIX;
typedef void ** MATRIX;
#define typeDouble 1
#define typeInt 2
#define typeChar 3
#define typeFloat 4

#define MEM_ERROR -3000
#define INVALID_TYPE -4000
#define INVALID_ROWCOL -5000

/* Define NULL pointer value */

#ifndef NULL
#define NULL    ((void *)0)
#endif

int CheckError(void);

double *NewDVector(int Num);
int *NewIVector(int Num);
unsigned long *NewLVector(int Num); //changed here
char *NewCVector(int Num);
void *NewVector(int size);
void 	FreeVector(void *vector);

MATRIX NewMatrix(int type,int NumRow,int NumCol);
void FreeMatrix(MATRIX m);
void AddRowToMatrix(MATRIX *pm,int AddNumRow);
int GetMatrixNumRow(MATRIX m);
int GetMatrixNumCol(MATRIX m);
void *Matrix2FArray(MATRIX m);
void *Matrix2CArray(MATRIX m);
MATRIX CArray2Matrix(void *cm);
MATRIX FArray2Matrix(void *cm);
MATRIX FArray2Matrix2(void *cm,int row,int col,int type);
void FreeArray(void *cm);

/*Double Matrix*/
DMATRIX NewDMatrix(int NumRow,int NumCol);
void FreeDMatrix(DMATRIX m);
void AddRowToDMatrix(DMATRIX *pm,int AddNumRow);
int GetDMatrixNumRow(DMATRIX m);
int GetDMatrixNumCol(DMATRIX m);
double * DMatrix2CArray(DMATRIX m);
double * DMatrix2FArray(DMATRIX m);
DMATRIX CArray2DMatrix(double *cm);
DMATRIX FArray2DMatrix(double *cm);
DMATRIX FArray2DMatrix2(double *cm,int row,int col);
void FreeDArray(double *cm);
FMATRIX DMatrix2FMatrix(DMATRIX m);


/*int Matrix*/
IMATRIX NewIMatrix(int NumRow,int NumCol);
void FreeIMatrix(IMATRIX m);
void AddRowToIIMatrix(IMATRIX *m,int AddNumRow);
int GetIMatrixNumRow(IMATRIX m);
int GetIMatrixNumCol(IMATRIX m);

/*char Matrix*/
CMATRIX NewCMatrix(int NumRow,int NumCol);
void FreeCMatrix(CMATRIX m);
void AddRowToICMatrix(CMATRIX *m,int AddNumRow);
int GetCMatrixNumRow(CMATRIX m);
int GetCMatrixNumCol(CMATRIX m);

/*float Matrix*/
FMATRIX NewFMatrix(int NumRow,int NumCol);
void FreeFMatrix(FMATRIX m);
void AddRowToIFMatrix(FMATRIX *m,int AddNumRow);
int GetFMatrixNumRow(FMATRIX m);
int GetFMatrixNumCol(FMATRIX m);
float * FMatrix2FArray(FMATRIX m);
void FreeFArray(float *cm);

