#include "doe_Matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


int GetTypeSize(int type);
void InitializeMatrix(MATRIX m, int StartRow,int EndRow);

static int MatrixErrorCode=0;

int CheckMatrixError(void)
{
return(MatrixErrorCode);
}
/*
MATRIX is a void-type pointer for a void-type pointer array, each
of which is a head pointer of a void-type array saving the matrix content.
The first three elements are used for:
row,col,and the type of each matrix element.

All types of MATRIX is initialized to 0
*/

int *NewIVector(int Num)
{
	int *m,i;
	if(Num<=0) return(0);
    m=(int *)malloc((unsigned) Num*sizeof(int));
	for(i=0;i<Num;i++) m[i]=0;
	return(m);
}

double *NewDVector(int Num)
{
	double *m;
	int i;
	if(Num<=0) return(0);
    m=(double *)malloc((unsigned) Num*sizeof(double));
	for(i=0;i<Num;i++) m[i]=0;
	return(m);
}

unsigned long *NewLVector(int Num) //changed here
{
	unsigned long *m;
	int i;
	if(Num<=0) return(0);
    m=(unsigned long *)malloc((unsigned) Num*sizeof(unsigned long));
	for(i=0;i<Num;i++) m[i]=0;
	return(m);
}

char *NewCVector(int Num)
{
	char *m;
	int i;
	if(Num<=0) return(0);
    m=(char *)malloc((unsigned) Num*sizeof(char));
	for(i=0;i<Num;i++) m[i]=0;
	return(m);
}

void *NewVector(int size)
{
	return(malloc(size));
}

void FreeVector(void *m)
{
	free(m);
}


MATRIX NewMatrix(int type,int NumRow,int NumCol)
{
  int i,ElementSize;
  void **m;
  ElementSize=GetTypeSize(type);
  if(NumRow<=0||NumCol<=0)
    {
    MatrixErrorCode=INVALID_ROWCOL;
    return(0);
    }
  if((ElementSize=GetTypeSize(type))==INVALID_TYPE)
    {
    MatrixErrorCode=INVALID_TYPE;
    return(0);
    }
  m=(void **)malloc((unsigned) (NumRow+3)*sizeof(void*));
  if (!m) {MatrixErrorCode=MEM_ERROR; return (0);}
  m[0]=(void *)NumRow;
  m[1]=(void *)NumCol;
  m[2]=(void *)type;
  m=&m[3];
  if (!m) {MatrixErrorCode=MEM_ERROR; return (0);}
  for(i=0;i<NumRow;i++)
        {
         m[i]=malloc((unsigned) NumCol*ElementSize);
         if(!m[i]) {MatrixErrorCode=MEM_ERROR; return(0);}
        }
  InitializeMatrix(m,0,NumRow);
  return(m);
}

/* get c-style array. All the elements are arranged sequenically
in c-style in the memory*/
void *Matrix2CArray(MATRIX m)
{
void *cm;
int row,col,type,ElementSize,i,j;
if(m==0) return(0);
row=(uintptr_t)m[-3];//row=(int)m[-3];
col=(uintptr_t)m[-2];//col=(int)m[-2];
type=(uintptr_t)m[-1];//type=(int)m[-1];
col=GetMatrixNumCol(m);
ElementSize=GetTypeSize(type);
cm=malloc((unsigned)row*col*ElementSize+3*sizeof(int));
if(!cm) {MatrixErrorCode=MEM_ERROR; return (0);}
((int *)cm)[0]=row;
((int *)cm)[1]=col;
((int *)cm)[2]=type;
cm=&(((int *)cm)[3]);
switch(type)
 {
  case typeDouble:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          *((double *)cm+i*col+j)=((double **)m)[i][j];
    break;
  case typeFloat:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          *((float *)cm+i*col+j)=((float **)m)[i][j];
    break;
  case typeInt:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          *((int *)cm+i*col+j)=((int **)m)[i][j];
    break;
  case typeChar:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          *((char *)cm+i*col+j)=((char **)m)[i][j];
  }
return(cm);
}

/* get Fortran-style array.
 All the elements are arranged sequenically	in Fortran-style
 in the memory*/
void *Matrix2FArray(MATRIX m)
{
void *cm;
int row,col,type,ElementSize,i,j;
if(m==0) return(0);
row=(uintptr_t)m[-3];//row=(int)m[-3];
col=(uintptr_t)m[-2];//col=(int)m[-2];
type=(uintptr_t)m[-1];//type=(int)m[-1];
ElementSize=GetTypeSize(type);
cm=malloc((unsigned)row*col*ElementSize+3*sizeof(int));
if(!cm) {MatrixErrorCode=MEM_ERROR; return (0);}
((int *)cm)[0]=row;
((int *)cm)[1]=col;
((int *)cm)[2]=type;
cm=&(((int *)cm)[3]);
switch(type)
 {
  case typeDouble:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          *((double *)cm+j*row+i)=((double **)m)[i][j];
    break;
  case typeFloat:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          *((float *)cm+j*row+i)=((float **)m)[i][j];
    break;
  case typeInt:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          *((int *)cm+j*row+i)=((int **)m)[i][j];
    break;
  case typeChar:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          *((char *)cm+j*row+i)=((char **)m)[i][j];
  }
return(cm);
}

/* change fortan-styled array to MATRIX.
*/
MATRIX FArray2Matrix(void *cm)
{
int row,col,type,i,j;
MATRIX m;
if(cm==0) return(0);
row=((int *)cm)[-3];
col=((int *)cm)[-2];
type=((int *)cm)[-1];
m=NewMatrix(type,row,col);
if(!m) {MatrixErrorCode=MEM_ERROR; return (0);}
switch(type)
 {
  case typeDouble:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((double **)m)[i][j]=*((double *)cm+j*row+i);
    break;
  case typeFloat:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((float **)m)[i][j]=*((float *)cm+j*row+i);
    break;
  case typeInt:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((int **)m)[i][j]=*((int *)cm+j*row+i);
    break;
  case typeChar:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((char **)m)[i][j]=*((char *)cm+j*row+i);
  }
return(m);
}

DMATRIX FArray2DMatrix2(double *cm,int row,int col)
{
	return((DMATRIX)FArray2Matrix2(cm,row,col,typeDouble));
}

MATRIX FArray2Matrix2(void *cm,int row,int col,int type)
{
int i,j;
MATRIX m;
if(cm==0) return(0);
m=NewMatrix(type,row,col);
if(!m) {MatrixErrorCode=MEM_ERROR; return (0);}
switch(type)
 {
  case typeDouble:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((double **)m)[i][j]=*((double *)cm+j*row+i);
    break;
  case typeFloat:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((float **)m)[i][j]=*((float *)cm+j*row+i);
    break;
  case typeInt:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((int **)m)[i][j]=*((int *)cm+j*row+i);
    break;
  case typeChar:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((char **)m)[i][j]=*((char *)cm+j*row+i);
  }
return(m);
}

/* change C-styled array to MATRIX.
*/
MATRIX CArray2Matrix(void *cm)
{
int row,col,type,/*ElementSize,*/i,j;
MATRIX m;
if(cm==0) return(0);
row=((int *)cm)[-3];
col=((int *)cm)[-2];
type=((int *)cm)[-1];
//ElementSize=GetTypeSize(type);
m=NewMatrix(type,row,col);
if(!m) {MatrixErrorCode=MEM_ERROR; return (0);}
switch(type)
 {
  case typeDouble:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((double **)m)[i][j]=*((double *)cm+i*col+j);
    break;
  case typeFloat:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((float **)m)[i][j]=*((float *)cm+i*col+j);
    break;
  case typeInt:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((int **)m)[i][j]=*((int *)cm+i*col+j);
    break;
  case typeChar:
    for(i=0;i<row;i++)
       for(j=0;j<col;j++)
          ((char **)m)[i][j]=*((char *)cm+i*col+j);
  }
return(m);
}


void FreeMatrix(MATRIX m)
{
int i,row;
if(m==0)return;
row=(uintptr_t)m[-3];//row=(int)m[-3];
for(i=row-1;i>=0;i--) free(m[i]);
free((void*)(&m[-3]));
}

void FreeArray(void *cm)
{
free((void*)(&(((int *)cm)[-3])));
}

void AddRowToMatrix(MATRIX *pm, int AddNumRow)
{
  int i,OldNumRow,NewNumRow,NumCol,type,ElementSize;
  MATRIX m;
  if(AddNumRow<=0) return;
  m=*pm;
  if(m==0) return;
  m=&m[-3];
  OldNumRow=(uintptr_t)m[0];//OldNumRow=(int)m[0];
  NumCol=(uintptr_t)m[1];//NumCol=(int)m[1];
  type = (uintptr_t)m[2];//type=(int)m[2];
  NewNumRow=OldNumRow+AddNumRow;
  m[0]=(void *) NewNumRow;
  ElementSize=GetTypeSize(type);
  m = (void **)realloc((void *)(&m[0]),
		  (unsigned) (NewNumRow+3)*sizeof(void*));
  if (!m)
     {MatrixErrorCode=MEM_ERROR;
      *pm=0;
      return;
     }
  m=&m[3];
  for(i=OldNumRow; i<NewNumRow; i++)
    {
     m[i] = malloc((unsigned) (NumCol*ElementSize));
     if (!m[i])
       {MatrixErrorCode=MEM_ERROR;
        *pm=0;
        return;
       }
    }
  InitializeMatrix(m,OldNumRow,NewNumRow-1);
  *pm=m;
}

int GetTypeSize(int type)
{
  switch (type)
    {
    case typeDouble:
      return(sizeof(double));
    case typeFloat:
      return(sizeof(float));
    case typeChar:
      return(sizeof(char));
    case typeInt:
      return(sizeof(int));
    default:
      return(INVALID_TYPE);
    }
}

void InitializeMatrix(MATRIX m,int StartRow,int EndRow)
/*Note: the first row is 0*/
{
int i,j,row,col,type;
row=(uintptr_t)m[-3];//row=(int)m[-3];
col=(uintptr_t)m[-2];//col=(int)m[-2];
type=(uintptr_t)m[-1];//type=(int)m[-1];
StartRow=StartRow>=0?StartRow:0;
EndRow=EndRow<row?EndRow:row-1;

switch(type)
 {
  case typeDouble:
    for(i=StartRow;i<=EndRow;i++)
      for(j=0;j<col;j++)
        ((double **)m)[i][j]=0;
    break;
  case typeFloat:
    for(i=StartRow;i<=EndRow;i++)
      for(j=0;j<col;j++)
        ((float **)m)[i][j]=(float)0;
    break;
  case typeInt:
    for(i=StartRow;i<=EndRow;i++)
      for(j=0;j<col;j++)
        ((int **)m)[i][j]=0;
    break;
  case typeChar:
    for(i=StartRow;i<=EndRow;i++)
      for(j=0;j<col;j++)
        ((char **)m)[i][j]=0;
  }
}

int GetMatrixNumRow(MATRIX m)
{
if(m==0) return 0;
return((uintptr_t)(m[-3]));//return((int)(m[-3]));
}

int GetMatrixNumCol(MATRIX m)
{
if(m==0) return 0;
return((uintptr_t)(m[-2]));//return((int)(m[-2]));
}

/*Double Matrix*/
DMATRIX NewDMatrix(int NumRow,int NumCol)
{
return((DMATRIX)(NewMatrix(typeDouble,NumRow,NumCol)));
}

void FreeDMatrix(DMATRIX m)
{
FreeMatrix((MATRIX)m);
}

void AddRowToDMatrix(DMATRIX *pm, int AddNumRow)
{
AddRowToMatrix((MATRIX *) pm,AddNumRow);
}

int GetDMatrixNumRow(DMATRIX m)
{
return(GetMatrixNumRow((MATRIX)(m)));
}

int GetDMatrixNumCol(DMATRIX m)
{
return(GetMatrixNumCol((MATRIX)(m)));
}

double * DMatrix2CArray(DMATRIX m)
{
return((double *)Matrix2CArray((MATRIX)m));
}

double * DMatrix2FArray(DMATRIX m)
{
return((double *)Matrix2FArray((MATRIX)m));
}

DMATRIX CArray2DMatrix(double *cm)
{
return((DMATRIX)CArray2Matrix((void *)cm));
}

DMATRIX FArray2DMatrix(double *cm)
{
return((DMATRIX)FArray2Matrix((void *)cm));
}

void FreeDArray(double *cm)
{
  FreeArray((void *)cm);
}

FMATRIX DMatrix2FMatrix(DMATRIX m)
{
int row,col,i,j;
FMATRIX nm;
row=GetDMatrixNumRow(m);
col=GetDMatrixNumCol(m);
nm=NewFMatrix(row,col);
for(i=0;i<row;i++)
  for(j=0;j<col;j++)
    nm[i][j]=(float)m[i][j];
return(nm);
}

/*int Matrix*/
IMATRIX NewIMatrix(int NumRow,int NumCol)
{
return((IMATRIX)NewMatrix(typeInt,NumRow,NumCol));
}

void FreeIMatrix(IMATRIX m)
{
FreeMatrix((MATRIX)m);
}

void AddRowToIMatrix(IMATRIX *pm, int AddNumRow)
{
AddRowToMatrix((MATRIX *) pm,AddNumRow);
}

int GetIMatrixNumRow(IMATRIX m)
{
return(GetMatrixNumRow((MATRIX)(m)));
}

int GetIMatrixNumCol(IMATRIX m)
{
return(GetMatrixNumCol((MATRIX)m));
}

/*char Matrix*/
CMATRIX NewCMatrix(int NumRow,int NumCol)
{
return((CMATRIX)NewMatrix(typeChar,NumRow,NumCol));
}

void FreeCMatrix(CMATRIX m)
{
FreeMatrix((MATRIX)m);
}

void AddRowToCMatrix(CMATRIX *pm, int AddNumRow)
{
AddRowToMatrix((MATRIX *) pm,AddNumRow);
}

int GetCMatrixNumRow(CMATRIX m)
{
return(GetMatrixNumRow((MATRIX)m));
}

int GetCMatrixNumCol(CMATRIX m)
{
return(GetMatrixNumCol((MATRIX)m));
}

/*float Matrix*/
FMATRIX NewFMatrix(int NumRow,int NumCol)
{
return((FMATRIX)NewMatrix(typeFloat,NumRow,NumCol));
}

void FreeFMatrix(FMATRIX m)
{
FreeMatrix((MATRIX)m);
}

void AddRowToFMatrix(FMATRIX *pm, int AddNumRow)
{
AddRowToMatrix((MATRIX *) pm,AddNumRow);
}

int GetFMatrixNumRow(FMATRIX m)
{
return(GetMatrixNumRow((MATRIX)(m)));
}

int GetFMatrixNumCol(FMATRIX m)
{
return(GetMatrixNumCol((MATRIX)(m)));
}

float * FMatrix2FArray(FMATRIX m)
{
return((float *)Matrix2FArray((MATRIX)m));
}

void FreeFArray(float *cm)
{
  FreeArray((void *)cm);
}
