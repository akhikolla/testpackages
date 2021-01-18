#include "doe_Index.h"
#include "doe_DesignUtil.h"

#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
void indexx(unsigned  int n, double arr[], unsigned  int indx[]);
void indexx_int(unsigned  int n, int arr[], unsigned  int indx[]);

/* array arr[0..n-1], index[0..n-1] (values from 0 to n-1)*/
void indexx1(unsigned  int n,double arr[], unsigned  int indx[])
{
	double *arr1;
	unsigned  int *indx1,i;

	arr1=arr-1;
	indx1=indx-1;
	indexx(n,arr1,indx1);
	for(i=0;i<n;i++) indx[i]--;
}

/*Indexes an array arr[1..n], i.e., outputs the array indx[1..n] such that arr[indx[j]] is
in ascending order for j = 1, 2, . . . ,N. The input quantities n and arr are not changed.*/
void indexx(unsigned  int n, double arr[], unsigned  int indx[])
{
	unsigned  int i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	double a;

	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=1;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[l]]) {
				SWAP(indx[l+1],indx[l])
			}
			i=l+1;
			j=ir;
			indxt=indx[l];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			//if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}


/* array arr[0..n-1], index[0..n-1] (values from 0 to n-1)*/
void indexx2(unsigned  int n,int arr[], unsigned  int indx[])
{
	int *arr1;
	unsigned  int *indx1,i;

	arr1=arr-1;
	indx1=indx-1;
	indexx_int(n,arr1,indx1);
	for(i=0;i<n;i++) indx[i]--;
}

/*Indexes an array arr[1..n], i.e., outputs the array indx[1..n] such that arr[indx[j]] is
in ascending order for j = 1, 2, . . . ,N. The input quantities n and arr are not changed.*/
void indexx_int(unsigned  int n, int arr[], unsigned  int indx[])
{
	unsigned  int i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	int a;

	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=1;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[l]]) {
				SWAP(indx[l+1],indx[l])
			}
			i=l+1;
			j=ir;
			indxt=indx[l];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			//if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
