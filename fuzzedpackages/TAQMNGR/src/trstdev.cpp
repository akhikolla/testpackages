#include <taq.h>


void trstdev(vector<double> vec, vector<double> ones, int size, double &trMean, double &StdDev)
{
	double sum = 0, sumtr = 0, sumq = 0, ntr = 0;
	double var, media, mediaq, comV = 0, comO = 0, *comVp, *comOp;
	comVp = &comV;
	comOp = &comO;
	
	
	for(int k = 0; k < size; k++)
	{
		*comVp = vec[k];
		*comOp = ones[k];
		sumtr += (*comVp * *comOp);
		ntr += *comOp;
		sum += *comVp;
		sumq += (*comVp * *comVp);
	}
	
	media = sum/size;
	mediaq = sumq/size;
	trMean = sumtr/ntr;
	if(size == 1){
	  var = 0;
	}
	else
	{
	  var = (mediaq - (media * media)) * (size/(size - 1));
	}
	if(var < 0)
	{
		var = 0;
	}
	StdDev = sqrt(var);
}
