#include<cstdlib>
#include<algorithm>
#include<cmath>
#include<string>
#include<ctime>
#include<Rcpp.h>

using namespace Rcpp;

double Summation(const NumericVector& x, const NumericVector& sigma, const int& LowI, const int& UppI)
{
	// LowI should be by default equal to 0
	// UppI should be by default equal to SampleSize - 1
	double SumRes = x[LowI]/(sigma[LowI]*sigma[LowI]);
	double SumSigma = 1/(sigma[LowI]*sigma[LowI]);
	for (int i = LowI + 1; i <= UppI; i++)
	{
		SumRes += x[i]/(sigma[i]*sigma[i]);
		SumSigma += 1/(sigma[i]*sigma[i]);
	}
	return SumRes/SumSigma;
}

double LogLikelihood(const NumericVector& y, const NumericVector& sigma, const int& LowI, const int& UppI) // Calculate the LogLiklihood for points with indices between LowI and UppI
{
	// Calculate the average
	double MeanY = Summation(y, sigma, LowI, UppI); // division by the number of elements is included
	double LogLikhood = 0;
	for (int i = LowI; i <= UppI; i++)
	{
		LogLikhood += (y[i] - MeanY)*(y[i] - MeanY) / (sigma[i] * sigma[i]);
	}
	return LogLikhood;
}

void BinaryConfig(unsigned long long int c, int*& Config, int& l, const int& Shift, const int& start)
{
	unsigned long long int residu = c;
	int counter = 0; l = 0;
	while (residu > 1)
	{
		if (residu % 2 == 1)
		{
			Config[l+start] = counter+Shift;
			l++;
		}

		residu /= 2;
		counter++;
	}
	if (residu == 1)
	{
		Config[l+start] = counter+Shift;
		l++;
	}
}
void RankUpdate(IntegerVector& Lower, IntegerVector& Upper, const int* InqPosi, const int& l, const int& n)
{
	for (int i = 0; i <= InqPosi[0]; i++)
	{
		Lower[i] = 0;
		if (InqPosi[0] > Upper[i]) {
			Upper[i] = InqPosi[0];
		}
	}
	int j = 0;
	
	while (j <= (l - 2))
	{
		for (int i = InqPosi[j] + 1; i <= InqPosi[j + 1]; i++)
		{
			if (InqPosi[j] + 1 < Lower[i]) {
				Lower[i] = InqPosi[j] + 1;
			}
			if (InqPosi[j + 1] > Upper[i]) {
				Upper[i] = InqPosi[j + 1];
			}
		}
		j++;
	}
	for (int i = InqPosi[l - 1] + 1; i < n; i++)
	{
		if (InqPosi[l - 1] + 1 < Lower[i]) {
			Lower[i] = InqPosi[l - 1] + 1;
		}
		Upper[i] = n - 1;
	}
}
unsigned long long int binomialCoeff(int n, int k)
{
	if (k>n) return 0;
	unsigned long long int res = 1;

	// Since C(n, k) = C(n, n-k)
	if (k > n - k){
		k = n - k;
	}

	// Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
	for (int i = 0; i < k; ++i)
	{
		res *= (n - i);
		res /= (i + 1);
	}

	return res;
}
bool PAVACheck(const NumericVector& y_temp, const NumericVector& sigma_temp, const int& l, const int* InqPosi, const int& n)
{
	if(l==0) return false;
	double MeanBlock1 = Summation(y_temp,sigma_temp,0,InqPosi[0]);
	double MeanBlock2;
	bool CheckPAVA = false;
	int k = 0;
	while (k <= (l - 2))
	{
		MeanBlock2 = Summation(y_temp,sigma_temp,InqPosi[k] + 1,InqPosi[k + 1]);
		if(MeanBlock1>MeanBlock2)
		{
			CheckPAVA = true;
		}
		else
		{
			MeanBlock1 = MeanBlock2;
		}
		k++;
	}
	if(CheckPAVA == false)
	{
		MeanBlock2 = Summation(y_temp,sigma_temp,InqPosi[l - 1] + 1,n-1);
		if(MeanBlock1>MeanBlock2)
		{
			CheckPAVA = true;
		}
	}
	return CheckPAVA;
}
void CorrectPermutationsF(const NumericVector& y, const NumericVector& sigma, const NumericVector& crit, IntegerVector& Lower, IntegerVector& Upper, const int* InqPosi, const int& l, const int& n, const bool& EqSigma)
{
	NumericVector y_temp(n);
	NumericVector sigma_temp(n);
	IntegerVector Lower_temp(n);
	IntegerVector Upper_temp(n);
	
	double LR=0;
	int k = 0;
	bool CheckPAVA = false;
	for(int i=1; i<=n-1; i++)
	{
		int j;
		for(j=1; j<=n-i; j++)
		{
			// Set things into order
			for (int u = 0; u<n; u++)
			{
				Lower_temp[u] = u;
				Upper_temp[u] = u;
				y_temp[u] = y[u];
				sigma_temp[u] = sigma[u];
			}
			// Apply the permutation
			for(int s=j; s<=j+i-1; s++)
			{
	    		y_temp[s] = y[s-1];
	  			sigma_temp[s] = sigma[s-1];
	  			Lower_temp[s] = s-1;
	  			Upper_temp[s] = s-1;
			}
			y_temp[j-1] = y[j+i-1];
			sigma_temp[j-1] = sigma[j+i-1];
			Lower_temp[j-1] = j+i-1;
			Upper_temp[j-1] = j+i-1;
			
			// Test if a PAVA is needed. If yes then break.
			// Could improve this part by using the means in the calculation of the LR later
			CheckPAVA = PAVACheck(y_temp, sigma_temp, l, InqPosi, n);
			if(CheckPAVA == true)
			{
				continue;
			}
			// Calculate the LR corresponding to the permuted vector
			LR = LogLikelihood(y_temp,sigma_temp,0,InqPosi[0]);
			k = 0;
			while (k <= (l - 2))
			{
				LR += LogLikelihood(y_temp,sigma_temp,InqPosi[k] + 1,InqPosi[k + 1]);
				k++;
			}
			LR += LogLikelihood(y_temp,sigma_temp,InqPosi[l - 1] + 1,n-1);		
			if(LR < crit[l])
			{
				RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);
				// Adjust the CI for the permuted centers
				for(int s=j; s<=j+i-1; s++)
				{
		  			Lower[s-1] = fmin(Lower[s-1], Lower_temp[s]);	
		  			Upper[s-1] = fmax(Upper[s-1], Upper_temp[s]);
				}
				Lower[j+i-1] = fmin(Lower[j+i-1], Lower_temp[j-1]);
				Upper[j+i-1] = fmax(Upper[j+i-1], Upper_temp[j-1]);
				
			}
			else
			{
				if(EqSigma) break;
			}
			
		}
		if(j<n-i && EqSigma) break;
	}

}

void CorrectPermutationsB(const NumericVector& y, const NumericVector& sigma, const NumericVector& crit, IntegerVector& Lower, IntegerVector& Upper, const int* InqPosi, const int& l, const int& n, const bool& EqSigma)
{
	NumericVector y_temp(n);
	NumericVector sigma_temp(n);
	IntegerVector Lower_temp(n);
	IntegerVector Upper_temp(n);
	
	double LR=0;
	int k = 0;
	bool CheckPAVA = false;
	for(int i=1; i<=n-1; i++)
	{
		int j;
		for(j=1; j<=n-i; j++)
		{
			// Set things into order
			for (int u = 0; u<n; u++)
			{
				Lower_temp[u] = u;
				Upper_temp[u] = u;
				y_temp[u] = y[u];
				sigma_temp[u] = sigma[u];
			}
			// Apply the permutation
			for(int s=j; s<=j+i-1; s++)
			{
	    		y_temp[s-1] = y[s];
	  			sigma_temp[s-1] = sigma[s];
	  			Lower_temp[s-1] = s;
	  			Upper_temp[s-1] = s;
			}
			y_temp[j+i-1] = y[j-1];
			sigma_temp[j+i-1] = sigma[j-1];
			Lower_temp[j+i-1] = j-1;
			Upper_temp[j+i-1] = j-1;
			
			// Test if a PAVA is needed. If yes then break.
			// Could improve this part by using the means in the calculation of the LR later
			CheckPAVA = PAVACheck( y_temp, sigma_temp, l, InqPosi, n);
			if(CheckPAVA == true)
			{
				continue;
			}
			// Calculate the LR corresponding to the permuted vector
			LR = LogLikelihood(y_temp,sigma_temp,0,InqPosi[0]);
			k = 0;
			while (k <= (l - 2))
			{
				LR += LogLikelihood(y_temp,sigma_temp,InqPosi[k] + 1,InqPosi[k + 1]);
				k++;
			}
			LR += LogLikelihood(y_temp,sigma_temp,InqPosi[l - 1] + 1,n-1);		
			
			if(LR < crit[l])
			{
				RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);
				// Adjust the CI for the permuted centers
				for(int s=j; s<=j+i-1; s++)
				{
		  			Lower[s] = fmin(Lower[s], Lower_temp[s-1]);	
		  			Upper[s] = fmax(Upper[s], Upper_temp[s-1]);
				}
				Lower[j-1] = fmin(Lower[j-1], Lower_temp[j+i-1]);
				Upper[j-1] = fmax(Upper[j-1], Upper_temp[j+i-1]);
				/*if(Upper[3]==8) {
					Rcout<<"\n I = "<<i<<", J = "<<j<<"\n InqPosi = (";
					for(int ss = 0; ss<l; ss++)
						Rcout<<InqPosi[ss]<<",";
					Rcout<<")\n";
				}*/
			}
			else
			{
				if(EqSigma) break;
			}
		}
		if(j<n-i && EqSigma) break;
	}
}
void UnrankCombin(int*& S, unsigned long long int m, int k, unsigned long long int**& CnkMat)
{
	int i = k - 1;
	while (i >= 0)
	{
		int l = i;
		while (CnkMat[l][i + 1] <= m)
		{
			l++;
		}
		S[i]= l - 1;
		m = m - CnkMat[l - 1][i + 1];
		i--;
	}

}
// An auxiliary function which does not permute anything. It is expected that the y's are not ordered and that a PAVA check has to be made.
void PartitioningRankingLevel(const NumericVector& y, const NumericVector& sigma, IntegerVector& Lower, IntegerVector& Upper, const NumericVector& crit, unsigned long long int**& CnkMat, const int& n, const bool& trace)
{
	bool CheckPAVA = false;
	double Likelihood0 = 0;
	// Calculate the Likelihood matrix of the blocks.
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
		}
	}
	
	int* InqPosi = new int[n];
	//if(trace == true) {
	//		Rcout<<"Processed levels:";
	//}
	for (int l = 1; l <= n - 2; l++)
	{
		//			int l = 3; int i = 8;
	//	if(trace == true) Rcout<<l<<".";
		unsigned long long int m = CnkMat[n - 1][l];
		for (unsigned long long int c = 0; c<m; c++)
		{
			UnrankCombin(InqPosi, c, l, CnkMat);
			// Test if a PAVA is needed. If yes then break.
			// Could improve this part by using the means in the calculation of the LR later
			CheckPAVA = PAVACheck(y, sigma, l, InqPosi, n);
			if(CheckPAVA == true)
			{
				continue;
			}
			// Inside each configuration, calculate each group's share in the likelihood
			Likelihood0 = LikelihoodMat[0][InqPosi[0]];
			int j = 0;
			while (j <= (l - 2))
			{
				Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
				j++;
			}
			Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];
			// Update the ranking
			if (Likelihood0<crit[l])
			{
				RankUpdate(Lower, Upper, InqPosi, l, n);				
			}
		}
		
	}
	delete[] InqPosi;

	/*NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}*/
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] LikelihoodMat[i];
	}
	delete[] LikelihoodMat;
	
	//return CIs;
}


// [[Rcpp::export]]
NumericMatrix PartitioningRankingLevelEqSig(NumericVector y, NumericVector sigma, NumericVector crit, int n, bool trace)
{
	
	// Calculate the Likelihood matrix of the blocks.
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
		}
	}
	// Calculate the C_n^k matrix in order to use it in the UnrankCombin function
	unsigned long long int ** CnkMat = new unsigned long long int*[n];
	for (int i = 0; i<n; i++)
	{
		CnkMat[i] = new unsigned long long int[n];
		CnkMat[i][i] = 1;
		for (int j = 0; j<i; j++)
		{
			CnkMat[i][j] = binomialCoeff(i, j);
			CnkMat[j][i] = 0;
		}
	}

	// Definition and initialization of the vector of ranks
	IntegerVector Lower(n);
	IntegerVector Upper(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
	}
	// Test the upper level with all equalities
	double Likelihood0 = LikelihoodMat[0][n - 1];
	if (Likelihood0<crit[0])
	{
		for (int i = 0; i<n; i++)
		{
			Lower[i] = 0;
			Upper[i] = n - 1;
		}
		if(trace == true) {
			Rcout << "Process ended with trivial confidence intervals.\n";
		}
	}
	else
	{
		int* InqPosi = new int[n];
		if(trace == true) {
				Rcout<<"Processed levels:";
		}
		for (int l = 1; l <= n - 2; l++)
		{
			//			int l = 3; int i = 8;
			if(trace == true) Rcout<<l<<".";
			unsigned long long int m = CnkMat[n - 1][l];
			for (unsigned long long int c = 0; c<m; c++)
			{
				UnrankCombin(InqPosi, c, l, CnkMat);
				// Inside each configuration, calculate each group's share in the likelihood
				Likelihood0 = LikelihoodMat[0][InqPosi[0]];
				int j = 0;
				while (j <= (l - 2))
				{
					Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
					j++;
				}
				Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];
				// Update the ranking
				if (Likelihood0<crit[l])
				{
					RankUpdate(Lower, Upper, InqPosi, l, n);
					CorrectPermutationsF(y, sigma, crit, Lower, Upper, InqPosi, l, n, true);
					CorrectPermutationsB(y, sigma, crit, Lower, Upper, InqPosi, l, n, true);
					
				}
				//	CorrectPermutationsF(y, sigma, crit, Lower, Upper, InqPosi, l, n, EqSigma, NbOfPermut);
				//	CorrectPermutationsB(y, sigma, crit, Lower, Upper, InqPosi, l, n, EqSigma);

			}
			
		}
		delete[] InqPosi;
	}
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] LikelihoodMat[i];
	}
	delete[] LikelihoodMat;
	for (int i = 0; i<n; i++)
	{
		delete[] CnkMat[i];
	}
	delete[] CnkMat;
	return CIs;
}


double PartitionCoverage(int***& ResCIsMat, const IntegerVector& Lower_temp, const IntegerVector& Upper_temp, const int& n, const int& MM, const int& a)
{
	int coveragePartition = MM;
	bool coverage;
	for(int j = 0; j<MM; j++)
	{
		coverage = true;
		for(int i = 0; i<n; i++)
		{
			//if(InqPosi[0]==2 && InqPosi[1]==4 && InqPosi[2]==8 && l==3 && trace == true)Rcout<<"["<<ResCIs(i,0)<<","<<ResCIs(i,1)<<"] contains ["<<Lower_temp[i]<<","<<Upper_temp[i]<<"]\n";
			if(ResCIsMat[a][i][2*j]>Lower_temp[i] || ResCIsMat[a][i][2*j+1]<Upper_temp[i])
			{
				coverage = false;
				break;
			}
		}
		//Rcout<<"....\n";
		if(coverage == false)
		{
			coveragePartition--;			
		}
	}
	return 1.0*coveragePartition/(MM*1.0);
}
void PartitioningRankingGeneralProcInit(int***& ResCIsMat, int***& ResCIsGridMat, int*& AlphaRescaled, IntegerVector& Lower, IntegerVector& Upper, const IntegerVector& EmpOrderInit, unsigned long long int**& CnkMat, const NumericMatrix& crit, const int& n, const int& MM, const int& gridSize, const double& alpha, const bool& trace)
{
	double* coveragePartition = new double[gridSize];
	int* InqPosi = new int[n];
	IntegerVector Lower_temp(n), Upper_temp(n), EmpOrder = seq(0,n-1);
	bool coverage;
	double minCoverage;
	int alphaRescaled;
	for (int l = 1; l <= n - 2; l++)
	{
		unsigned long long int m = CnkMat[n - 1][l];
		for (unsigned long long int c = 0; c<m; c++)
		{
			UnrankCombin(InqPosi, c, l, CnkMat);
			Lower_temp = clone(EmpOrder);
			Upper_temp = clone(EmpOrder);
			RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);

			coveragePartition[0] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, 0)-(1-alpha));
			minCoverage = coveragePartition[0]; alphaRescaled = 0;
			for(int a = 1; a<gridSize; a++)
			{
				coveragePartition[a] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, a)-(1-alpha));
				if(coveragePartition[a]<minCoverage)
				{
					alphaRescaled = a;
					minCoverage = coveragePartition[a];
				}
			}
			AlphaRescaled[(l-1)*CnkMat[n - 1][l-1] + c] = alphaRescaled;
			coverage = true;
			for(int i = 0; i<n; i++)
			{
				if(ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][0]>Lower_temp[i] || ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][1]<Upper_temp[i])
				{
					coverage = false;
					break;
				}
			}
			if(coverage == true)
			{
				for(int s = 0 ; s<n; s++)
				{
					
					Lower[s] = fmin(Lower[s], Lower_temp[s]);
					Upper[s] = fmax(Upper[s], Upper_temp[s]);
					
				}
			}
		}
		
	}
	delete[] InqPosi;
	delete[] coveragePartition;
}
void PartitioningRankingGeneralProc(int***& ResCIsMat, int***& ResCIsGridMat, int*& AlphaRescaled, IntegerVector& Lower, IntegerVector& Upper, const IntegerVector& EmpOrderInit, unsigned long long int**& CnkMat, const NumericMatrix& crit, const int& n, const int& MM, const int& gridSize, const double& alpha, const bool& trace)
{
	int* InqPosi = new int[n];
	IntegerVector Lower_temp(n), Upper_temp(n), EmpOrder = seq(0,n-1);
	bool coverage;
	int alphaRescaled;
	for (int l = 1; l <= n - 2; l++)
	{
		unsigned long long int m = CnkMat[n - 1][l];
		for (unsigned long long int c = 0; c<m; c++)
		{
			UnrankCombin(InqPosi, c, l, CnkMat);
			Lower_temp = clone(EmpOrder);
			Upper_temp = clone(EmpOrder);
			RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);

			alphaRescaled = AlphaRescaled[(l-1)*CnkMat[n - 1][l-1] + c];
			coverage = true;
			for(int i = 0; i<n; i++)
			{
				if(ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][0]>Lower_temp[i] || ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][1]<Upper_temp[i])
				{
					coverage = false;
					break;
				}
			}
			if(coverage == true)
			{
				for(int s = 0 ; s<n; s++)
				{
					
					Lower[s] = fmin(Lower[s], Lower_temp[s]);
					Upper[s] = fmax(Upper[s], Upper_temp[s]);
					
				}
			}
		}
		
	}
	delete[] InqPosi;
}

void RescalTopHpo(int***& ResCIsMat, int***& ResCIsGridMat, IntegerVector& Lower, IntegerVector& Upper, const IntegerVector& EmpOrderInit, const int& n, const int& MM, const int& gridSize, const double& alpha)
{
	double* coveragePartition = new double[gridSize];
	IntegerVector Lower_temp(n), Upper_temp(n);
	bool coverage;
	double minCoverage;
	int alphaRescaled;

	for(int i = 0; i < n; i++)
	{
		Lower_temp[i] = 0;
		Upper_temp[i] = n-1;
	}
	
	
			coveragePartition[0] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, 0)-(1-alpha));
			minCoverage = coveragePartition[0]; alphaRescaled = 0;
			for(int a = 1; a<gridSize; a++)
			{
				coveragePartition[a] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, a)-(1-alpha));
				if(coveragePartition[a]<minCoverage)
				{
					alphaRescaled = a;
					minCoverage = coveragePartition[a];
				}
			}
			
			coverage = true;
			for(int i = 0; i<n; i++)
			{
				if(ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][0]>Lower_temp[i] || ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][1]<Upper_temp[i])
				{
					coverage = false;
					break;
				}
			}
			if(coverage == true)
			{
				for(int s = 0 ; s<n; s++)
				{
					
					Lower[s] = 0;
					Upper[s] = n-1;
					
				}
			}
		
		
	
	//delete[] InqPosi;
	delete[] coveragePartition;
}

// [[Rcpp::export]]
NumericMatrix PartitioningRankingLevelEqSigRescaled(NumericVector y, NumericVector sigma, NumericMatrix crit, NumericMatrix SampleCoverage, int MM, int n, int NbOfPermut, double alpha, int gridSize, bool trace)
{
	int powN2 = 1;
	for(int i = 0; i<n; i++)
	{
		powN2 *=2; 
	}
	//Rcout<<"powN2 = "<<powN2<<"\n";
	//Rcout<<"\n Start : gridSize = "<<gridSize<<", MM = "<<MM<<"\n";
	if(trace == true) Rcout<<"Initializing variables...\n";
	int*** ResCIsMat = new int**[gridSize]; // alpha values range in 0.05, ..., 0.4
	int*** ResCIsMat_temp = new int**[gridSize];
	int*** ResCIsGridMat = new int**[gridSize];
	int* AlphaRescaled = new int[powN2];
	for(int i = 0; i<gridSize; i++)
	{
		ResCIsMat[i] = new int*[n];
		ResCIsMat_temp[i] = new int*[n];
		ResCIsGridMat[i] = new int*[n];
		for(int j = 0; j<n; j++)
		{
			ResCIsMat[i][j] = new int[2*MM];
			ResCIsMat_temp[i][j] = new int[2*MM];
			ResCIsGridMat[i][j] = new int[2];
		}
	}
	
	NumericMatrix ResCIs(n,2);
	if(trace == true) Rcout<<"Calculating SCIs for the ranks of a sample.\n";
	for(int a = 0; a<gridSize; a++)
	{
		for(int s = 0; s<MM; s++)
		{					
			// Caclulate simultaneous CIs for ranks using sample i
			ResCIs = PartitioningRankingLevelEqSig(SampleCoverage.column(s), sigma, crit.column(a), n, false);
			for(int i = 0; i<n; i++)
			{
				ResCIsMat[a][i][2*s] = ResCIs(i,0)-1; // Even for Lower CIs
				ResCIsMat[a][i][2*s+1] = ResCIs(i,1)-1; // Odd for Upper CIs	
			}
		}
		ResCIs = PartitioningRankingLevelEqSig(y, sigma, crit.column(a), n, false);
		for(int i = 0; i<n; i++)
		{
			ResCIsGridMat[a][i][0] = ResCIs(i,0)-1;
			ResCIsGridMat[a][i][1] = ResCIs(i,1)-1;
		}
		
	}
		
	IntegerVector RandPermutInit, RandPermut;
	RandPermutInit = seq(0,n-1);
	// Calculate the C_n^k matrix in order to use it in the UnrankCombin function
	unsigned long long int ** CnkMat = new unsigned long long int*[n];
	for (int i = 0; i<n; i++)
	{
		CnkMat[i] = new unsigned long long int[n];
		CnkMat[i][i] = 1;
		for (int j = 0; j<i; j++)
		{
			CnkMat[i][j] = binomialCoeff(i, j);
			CnkMat[j][i] = 0;
		}
	}

	// Definition and initialization of the vector of ranks
	IntegerVector Lower(n), Lower_temp(n), EmpOrderInit(n);
	IntegerVector Upper(n), Upper_temp(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
		EmpOrderInit[i] = i;
	}
	// Calculate first the adjustment for the hypothesis mu_1=...=mu_n
	RescalTopHpo(ResCIsMat, ResCIsGridMat, Lower, Upper, EmpOrderInit, n, MM, gridSize, alpha);
	
	// Calculate CIs for ranks using the general proc without any permutations (correctly ordered hypotheses)
	if(trace == true) Rcout<<"Perform the correctly ordered hypotheses...\n";
	PartitioningRankingGeneralProcInit(ResCIsMat, ResCIsGridMat, AlphaRescaled, Lower, Upper, EmpOrderInit, CnkMat, crit, n, MM, gridSize, alpha, trace);
	bool CheckTrivialCIs = true;
	for(int i = 0; i<n; i++)
	{
		if(Lower[i] != 0 || Upper[i] != n-1)
		{
			CheckTrivialCIs = false;
			break;
		}
	}
	if(CheckTrivialCIs == false){
	if(trace == true) {Rcout<<"\n Forward permutations ("<<n-1<<" permutations)\n";}
	//### part1: permute (i,i+1,i+2,...)
	for(int I = 1; I<=n-1; I++)//#row index
	{
		if(trace == true) {Rcout<<I<<".";}
		for(int J = 1; J <= n-I; J++)//#column index
		{
			CheckTrivialCIs = true;
			for(int i = 0; i<n; i++)
			{
				if(Lower[i] != 0 || Upper[i] != n-1)
				{
					CheckTrivialCIs = false;
					break;
				}
			}
			if(CheckTrivialCIs == true) break;
			for(int s = 0; s<n; s++)
			{
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}			
			}
			
			//# Permute
			for(int s = J; s<=J+I-1; s++)
			{
				Lower_temp[s] = EmpOrderInit[s-1]; 
				Upper_temp[s] = EmpOrderInit[s-1];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s-1][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s-1][2*j+1];
					}
				}			
			}
			Lower_temp[J-1] = EmpOrderInit[J+I-1];	
			Upper_temp[J-1] = EmpOrderInit[J+I-1];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][J-1][2*j] = ResCIsMat[i][J+I-1][2*j];
					ResCIsMat_temp[i][J-1][2*j+1] = ResCIsMat[i][J+I-1][2*j+1];
				}
			}
		
			//# Calculate the SCI for ranks
			PartitioningRankingGeneralProc(ResCIsMat_temp, ResCIsGridMat, AlphaRescaled, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);
			
			//# Permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
				Lower[s-1] = fmin(Lower[s-1], Lower_temp[s]);
				Upper[s-1] = fmax(Upper[s-1], Upper_temp[s]);
			}
			Lower[J+I-1] = fmin(Lower[J+I-1], Lower_temp[J-1]);
			Upper[J+I-1] = fmax(Upper[J+I-1], Upper_temp[J-1]);
		}
	}
	
	if(trace == true) {Rcout<<"\n Applying backward permutations.\n";}
	int I;
	//### Part 2: permute (i+4,i+3,..)
	for(int J = 1; J<=n-1; J++)//#column index
	{
		if(trace == true) {Rcout<<J<<".";}
		I = 2;
		while(I<=(n-J))//#row index
		{
			CheckTrivialCIs = true;
			for(int i = 0; i<n; i++)
			{
				if(Lower[i] != 0 || Upper[i] != n-1)
				{
					CheckTrivialCIs = false;
					break;
				}
			}
			if(CheckTrivialCIs == true) break;

			for(int s = 0; s<n; s++)
			{
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}
			}
			//# Permute
			for(int s = J; s<= J+I-1; s++)
			{
			  	Lower_temp[s-1] = EmpOrderInit[s];
				Upper_temp[s-1] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s-1][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s-1][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}
			}
			Lower_temp[J+I-1] = EmpOrderInit[J-1];
			Upper_temp[J+I-1] = EmpOrderInit[J-1];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][J+I-1][2*j] = ResCIsMat[i][J-1][2*j];
					ResCIsMat_temp[i][J+I-1][2*j+1] = ResCIsMat[i][J-1][2*j+1];
				}
			}
						
			PartitioningRankingGeneralProc(ResCIsMat_temp, ResCIsGridMat, AlphaRescaled, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);
			
			//# permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
			  Lower[s] = fmin(Lower[s], Lower_temp[s-1]);
			  Upper[s] = fmax(Upper[s], Upper_temp[s-1]);
			}
			Lower[J-1] = fmin(Lower[J-1], Lower_temp[J+I-1]);
			Upper[J-1] = fmax(Upper[J-1], Upper_temp[J+I-1]);
		 	I = I+1;
		}
	}
	
	if(trace == true) {Rcout<<"\n Applying "<<NbOfPermut<<" Random permutations.\n";}
	for(int I = 1; I<=NbOfPermut; I++)//#row index
	{
		CheckTrivialCIs = true;
		for(int i = 0; i<n; i++)
		{
			if(Lower[i] != 0 || Upper[i] != n-1)
			{
				CheckTrivialCIs = false;
				break;
			}
		}
		if(CheckTrivialCIs == true) break;

		if((trace == true) && (I % 500000 == 0)) {Rcout<<I<<".";}
		//# Permute
		RandPermut = sample(RandPermutInit,n);//clone(RandPermutInit);
		//std::random_shuffle(RandPermut.begin(), RandPermut.end());
		for(int s = 0; s<n; s++)
		{
			Lower_temp[s] = EmpOrderInit[RandPermut[s]];
			Upper_temp[s] = EmpOrderInit[RandPermut[s]];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][RandPermut[s]][2*j];
					ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][RandPermut[s]][2*j+1];
				}
			}
		}
		PartitioningRankingGeneralProc(ResCIsMat_temp, ResCIsGridMat, AlphaRescaled, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);

		//# Permute the ranks
		for(int s = 0; s<n; s++)
		{
			Lower[RandPermut[s]] = fmin(Lower[RandPermut[s]], Lower_temp[s]);
			Upper[RandPermut[s]] = fmax(Upper[RandPermut[s]], Upper_temp[s]);
		}
	}
	}
	Rcout<<"\n";
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] CnkMat[i];
	}
	delete[] CnkMat;
	for(int a = 0; a<gridSize; a++)
	{
		for(int i = 0; i<n; i++)
		{
			delete[] ResCIsMat[a][i];
			delete[] ResCIsMat_temp[a][i];
			delete[] ResCIsGridMat[a][i];
		}
		delete[] ResCIsMat[a];
		delete[] ResCIsMat_temp[a];
		delete[] ResCIsGridMat[a];
	}
	delete[] ResCIsMat;
	delete[] ResCIsMat_temp;
	delete[] ResCIsGridMat;
	delete[] AlphaRescaled;

	return CIs;
}


// [[Rcpp::export]]
NumericMatrix PartitioningRankingLevelUneqSig(NumericVector y, NumericVector sigma, NumericVector crit, int n, bool trace, const int NbOfPermut, bool SwapPerm)
{
	NumericVector y_temp(n), sigma_temp(n);
	IntegerVector RandPermutInit, RandPermut;
	RandPermutInit = seq(0,n-1);
	// Calculate the C_n^k matrix in order to use it in the UnrankCombin function
	unsigned long long int ** CnkMat = new unsigned long long int*[n];
	for (int i = 0; i<n; i++)
	{
		CnkMat[i] = new unsigned long long int[n];
		CnkMat[i][i] = 1;
		for (int j = 0; j<i; j++)
		{
			CnkMat[i][j] = binomialCoeff(i, j);
			CnkMat[j][i] = 0;
		}
	}

	// Definition and initialization of the vector of ranks
	IntegerVector Lower(n), Lower_temp(n), EmpOrderInit(n);
	IntegerVector Upper(n), Upper_temp(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
		EmpOrderInit[i] = i;
	}
	if(trace == true) {Rcout<<"\n Start without permutations...\n";}
	// Perform a partitioning without any permutation
	PartitioningRankingLevel(y, sigma, Lower, Upper, crit, CnkMat, n, trace);
	
	if(trace == true) {Rcout<<"\n Forward permutations ("<<n-1<<" permutations)\n";}
	//### part1: permute (i,i+1,i+2,...)
	for(int I = 1; I<=n-1; I++)//#row index
	{
		if(trace == true) {Rcout<<I<<".";}
		for(int J = 1; J <= n-I; J++)//#column index
		{
			//y_temp = y; sigma_temp = sigma; POINTER COPYING !!!
			//EmpOrder_temp = EmpOrder;
			for(int s = 0; s<n; s++)
			{
				y_temp[s] = y[s]; 
				sigma_temp[s] = sigma[s];
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
			}
			
			//# Permute
			for(int s = J; s<=J+I-1; s++)
			{
				y_temp[s] = y[s-1];
				sigma_temp[s] = sigma[s-1];
				Lower_temp[s] = EmpOrderInit[s-1]; 
				Upper_temp[s] = EmpOrderInit[s-1];
			}
			y_temp[J-1] = y[J+I-1];
			sigma_temp[J-1] = sigma[J+I-1];
			Lower_temp[J-1] = EmpOrderInit[J+I-1];	
			Upper_temp[J-1] = EmpOrderInit[J+I-1];
			
			//# Calculate the SCI for ranks
			PartitioningRankingLevel(y_temp, sigma_temp, Lower_temp, Upper_temp, crit, CnkMat, n, trace);
			
			
			//# Permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
				Lower[s-1] = fmin(Lower[s-1], Lower_temp[s]);
				Upper[s-1] = fmax(Upper[s-1], Upper_temp[s]);
			}
			Lower[J+I-1] = fmin(Lower[J+I-1], Lower_temp[J-1]);
			Upper[J+I-1] = fmax(Upper[J+I-1], Upper_temp[J-1]);
			//if(trace == true) {Rcout<<"\n";}
			
			//if(sum(Lower==rep(1,n) & Upper==rep(n,n)) == n) return(list(Lower=Lower,Upper=Upper))
		}
	}
	
	if(trace == true) {Rcout<<"\n Applying backward permutations.\n";}
	int I;
	//### Part 2: permute (i+4,i+3,..)
	for(int J = 1; J<=n-1; J++)//#column index
	{
		if(trace == true) {Rcout<<J<<".";}
		I = 2;
		while(I<=(n-J))//#row index
		{
			//y_temp = y; sigma_temp = sigma; EmpOrder_temp = EmpOrder;
			for(int s = 0; s<n; s++)
			{
				y_temp[s] = y[s]; 
				sigma_temp[s] = sigma[s];
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
			}
			//# Permute
			for(int s = J; s<= J+I-1; s++)
			{
			  	y_temp[s-1] = y[s];
			  	sigma_temp[s-1] = sigma[s];
			  	Lower_temp[s-1] = EmpOrderInit[s];
				Upper_temp[s-1] = EmpOrderInit[s];
			}
			y_temp[J+I-1] = y[J-1];
			sigma_temp[J+I-1] = sigma[J-1];
			Lower_temp[J+I-1] = EmpOrderInit[J-1];
			Upper_temp[J+I-1] = EmpOrderInit[J-1];
						
			//# Calculate the SCI for ranks
			//#res = ApproximatePartitionWrongOrder(y_temp, sigma_temp, EmpOrder_temp,critFun, Intercept)
			PartitioningRankingLevel(y_temp, sigma_temp, Lower_temp, Upper_temp, crit, CnkMat, n, trace);
			
			
			//# permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
			  Lower[s] = fmin(Lower[s], Lower_temp[s-1]);
			  Upper[s] = fmax(Upper[s], Upper_temp[s-1]);
			}
			Lower[J-1] = fmin(Lower[J-1], Lower_temp[J+I-1]);
			Upper[J-1] = fmax(Upper[J-1], Upper_temp[J+I-1]);
		 	I = I+1;
			//if(sum(Lower==rep(1,n) & Upper==rep(n,n)) == n) return(list(Lower=Lower,Upper=Upper))
			//if(Lower[10] == 1) Rcout<<"I="<<I<<", J="<<J<<"\n";
			//if(trace == true) {Rcout<<"\n";}
		}
	}
	
	if(SwapPerm == true)
	{
		if(trace == true) {Rcout<<"\n Apply "<<n<<" swap permutations \n";}// They do help some times
		for(int I = 1; I<=n-1; I++)
		{
			if(trace == true) {Rcout<<I<<".";}
			for(int J = I+1; J<=n; J++)
			{
				// Swap mu_i with mu_j
				y_temp[I-1] = y[J-1];
				y_temp[J-1] = y[I-1];
				sigma_temp[I-1] = sigma[J-1];
				sigma_temp[J-1] = sigma[I-1];
				Lower_temp[I-1] = EmpOrderInit[J-1];
				Lower_temp[J-1] = EmpOrderInit[I-1];
				Upper_temp[I-1] = EmpOrderInit[J-1];
				Upper_temp[J-1] = EmpOrderInit[I-1];
				
				//Rcout<<"\n (I,J)=("<<I<<","<<J<<").";
				PartitioningRankingLevel(y_temp, sigma_temp, Lower_temp, Upper_temp, crit, CnkMat, n, trace);
				
				Lower[I-1] = fmin(Lower[I-1], Lower_temp[J-1]);
				Lower[J-1] = fmin(Lower[J-1], Lower_temp[I-1]);
				Upper[I-1] = fmax(Upper[I-1], Upper_temp[J-1]);
				Upper[J-1] = fmax(Upper[J-1], Upper_temp[I-1]);
				
				// Swap back
				y_temp[I-1] = y[I-1];
				y_temp[J-1] = y[J-1];
				sigma_temp[I-1] = sigma[I-1];
				sigma_temp[J-1] = sigma[J-1];
				Lower_temp[I-1] = EmpOrderInit[I-1];
				Lower_temp[J-1] = EmpOrderInit[J-1];
				Upper_temp[I-1] = EmpOrderInit[I-1];
				Upper_temp[J-1] = EmpOrderInit[J-1];
			}
		}
	}

	if(trace == true) {Rcout<<"\n Apply "<<NbOfPermut<<" random permutations \n";}
	
	for(int I = 1; I<=NbOfPermut; I++)//#row index
	{
		if((trace == true) && (I % 100 == 0)) {Rcout<<I<<".";}
					
		//# Permute
		RandPermut = sample(RandPermutInit,n);//clone(RandPermutInit);
		//std::random_shuffle(RandPermut.begin(), RandPermut.end(), randWrapper);
		for(int s = 0; s<n; s++)
		{
			y_temp[s] = y[RandPermut[s]]; 
		//	Rcout<<y_temp[s]<<",";
			sigma_temp[s] = sigma[RandPermut[s]];
			Lower_temp[s] = EmpOrderInit[RandPermut[s]];
			Upper_temp[s] = EmpOrderInit[RandPermut[s]];
		}
		//Rcout<<"\n";
		PartitioningRankingLevel(y_temp, sigma_temp, Lower_temp, Upper_temp, crit, CnkMat, n, trace);

		//# Permute the ranks
		for(int s = 0; s<n; s++)
		{
			Lower[RandPermut[s]] = fmin(Lower[RandPermut[s]], Lower_temp[s]);
			Upper[RandPermut[s]] = fmax(Upper[RandPermut[s]], Upper_temp[s]);
		}
		//if(trace == true) {Rcout<<"\n";}
	}
	

	
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] CnkMat[i];
	}
	delete[] CnkMat;
	return CIs;
}

// [[Rcpp::export]]
NumericMatrix PartitioningRankingBlockCorrectOrder(NumericVector y, NumericVector sigma, NumericVector crit, NumericVector MinBlock, NumericVector MaxBlock, IntegerVector Lower, IntegerVector Upper, int n, bool trace)
{
	// Calculate the Likelihood matrix of the blocks
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
		}
	}
	// Calculate the vector of powers to 2.
	unsigned long long int* PowToN = new unsigned long long int[n];
	PowToN[0] = 1;
	for (int i = 1; i < n; i++)
	{
		PowToN[i] = 2 * PowToN[i-1];
	}

	// Test the upper level with all equalities
	double Likelihood0 = LikelihoodMat[0][n - 1];
	if (Likelihood0<crit[0])
	{
		for (int i = 0; i<n; i++)
		{
			Lower[i] = 0;
			Upper[i] = n - 1;
		}
		if(trace == true){
			Rcout << "Process ended with trivial confidence intervals.\n";
		} 
	}
	else
	{
		int ConfigBase[2]; // It points out to the extremeties of the tested block.
		int* InqPosi = new int[n - 1];
		int* ConfigCompRight = new int[n - 1]; 
		int ConfigBaseLen = 2;
		int BlockConfig[2];
		for (int i = 0; i < n - 1; i++)
		{
			if(trace == true) {
				Rcout << i << ".";
			}
			if (MaxBlock[i] >= (MinBlock[i] + 1))
			{
				int k = MaxBlock[i];
				while (k >= MinBlock[i] + 1)
				{
					if (i > 0)
					{
						ConfigBase[0] = i - 1;
						ConfigBase[1] = k + i;
						ConfigBaseLen = 2;
					}
					if (i == 0)
					{
						ConfigBase[0] = k;
						ConfigBaseLen = 1;
					}
					if (k + i >= n-1)
					{
						ConfigBase[0] = i - 1;
						ConfigBaseLen = 1;
					}
					if (i == 0 && k + i >= n-1) {
						ConfigBaseLen = 0;
					}
					ConfigCompRight[0] = ConfigBase[0];
					ConfigCompRight[1] = ConfigBase[1];
					BlockConfig[0] = ConfigBase[0];
					BlockConfig[1] = ConfigBase[1];
					if (BlockConfig[1] > n - 1) {
						BlockConfig[1] = n - 1;
					}

					/********************************************************/
					/************ First Case: Left configuration ***********/
					/********************************************************/

					if (k + i >= n - 2) // substract 1 for subscribts diff
					{
						BlockConfig[1] = n - 1;
						// In this case, there is only one side where Partitioing must be done; the left side.
						// The partitions must cover mu_0 to mu_{i-2}
						unsigned long long int m = 2;
						unsigned long long int c = 1;
						if (i > 2)
						{
							int j = 0, l;
							m = PowToN[i-1];
							// Test the top hypothesis corresponding to mu_1=...=mu_{i-1}<mu_i=..=mu_n
							InqPosi[0] = ConfigBase[0];
							l = 1;
							if (ConfigBaseLen == 2) 
							{
								InqPosi[1] = ConfigBase[1];
								l++;
							}
							// Inside each configuration, calculate each group's share in the likelihood
							Likelihood0 = LikelihoodMat[0][InqPosi[0]];
							j = 0;
							while (j <= (l - 2))
							{
								Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
								j++;
							}
							Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];
							if (Likelihood0 < crit[l])
							{
								
								RankUpdate(Lower, Upper, InqPosi, l, n);
								
								break;
							}
							// Check the significance of the current block to the actual ranking. If not significant, then there is no need to start the partitioning. Smaller blocks will also be the same.
							//if (Lower[BlockConfig[1]]<=BlockConfig[0] + 1 && Upper[BlockConfig[0] + 1]>=BlockConfig[1]) break;
							// Check the block itself without any additions to the left or the right
							for (c = 0; c<m ; c++)
							{
								BinaryConfig(c, InqPosi, l, 0, 0);
								// Add ConfigBase to InqPosi. Keep in mind that we only use the l+ConfigBaseLen first elements.
								InqPosi[l] = ConfigBase[0];
								InqPosi[l + 1] = ConfigBase[1];
								// Update the length of InqPosi
								l += ConfigBaseLen; // If ConfigBase is empty, the previous attributions will not have effect on the configuration since its length is controled by l.


								/*********************************/
								// Inside each configuration, calculate each group's share in the likelihood
								Likelihood0 = LikelihoodMat[0][InqPosi[0]];
								j = 0;
								while (j <= (l - 2))
								{
									Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
									j++;
								}
								Likelihood0 += LikelihoodMat[InqPosi[l - 1] + 1][n - 1];

								// Update the ranking
								if (Likelihood0 < crit[l])
								{
									RankUpdate(Lower, Upper, InqPosi, l, n);
									break;// The block was accepted once, and no need to check others.
								}
							}
							if (c < m) break;
						}
					}



					/***************************************************************************/
					/************ Second Case: Right and possibly Left configuration ***********/
					/***************************************************************************/

					else
					{
						unsigned long long int m = 2;
						unsigned long long int c = 1; // I need the counter to check if I was able to find a hypothesis that was not rejected.
						//unsigned long long int c1 = 1;
						if (n - k - i > 2)
						{
							
							// Test the top level in this sub-partition
							int j = 0, l;
							m = PowToN[n - k - i - 2]; // The right part contains n-i-k-1 centers.
							// Check the significance of the current block to the actual ranking. If not significant, then there is no need to start the partitioning. Smaller blocks will also be the same.
							//if (Lower[BlockConfig[1]]<=BlockConfig[0] + 1 && Upper[BlockConfig[0] + 1]>=BlockConfig[1]) break;
							for (c = 0; c < m ; c++)
							{
								BinaryConfig(c, ConfigCompRight, l, k + i + 1, ConfigBaseLen);// ConfigBase if existed, is already included in ConfigCompRight.
								// Update the length of ConfigCompRight.
								l += ConfigBaseLen;
								/*********************************/
								// Inside each configuration, calculate each group's share in the likelihood
								Likelihood0 = LikelihoodMat[0][ConfigCompRight[0]];
								j = 0;
								while (j <= (l - 2))
								{
									Likelihood0 += LikelihoodMat[ConfigCompRight[j] + 1][ConfigCompRight[j + 1]];
									j++;
								}
								Likelihood0 += LikelihoodMat[ConfigCompRight[l - 1] + 1][n - 1];

								// Update the ranking
								if (Likelihood0 < crit[l])
								{
									RankUpdate(Lower, Upper, ConfigCompRight, l, n);
									break;// The block was accepted once, and no need to check others.
								}


								// Add a configuration to the left of ConfigBase in case there is suffcient place.
								/***************Third case: a configuration to the left and a configuration to the right*******************/

								// The top hypothesis in this sub-partitioning was just tested.
								unsigned long long int cc = 1;
								//unsigned long long int cc1 = 1;
								unsigned long long int mm = 2; // Initial values here are necessary to cancel the break hereafter in case, we do not enter the loop over cc.
								if (i > 1)
								{
									int ll;
									mm = PowToN[i - 1];
									// Notice that the upper hypothesis in this sub-partitioning corresponds to testing the block itself. This is already done in the initial CIs.
									for (cc = 0; cc<mm ; cc++) // If the left side does not contain any elements, this loop wont start.
									{

										BinaryConfig(cc, InqPosi, ll, 0, 0);
										// Add now the configuration in the right side
										for (int iter = 0; iter < l; iter++) InqPosi[iter + ll] = ConfigCompRight[iter];
										// Update the length of InqPosi
										ll += l;

										/*********************************/
										// Inside each configuration, calculate each group's share in the likelihood
										Likelihood0 = LikelihoodMat[0][InqPosi[0]];
										j = 0;
										while (j <= (ll - 2))
										{
											Likelihood0 += LikelihoodMat[InqPosi[j] + 1][InqPosi[j + 1]];
											j++;
										}
										Likelihood0 += LikelihoodMat[InqPosi[ll - 1] + 1][n - 1];

										// Update the ranking
										if (Likelihood0 < crit[ll])
										{
											RankUpdate(Lower, Upper, InqPosi, ll, n);
											break;// The block was accepted once, and no need to check others.
										}
									}
									if (cc < mm) break;
								}

							}
							if (c < m) break;
						}
					}
					k--;
				}

			}

		}

		

		delete[] InqPosi;
		delete[] ConfigCompRight;

	}

	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] LikelihoodMat[i];
	}
	delete[] LikelihoodMat;
	delete[] PowToN;
	
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	return CIs;
}


// An internal function for an instance of observations.
NumericMatrix OnlyBlockRanking_instance(const NumericVector& y, const NumericVector& sigma, const NumericVector& crit, IntegerVector& Lower, IntegerVector& Upper,  int n)
{
	// Detect equal-sigma for faster implementation
	//bool EqSigma = true;
	/*for(int i=0;i<n-1;i++)
	{
		if(sigma[i] != sigma[i+1])
		{
			EqSigma = false;
			break;
		}
	}*/
	bool CheckPAVA = false;
	// Calculate the Likelihood matrix of the blocks.
	double** LikelihoodMat = new double*[n];
	for (int i = 0; i<n; i++)
	{
		LikelihoodMat[i] = new double[n];
		for (int j = i; j<n; j++)
		{
			LikelihoodMat[i][j] = LogLikelihood(y, sigma, i, j);
		}
	}
	
	// Definition and initialization of the vector of ranks
	/*IntegerVector Lower(n);
	IntegerVector Upper(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
	}*/
	// Test the upper level with all equalities
	double Likelihood0 = LikelihoodMat[0][n - 1];

		int* InqPosi = new int[n];
		int l,Lleft,Lright;
		for (int i = 0; i <= n - 2; i++)
		{
			Lleft = 0;
			if(i>0) 
			{
				for(int s = 0; s<i-1; s++)
				{
					InqPosi[Lleft] = s;
					Lleft++;
				}
				InqPosi[Lleft] = i-1;
				Lleft++;	
			}
			 //Rcout<<i<<".";
			for (int j = i+1; j<n; j++)
			{
				// NOTE: When i = 0 and j = n-1, InqPosi is empty. It corresponds to mu_1=...=mu_n. This is already tested and will not be repeated here again.
				if(i ==0 && j == n-1) continue;
				Lright = 0;
				if(j<n-1) 
				{
					InqPosi[Lleft+Lright]=j;
					Lright++;
					for(int s = j; s<=n-3; s++)
					{
						InqPosi[Lleft+Lright] = s+1;
						Lright++;
					}
				}
				
				l = Lleft + Lright;
				//Rcout<<"\n l="<<l<<". InqPosi=("<<InqPosi[0]<<",...,"<<InqPosi[l-1]<<"\n";
				/*Rcout<<"\n InqPosi = (";
				for(int ss = 0; ss<l; ss++) Rcout<<InqPosi[ss]<<",";
				Rcout<<"\n";*/
				CheckPAVA = PAVACheck(y, sigma, l, InqPosi, n);
				if(CheckPAVA == true)
				{
					continue;
				}
				// Inside each configuration, calculate each group's share in the likelihood
				Likelihood0 = LikelihoodMat[i][j];
				// Update the ranking
				if (Likelihood0<crit[l])
				{
					RankUpdate(Lower, Upper, InqPosi, l, n);
				}
			}
		}
		delete[] InqPosi;
	
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] LikelihoodMat[i];
	}
	delete[] LikelihoodMat;
	
	return CIs;
}

// ----------------------------------------------------------
// .......... To do Function ................................
// NumericMatrix OnlyBlockRanking(NumericVector y, NumericVector sigma, NumericVector crit, int n, bool trace, const int NbOfPermut)
// ---------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix OnlyBlockRanking(NumericVector y, NumericVector sigma, NumericVector crit, int n, bool trace, const int NbOfPermut, bool SwapPerm)
{
	NumericVector y_temp(n), sigma_temp(n);
	IntegerVector RandPermutInit, RandPermut;
	RandPermutInit = seq(0,n-1);
	// Calculate the C_n^k matrix in order to use it in the UnrankCombin function
	

	// Definition and initialization of the vector of ranks
	IntegerVector Lower(n), Lower_temp(n), EmpOrderInit(n);
	IntegerVector Upper(n), Upper_temp(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
		EmpOrderInit[i] = i;
	}
	if(trace == true) {Rcout<<"\n Start without permutations...\n";}
	// Perform a partitioning without any permutation
	//OnlyBlockRanking_instance(y, sigma, crit, Lower, Upper, n);
	
	if(trace == true) {Rcout<<"\n Forward permutations ("<<n-1<<" permutations)\n";}
	//### part1: permute (i,i+1,i+2,...)
	for(int I = 1; I<=n-1; I++)//#row index
	{
		if(trace == true) {Rcout<<I<<".";}
		for(int J = 1; J <= n-I; J++)//#column index
		{
			//y_temp = y; sigma_temp = sigma; POINTER COPYING !!!
			//EmpOrder_temp = EmpOrder;
			for(int s = 0; s<n; s++)
			{
				y_temp[s] = y[s]; 
				sigma_temp[s] = sigma[s];
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
			}
			
			//# Permute
			for(int s = J; s<=J+I-1; s++)
			{
				y_temp[s] = y[s-1];
				sigma_temp[s] = sigma[s-1];
				Lower_temp[s] = EmpOrderInit[s-1]; 
				Upper_temp[s] = EmpOrderInit[s-1];
			}
			y_temp[J-1] = y[J+I-1];
			sigma_temp[J-1] = sigma[J+I-1];
			Lower_temp[J-1] = EmpOrderInit[J+I-1];	
			Upper_temp[J-1] = EmpOrderInit[J+I-1];
			
			//# Calculate the SCI for ranks
			OnlyBlockRanking_instance(y_temp, sigma_temp, crit, Lower_temp, Upper_temp, n);
			
			
			//# Permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
				Lower[s-1] = fmin(Lower[s-1], Lower_temp[s]);
				Upper[s-1] = fmax(Upper[s-1], Upper_temp[s]);
			}
			Lower[J+I-1] = fmin(Lower[J+I-1], Lower_temp[J-1]);
			Upper[J+I-1] = fmax(Upper[J+I-1], Upper_temp[J-1]);
			//if(trace == true) {Rcout<<"\n";}
			
			//if(sum(Lower==rep(1,n) & Upper==rep(n,n)) == n) return(list(Lower=Lower,Upper=Upper))
		}
	}
	
	if(trace == true) {Rcout<<"\n Applying backward permutations.\n";}
	int I;
	//### Part 2: permute (i+4,i+3,..)
	for(int J = 1; J<=n-1; J++)//#column index
	{
		if(trace == true) {Rcout<<J<<".";}
		I = 2;
		while(I<=(n-J))//#row index
		{
			//y_temp = y; sigma_temp = sigma; EmpOrder_temp = EmpOrder;
			for(int s = 0; s<n; s++)
			{
				y_temp[s] = y[s]; 
				sigma_temp[s] = sigma[s];
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
			}
			//# Permute
			for(int s = J; s<= J+I-1; s++)
			{
			  	y_temp[s-1] = y[s];
			  	sigma_temp[s-1] = sigma[s];
			  	Lower_temp[s-1] = EmpOrderInit[s];
				Upper_temp[s-1] = EmpOrderInit[s];
			}
			y_temp[J+I-1] = y[J-1];
			sigma_temp[J+I-1] = sigma[J-1];
			Lower_temp[J+I-1] = EmpOrderInit[J-1];
			Upper_temp[J+I-1] = EmpOrderInit[J-1];
						
			//# Calculate the SCI for ranks
			//#res = ApproximatePartitionWrongOrder(y_temp, sigma_temp, EmpOrder_temp,critFun, Intercept)
			OnlyBlockRanking_instance(y_temp, sigma_temp, crit, Lower_temp, Upper_temp, n);
			
			
			//# permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
			  Lower[s] = fmin(Lower[s], Lower_temp[s-1]);
			  Upper[s] = fmax(Upper[s], Upper_temp[s-1]);
			}
			Lower[J-1] = fmin(Lower[J-1], Lower_temp[J+I-1]);
			Upper[J-1] = fmax(Upper[J-1], Upper_temp[J+I-1]);
		 	I = I+1;
			
		}
	}
	
	if(SwapPerm == true)
	{
		if(trace == true) {Rcout<<"\n Apply "<<n<<" swap permutations \n";}// They do help some times
		for(int I = 0; I<n-1; I++)
		{
			if(trace == true) {Rcout<<I<<".";}
			for(int J = I+1; J<n; J++)
			{
				// Swap mu_i with mu_j
				y_temp[I] = y[J];
				y_temp[J] = y[I];
				sigma_temp[I] = sigma[J];
				sigma_temp[J] = sigma[I];
				Lower_temp[I] = EmpOrderInit[J];
				Lower_temp[J] = EmpOrderInit[I];	
				Upper_temp[I] = EmpOrderInit[J];
				Upper_temp[J] = EmpOrderInit[I];
				
				//Rcout<<"\n (I,J)=("<<I<<","<<J<<").";
				OnlyBlockRanking_instance(y_temp, sigma_temp, crit, Lower_temp, Upper_temp, n);
				
				Lower[I] = fmin(Lower[I], Lower_temp[J]);
				Lower[J] = fmin(Lower[J], Lower_temp[I]);
				Upper[I] = fmax(Upper[I], Upper_temp[J]);
				Upper[J] = fmax(Upper[J], Upper_temp[I]);
				
				// Swap back
				y_temp[I] = y[I];
				y_temp[J] = y[J];
				sigma_temp[I] = sigma[I];
				sigma_temp[J] = sigma[J];
				Lower_temp[I] = EmpOrderInit[I];
				Lower_temp[J] = EmpOrderInit[J];
				Upper_temp[I] = EmpOrderInit[I];
				Upper_temp[J] = EmpOrderInit[J];
			}
		}
	}

	if(trace == true) {Rcout<<"\n Apply "<<NbOfPermut<<" random permutations \n";}
	
	for(int I = 1; I<=NbOfPermut; I++)//#row index
	{
		if((trace == true) && (I % 100 == 0)) {Rcout<<I<<".";}
					
		//# Permute
		RandPermut = sample(RandPermutInit,n);//clone(RandPermutInit);
		//std::random_shuffle(RandPermut.begin(), RandPermut.end());
		for(int s = 0; s<n; s++)
		{
			y_temp[s] = y[RandPermut[s]]; 
		//	Rcout<<y_temp[s]<<",";
			sigma_temp[s] = sigma[RandPermut[s]];
			Lower_temp[s] = EmpOrderInit[RandPermut[s]];
			Upper_temp[s] = EmpOrderInit[RandPermut[s]];
		}
		//Rcout<<"\n";
		OnlyBlockRanking_instance(y_temp, sigma_temp, crit, Lower_temp, Upper_temp, n);

		//# Permute the ranks
		for(int s = 0; s<n; s++)
		{
			Lower[RandPermut[s]] = fmin(Lower[RandPermut[s]], Lower_temp[s]);
			Upper[RandPermut[s]] = fmax(Upper[RandPermut[s]], Upper_temp[s]);
		}
		//if(trace == true) {Rcout<<"\n";}
	}
	

	
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	return CIs;
}


// ----------------------------------------------------------
// ...........................................
//................ Bracketing  ...................
// ...........................................
// ----------------------------------------------------------
/*double critFun(const int& x, const double& Slop)
{
	return Slop*x;
}*/
void WhichBounds(const NumericVector& y, const int& I, const int& J, int& minInd, int& maxInd)
{ // Determine the position of the min and max values
	minInd = 1; maxInd = J-I+1;
	int minY = y[I], maxY = y[J];
	for(int i = 1; i<J-I+1; i++)
	{
		if(y[i+I]<minY)
		{
			minY = y[i+I];
			minInd = i+1;
		}
		if(y[i+I]>maxY)
		{
			maxY = y[i+I];
			maxInd = i+1;
		}
	}
}

void IndividContribs(const NumericVector& y_temp, const NumericVector& sigma_temp, double**& LogL, const int& K, const int& L, const double& Binf, const double& Bsup, double**& IndividContribBlock, double***& AverageBlock, const double& Slop, const double& Intercept, const int& n)
{
	// The contribution of a block of means should be taken only when it is negative.
	// If it is positive, then a better option is mu_i<...<mu_j which contributes by zero.
	// However, this is only valid if y_i<...<y_j which is not the case with partially correctly ordered hypotheses and needs to be checked using Order_ij.
	// Here the block is mu_{K}, mu_{K+1}, ..., mu_{L}. The contribution matrixes count between 1 and L-K+1
	double MinContribBlock;
	int indMinContrib;//, minYBlock_ind, maxYBlock_ind;
	bool Order_ij;
	// Initialize the matrices of averages and contributions if previously used.
	for(int i = 1; i<= L-K+1; i++)
	{
		for(int j = 2; j<= L-K+1; j++)
		{
			IndividContribBlock[i][j] = LogL[i+K-1][j+K-1] - (j-i)*Slop;// We start by mu_i=...=mu_j. 
			AverageBlock[i][j][1] = Summation(y_temp, sigma_temp, K-1+i, K-1+j);
			AverageBlock[i][j][2] = AverageBlock[i][j][1];
		}
		AverageBlock[i][i][1] = y_temp[i+K-1];
		AverageBlock[i][i][2] = y_temp[i+K-1];
	}
	//### If the Block has several sublocks, I need to keep track of the smallest average and the largest one !
	if(L-K+1>=2)
	{
		for(int j = 2; j<=(L-K+1); j++)
		{
			for(int i = j-1; i>=1; i--)
			{
			 	
				
				indMinContrib = 0; MinContribBlock = IndividContribBlock[i][j];
				for(int s = 1; s<=(j-i); s++)
				{
				  if(AverageBlock[i][i+s-1][2] <= AverageBlock[i+s][j][1] && AverageBlock[i][i+s-1][1]>=Binf && AverageBlock[i+s][j][2]<=Bsup)// # order is respected and no PAVA
				  {
					
					if(IndividContribBlock[i][i+s-1]+IndividContribBlock[i+s][j] < IndividContribBlock[i][j])
					{
					  //IndividContribBlock[i][j] = MinContribBlock;
					  IndividContribBlock[i][j] = IndividContribBlock[i][i+s-1]+IndividContribBlock[i+s][j];
					  AverageBlock[i][j][1] = AverageBlock[i][i+s-1][1];
					  AverageBlock[i][j][2] = AverageBlock[i+s][j][2];
					  indMinContrib = s;
					}
				  }
				}
				
				if(indMinContrib == 0)
				{ 
					// # block contributes as mu_i<...<mu_j OR mu_i=...=mu_j
					// If y_i<...<y_j, then we have to choose between one of the two options. Otherwise, it is only the second one.
					if(IndividContribBlock[i][j] < 0)
					{
						// Here the configuration mu_i=...=mu_j contributes better anyway. Do not do any thing
						continue;
					}
					Order_ij = true; // y_i<...<y_j
					for(int s = 1; s<=(j-i); s++)
					{
						if(y_temp[i+s+K-2]>y_temp[i+s+K-1])
						{
							Order_ij = false;
						}
					}
					if(Order_ij == true && IndividContribBlock[i][j] > 0)
					{
						// We have to choose. If mu_i=...=mu_j contributes positively, then mu_i<...<mu_j is the choice.
						AverageBlock[i][j][1] = y_temp[i+K-1];
						AverageBlock[i][j][2] = y_temp[j+K-1];
						IndividContribBlock[i][j] = 0.0;
					}
					
				}
				
			}
		}
	}
}
 


void ApproximatePartitionWrongOrderFast(const NumericVector& y_temp, const NumericVector& sigma_temp, const IntegerVector& EmpOrder_temp, int* &Lower, int*& Upper, double**& LogL, double***& AverageBlock, double***& AverageBlock_temp, double**& IndividContribBlock, double**& IndividContribBlock_temp, const double& Slop, const double& Intercept, const int& n, const double& maxY, const double& minY)
{
	double Contrib, ContribR, UpAv, LowAv, AvCurrBlock;
	
	//#### If the Block has several sublocks, I need to keep track of the smallest average and the largest one !
	for(int j = 2 ; j<=n; j++)
	{
		for(int i = j-1; i>=1; i--)
		{
			LogL[i][j] = LogLikelihood(y_temp, sigma_temp, i, j);
		}
	}
	
	//res = IndividContribs(y_temp,sigma_temp,LogL, critFun,1,n,minY,maxY,TRUE)
	IndividContribs(y_temp, sigma_temp, LogL, 1, n, minY, maxY, IndividContribBlock, AverageBlock, Slop, Intercept, n);

	//# Test the blocks by adding hypotheses to the left and to the right of it
	for(int i = 1; i<=n; i++)
	{
		Lower[i] = EmpOrder_temp[i];
		Upper[i] = EmpOrder_temp[i];
	}
	//# Treat the case of mu_1=...=mu_j + hypothesis
	for(int j = n-1; j>=2; j--)
	{ 
		AvCurrBlock = Summation(y_temp, sigma_temp, 1, j);
		if(AvCurrBlock<=AverageBlock[j+1][n][1])
		{
			if(LogL[1][j] - (j-1)*Slop + IndividContribBlock[j+1][n] - Intercept < 0)// # Block [1,j] is accepted. Update the CIs.
			{
			  	for(int s = 1; s<=j; s++)
			  	{
					Lower[s] = 1;
			    	//Upper[s] = fmax(Upper[s], j);
			    	if(Upper[s]<j) Upper[s] = j;
				}
				break;
			}
		}
		else //# The minimally contributing block does not respect the ordering. We need to find another one with slightly higher contrib but respects the order
		{
		    IndividContribs(y_temp,sigma_temp,LogL,j+1,n,AvCurrBlock,maxY,IndividContribBlock_temp, AverageBlock_temp, Slop, Intercept, n);
		    LowAv = AverageBlock_temp[1][n-j][1]; Contrib = IndividContribBlock_temp[1][n-j];
		    if(AvCurrBlock < LowAv)
		    {	  
			  if(LogL[1][j] - (j-1)*Slop + Contrib - Intercept < 0) //# Block [1,j] is accepted. Update the CIs.
			  {
			  	for(int s = 1; s<=j; s++)
			  	{
					Lower[s] = 1;
			    	//Upper[s] = fmax(Upper[s], j);
			    	if(Upper[s]<j) Upper[s] = j;
				}
			    break;
			  }
		    }
		}
	}
	//# Treat the case of hypothesis + mu_i=...=mu_n
	for(int i = 2; i <= n-1; i++)
	{
		AvCurrBlock = Summation(y_temp, sigma_temp, i, n);
		if(AverageBlock[1][i-1][2] <= AvCurrBlock)
		{
			if(LogL[i][n] - (n-i)*Slop + IndividContribBlock[1][i-1] - Intercept < 0)// # Block [i,n] is accepted. Update the CIs.
			{
				for(int s = i; s <= n; s++)
				{
					//Lower[s] = fmin(Lower[s],i);
					if(Lower[s]>i) Lower[s] = i;
			  		Upper[s] = n;
				}
			  
			  break;
			}
		}
		else
		{
		    IndividContribs(y_temp,sigma_temp,LogL,1,i-1,minY,AvCurrBlock,IndividContribBlock_temp, AverageBlock_temp, Slop, Intercept,n);
		    UpAv = AverageBlock_temp[1][i-1][2]; Contrib = IndividContribBlock_temp[1][i-1];
		    if(AvCurrBlock > UpAv)
		    {
				if(LogL[i][n] - (n-i)*Slop + Contrib - Intercept < 0)// # Block [i,n] is accepted. Update the CIs.
				{
				    for(int s = i; s <= n; s++)
					{
						//Lower[s] = fmin(Lower[s],i);
						if(Lower[s]>i) Lower[s] = i;
				  		Upper[s] = n;
					}
				    break;
				}
			
		    }
		}
	}
	//# Treat the case of hypothesis + mu_i=...=mu_j + hypothesis
	for(int i = 2; i<=n-2; i++)
	{
		for(int j = n-1; j >= i+1; j--)
		{
			AvCurrBlock = Summation(y_temp, sigma_temp, i, j);
			if(AverageBlock[1][i-1][2] <= AvCurrBlock && AvCurrBlock <= AverageBlock[j+1][n][1])
			{
				if(LogL[i][j] - (j-i)*Slop + IndividContribBlock[1][i-1] + IndividContribBlock[j+1][n] - Intercept < 0)// # Block [i,j] is accepted. Update the CIs.
				{
					for(int s = i; s <= j; s++)
					{
						//Lower[s] = fmin(Lower[s],i);
						if(Lower[s]>i) Lower[s] = i;
				  		//Upper[s] = fmax(Upper[s],j);
				  		if(Upper[s]<j) Upper[s] = j;
					}
					break;
				}
			}
			else
			{
			    IndividContribs(y_temp,sigma_temp,LogL,1,i-1,minY,AvCurrBlock,IndividContribBlock_temp, AverageBlock_temp, Slop, Intercept, n);
			    UpAv = AverageBlock_temp[1][i-1][2]; Contrib = IndividContribBlock_temp[1][i-1];
			    IndividContribs(y_temp,sigma_temp,LogL,j+1,n,AvCurrBlock,maxY,IndividContribBlock_temp, AverageBlock_temp, Slop, Intercept, n);
			    LowAv = AverageBlock_temp[1][n-j][1]; ContribR = IndividContribBlock_temp[1][n-j] ;
			    if(AvCurrBlock > UpAv && AvCurrBlock < LowAv)//# THE CONDITION MUST BE ADAPTED TO THE MINIMUM CONTRIB BLOCK if it is not possible to bind with the left part, then move on !
			    {
					if(LogL[i][j] - (j-i)*Slop + ContribR  + Contrib - Intercept < 0) //# Block [i,j] is accepted. Update the CIs.
					{
						for(int s = i; s <= j; s++)
						{
							//Lower[s] = fmin(Lower[s],i);
							if(Lower[s]>i) Lower[s] = i;
				  			//Upper[s] = fmax(Upper[s], j);
				  			if(Upper[s]<j) Upper[s] = j;
						}
					  	break;
					}
			    }
			}
			
		}
	}
}


// [[Rcpp::export]]
NumericMatrix ApproximatePartitionPermutations(NumericVector yInit, NumericVector sigmaInit, IntegerVector LowerInit, IntegerVector UpperInit, int n, double Slop, double Intercept, double minY, double maxY, bool trace, const bool SwapPerm, const int NbOfPermut)
{
	IntegerVector Lower(n+1), Upper(n+1), RandPermutInit, RandPermut;
	RandPermutInit = seq(1,n);
	
	for(int i = 1; i<=n; i++)
	{
		Lower[i] = LowerInit[i-1];
		Upper[i] = UpperInit[i-1];
	}
	// Define all matrices once and for all and avoid all allocations and deletion.
 	IntegerVector EmpOrder = seq(0,n);
 	NumericVector y(n+1), y_temp(n+1), sigma(n+1), sigma_temp(n+1); 
	IntegerVector EmpOrder_temp(n+1);
 	// Rescale the index from 1 in R to 1 in C++ by ignoring the first element.
	for(int i = 1; i<=n; i++)
	{
		y[i] = yInit[i-1];
		sigma[i] = sigmaInit[i-1];
	}
	
	int* Lower_temp = new int[n+1]; int* Upper_temp = new int[n+1];
	// I AM STILL USING THE INDEX + 1
	double** IndividContribBlock = new double*[n+1];
	double** IndividContribBlock_temp = new double*[n+1];
	double*** AverageBlock = new double**[n+1]; // array(dim = c(L-K+1,L-K+1,2)) # a block with a single observation has a mean equal to it.
	double*** AverageBlock_temp = new double**[n+1]; // array(dim = c(L-K+1,L-K+1,2)) # a block with a single observation has a mean equal to it.
	for(int i = 0; i<=n; i++)
	{
		AverageBlock[i] = new double*[n+1];
		AverageBlock_temp[i] = new double*[n+1];
		IndividContribBlock[i] = new double[n+1];
		IndividContribBlock_temp[i] = new double[n+1];
		for(int j = 0; j<=n; j++)
		{
			AverageBlock[i][j] = new double[3];
			AverageBlock_temp[i][j] = new double[3];
			AverageBlock[i][j][0] = 0; AverageBlock[i][j][1] = 0; AverageBlock[i][j][2] = 0;
			AverageBlock_temp[i][j][0] = 0; AverageBlock_temp[i][j][1] = 0; AverageBlock_temp[i][j][2] = 0;
		}
	}
	
	
	//# Calculate the matrix of Log-likelihood contributions
	double** LogL = new double*[n+1];
	for(int i = 0; i<=n; i++)
	{
		LogL[i] = new double[n+1];
		for(int j = 0; j<=n; j++)
		{
			LogL[i][j] = 0;
		}
	}
	
	if(trace == true) {Rcout<<"\n Forward permutations ("<<n-1<<" permutations)\n";}
	//### part1: permute (i,i+1,i+2,...)
	for(int I = 1; I<=n-1; I++)//#row index
	{
		if(trace == true) {Rcout<<I<<".";}
		for(int J = 1; J <= n-I; J++)//#column index
		{
			//y_temp = y; sigma_temp = sigma; POINTER COPYING !!!
			//EmpOrder_temp = EmpOrder;
			for(int s = 1; s<=n; s++)
			{
				y_temp[s] = y[s]; 
				sigma_temp[s] = sigma[s];
				EmpOrder_temp[s] = EmpOrder[s];
			}
			
			//# Permute
			for(int s = J; s<=J+I-1; s++)
			{
				y_temp[s+1] = y[s];
				sigma_temp[s+1] = sigma[s];
				EmpOrder_temp[s+1] = EmpOrder[s]; 
			}
			y_temp[J] = y[J+I];
			sigma_temp[J] = sigma[J+I];
			EmpOrder_temp[J] = EmpOrder[J+I];	
			
			//# Calculate the SCI for ranks
			//#res = ApproximatePartitionWrongOrder(y_temp, sigma_temp, EmpOrder_temp,critFun, Intercept)
			ApproximatePartitionWrongOrderFast(y_temp, sigma_temp, EmpOrder_temp, Lower_temp, Upper_temp, LogL, AverageBlock, AverageBlock_temp, IndividContribBlock, IndividContribBlock_temp, Slop, Intercept, n, maxY, minY);
		
			//# Permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
				Lower[s] = fmin(Lower[s], Lower_temp[s+1]);
				Upper[s] = fmax(Upper[s], Upper_temp[s+1]);
			}
			Lower[J+I] = fmin(Lower[J+I], Lower_temp[J]);
			Upper[J+I] = fmax(Upper[J+I], Upper_temp[J]);
			
			//if(sum(Lower==rep(1,n) & Upper==rep(n,n)) == n) return(list(Lower=Lower,Upper=Upper))
		}
	}
	
	if(trace == true) {Rcout<<"\n Applying backward permutations.\n";}
	int I;
	//### Part 2: permute (i+4,i+3,..)
	for(int J = 1; J<=n-1; J++)//#column index
	{
		if(trace == true) {Rcout<<J<<".";}
		I = 2;
		while(I<=(n-J))//#row index
		{
			//y_temp = y; sigma_temp = sigma; EmpOrder_temp = EmpOrder;
			for(int s = 1; s<=n; s++)
			{
				y_temp[s] = y[s]; 
				sigma_temp[s] = sigma[s];
				EmpOrder_temp[s] = EmpOrder[s];
			}
			//# Permute
			for(int s = J; s<= J+I-1; s++)
			{
			  y_temp[s] = y[s+1];
			  sigma_temp[s] = sigma[s+1];
			  EmpOrder_temp[s] = EmpOrder[s+1];
			}
			y_temp[J+I] = y[J];
			sigma_temp[J+I] = sigma[J];
			EmpOrder_temp[J+I] = EmpOrder[J];
			
			//# Calculate the SCI for ranks
			//#res = ApproximatePartitionWrongOrder(y_temp, sigma_temp, EmpOrder_temp,critFun, Intercept)
			ApproximatePartitionWrongOrderFast(y_temp, sigma_temp, EmpOrder_temp, Lower_temp, Upper_temp, LogL, AverageBlock, AverageBlock_temp, IndividContribBlock, IndividContribBlock_temp, Slop, Intercept, n, maxY, minY);
			//# permute the ranks
			
			for(int s = J; s<=J+I-1; s++)
			{
			  Lower[s+1] = fmin(Lower[s+1], Lower_temp[s]);
			  Upper[s+1] = fmax(Upper[s+1], Upper_temp[s]);
			}
			Lower[J] = fmin(Lower[J], Lower_temp[J+I]);
			Upper[J] = fmax(Upper[J], Upper_temp[J+I]);
		 	I = I+1;
			//if(sum(Lower==rep(1,n) & Upper==rep(n,n)) == n) return(list(Lower=Lower,Upper=Upper))
			//if(Lower[10] == 1) Rcout<<"I="<<I<<", J="<<J<<"\n";
		}
	}

	if(SwapPerm == true) 
	{
		// They do help some times
		if(trace == true) {Rcout<<"\n Apply "<<n<<" swap permutations \n";}
		for(int I = 1; I<=n-1; I++)
		{
			if(trace == true) {Rcout<<I<<".";}
			for(int J = I+1; J<=n; J++)
			{
				// Swap mu_i with mu_j
				y_temp[I] = y[J];
				y_temp[J] = y[I];
				sigma_temp[I] = sigma[J];
				sigma_temp[J] = sigma[I];
				EmpOrder_temp[I] = EmpOrder[J];
				EmpOrder_temp[J] = EmpOrder[I];
				//Rcout<<"\n (I,J)=("<<I<<","<<J<<").";
				ApproximatePartitionWrongOrderFast(y_temp, sigma_temp, EmpOrder_temp, Lower_temp, Upper_temp, LogL, AverageBlock, AverageBlock_temp, IndividContribBlock, IndividContribBlock_temp, Slop, Intercept, n, maxY, minY);
				
				Lower[I] = fmin(Lower[I], Lower_temp[J]);
				Lower[J] = fmin(Lower[J], Lower_temp[I]);
				Upper[I] = fmax(Upper[I], Upper_temp[J]);
				Upper[J] = fmax(Upper[J], Upper_temp[I]);
				
				// Swap back
				y_temp[I] = y[I];
				y_temp[J] = y[J];
				sigma_temp[I] = sigma[I];
				sigma_temp[J] = sigma[J];
				EmpOrder_temp[I] = EmpOrder[I];
				EmpOrder_temp[J] = EmpOrder[J];
			}
		}
	}
	if(trace == true) {Rcout<<"\n Apply "<<NbOfPermut<<" random permutations \n";}
	
	for(int I = 1; I<=NbOfPermut; I++)//#row index
	{
		if((trace == true) && (I % 10000 == 0)) {Rcout<<I<<".";}
					
		//# Permute
		RandPermut = sample(n,n);//clone(RandPermutInit);
		//std::random_shuffle(RandPermut.begin(), RandPermut.end());
		for(int s = 1; s<=n; s++)
		{
			y_temp[s] = y[RandPermut[s-1]]; 
			sigma_temp[s] = sigma[RandPermut[s-1]];
			EmpOrder_temp[s] = EmpOrder[RandPermut[s-1]];
		}
		ApproximatePartitionWrongOrderFast(y_temp, sigma_temp, EmpOrder_temp, Lower_temp, Upper_temp, LogL, AverageBlock, AverageBlock_temp, IndividContribBlock, IndividContribBlock_temp, Slop, Intercept, n, maxY, minY);

		//# Permute the ranks
		for(int s = 1; s<=n; s++)
		{
			Lower[RandPermut[s-1]] = fmin(Lower[RandPermut[s-1]], Lower_temp[s]);
			Upper[RandPermut[s-1]] = fmax(Upper[RandPermut[s-1]], Upper_temp[s]);
		}
		
	}
	
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i+1];
		CIs[n+i] = Upper[i+1];
	}
	
	// Free some space.
	for (int i = 0; i<=n; i++)
	{
		delete[] LogL[i];
		delete[] IndividContribBlock[i];
		delete[] IndividContribBlock_temp[i];
		for(int j = 0; j<=n; j++)
		{
			delete[] AverageBlock[i][j];
			delete[] AverageBlock_temp[i][j];
		}
		delete[] AverageBlock[i];
		delete[] AverageBlock_temp[i];
	}
	delete[] LogL;
	delete[] IndividContribBlock;
	delete[] IndividContribBlock_temp;
	delete[] AverageBlock;
	delete[] AverageBlock_temp;
	delete[] Lower_temp;
	delete[] Upper_temp;
	
	return CIs;

}




// This function is the same as the one written in the R code, but I'm using it here only to make things run faster.
NumericMatrix TukeyProc(const NumericVector& y, const NumericVector& sigma, const double& qq, const int& n)
{
  	double stat;
  	int countPos, countNeg;
  	NumericMatrix ResCIs(n,2);
  	for(int j = 0; j < n; j++)
  	{
  		countPos = 0; 
  		countNeg = 0;
		for(int i = 0; i < n; i++)
		{
			stat = (y[j]-y[i])/sqrt(sigma[j]*sigma[j]+sigma[i]*sigma[i]);
			//if(j == 3) Rcout<<y[j]<<"-"<<y[i]<<", sigmas = "<<sigma[j]<<","<<sigma[i]<<". Stat = "<<stat<<",";
			if(stat>qq)
			{
				countPos++;
			}
			if(stat<(-qq))
			{
				countNeg++;
			}
		}
		//Rcout<<"\n";
		ResCIs(j,0)=1+countPos;
		ResCIs(j,1)=n-countNeg;
  	}
	return ResCIs;
}

// [[Rcpp::export]]
NumericMatrix TukeyRankingLevelEqSigRescaled(NumericVector y, NumericVector sigma, NumericMatrix crit, NumericMatrix SampleCoverage, int MM, int n, int NbOfPermut, double alpha, int gridSize, bool trace)
{
	int powN2 = 1;
	for(int i = 0; i<n; i++)
	{
		powN2 *=2; 
	}
	if(trace == true) Rcout<<"Initializing variables...\n";
	int*** ResCIsMat = new int**[gridSize]; // alpha values range in 0.05, ..., 0.4
	int*** ResCIsMat_temp = new int**[gridSize];
	int*** ResCIsGridMat = new int**[gridSize];
	int* AlphaRescaled = new int[powN2];
	for(int i = 0; i<gridSize; i++)
	{
		ResCIsMat[i] = new int*[n];
		ResCIsMat_temp[i] = new int*[n];
		ResCIsGridMat[i] = new int*[n];
		for(int j = 0; j<n; j++)
		{
			ResCIsMat[i][j] = new int[2*MM];
			ResCIsMat_temp[i][j] = new int[2*MM];
			ResCIsGridMat[i][j] = new int[2];
		}
	}
	NumericMatrix ResCIs(n,2);
	if(trace == true) Rcout<<"Calculating SCIs for the ranks of a sample.\n";
	for(int a = 0; a<gridSize; a++)
	{
		for(int s = 0; s<MM; s++)
		{					
			// Caclulate simultaneous CIs for ranks using sample i
			ResCIs = TukeyProc(SampleCoverage.column(s), sigma, crit[a], n);
			for(int i = 0; i<n; i++)
			{
				ResCIsMat[a][i][2*s] = ResCIs(i,0)-1; // Even for Lower CIs
				ResCIsMat[a][i][2*s+1] = ResCIs(i,1)-1; // Odd for Upper CIs	
			}
		}
		ResCIs = TukeyProc(y, sigma, crit[a], n);

		for(int i = 0; i<n; i++)
		{
			ResCIsGridMat[a][i][0] = ResCIs(i,0)-1;
			ResCIsGridMat[a][i][1] = ResCIs(i,1)-1;
		}
	}

	//NumericVector y_temp(n), sigma_temp(n);
	IntegerVector RandPermutInit, RandPermut;
	RandPermutInit = seq(0,n-1);
	// Calculate the C_n^k matrix in order to use it in the UnrankCombin function
	unsigned long long int ** CnkMat = new unsigned long long int*[n];
	for (int i = 0; i<n; i++)
	{
		CnkMat[i] = new unsigned long long int[n];
		CnkMat[i][i] = 1;
		for (int j = 0; j<i; j++)
		{
			CnkMat[i][j] = binomialCoeff(i, j);
			CnkMat[j][i] = 0;
		}
	}

	// Definition and initialization of the vector of ranks
	IntegerVector Lower(n), Lower_temp(n), EmpOrderInit(n);
	IntegerVector Upper(n), Upper_temp(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
		EmpOrderInit[i] = i;
	}
	// Calculate first the adjustment for the hypothesis mu_1=...=mu_n
	RescalTopHpo(ResCIsMat, ResCIsGridMat, Lower, Upper, EmpOrderInit, n, MM, gridSize, alpha);
	
	// Calculate CIs for ranks using the general proc without any permutations (correctly ordered hypotheses)
	if(trace == true) Rcout<<"Perform the correctly ordered hypotheses...\n";
	PartitioningRankingGeneralProcInit(ResCIsMat, ResCIsGridMat, AlphaRescaled, Lower, Upper, EmpOrderInit, CnkMat, crit, n, MM, gridSize, alpha, trace);
	
	bool CheckTrivialCIs = true;
	if(trace == true) {Rcout<<"\n Forward permutations ("<<n-1<<" sets of permutations)\n";}
	//### part1: permute (i,i+1,i+2,...)
	for(int I = 1; I<=n-1; I++)//#row index
	{
		if(trace == true) {Rcout<<I<<".";}
		for(int J = 1; J <= n-I; J++)//#column index
		{
			CheckTrivialCIs = true;
			for(int i = 0; i<n; i++)
			{
				if(Lower[i] != 0 && Upper[i] != n-1)
				{
					CheckTrivialCIs = false;
					break;
				}
			}
			if(CheckTrivialCIs == true) break;

			for(int s = 0; s<n; s++)
			{
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}			
			}
			
			//# Permute
			for(int s = J; s<=J+I-1; s++)
			{
				Lower_temp[s] = EmpOrderInit[s-1]; 
				Upper_temp[s] = EmpOrderInit[s-1];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s-1][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s-1][2*j+1];
					}
				}			
			}
			Lower_temp[J-1] = EmpOrderInit[J+I-1];	
			Upper_temp[J-1] = EmpOrderInit[J+I-1];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][J-1][2*j] = ResCIsMat[i][J+I-1][2*j];
					ResCIsMat_temp[i][J-1][2*j+1] = ResCIsMat[i][J+I-1][2*j+1];
				}
			}
		
			//# Calculate the SCI for ranks
			PartitioningRankingGeneralProc(ResCIsMat_temp, ResCIsGridMat, AlphaRescaled, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);
			
			//# Permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
				Lower[s-1] = fmin(Lower[s-1], Lower_temp[s]);
				Upper[s-1] = fmax(Upper[s-1], Upper_temp[s]);
			}
			Lower[J+I-1] = fmin(Lower[J+I-1], Lower_temp[J-1]);
			Upper[J+I-1] = fmax(Upper[J+I-1], Upper_temp[J-1]);
		}
	}
	
	if(trace == true) {Rcout<<"\n Backward permutations ("<<n-1<<" sets of permutations)\n";}
	int I;
	//### Part 2: permute (i+4,i+3,..)
	for(int J = 1; J<=n-1; J++)//#column index
	{
		if(trace == true) {Rcout<<J<<".";}
		I = 2;
		while(I<=(n-J))//#row index
		{
			CheckTrivialCIs = true;
			for(int i = 0; i<n; i++)
			{
				if(Lower[i] != 0 && Upper[i] != n-1)
				{
					CheckTrivialCIs = false;
					break;
				}
			}
			if(CheckTrivialCIs == true) break;

			for(int s = 0; s<n; s++)
			{
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}
			}
			//# Permute
			for(int s = J; s<= J+I-1; s++)
			{
			  	Lower_temp[s-1] = EmpOrderInit[s];
				Upper_temp[s-1] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s-1][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s-1][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}
			}
			Lower_temp[J+I-1] = EmpOrderInit[J-1];
			Upper_temp[J+I-1] = EmpOrderInit[J-1];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][J+I-1][2*j] = ResCIsMat[i][J-1][2*j];
					ResCIsMat_temp[i][J+I-1][2*j+1] = ResCIsMat[i][J-1][2*j+1];
				}
			}
						
			PartitioningRankingGeneralProc(ResCIsMat_temp, ResCIsGridMat, AlphaRescaled, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);
			
			//# permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
			  Lower[s] = fmin(Lower[s], Lower_temp[s-1]);
			  Upper[s] = fmax(Upper[s], Upper_temp[s-1]);
			}
			Lower[J-1] = fmin(Lower[J-1], Lower_temp[J+I-1]);
			Upper[J-1] = fmax(Upper[J-1], Upper_temp[J+I-1]);
		 	I = I+1;
		}
	}
	
	if(trace == true) {Rcout<<"\n Applying "<<NbOfPermut<<" Random permutations.\n";}
	
	for(int I = 1; I<=NbOfPermut; I++)//#row index
	{
		CheckTrivialCIs = true;
		for(int i = 0; i<n; i++)
		{
			if(Lower[i] != 0 && Upper[i] != n-1)
			{
				CheckTrivialCIs = false;
				break;
			}
		}
		if(CheckTrivialCIs == true) break;

		if((trace == true) && (I % 500000 == 0)) {Rcout<<I<<".";}
		//# Permute
		RandPermut = sample(RandPermutInit,n);//clone(RandPermutInit);
		//std::random_shuffle(RandPermut.begin(), RandPermut.end());
		for(int s = 0; s<n; s++)
		{
			Lower_temp[s] = EmpOrderInit[RandPermut[s]];
			Upper_temp[s] = EmpOrderInit[RandPermut[s]];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][RandPermut[s]][2*j];
					ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][RandPermut[s]][2*j+1];
				}
			}
		}
		PartitioningRankingGeneralProc(ResCIsMat_temp, ResCIsGridMat, AlphaRescaled, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);

		//# Permute the ranks
		for(int s = 0; s<n; s++)
		{
			Lower[RandPermut[s]] = fmin(Lower[RandPermut[s]], Lower_temp[s]);
			Upper[RandPermut[s]] = fmax(Upper[RandPermut[s]], Upper_temp[s]);
		}
	}
	
	Rcout<<"\n";
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] CnkMat[i];
	}
	delete[] CnkMat;
	for(int a = 0; a<gridSize; a++)
	{
		for(int i = 0; i<n; i++)
		{
			delete[] ResCIsMat[a][i];
			delete[] ResCIsMat_temp[a][i];
			delete[] ResCIsGridMat[a][i];
		}
		delete[] ResCIsMat[a];
		delete[] ResCIsMat_temp[a];
		delete[] ResCIsGridMat[a];
	}
	delete[] ResCIsMat;
	delete[] ResCIsMat_temp;
	delete[] ResCIsGridMat;
	delete[] AlphaRescaled;
	return CIs;
}
void PartitioningRankingGeneralProcTuk(int***& ResCIsMat, int***& ResCIsGridMat, IntegerVector& Lower, IntegerVector& Upper, const IntegerVector& EmpOrderInit, unsigned long long int**& CnkMat, const NumericMatrix& crit, const int& n, const int& MM, const int& gridSize, const double& alpha, const bool& trace)
{
	double* coveragePartition = new double[gridSize];
	int* InqPosi = new int[n];
	IntegerVector Lower_temp(n), Upper_temp(n), EmpOrder = seq(0,n-1);
	bool coverage;
	double minCoverage;
	int alphaRescaled;
	for (int l = 1; l <= n - 2; l++)
	{
		unsigned long long int m = CnkMat[n - 1][l];
		for (unsigned long long int c = 0; c<m; c++)
		{
			UnrankCombin(InqPosi, c, l, CnkMat);
			Lower_temp = clone(EmpOrder);
			Upper_temp = clone(EmpOrder);
			RankUpdate(Lower_temp, Upper_temp, InqPosi, l, n);

			coveragePartition[0] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, 0)-(1-alpha));
			minCoverage = coveragePartition[0]; alphaRescaled = 0;
			for(int a = 1; a<gridSize; a++)
			{
				coveragePartition[a] = fabs(PartitionCoverage(ResCIsMat, Lower_temp, Upper_temp, n, MM, a)-(1-alpha));
				if(coveragePartition[a]<minCoverage)
				{
					alphaRescaled = a;
					minCoverage = coveragePartition[a];
				}
			}
			//AlphaRescaled[(l-1)*CnkMat[n - 1][l-1] + c] = alphaRescaled;
			coverage = true;
			for(int i = 0; i<n; i++)
			{
				if(ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][0]>Lower_temp[i] || ResCIsGridMat[alphaRescaled][EmpOrderInit[i]][1]<Upper_temp[i])
				{
					coverage = false;
					break;
				}
			}
			if(coverage == true)
			{
				for(int s = 0 ; s<n; s++)
				{
					
					Lower[s] = fmin(Lower[s], Lower_temp[s]);
					Upper[s] = fmax(Upper[s], Upper_temp[s]);
					
				}
			}
		}
		
	}
	delete[] InqPosi;
	delete[] coveragePartition;
}

// [[Rcpp::export]]
NumericMatrix TukeyRankingLevelUneqSigRescaled(NumericVector y, NumericVector sigma, NumericMatrix crit, NumericMatrix SampleCoverage, int MM, int n, int NbOfPermut, double alpha, int gridSize, bool trace)
{
	/*int powN2 = 1;
	for(int i = 0; i<n; i++)
	{
		powN2 *=2; 
	}*/
	if(trace == true) Rcout<<"Initializing variables...\n";
	int*** ResCIsMat = new int**[gridSize]; // alpha values range in 0.05, ..., 0.4
	int*** ResCIsMat_temp = new int**[gridSize];
	int*** ResCIsGridMat = new int**[gridSize];
	//int* AlphaRescaled = new int[powN2];
	for(int i = 0; i<gridSize; i++)
	{
		ResCIsMat[i] = new int*[n];
		ResCIsMat_temp[i] = new int*[n];
		ResCIsGridMat[i] = new int*[n];
		for(int j = 0; j<n; j++)
		{
			ResCIsMat[i][j] = new int[2*MM];
			ResCIsMat_temp[i][j] = new int[2*MM];
			ResCIsGridMat[i][j] = new int[2];
		}
	}
	NumericMatrix ResCIs(n,2);
	if(trace == true) Rcout<<"Calculating SCIs for the ranks of a sample.\n";
	for(int a = 0; a<gridSize; a++)
	{
		for(int s = 0; s<MM; s++)
		{					
			// Caclulate simultaneous CIs for ranks using sample i
			ResCIs = TukeyProc(SampleCoverage.column(s), sigma, crit[a], n);
			for(int i = 0; i<n; i++)
			{
				ResCIsMat[a][i][2*s] = ResCIs(i,0)-1; // Even for Lower CIs
				ResCIsMat[a][i][2*s+1] = ResCIs(i,1)-1; // Odd for Upper CIs	
			}
		}
		ResCIs = TukeyProc(y, sigma, crit[a], n);

		for(int i = 0; i<n; i++)
		{
			ResCIsGridMat[a][i][0] = ResCIs(i,0)-1;
			ResCIsGridMat[a][i][1] = ResCIs(i,1)-1;
		}
	}

	//NumericVector y_temp(n), sigma_temp(n);
	IntegerVector RandPermutInit, RandPermut;
	RandPermutInit = seq(0,n-1);
	// Calculate the C_n^k matrix in order to use it in the UnrankCombin function
	unsigned long long int ** CnkMat = new unsigned long long int*[n];
	for (int i = 0; i<n; i++)
	{
		CnkMat[i] = new unsigned long long int[n];
		CnkMat[i][i] = 1;
		for (int j = 0; j<i; j++)
		{
			CnkMat[i][j] = binomialCoeff(i, j);
			CnkMat[j][i] = 0;
		}
	}

	// Definition and initialization of the vector of ranks
	IntegerVector Lower(n), Lower_temp(n), EmpOrderInit(n);
	IntegerVector Upper(n), Upper_temp(n);
	for (int i = 0; i<n; i++)
	{
		Lower[i] = i;
		Upper[i] = i;
		EmpOrderInit[i] = i;
	}
	// Calculate first the adjustment for the hypothesis mu_1=...=mu_n
	RescalTopHpo(ResCIsMat, ResCIsGridMat, Lower, Upper, EmpOrderInit, n, MM, gridSize, alpha);
	
	// Calculate CIs for ranks using the general proc without any permutations (correctly ordered hypotheses)
	if(trace == true) Rcout<<"Perform the correctly ordered hypotheses...\n";
	PartitioningRankingGeneralProcTuk(ResCIsMat, ResCIsGridMat, Lower, Upper, EmpOrderInit, CnkMat, crit, n, MM, gridSize, alpha, trace);
	
	bool CheckTrivialCIs = true;
	if(trace == true) {Rcout<<"\n Forward permutations ("<<n-1<<" sets of permutations)\n";}
	//### part1: permute (i,i+1,i+2,...)
	for(int I = 1; I<=n-1; I++)//#row index
	{
		if(trace == true) {Rcout<<I<<".";}
		for(int J = 1; J <= n-I; J++)//#column index
		{
			CheckTrivialCIs = true;
			for(int i = 0; i<n; i++)
			{
				if(Lower[i] != 0 && Upper[i] != n-1)
				{
					CheckTrivialCIs = false;
					break;
				}
			}
			if(CheckTrivialCIs == true) break;

			for(int s = 0; s<n; s++)
			{
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}			
			}
			
			//# Permute
			for(int s = J; s<=J+I-1; s++)
			{
				Lower_temp[s] = EmpOrderInit[s-1]; 
				Upper_temp[s] = EmpOrderInit[s-1];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s-1][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s-1][2*j+1];
					}
				}			
			}
			Lower_temp[J-1] = EmpOrderInit[J+I-1];	
			Upper_temp[J-1] = EmpOrderInit[J+I-1];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][J-1][2*j] = ResCIsMat[i][J+I-1][2*j];
					ResCIsMat_temp[i][J-1][2*j+1] = ResCIsMat[i][J+I-1][2*j+1];
				}
			}
		
			//# Calculate the SCI for ranks
			PartitioningRankingGeneralProcTuk(ResCIsMat_temp, ResCIsGridMat, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);
			
			//# Permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
				Lower[s-1] = fmin(Lower[s-1], Lower_temp[s]);
				Upper[s-1] = fmax(Upper[s-1], Upper_temp[s]);
			}
			Lower[J+I-1] = fmin(Lower[J+I-1], Lower_temp[J-1]);
			Upper[J+I-1] = fmax(Upper[J+I-1], Upper_temp[J-1]);
		}
	}
	
	if(trace == true) {Rcout<<"\n Backward permutations ("<<n-1<<" sets of permutations)\n";}
	int I;
	//### Part 2: permute (i+4,i+3,..)
	for(int J = 1; J<=n-1; J++)//#column index
	{
		if(trace == true) {Rcout<<J<<".";}
		I = 2;
		while(I<=(n-J))//#row index
		{
			CheckTrivialCIs = true;
			for(int i = 0; i<n; i++)
			{
				if(Lower[i] != 0 && Upper[i] != n-1)
				{
					CheckTrivialCIs = false;
					break;
				}
			}
			if(CheckTrivialCIs == true) break;

			for(int s = 0; s<n; s++)
			{
				Lower_temp[s] = EmpOrderInit[s];
				Upper_temp[s] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}
			}
			//# Permute
			for(int s = J; s<= J+I-1; s++)
			{
			  	Lower_temp[s-1] = EmpOrderInit[s];
				Upper_temp[s-1] = EmpOrderInit[s];
				for(int i = 0; i<gridSize; i++)
				{
					for(int j = 0; j<MM; j++)
					{
						ResCIsMat_temp[i][s-1][2*j] = ResCIsMat[i][s][2*j];
						ResCIsMat_temp[i][s-1][2*j+1] = ResCIsMat[i][s][2*j+1];
					}
				}
			}
			Lower_temp[J+I-1] = EmpOrderInit[J-1];
			Upper_temp[J+I-1] = EmpOrderInit[J-1];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][J+I-1][2*j] = ResCIsMat[i][J-1][2*j];
					ResCIsMat_temp[i][J+I-1][2*j+1] = ResCIsMat[i][J-1][2*j+1];
				}
			}
						
			PartitioningRankingGeneralProcTuk(ResCIsMat_temp, ResCIsGridMat, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);
			
			//# permute the ranks
			for(int s = J; s<=J+I-1; s++)
			{
			  Lower[s] = fmin(Lower[s], Lower_temp[s-1]);
			  Upper[s] = fmax(Upper[s], Upper_temp[s-1]);
			}
			Lower[J-1] = fmin(Lower[J-1], Lower_temp[J+I-1]);
			Upper[J-1] = fmax(Upper[J-1], Upper_temp[J+I-1]);
		 	I = I+1;
		}
	}
	
	if(trace == true) {Rcout<<"\n Applying "<<NbOfPermut<<" Random permutations.\n";}
	
	for(int I = 1; I<=NbOfPermut; I++)//#row index
	{
		CheckTrivialCIs = true;
		for(int i = 0; i<n; i++)
		{
			if(Lower[i] != 0 && Upper[i] != n-1)
			{
				CheckTrivialCIs = false;
				break;
			}
		}
		if(CheckTrivialCIs == true) break;

		if((trace == true) && (I % 500000 == 0)) {Rcout<<I<<".";}
		//# Permute
		RandPermut = sample(RandPermutInit,n);//clone(RandPermutInit);
		//std::random_shuffle(RandPermut.begin(), RandPermut.end());
		for(int s = 0; s<n; s++)
		{
			Lower_temp[s] = EmpOrderInit[RandPermut[s]];
			Upper_temp[s] = EmpOrderInit[RandPermut[s]];
			for(int i = 0; i<gridSize; i++)
			{
				for(int j = 0; j<MM; j++)
				{
					ResCIsMat_temp[i][s][2*j] = ResCIsMat[i][RandPermut[s]][2*j];
					ResCIsMat_temp[i][s][2*j+1] = ResCIsMat[i][RandPermut[s]][2*j+1];
				}
			}
		}
		PartitioningRankingGeneralProcTuk(ResCIsMat_temp, ResCIsGridMat, Lower_temp, Upper_temp, clone(Lower_temp), CnkMat, crit, n, MM, gridSize, alpha, trace);

		//# Permute the ranks
		for(int s = 0; s<n; s++)
		{
			Lower[RandPermut[s]] = fmin(Lower[RandPermut[s]], Lower_temp[s]);
			Upper[RandPermut[s]] = fmax(Upper[RandPermut[s]], Upper_temp[s]);
		}
	}
	
	Rcout<<"\n";
	NumericMatrix CIs(n,2);
	for(int i = 0; i<n; i++)
	{
		CIs[i] = Lower[i]+1;
		CIs[n+i] = Upper[i]+1;
	}
	
	// Free some space.
	for (int i = 0; i<n; i++)
	{
		delete[] CnkMat[i];
	}
	delete[] CnkMat;
	for(int a = 0; a<gridSize; a++)
	{
		for(int i = 0; i<n; i++)
		{
			delete[] ResCIsMat[a][i];
			delete[] ResCIsMat_temp[a][i];
			delete[] ResCIsGridMat[a][i];
		}
		delete[] ResCIsMat[a];
		delete[] ResCIsMat_temp[a];
		delete[] ResCIsGridMat[a];
	}
	delete[] ResCIsMat;
	delete[] ResCIsMat_temp;
	delete[] ResCIsGridMat;
	return CIs;
}

