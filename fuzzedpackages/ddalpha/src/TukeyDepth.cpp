/*
  File:             TukeyDepth.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  28.02.2013
  Last revised:     13.11.2015
  
  Computation of the random Tukey data depth.

  For a description of the algorithm, see:
    Cuesta-Albertos, J. A. and Nieto-Reyes, A. (2008). The random Tukey depth. Computational Statistics & Data Analysis 52, 11 (July 2008), 4979-4988.
    Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world data with the DDalpha-procedure. Mimeo.
*/

#include "stdafx.h"

static int CompareAsc(OrderRec x, OrderRec y)
{
	return (x.value < y.value);
}

static int CompareDec(OrderRec x, OrderRec y)
{
	return (x.value > y.value);
}

void GetPrjDepths(double* projection, int n, TVariables& cardinalities, unsigned curClass, TVariables *prjDepths){
	//Collecting basic statistics
	int beginIndex = 0;
	for (unsigned i = 0; i < cardinalities.size(); i++){
		if (i >= curClass){break;}
		beginIndex += cardinalities[i];
	}
	int endIndex = beginIndex + cardinalities[curClass] - 1;

	//Preparing structures
	vector<OrderRec> prjSort(n);
	for (int i = 0; i < n; i++){
		prjSort[i].order = i;
		prjSort[i].value = projection[i];
	}
	
	//Calculating projection depths
	TVariables depthsForwards(n);
	TVariables depthsBackwards(n);
	//Forwards
	sort(prjSort.begin(), prjSort.end(), CompareAsc);
	int curDepth = 0;
	for (int i = 0; i < n; i++){
		if ((prjSort[i].order >= beginIndex) && (prjSort[i].order <= endIndex)){curDepth++;}
		depthsForwards[prjSort[i].order] = curDepth;
	}
	//Backwards
	sort(prjSort.begin(), prjSort.end(), CompareDec);
	curDepth = 0;
	for (int i = 0; i < n; i++){
		if ((prjSort[i].order >= beginIndex) && (prjSort[i].order <= endIndex)){curDepth++;}
		depthsBackwards[prjSort[i].order] = curDepth;
	}
	//Merge
	for (int i = 0; i < n; i++){
		if (depthsForwards[i] < depthsBackwards[i]){
			(*prjDepths)[i] = depthsForwards[i];
		}else{
			(*prjDepths)[i] = depthsBackwards[i];
		}
	}
}

inline void GetPtPrjDepths(double* projection, int n, double point, TVariables& cardinalities, double* ptPrjDepths){
	int q = cardinalities.size();
	for (int i = 0; i < q; i++){
		int beginIndex = 0;
		for (int j = 0; j < q; j++){
			if (j >= i){break;}
			beginIndex += cardinalities[j];
		}
		int endIndex = beginIndex + cardinalities[i];
		int nPtsBelow = 0;
		int nPtsAbove = 0;
		for (int j = beginIndex; j < endIndex; j++){
			if (projection[j] <= point){nPtsBelow++;}
			if (projection[j] >= point){nPtsAbove++;}
		}
		ptPrjDepths[i] = (nPtsBelow <= nPtsAbove)?(double)nPtsBelow:(double)nPtsAbove;
	}
}

//Indexing from zero
void GetDSpace(TDMatrix points, int n, int d, TVariables& cardinalities, int k, bool atOnce, TDMatrix dSpace, TDMatrix directions, TDMatrix projections){
	//1. Collecting basic statistics
	int q = cardinalities.size();

	if (!atOnce){
		TDMatrix ptPrjDepths = newM(k, q);
		for (int i = 0; i < n; i++){
/*			TMatrix dir(k, TPoint(d));
			TMatrix proj(k, TPoint(q));*/
			GetDepths(points[i], points, n, d, cardinalities, k, false, directions, projections, dSpace[i], ptPrjDepths);
		}
		deleteM(ptPrjDepths);
		return;
	}
	GetDirections(directions, k, d);
	GetProjections(points, n, d, directions, k, projections);
	//2. Calculate projection depths
	vector<vector<TVariables> > prjDepths(k, vector<TVariables>(q, TVariables(n)));
	for (int i = 0; i < k; i++){
		for (int j = 0; j < q; j++){
			GetPrjDepths(projections[i], n, cardinalities, j, &prjDepths[i][j]);
		}
	}
	//3. Merge depths
	for (int i = 0; i < n; i++){
		for (int j = 0; j < q; j++){
			dSpace[i][j] = cardinalities[j] + 1;
		}
	}
	for (int i = 0; i < k; i++){
		for (int j = 0; j < q; j++){
			for (int l = 0; l < n; l++){
				if (prjDepths[i][j][l] < dSpace[l][j]){
					dSpace[l][j] = prjDepths[i][j][l];
				}
			}
		}
	}
	for (int i = 0; i < q; i++){
		for (int j = 0; j < n; j++){
			dSpace[j][i] /= cardinalities[i];
		}
	}
}

void GetDepths(double* point, TDMatrix points, int n, int d, 
	TVariables& cardinalities, int k, bool atOnce, 
	TDMatrix directions, TDMatrix projections, double* depths,
	TDMatrix ptPrjDepths /*accu, k*q */){
	//1. Collecting basic statistics
	int q = cardinalities.size();
	if (!atOnce){
		GetDirections(directions, k, d);
		GetProjections(points, n, d, directions, k, projections);
	}	
	//2. Calculate projection depths
	TPoint pointProjections(k);
	for (int i = 0; i < k; i++){
		double curPrj = 0;
		for (int j = 0; j < d; j++){
			curPrj += point[j]*directions[i][j];
		}
		pointProjections[i] = curPrj;
	}
	for (int i = 0; i < k; i++){
		GetPtPrjDepths(projections[i], n, pointProjections[i], cardinalities, ptPrjDepths[i]);
	}
	//3. Merge depths
	for (int i = 0; i < q; i++){
		depths[i] = cardinalities[i] + 1;
	}
	for (int i = 0; i < k; i++){
		for (int j = 0; j < q; j++){
			if (ptPrjDepths[i][j] < depths[j]){
				depths[j] = ptPrjDepths[i][j];
			}
		}
	}
	for (int i = 0; i < q; i++){
		depths[i] /= (double)cardinalities[i];
	}
}
