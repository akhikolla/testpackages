/*
File:             DKnn.cpp
Created by:       Oleksii Pokotylo
First published:
Last revised:

The realization of the Depth-based KNN classifier of Paindaveine and Van Bever (2015).
*/

#include "stdafx.h"

// compare UPoints descending
static int Compare(UPoint p1, UPoint p2){
	return (p1.value > p2.value);
}


void CountDepths(TDMatrix learnpoints, int* learnlabels, int nlearn, int d, TDMatrix checkpoints, int ncheck,
	int depthType, vector<UPoint>& depths, double* tempDepths, 
	TVariables car, TDMatrix dirs, TDMatrix prjs, TDMatrix ptPrjDepths, int ndir // halfspace
	){
	if (depthType == 1){
		for (int i = 0; i < ncheck; i++){
			GetDepths(checkpoints[i], learnpoints, nlearn, d, car,
				ndir, i != 0, dirs, prjs, &(depths[i].value), ptPrjDepths);
			depths[i].pattern = learnlabels[i];
		}
		return;
	}

	if (depthType == 2)
		MahalanobisDepth(learnpoints, checkpoints, d, nlearn, ncheck, 1, tempDepths);
	if (depthType == 3)
		SimplicialDepthsApx(learnpoints, checkpoints, d, nlearn, ncheck, choose(nlearn, d)*0.05, tempDepths);

	for (int i = 0; i < ncheck; i++){
		depths[i].value = tempDepths[i];
		depths[i].pattern = learnlabels[i];
	}
}



/*parameter cv: true - return classes for each k, false - return classes for kMax only*/
void knnGetClasses(TDMatrix learnpoints, int* learnlabels, int nlearn, int d,
	int nClasses, TDMatrix checkpoints, int ncheck, int kMax, bool cv,
	int depthType,
	int* classes /*for cv matrix [ncheck*kMax], for classification vector [ncheck]*/
	){
	// create the data structure for the reflected points
	double* arr = new double[nlearn*d];
	TDMatrix reflected = new double*[nlearn * 2];
	for (int i = 0; i < nlearn; i++){
		reflected[2 * i] = learnpoints[i];
		reflected[2 * i + 1] = arr + i*d;
	}

	vector<UPoint> depths(nlearn);
	double* tempDepths = new double[nlearn];
	int ndir = 1000;
	TVariables car(1, nlearn * 2);
	TDMatrix dirs; if (depthType == 1) dirs = newM(ndir, d);
	TDMatrix prjs; if (depthType == 1) prjs = newM(ndir, nlearn * 2);
	TDMatrix ptPrjDepths; if (depthType == 1) ptPrjDepths = newM(ndir, 1);

	for (int p = 0; p < ncheck; p++){
		double* point = checkpoints[p];
		// reflect
		for (int i = 0; i < nlearn; i++){
			for (int j = 0; j < d; j++){
				reflected[2 * i + 1][j] = 2 * point[j] - reflected[2 * i][j];
			}
		}

		// count depths of learn in reflected
		CountDepths(reflected, learnlabels, nlearn * 2, d, learnpoints, nlearn,
			depthType, depths, tempDepths,
			car, dirs, prjs, ptPrjDepths, ndir // halfspace
			);
	/*	for (int i = 0; i < nlearn; i++){
			GetDepths(learnpoints[i], reflected, nlearn * 2, d, car,
				ndir, i != 0, dirs, prjs, &(depths[i].value), ptPrjDepths);
			depths[i].pattern = learnlabels[i];
		}*/

		sort(depths.begin(), depths.end(), Compare);

		TVariables counts(nClasses+1, 0); 
		for (int i = 1; i <= nClasses; i++) counts[i] = 0;

		int prevmax = 0, prevclass = -1;
		for (int k = 1; k <= kMax; k++){
			counts[depths[k - 1].pattern]++;
			int max = -1, clmax = -1;
			for (int cl = 1; cl <= nClasses; cl++){
				if (counts[cl]>max){
					max = counts[cl];
					clmax = cl;
				} else if (max == counts[cl] && max == prevmax){
					// if the same number of neighbors, use the prev rule
					clmax = prevclass;
				}
			}
			if (cv) classes[p*kMax + k - 1] = clmax;
			prevclass = clmax;
			prevmax = max;
		}
		if (!cv)
			classes[p] = prevclass;
	}

	delete[] tempDepths;
	if (depthType == 1){
		deleteM(dirs);
		deleteM(prjs);
		deleteM(ptPrjDepths);
	}

	delete[] reflected;
	delete[] arr;
}

int DKnnCv(TDMatrix points, int n, int d, int* labels, int kMax, int depthType, int chunkNumber){
	set<int> unique_labels; unique_labels.insert(labels, labels + n - 1);
	int nClasses = unique_labels.size();

	int chunkSize = ceil((double)n / chunkNumber);

	TDMatrix learnpoints = new double*[n - chunkSize + 1];
	TDMatrix checkpoints = new double*[chunkSize];
	int* learnlabels = new int[n - chunkSize + 1];
	int* checklabels = new int[chunkSize];
	int* testlabels = new int[n];
	int* classes = new int[n*kMax];

	int chunk = 0;
	for (int j = 0, l = 0, c = 0; j < n; j++){
		if (j%chunkNumber){
			learnpoints[l] = points[j];
			learnlabels[l++] = labels[j];
		}
		else{
			checkpoints[c] = points[j];
			checklabels[c++] = labels[j];
		}
	}

	bool bigch = true;
	int hadObjects = 0;
	for (; chunk < chunkNumber; chunk++){
		if (chunk > 0){
			if (bigch && (chunkNumber)*(chunkSize - 1) + chunk == n){
				bigch = false;
				chunkSize--;
				learnpoints[n - chunkSize - 1] = points[n - 1];
				learnlabels[n - chunkSize - 1] = labels[n - 1];
			}

			for (int i = 0; i < chunkSize; i++){
				checkpoints[i] = learnpoints[(chunkNumber - 1)*i + chunk - 1];
				checklabels[i] = learnlabels[(chunkNumber - 1)*i + chunk - 1];
				learnpoints[(chunkNumber - 1)*i + chunk - 1] = points[chunkNumber*i + chunk - 1];
				learnlabels[(chunkNumber - 1)*i + chunk - 1] = labels[chunkNumber*i + chunk - 1];
			}
		}

		knnGetClasses(learnpoints, learnlabels, n - chunkSize, d, nClasses, 
			checkpoints, chunkSize, kMax, true, depthType, classes + hadObjects*kMax);
		//store checklabels
		memcpy(testlabels + hadObjects, checklabels, chunkSize*sizeof(int));
		hadObjects += chunkSize;
	}

	// run over k, count errors
	int kopt = 1, opterr = n;
	for (int k = 1; k <= kMax; k++){
		int err = 0;
		for (int p = 0; p < n; p++){
			if (classes[p*kMax + k - 1] != testlabels[p])
				err++;
		}
		if (err < opterr){
			opterr = err;
			kopt = k;
		}
	}

	delete[] learnpoints;
	delete[] checkpoints;
	delete[] learnlabels;
	delete[] checklabels;
	delete[] testlabels;
	delete[] classes;
	return kopt;
}

void DKnnClassify(TDMatrix points, int n, int d, int* labels, TDMatrix objects, int nobjects, int k, int depthType, int* classes){
	set<int> unique_lasbels; unique_lasbels.insert(labels, labels + n - 1);
	int nClasses = unique_lasbels.size();

	knnGetClasses(points, labels, n, d,
		nClasses, objects, nobjects, k, false,
		depthType,
		classes);
}


