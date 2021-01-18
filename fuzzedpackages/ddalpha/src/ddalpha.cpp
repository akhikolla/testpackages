/*
  File:             ddalpha.cpp
  Created by:       Pavlo Mozharovskyi, Oleksii Pokotylo
  First published:  28.02.2013
  Last revised:     20.02.2019

  Defines the exported functions for the 'ddalpha'-package.

  For a description of the algorithm, see:
    Lange, T., Mosler, K. and Mozharovskyi, P. (2012). Fast nonparametric classification based on data depth. Statistical Papers.
    Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world data with the DDalpha-procedure. Mimeo.
*/

#include "stdafx.h"

#define EOF (-1)

#ifdef __cplusplus
extern "C" {
#endif

boost::random::rand48 rEngine;
boost::random::normal_distribution<double> normDist;

void Sum(double *a, double *b, double *res){
	res[0] = a[0] + b[0];
}

void setSeed(int random_seed){
	if (random_seed != 0) {
		setseed(random_seed);
		rEngine.seed(random_seed);
	}
	else {
		setseed(time(NULL));
		rEngine.seed(time(NULL));
	}
}

void IsInConvexes(double *points, int *dimension, int *cardinalities, int *numClasses, double *objects, int *numObjects, int *seed, int *isInConvexes){
	setSeed(*seed);
	int numPoints = 0;for (int i = 0; i < numClasses[0]; i++){numPoints += cardinalities[i];}
	TMatrix x(numPoints);
	for (int i = 0; i < numPoints; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TMatrix o(numObjects[0]);
	for (int i = 0; i < numObjects[0]; i++){o[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numObjects[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			o[i][j] = objects[i * dimension[0] + j];
		}
	}
	TVariables cars(numClasses[0]);
	for (int i = 0; i < numClasses[0]; i++){
		cars[i] = cardinalities[i];
	}
	TIntMatrix answers(o.size());
	int error = 0;
	InConvexes(x, cars, o, error, &answers);
	for (int i = 0; i < numObjects[0]; i++)
    for (int j = 0; j < numClasses[0]; j++){
  		isInConvexes[numClasses[0]*i+j] = answers[i][j];
  	}
}

void ZDepth(double *points, double *objects, int *numPoints, int *numObjects, int *dimension, int *seed, double *depths){
	setSeed(*seed);
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TPoint means(*dimension); TPoint sds(*dimension);
	GetMeansSds(x, &means, &sds);
	Standardize(x, means, sds);
	TMatrix z(numObjects[0]);
	for (int i = 0; i < numObjects[0]; i++){z[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numObjects[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			z[i][j] = objects[i * dimension[0] + j];
		}
		Standardize(z[i], means, sds);
		int error;
		depths[i] = ZonoidDepth(x, z[i], error);
	}
}

void HDepth(double *points, double *objects, int *numObjects, int *dimension, 
	int *cardinalities, int *numClasses, double *directions, double *projections, 
	int *k, int *sameDirs, int *seed, double *depths){
	setSeed(*seed);
	int numPoints = 0;for (int i = 0; i < numClasses[0]; i++){numPoints += cardinalities[i];}
	TDMatrix x = asMatrix(points, numPoints, dimension[0]);
	TDMatrix z = asMatrix(objects, numObjects[0], dimension[0]);

	TVariables cars(numClasses[0]);
	for (int i = 0; i < numClasses[0]; i++){
		cars[i] = cardinalities[i];
	}

	TDMatrix dirs = asMatrix(directions, k[0], *dimension);
	TDMatrix prjs = asMatrix(projections,k[0], numPoints);
	TDMatrix ptPrjDepths = newM(*k, *numClasses);
	for (int i = 0; i < numObjects[0]; i++){
		GetDepths(z[i], x, numPoints, *dimension, cars, 
			k[0], i == 0 ? 0 : sameDirs[0] /*at the first step fill the matrices*/, 
			dirs, prjs, depths + i * numClasses[0], ptPrjDepths);
	/*	for (int j = 0; j < numClasses[0]; j++){
			depths[i * numClasses[0] + j] = dps[j];
		}*/
	}
	deleteM(ptPrjDepths);

/*	if (*sameDirs){
		for (int i = 0; i < k[0] * dimension[0]; i++){
			directions[i] = dirs[i / dimension[0]][i%dimension[0]];
		}
		for (int i = 0; i < k[0] * numPoints; i++){
			projections[i] = prjs[i / numPoints][i%numPoints];
		}
		}
	deleteM(dirs);
	deleteM(prjs);
*/

	delete[] x;
	delete[] z;
	delete[] dirs;
	delete[] prjs;
}

void HDSpace(double *points, int *dimension, int *cardinalities, int *numClasses,
	int *k, int *sameDirs, int *seed, double *dSpace, double *directions, double *projections){
	setSeed(*seed);
	int numPoints = 0;for (int i = 0; i < numClasses[0]; i++){numPoints += cardinalities[i];}
	TDMatrix x = asMatrix(points, numPoints, *dimension);

	TVariables cars(numClasses[0]);
	for (int i = 0; i < numClasses[0]; i++){
		cars[i] = cardinalities[i];
	}
	TDMatrix dsps = asMatrix(dSpace, numPoints, *numClasses);
	TDMatrix dirs = asMatrix(directions, k[0], (*dimension));
	TDMatrix prjs = asMatrix(projections, k[0], (numPoints));
	GetDSpace(x, numPoints, *dimension, cars, k[0], sameDirs[0], dsps, dirs, prjs);
/*	for (int i = 0; i < numPoints*numClasses[0]; i++){
		dSpace[i] = dsps[i/numClasses[0]][i%numClasses[0]];
	}
*/
	/*	if (*sameDirs){
	for (int i = 0; i < k[0] * dimension[0]; i++){
	directions[i] = dirs[i / dimension[0]][i%dimension[0]];
	}
	for (int i = 0; i < k[0] * numPoints; i++){
	projections[i] = prjs[i / numPoints][i%numPoints];
	}
	}
	deleteM(dirs);
	deleteM(prjs);
	*/

	delete[] x;
	delete[] dsps;
	delete[] dirs;
	delete[] prjs;
}

void HDepthSpaceEx(double *points, double *objects, int *cardinalities, int *numClasses, int *numObjects,
	int *dimension, int *algNo, double *depths){
	
	double(*func)(double *z, double **xx, int n, int d);
	switch ((HDalgs)*algNo)
	{
	case recursive:
		func = &HD_Rec; break;
	case plane:
		func = &HD_Comb2; break;
	case line:
		func = &HD_Comb; break;
	default:
		func = 0; break;
	}

	TDMatrix x = asMatrix(objects, *numObjects, *dimension);
	int classBegin = 0;

	if (func)
	for (int c = 0; c < *numClasses; c++){
		TDMatrix X = asMatrix(points+classBegin, cardinalities[c], *dimension);
	//	printMatrix(X, cardinalities[c], *dimension);
		for (int i = 0; i < *numObjects; i++){
			depths[c * (*numObjects) + i] = func(x[i], X, cardinalities[c], *dimension);
		}
		classBegin += cardinalities[c]* *dimension;
		delete[] X;
	}
	delete[] x;
}

void HDepthEx(double *points, double *objects, int *numPoints, int *numObjects, int *dimension, int *algNo, double *depths){

	double(*func)(double *z, double **xx, int n, int d);
	switch ((HDalgs)*algNo)
	{
	case recursive:
		func = &HD_Rec; break;
	case plane:
		func = &HD_Comb2; break;
	case line:
		func = &HD_Comb; break;
	default:
		func = 0; break;
	}

	TDMatrix X = asMatrix(points, *numPoints, *dimension);
	TDMatrix x = asMatrix(objects, *numObjects, *dimension);

	if (func)
	for (int i = 0; i < *numObjects; i++){
		depths[i] = func(x[i], X, *numPoints, *dimension);
	}
	delete[] X;
	delete[] x;
}

void MahalanobisDepth(double *points, double *objects, int *numPoints, int *numObjects, int *dimension, double* MCD, double *depths){
	TDMatrix X = asMatrix(points, *numPoints, *dimension);
	TDMatrix x = asMatrix(objects, *numObjects, *dimension);

	MahalanobisDepth(X, x, *dimension, *numPoints, *numObjects, *MCD, depths);

	delete[] X;
	delete[] x;
}

void OjaDepth(double *points, double *objects, int *numPoints, int *numObjects, int *dimension, int *seed, int* exact, int *k, int *useCov, double *covEst, double *depths){
	setSeed(*seed);
	TDMatrix X = asMatrix(points, *numPoints, *dimension);
	TDMatrix x = asMatrix(objects, *numObjects, *dimension);
	TDMatrix cov = asMatrix(covEst, *dimension, *dimension);

	if (*exact)
		OjaDepthsEx(X, x, *dimension, *numPoints, *numObjects, *useCov, cov, depths);
	else{
		long long K = ((long long)2000000000)*k[0] + k[1];
		OjaDepthsApx(X, x, *dimension, *numPoints, *numObjects, K, *useCov, cov, depths);
	}
	delete[] X;
	delete[] x;
	delete[] cov;
}

void SimplicialDepth(double *points, double *objects, int *numPoints, int *numObjects, int *dimension, int *seed, int* exact, int *k, double *depths){
	setSeed(*seed);
	TDMatrix X = asMatrix(points, *numPoints, *dimension);
	TDMatrix x = asMatrix(objects, *numObjects, *dimension);

	if (*dimension == 2)
		SimplicialDepths2(X, x, *numPoints, *numObjects, depths);
	else if (*exact)
		SimplicialDepthsEx(X, x, *dimension, *numPoints, *numObjects, depths);
	else {
		long long K = ((long long)2000000000)*k[0] + k[1];
		SimplicialDepthsApx(X, x, *dimension, *numPoints, *numObjects, K, depths);
	}
	delete[] X;
	delete[] x;
}

void AlphaLearn(double *points, int *numPoints, int *dimension, int *cardinalities, int *degree, int *minFeatures, double *ray){
	

	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables y(numPoints[0]);
	for (int i = 0; i < cardinalities[0]; i++){y[i] = 1;}
	for (int i = cardinalities[0]; i < numPoints[0]; i++){y[i] = -1;}
	TMatrix _x;
	ExtendWithProducts(x, degree[0], &_x);
	TPoint direction;
	OUT_ALPHA = true;
	Learn(_x, y, minFeatures[0], &direction);
	ray[0] = degree[0];
	for (int i = 0; i < direction.size(); i++){
		ray[i + 1] = direction[i];
	}
}

void AlphaLearnCV(double *points, int *numPoints, int *dimension, int *cardinalities, int *upToPower, int *numFolds, int *minFeatures, int *debug, double *ray){
	TMatrix x(numPoints[0],TPoint(dimension[0]));
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables y(numPoints[0]);
	for (int i = 0; i < cardinalities[0]; i++){y[i] = 1;}
	for (int i = cardinalities[0]; i < numPoints[0]; i++){y[i] = -1;}
	TPoint direction; unsigned int power;
	OUT_ALPHA = (debug[0])!=0;
	LearnCV(x, y, minFeatures[0], upToPower[0], numFolds[0], &direction, &power);
	ray[0] = power;
	for (int i = 0; i < direction.size(); i++){
		ray[i + 1] = direction[i];
	}
}

void AlphaClassify(double *points, int *numPoints, int *dimension, int *degree, double *ray, int *output){
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TMatrix _x;
	ExtendWithProducts(x, degree[0], &_x);
	TPoint direction(_x[0].size());
	for (int i = 0; i < _x[0].size(); i++){
		direction[i] = ray[i + 1];
	}
	TVariables y;
	Classify(_x, direction, &y);
	for (int i = 0; i < numPoints[0]; i++){
		output[i] = y[i];
	}
}

void KnnAffInvLearnJK(double *points, int *dimension, int *cardinalities, int *maxk, int *k){
	int numPoints = cardinalities[0] + cardinalities[1];
	TMatrix x(numPoints);
	for (int i = 0; i < numPoints; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables cars(2);cars[0] = cardinalities[0];cars[1] = cardinalities[1];
	k[0] = GetK_JK_Binary(x, cars, maxk[0]);
}

void KnnAffInvClassify(double *objects, int *numObjects, double *points, int *dimension, int *cardinalities, int *k, int *output){
	int numPoints = cardinalities[0] + cardinalities[1];
	TMatrix x(numPoints);
	for (int i = 0; i < numPoints; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables cars(2);cars[0] = cardinalities[0];cars[1] = cardinalities[1];
	TMatrix y(numObjects[0]);
	for (int i = 0; i < numObjects[0]; i++){y[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numObjects[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			y[i][j] = objects[i * dimension[0] + j];
		}
	}
	TVariables result;
	Knn_Classify_Binary(y, x, cars, k[0], &result);
	for (int i = 0; i < numObjects[0]; i++){
		output[i] = result[i];
	}
}

void KnnLearnJK(double *points, int *labels, int *numPoints, int *dimension, int *kmax, int *distType, int *k){
	TMatrix x(numPoints[0]);TVariables y(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){
		x[i] = TPoint(dimension[0]);
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
		y[i] = labels[i];
	}
	k[0] = KnnCv(x, y, kmax[0], distType[0], 0);
}

void KnnClassify(double *objects, int *numObjects, double *points, int *labels, int *numPoints, int *dimension, int *k, int *distType, int *output){
	TMatrix x(numPoints[0]);TVariables y(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){
		x[i] = TPoint(dimension[0]);
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
		y[i] = labels[i];
	}
	TMatrix z(numObjects[0]);
	for (int i = 0; i < numObjects[0]; i++){
		z[i] = TPoint(dimension[0]);
		for (int j = 0; j < dimension[0]; j++){
			z[i][j] = objects[i * dimension[0] + j];
		}
	}
	TVariables result;
	Knn(z, x, y, k[0], distType[0], &result);
	for (int i = 0; i < numObjects[0]; i++){
		output[i] = result[i];
	}
}

void DKnnLearnCv(double *points, int *labels, int *numPoints, int *dimension, int *kmax, int *depthType, int *k, int* chunkNumber, int *seed){
	setSeed(*seed);
	TDMatrix x = asMatrix(points, *numPoints, *dimension);
	*k = DKnnCv(x, *numPoints, *dimension, labels, *kmax, *depthType, *chunkNumber);
	delete[] x;
}

void DKnnClassify(double *objects, int *numObjects, double *points, int *labels, int *numPoints, int *dimension, int *k, int *depthType, int *seed, int *output){
	setSeed(*seed);
	TDMatrix x = asMatrix(points, *numPoints, *dimension);
	TDMatrix z = asMatrix(objects, *numObjects, *dimension);
	
	DKnnClassify(x, *numPoints, *dimension, labels, z, *numObjects, *k, *depthType, output);
	delete[] x;
	delete[] z;
}

void PolynomialLearnCV(double *points, int *numPoints, int *dimension, int *cardinalities, int *maxDegree, int *chunkNumber, int *seed, /*OUT*/ int *degree, /*OUT*/ int *axis, /*OUT*/ double *polynomial){
	setSeed(*seed);
	TDMatrix x = asMatrix(points, numPoints[0], dimension[0]);
	TVariables y(numPoints[0]);
	for (int i = 0; i < cardinalities[0]; i++){ y[i] = 1; }
	for (int i = cardinalities[0]; i < numPoints[0]; i++){ y[i] = -1; }
	
	TPoint pol = PolynomialLearnCV(x, cardinalities[0], cardinalities[1], *maxDegree, *chunkNumber, degree, axis);

	for (unsigned i = 0; i < pol.size(); i++){
		polynomial[i] = pol[i];
	}
	delete[] x;
}
/* everything implemented in R
void PolynomialClassify(double *points, int *numPoints, int *dimension, int *degree, double *ray, int *output){
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){ x[i] = TPoint(dimension[0]); }
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TMatrix _x;
	ExtendWithProducts(x, degree[0], &_x);
	TPoint direction(_x[0].size());
	for (unsigned i = 0; i < _x[0].size(); i++){
		direction[i] = ray[i + 1];
	}
	TVariables y;
	Classify(_x, direction, &y);
	for (int i = 0; i < numPoints[0]; i++){
		output[i] = y[i];
	}
}
*/

void ProjectionDepth(double *points, double *objects, int *numObjects,
					 int *dimension, int *cardinalities, int *numClasses,
					 double *directions, double *projections, int *k,
					 int *newDirs, int *seed, double *depths){
	setSeed(*seed);
	TVariables cars(numClasses[0]);
	int numPoints = 0;
	for (int i = 0; i < numClasses[0]; i++){
		numPoints += cardinalities[i];
		cars[i] = cardinalities[i];
	}
	TDMatrix x = asMatrix(points, numPoints, *dimension);
	TDMatrix z = asMatrix(objects, *numObjects, *dimension);

	TDMatrix dirs = asMatrix(directions, k[0], *dimension);
	TDMatrix prjs = asMatrix(projections, k[0], numPoints);
	TDMatrix _depths = asMatrix(depths, *numObjects, *numClasses);
	GetDepthsPrj(x, numPoints, *dimension, z, *numObjects, cars, 
		*k, *newDirs, _depths, dirs, prjs);
/*	for (int i = 0; i < numObjects[0]; i++){
		for (int j = 0; j < numClasses[0]; j++){
			depths[i * numClasses[0] + j] = _depths[i][j];
		}
	}
	if (newDirs[0]){
		for (int i = 0; i < k[0]*dimension[0]; i++){
			directions[i] = dirs[i/dimension[0]][i%dimension[0]];
		}
		for (int i = 0; i < k[0]*numPoints; i++){
			projections[i] = prjs[i/numPoints][i%numPoints];
		}
	}
	*/

	delete[] x;
	delete[] z;
	delete[] dirs;
	delete[] prjs;
	delete[] _depths;
}

void PotentialDepthsCount(double *points, int *numPoints, int *dimension, int *classes, int *cardinalities, double *testpoints, int *numTestPoints, int* kernelType, double *a, int* ignoreself, double *depths){
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){
		TPoint& curPoint = x[i];
		curPoint.resize(dimension[0]);
		for (int j = 0; j < dimension[0]; j++){
			curPoint[j] = points[i * dimension[0] + j];
		}
	}
	
	TMatrix xt(numTestPoints[0]);
	for (int i = 0; i < numTestPoints[0]; i++){
		TPoint& curPoint = xt[i];
		curPoint.resize(dimension[0]);
		for (int j = 0; j < dimension[0]; j++){
			curPoint[j] = testpoints[i * dimension[0] + j];
		}
	}

	TMatrix d(numTestPoints[0]);
	for (int i = 0; i < numTestPoints[0]; i++){
		d[i].resize(classes[0]);
	}
	TVariables car(classes[0]);
	for (int i = 0; i < classes[0]; i++){
		car[i] = cardinalities[i];
	}

	double (*Kernel) (TPoint& x, TPoint& y, double a) = 0;

	switch (*kernelType){
		case 1: Kernel = EDKernel; break;
		case 2: Kernel = GKernel; break;
		case 3: Kernel = EKernel; break;
		case 4: Kernel = TriangleKernel; break;
		case 5: Kernel = VarGKernel; break;
		default: throw "Unsupported kernel type";
	}
	
	PotentialDepths(x, car, xt, d, Kernel, *a, *ignoreself);

	for (int i = 0; i < numTestPoints[0]; i++){
		for (int j = 0; j < classes[0]; j++){
		//	depths[i * classes[0] + j] = d[i][j];
			depths[j * numTestPoints[0] + i] = d[i][j];
		}
	}
}

void BetaSkeletonDepth(double *points, double *objects, int *numPoints, int *numObjects, int *dimension, double* beta, int* distCode, double* p, double* sigma, double *depths){
  TDMatrix X = asMatrix(points, *numPoints, *dimension);
  TDMatrix x = asMatrix(objects, *numObjects, *dimension);
  TDMatrix s = asMatrix(sigma, *dimension, *dimension);

  LensDepth(X, x, *dimension, *numPoints, *numObjects, *beta, *distCode, *p, s, depths);
  
  delete[] X;
  delete[] x;
  delete[] s;
}

void SimplicialBandDepthF(double *objectsf, double *dataf, double *args, 
                          int *numObjects, int *numPoints, int *numArgs, 
                          int *dimension, int *modified, int *J, 
                          double *depths){
  // Structure the input data
  T3DMatrix x = as3DMatrix(objectsf, *numObjects, *numArgs, *dimension);
  T3DMatrix X = as3DMatrix(dataf, *numPoints, *numArgs, *dimension);
  // Delegate calculation of depths
  BandDepth(x, X, *numObjects, *numPoints, *numArgs, *dimension, 
            (bool)*modified, *J, depths);
  // Clean the memory
  for (int i = 0; i < *numPoints; i++){
    delete[] X[i];
  }
  delete[] X;
  for (int i = 0; i < *numObjects; i++){
    delete[] x[i];
  }
  delete[] x;
}

#ifdef __cplusplus
}
#endif
