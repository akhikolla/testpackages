/*
  File:             AlphaProcedure.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  28.02.2013
  Last revised:     13.11.2015

  Contains the modified alpha-procedure for the DDalpha-classifier.

  For a description of the algorithm, see:
    Lange, T., Mosler, K. and Mozharovskyi, P. (2012). Fast nonparametric classification based on data depth. Statistical Papers.
    Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world data with the DDalpha-procedure. Mimeo.
*/

#include "stdafx.h"

#define PI2 1.5707963267948966192313216916398

/* Definition of static variables*/
static TMatrix x;
static TVariables y;
static unsigned int numLess;
static unsigned int numMore;
static int difference;
static unsigned int n;
static unsigned int d;
static TVariables properties;
static Features features;
static TPoint curFeature;
static int numStartFeatures;

static int GetProducts(TPoint values, unsigned int power, TPoint *output){
	int d = values.size();
	if (power == 1){
		output->resize(d);
		for (int i = 0; i < d; i++){(*output)[i] = values[i];}
		return 0;
	}
	output->resize(0);
	TVariables indices(power);
	unsigned int counter = 0;
	while(indices[0] < d){
		output->push_back(1);
		for (unsigned int i = 0; i < power; i++){
			(*output)[counter] *= values[indices[i]];
		}
		counter++;
		int lastArray = power - 1;
		while(lastArray > 0 && indices[lastArray] == d - 1){lastArray--;}
		indices[lastArray]++;
		for (unsigned int i = lastArray; i < power; i++){
			indices[i] = indices[lastArray];
		}
		
	}
	return 0;
}

static int Compare(UPoint p1, UPoint p2){
	return (p1.value < p2.value);
}

static int UpdateCurFeature(){	
	double angle = -features[features.size() - 1].angle;
	unsigned int yAxisNumber = features[features.size() - 1].number;
	for (unsigned int i = 0; i < n; i++){
		curFeature[i] = curFeature[i]*cos(angle) - x[yAxisNumber][i]*sin(angle);
	}
	return 0;
}

static unsigned int DGetMinError(unsigned int yAxisNumber, Feature *yFeature){
	/* Calculate angle of each point to the xAxis and sort them */
	UPoints angles(n);
	for (unsigned int i = 0; i < n; i++){
		angles[i] = UPoint((x[yAxisNumber][i] == 0 && curFeature[i] == 0) ? 0 : y[i] , atan2(x[yAxisNumber][i], curFeature[i]));
	}
	sort(angles.begin(), angles.end(), Compare);
/*
	for (unsigned i = 0; i < n; i++){
		cout << (angles[i].pattern > 0 ? 1 : angles[i].pattern < 0 ? 0 : 3);
	}
	cout << endl;
*/
	/* Look for the optimal threshold */
	int leftDiff = 0; unsigned int optThreshold = 0; int maxCorr = 0; double nextValue = angles[0].value;
	for (unsigned i = 0; i < n - 1; i++){
		leftDiff += angles[i].pattern;
		if (angles[i + 1].value == nextValue){ continue; } nextValue = angles[i].value;
		int corr = max(numMore - leftDiff, numLess + leftDiff);
		//int corr = abs(leftDiff) + abs(difference - leftDiff); 
		if (corr > maxCorr){ maxCorr = corr; optThreshold = i; }
//		cout << i << " " << corr << "; ";
	}
//	cout << endl;

	/* Determine the angle of the separating direction */
	yFeature->angle = (angles[optThreshold].value + angles[optThreshold + 1].value) / 2. - PI2; yFeature->error = n - maxCorr; yFeature->number = yAxisNumber;
	return yFeature->error;
}

static unsigned int GetRay(TPoint *ray){
	ray->resize(d);
	double drivingAxis = 1;
	for (unsigned int i = features.size() - 1; i > 0; i--){
		(*ray)[features[i].number] = drivingAxis*sin(features[i].angle);
		drivingAxis = drivingAxis*cos(features[i].angle);
	}
	(*ray)[features[0].number] = drivingAxis;
	
	UPoints points(n);
	for (unsigned int i = 0; i < n; i++){
		points[i].pattern = y[i];
		for (unsigned int j = 0; j < d; j++){
			points[i].value += (*ray)[j]*x[j][i];
		}
#ifdef DEF_OUT_ALPHA
		if (OUT_ALPHA) Rcout << points[i].value << ", ";
#endif
	}
#ifdef DEF_OUT_ALPHA
	if (OUT_ALPHA) Rcout << endl;
#endif
	sort(points.begin(), points.end(), Compare);
	unsigned int numLeftLess = 0;
	unsigned int numLeftMore = 0;
	for (unsigned int i = 0; i < n; i++){
		if (points[i].value > 0){break;}
		if (points[i].pattern > 0){numLeftMore++;}else{numLeftLess++;}
	}
	unsigned int errorLeftLess = numLeftMore + numLess - numLeftLess;
	unsigned int errorLeftMore = numLeftLess + numMore - numLeftMore;
	if (errorLeftLess > errorLeftMore){
		for (unsigned int i = 0; i < d; i++){
			(*ray)[i] *= -1.;
		}
	}
#ifdef DEF_OUT_ALPHA
	if (OUT_ALPHA){
	Rcout << errorLeftLess << " " << errorLeftMore << " ";
	}
#endif
	return 0;
}

int Initialization(TMatrix input, TVariables output, unsigned int minFeatures){
	n = input.size(); if (n == 0){return -1;} // number of points
	if (output.size() != n){return -1;}
	d = input[0].size(); if (d == 0){return -1;} // space dimension		
	if (minFeatures == 0 || minFeatures > 2){return -1;}else{numStartFeatures = minFeatures;}
	
	/* Filling static variables x and y with input and output (transposing x) */
	x.resize(d);
	for (unsigned int i = 0; i < d; i++){
		x[i] = TPoint(n);
		for (unsigned int j = 0; j < n; j++){
			x[i][j] = input[j][i];
		}
	}
	y.resize(n); numLess = 0; numMore = 0;
	difference = 0;
	for (unsigned int i = 0; i < n; i++){
		y[i] = output[i];
		difference += y[i];
		if (y[i] > 0){numMore++;}else{numLess++;}
	}
	return 0;
}

int Alpha(TPoint *ray){
	/* 0. Subinitialization - clear auxiliary variables and empty and nonsignificant input axes */
	properties.resize(d); for (unsigned int i = 0; i < d; i++){properties[i] = i;} // initialize properties: all available
	features.clear();

	outMatrix(x);

	/* 1. Null-cycle */
	if (numStartFeatures == 2){ // start with two features?
		Feature optFeatureX;
		Feature optFeatureY;
		for (unsigned int i = 0; i < properties.size() - 1; i++){
			for (unsigned int j = i + 1; j < properties.size(); j++){
				/* Calculating minimal error on the plane of the i-th and the j-th properties */
				Feature tmpFeature;
				curFeature = x[properties[i]];
				unsigned int error = DGetMinError(properties[j], &tmpFeature);
#ifdef DEF_OUT_ALPHA
				if (OUT_ALPHA){
					Rcout << properties[i] << ", " << properties[j] << ", " << tmpFeature.angle << ", " << error << ", " << endl;
				}
#endif
				if (error < optFeatureY.error){optFeatureX.number = properties[i]; optFeatureY = tmpFeature;}
			}
		}
		features.push_back(optFeatureX);
		features.push_back(optFeatureY);
		for (unsigned int i = 0; i < properties.size(); i++){ // delete recently found X and Y properties
			if (properties[i] == optFeatureX.number){properties.erase(properties.begin() + i);}
			if (properties[i] == optFeatureY.number){properties.erase(properties.begin() + i);}
		}
		curFeature = x[features[0].number];
		UpdateCurFeature();
		outString("Feature 1:");
		outVector(curFeature);
	}

	/* 2. Main cycle */
	/* Searching an optimal feature space while empirical error rate decreases */	
	while(features[features.size() - 1].error > 0 && properties.size() > 0){
		Feature optFeature;
		for (unsigned int i = 0; i < properties.size(); i++){
			/* Calculating minimal error on the plane of the curFeature and the j-th properties */
			Feature tmpFeature;
			unsigned int error = DGetMinError(properties[i], &tmpFeature);
#ifdef DEF_OUT_ALPHA
			if (OUT_ALPHA){
				Rcout << properties[i] << ", " << tmpFeature.angle << ", " << error << ", " << endl;
			}
#endif
			if (error < optFeature.error){optFeature = tmpFeature;}
		}		
		if (optFeature.error < features[features.size() - 1].error){
			features.push_back(optFeature);
			for (unsigned int i = 0; i < properties.size(); i++){ // delete recently found property
				if (properties[i] == optFeature.number){properties.erase(properties.begin() + i);}
			}
			UpdateCurFeature();
			outString("Feature :");
			outVector(curFeature);
		}else{break;}
	}
	
	outString("Features:");
	outFeatures(features);

	/* Restoring the projection vector */
	GetRay(ray);
	return features[features.size() - 1].error;
}

int Classify(TMatrix input, TPoint weights, TVariables *output){
	/* 0. Initialization */
	unsigned int l = input.size(); if (l == 0){return -1;}
	unsigned int p = weights.size(); if (p == 0){return -1;} if (p > input[0].size()){return -1;}
	output->resize(l);
	/* 1. Classification of each point by comparison of their projections to 0 */
	for (unsigned int i = 0; i < l; i++){
		double curSum = 0;
		for (unsigned int j = 0; j < p; j++){curSum += weights[j]*input[i][j];}
		(*output)[i] = (curSum > 0) ? 1 : -1;
	}
	return 0;
}

int ExtendWithProducts(TMatrix input, unsigned int upToPower, TMatrix *output){
	unsigned int n = input.size();
	output->resize(n);
	for (unsigned int i = 0; i < n; i++){
		for(unsigned int j = 1; j <= upToPower; j++){
			TPoint extension;
			GetProducts(input[i], j, &extension);
			TPoint::iterator it;
			it = (*output)[i].end();
			(*output)[i].insert(it, extension.begin(), extension.end());
		}
	}
	return 0;
}

int Learn(TMatrix input, TVariables output, unsigned int minFeatures, TPoint *ray){
	if (Initialization(input, output, minFeatures) != 0){return -1;}
	return Alpha(ray);
}

int LearnCV(TMatrix input,  TVariables output, unsigned int minFeatures, unsigned int upToPower, unsigned int folds, TPoint *ray, unsigned int *power){
	bool oldOUT_ALPHA = OUT_ALPHA;
	OUT_ALPHA = false;
	unsigned int optDegree = 0;
	unsigned int optError = INT_MAX;
	unsigned int shortFolds = folds - 1;

	/* Get the optimal degree (outer cross-validation loop) */
	vector<TMatrix> spaceExtensions(upToPower);
	for (unsigned int i = 0; i < upToPower; i++){
		ExtendWithProducts(input, i + 1, &spaceExtensions[i]); // get the (i + 1)-th space extention
		Initialization(spaceExtensions[i], output, minFeatures); // initialize
		/* Prepare slider and start to cut data */
		unsigned sliderSize = (unsigned)ceil((double)n / folds); unsigned chSizeVal = n%folds - 1;
		TMatrix xSlider(sliderSize); TVariables ySlider(sliderSize);
		for (unsigned int j = 0; j < sliderSize; j++){
			xSlider[j] = TPoint(d);
			for (unsigned int k = 0; k < d; k++){xSlider[j][k] = x[k][j*shortFolds]; x[k].erase(x[k].begin() + j*shortFolds);}
			ySlider[j] = y[j*shortFolds]; y.erase(y.begin() + j*shortFolds);
			difference -= ySlider[j]; if (ySlider[j] > 0){numMore--;}else{numLess--;}			
		}
		n -= sliderSize;
		/* Cross-validation for the (i + 1)-th space extension (inner cross-validation loop) */
		unsigned int error = 0; TPoint p(0);
		double tmpXSliderVal; int tmpYSliderVal;
		for (unsigned int j = 0; j < folds; j++){
			/* Estimate the current cut */			
			Alpha(&p);
			TVariables res(0);
			Classify(xSlider, p, &res);
			for (unsigned int k = 0; k < sliderSize; k++){error += abs(res[k] - ySlider[k])/2;}
			/* Increment the pointer */
			if (j == shortFolds){break;}
			/* Replace the slider */
			if (j == chSizeVal){
				for (unsigned int l = 0; l < d; l++){x[l].push_back(xSlider[sliderSize - 1][l]);} y.push_back(ySlider[sliderSize - 1]);
				n++; difference += ySlider[sliderSize - 1]; if (ySlider[sliderSize - 1] > 0){numMore++;}else{numLess++;}
				sliderSize--; xSlider.erase(xSlider.begin() + sliderSize); ySlider.erase(ySlider.begin() + sliderSize);
//				for (unsigned int j = 0; j < d; j++){x[j].shrink_to_fit();} y.shrink_to_fit(); - IT IS TOO DANGEROUS
			}
			for (unsigned int k = 0; k < sliderSize; k++){
				for (unsigned int l = 0; l < d; l++){tmpXSliderVal = x[l][k*shortFolds + j]; x[l][k*shortFolds + j] = xSlider[k][l]; xSlider[k][l] = tmpXSliderVal;}
				difference += ySlider[k]; if (ySlider[k] > 0){numMore++;}else{numLess++;}
				tmpYSliderVal = y[k*shortFolds + j]; y[k*shortFolds + j] = ySlider[k]; ySlider[k] = tmpYSliderVal;
				difference -= ySlider[k]; if (ySlider[k] > 0){numMore--;}else{numLess--;}
			}
		}
		/* Check if we've got a better result */ 
		if (error < optError){optError = error; optDegree = i + 1; if (optError == 0){break;}}
	}

	OUT_ALPHA = oldOUT_ALPHA;
	/* Eventually get the classification ray */
	Initialization(spaceExtensions[optDegree - 1], output, minFeatures); // initialize
	power[0] = optDegree;
	return Alpha(ray);
}
