/*
File:             DKnn.h
Created by:       Oleksii Pokotylo
First published:  
Last revised:     

The realization of the Depth-based KNN classifier of Paindaveine and Van Bever (2015).
*/

int DKnnCv(TDMatrix points, int n, int d, int* labels, int kMax, int depthType, int chunkNumber);

void DKnnClassify(TDMatrix points, int n, int d, int* labels, TDMatrix objects, int nobjects, int k, int depthType, int* classes);
