/*
  File:             Polynomial.h
  Created by:       Oleksii Pokotylo
  First published:  07.05.2014
  Last revised:     07.05.2014

  Contains the polynomial classifier the DD-plot classification.

  For a description of the algorithm, see:
    Li, J., Cuesta-Albertos, J. A. and Liu, R. Y. (2012). DD-classifier: Nonparametric classification procedure based on
DD-plot, Journal of the American Statistical Association 107(498): 737 - 753.
*/

TPoint PolynomialLearnCV(TDMatrix input, int numClass1, int numClass2, int maxDegree, int chunkNumber, int *degree, int *axis);
