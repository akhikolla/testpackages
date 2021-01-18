/*
  File:             AlphaProcedure.h
  Created by:       Pavlo Mozharovskyi
  First published:  28.02.2013
  Last revised:     28.02.2013

  Contains the modified alpha-procedure for the DDalpha-classifier.

  For a description of the algorithm, see:
    Lange, T., Mosler, K. and Mozharovskyi, P. (2012). Fast nonparametric classification based on data depth. Statistical Papers.
    Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world data with the DDalpha-procedure. Mimeo.
*/

int ExtendWithProducts(TMatrix x, unsigned int upToPower, TMatrix *_x);
int Learn(TMatrix input, TVariables output, unsigned int minFeatures, TPoint *ray);
int LearnCV(TMatrix input, TVariables output, unsigned int minFeatures, unsigned int upToPower, unsigned int folds, TPoint *ray, unsigned int *power);
int Classify(TMatrix input, TPoint weights, TVariables *output);
