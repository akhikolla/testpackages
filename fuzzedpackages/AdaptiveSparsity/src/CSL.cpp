#include "CSL.h"

using namespace Rcpp;
using namespace std;
// using namespace arma doesn't seem to be recognized at least on Windows builds
// doing arma:: instead
using namespace arma;

////////////////////// OLS Estimation //////////////////////
arma::mat armaOLSestim (arma::mat x) {
  int m = x.n_rows;
	int n = x.n_cols;
  arma::mat mu = mean(x,0);
  arma::mat muMat = repmat(mu, m, 1);
  arma::mat LHat = zeros<arma::mat>(n,n);
	n = n-1;
	LHat(n, n) = 1/(sum(square(x.col(n) - mu(n)))/m);
	for (int i = n-1; i >= 0; i--) {
	    arma::mat X = x.cols(i+1, n) - muMat.cols(i+1, n);
	    arma::mat b = inv(trans(X)*X) * trans(X)*(x.col(i)-mu(i));
	    double s = sqrt(sum(square((x.col(i) - X*b - mu(i))))/m);
	    LHat.submat(i+1, i, n, i) = -b/2;
	    LHat(i,i) = 1/s;
	}
	return LHat;
}

SEXP olsInit(SEXP inputMatX){
	NumericMatrix data(inputMatX);
	arma::mat X(data.begin(), data.nrow(), data.ncol(), false);
	arma::mat olsL = armaOLSestim(X);
	return List::create(Named("OLS L") = olsL);
}


////////////////////// Adaptive Sparsity Algorithm //////////////////////
arma::mat LtoBeta (arma::mat LHat, arma::mat sigma) {
  int kNodes = LHat.n_rows;
  arma::mat betaMatrix = zeros(kNodes,kNodes);
	betaMatrix(kNodes-1, kNodes-1) = -1;
	for (int j = 0; j < (kNodes-1); j++) {
	    betaMatrix.col(j) = -LHat.col(j) * sigma(j);
	    betaMatrix(j,j) = -1;
	}
	return betaMatrix;
}

arma::mat betatoL (arma::mat betaMatrix, arma::mat sigma) {
  int kNodes = betaMatrix.n_rows;
	arma::mat LTemp = zeros(kNodes, kNodes);
	for (int j = 0; j <= (kNodes-1); j++) {
	    LTemp.col(j) = -betaMatrix.col(j)*1/sigma(j);
	}
	return LTemp;
}

double Q (arma::mat omegaTemp, arma::mat omegaHat, arma::mat w, int nSamples, arma::mat S, double epsilon) {
	int kNodes = omegaTemp.n_rows;
	double Q1a = 0;
	double Q1b = 0;
	double Q2 = 0;
    for (int i = 0; i <= (kNodes-1); i++) {
        Q1a = Q1a + nSamples*log(w(i));
        for (int j = 0; j <= (kNodes-1); j++) {
            Q1b = Q1b - omegaTemp(i,j)*S(i,j);
        }
        if (i != 0) {
            for (int j = 0; j <= (i-1); j++) {
                Q2 = Q2 + pow(omegaTemp(i,j),2)/(pow(omegaHat(i,j),2)+epsilon);
            }
        }
    }
    double Qout = -(Q1a + Q1b - Q2);
    return Qout;
}

arma::mat betaMatrixEstim (int nSamples, int kNodes, arma::mat betaMatrixTemp, arma::mat sigmaTemp, arma::mat omegaHat, double epsilon, arma::mat S) {
    arma::mat omegaTemp = omegaHat;
    for (int a = 1; a <= kNodes-1; a++) {
		// int mod = 20; // output counter every 20 iterations of beta estimation
		// int amod = a - floor(a/mod)*mod;
		// if (amod == 1) {
			// Rcout << a << endl;
		// }

        for (int b = 0; b <= a-1; b++) {
            double Q1c1 = 0;
            double Q1c0 = 0;
            Q1c1 = 2*S(a,a)/pow(sigmaTemp(b),2);
            arma::mat betaVector = betaMatrixTemp.col(b);
            betaVector(b) = 0;
            arma::mat tempQ1c0 = -2*(S(a,b) - trans(S.col(a))*betaVector)/(pow(sigmaTemp(b),2)) - Q1c1*betaMatrixTemp(a,b);
            Q1c0 = tempQ1c0(0,0);
            
	        double Q2c1part1 = 0;
	        double Q2c1part2 = 0;
	        double Q2c0part1 = 0;
	        double Q2c0part2 = 0;

	        for (int j = b; j <= a-1; j++) {
	            Q2c1part1 = Q2c1part1 + 2*pow(betaMatrixTemp(j,b),2) / ((pow(omegaHat(a,j),2)+epsilon) * pow(sigmaTemp(b),4));
	            Q2c0part1 = Q2c0part1 + betaMatrixTemp(j,b) / ((pow(omegaHat(a,j),2) + epsilon)* pow(sigmaTemp(b),2)) * omegaTemp(a,j);
	        }

	        if (a < kNodes-1) {
	            for (int i = a+1; i <= kNodes-1; i++) {
	                Q2c1part2 = Q2c1part2 +2*pow(betaMatrixTemp(i,b),2) / ((pow(omegaHat(i,a),2) + epsilon) * pow(sigmaTemp(b),4));
	                Q2c0part2 = Q2c0part2 + betaMatrixTemp(i,b)/((pow(omegaHat(i,a),2) + epsilon)*pow(sigmaTemp(b),2))*(omegaTemp(i,a));
	            }
	        }

	        double Q2c1 = Q2c1part1 + Q2c1part2;
	        double Q2c0 = 2*(Q2c0part1 + Q2c0part2) - Q2c1*betaMatrixTemp(a,b);
	        double c1 = Q1c1 + Q2c1;
	        double c0 = Q1c0 + Q2c0;
	        betaMatrixTemp(a,b) = -c0/c1;
	        arma::mat LTemp = betatoL(betaMatrixTemp, sigmaTemp);
	        omegaTemp = LTemp * trans(LTemp);
        }
    }
    return betaMatrixTemp;
}

arma::mat wEstim(int nSamples, int kNodes, arma::mat betaMatrixTemp, arma::mat w, arma::mat omegaHat, double epsilon, arma::mat S) {
    arma::mat wtemp = zeros(1, kNodes);
    for (int m = 0; m <= kNodes-1; m++) {
        double c0 = nSamples;
        double c1lkhd = 0;
        double c1prior = 0;
        double c2 = 0;
        for (int i = m; i <= kNodes-1; i++) {
            for (int j = 0; j <= kNodes-1; j++) {
                c1lkhd = c1lkhd - S(i,j)*betaMatrixTemp(i,m)*betaMatrixTemp(j,m);
            }
        }
        if (m != kNodes) {
            for (int j = m; j <= kNodes-2; j++) {
                for (int i = j+1; i <= kNodes-1; i++) {
                    c2 = c2 - 2*pow(betaMatrixTemp(i,m),2)*pow(betaMatrixTemp(j,m),2)/(pow(omegaHat(i,j),2) + epsilon);
                    double tempsum = 0;
                    for (int k = 0; k <= j; k++) {
                        if (k!= m) {
                            tempsum = tempsum + betaMatrixTemp(i,k)*betaMatrixTemp(j,k)*w(k);
                        }
                    }
                    c1prior = c1prior - 2/(pow(omegaHat(i,j),2)+epsilon)*tempsum*betaMatrixTemp(i,m)*betaMatrixTemp(j,m);
                }
            }
        }
        double c1 = c1lkhd + c1prior;
        double testIsReal = pow(c1,2) - 4*c2*c0;
        if (testIsReal >= 0) {
            if (c2  == 0) {
                double Root1 = -c0/c1;
                wtemp(m) = Root1;
            } else {
                double Root1 = (-c1 + sqrt(c1*c1 - 4*c0*c2))/(2*c2);
                double Root2 = (-c1 - sqrt(c1*c1 - 4*c0*c2))/(2*c2);
                if (Root1 >= 0 && (2*c2*Root1 + c1) < 0) {
                    wtemp(m) = Root1;
                } else if (Root2 >= 0 && (2*c2*Root2 + c1) < 0) {
                    wtemp(m) = Root2;
                } else {
					wtemp(m) = epsilon;
				}
            }
        }
    }
    return wtemp;
}

arma::mat runCSL(int it1, arma::mat LHat, int kNodes, int nSamples, double epsilon, arma::mat S, arma::mat ansQ) {
	arma::mat omegaHat = LHat * trans(LHat);
  arma::mat sigma = 1/LHat.diag(0);
  arma::mat w = pow(LHat.diag(0),2);
	arma::mat Beta = LtoBeta(LHat, sigma);
	
	double Qbegin;
	double Qend;
	double omegaDiff;
	double err_f;
	
	for (int iter = 0; iter < it1; iter++) {
		Rcout << "***** " << iter+1 << " *****" << endl;
		arma::mat omegaHatOld = omegaHat;
		Qbegin = Q(omegaHat, omegaHat, w, nSamples, S, epsilon);
		Rcout << Qbegin << endl;
		Beta = betaMatrixEstim (nSamples, kNodes, Beta, sigma, omegaHat, epsilon, S);
		w = wEstim(nSamples, kNodes, Beta, w, omegaHat, epsilon, S);
		sigma = sqrt(1/w);
		LHat = betatoL(Beta, sigma);
		arma::mat omegaTemp = LHat * trans(LHat);
		Qend = Q(omegaTemp, omegaHat, w, nSamples, S, epsilon);
		Rcout << Qend << endl;
		omegaHat = omegaTemp;
		
		err_f = sqrt(sum(sum(pow(abs(omegaHat - ansQ),2))));
		Rcout << err_f << endl;
		
		omegaDiff = max(max(abs(omegaHatOld - omegaHat)));
		Rcout << omegaDiff << endl;
		
		if (kNodes >= 150) {
			omegaHat.save("omegaHat.csv", csv_ascii);
		}
		
		if (omegaDiff < epsilon) {
			return omegaHat;
        }
	}
	
	return omegaHat;
}

SEXP CSL(SEXP it, SEXP init, SEXP k, SEXP n, SEXP eps, SEXP covMat, SEXP ansQMat){
	//iterations, LInit, kNodes, nSamples, epsilon, S
	NumericVector itVector(it);
	NumericMatrix initMatrix(init);
	NumericVector kVector(k);
	NumericVector nVector(n);
	NumericVector epsVector(eps);
	NumericMatrix covMatrix(covMat);
	NumericMatrix ansQMatrix(ansQMat);
	
	int iterations = itVector[0];
	arma::mat LInit(initMatrix.begin(), initMatrix.nrow(), initMatrix.ncol(), false);
	int kNodes = kVector[0];
	int nSamples = nVector[0];
	double epsilon = epsVector[0];
	arma::mat S(covMatrix.begin(), covMatrix.nrow(), covMatrix.ncol(), false);
	arma::mat ansQ(ansQMatrix.begin(), ansQMatrix.nrow(), ansQMatrix.ncol(), false);
	
	arma::mat omegaHat = runCSL(iterations, LInit, kNodes, nSamples, epsilon, S, ansQ);
	
	List res ;                                  
	res["omegaHat"] = omegaHat;
	return res ;
}
