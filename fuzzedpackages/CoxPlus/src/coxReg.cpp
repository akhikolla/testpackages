/*
coxReg.cpp
 * To do List
 * 1) Warning: May need to zero S1, S2, and riskset if the code is modified later on
 * 2) Singularity protection while inv or solve: maybe pivoting. read source code for inv and solve. Consequence, robust wald is not accurate.
 * Example, two id data used, results are different,all id are the same
 * 3) Double check potential errors (e.g., score residual split, grouping, how to choose delta_i, connection to EM)
 * 7) Remove NA in beta and I, which.sing= diag(v). Check singularity while inv or solve
 * 8) Template programming for transformation, Event to Integer vector
 * 9) Test statistics for frailty models and df computation (maybe not integer)
 * 10) For Rcpp and Armadillo, You cannot use multiple indices in a single [ ] expression!!!, use \[[^\]]+,.*?\] to find out all such wrong usage like [i,j]
 * 11) Armadillo can't save unsigned type like uvec
 * 12) Theta is hard to estimate as L is not sensitive to theta (too small in L)
 *

 */
#if !defined(DEBUG_MODE)
	#define ARMA_NO_DEBUG
#endif
//#define ARMA_BLAS_CAPITALS
#define ARMA_USE_BLAS
#define ARMA_USE_LAPACK
//#define ARMA_64BIT_WORD
// Use of BLAS and LAPACK and enabled by default, the above is for safety. Note that solve requires lapack
//#define ARMA_DONT_USE_LAPACK
//#define ARMA_DONT_USE_BLAS
#include <iostream>
#include <cstdio>
#include <RcppArmadillo.h>
#include <cmath>
#include <string>
#include <fstream>
#include <limits>
#include <assert.h>
#include <stdexcept>
//using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
class CoxReg {
	//Input, data sorted by strata, event, stop, start time, in descending order
	vec start, stop, event, weight;
	mat X; // p*n, column major
	Mat<int> strataBounds; // nStrata*2, starting from 0
	Mat<int> interMap; // support high order interactions, but might be slower
	umat fgroups; // consecutive group numbers starting from 0, plus the offset of that frailty in beta
	ivec thetagroups; // consecutive group numbers starting from 0
	uvec egroups; // consecutive group numbers starting from 0
	List par;
	int n; // number of observations
	int p; // number of covariates
	int nf; // number of frailty terms
	double eps;
	double tol;
	double tol_par;
	int maxIter, maxOutIter, thetaIterMax;
	int method; // 0-Breslow, 1-Efron, 2-Collective Cause, 3-Single Cause
	int verbose;
	int isREML;
	int isCutBeta;
	int isCutStep;
	int isUpdateDiggNum; // useful for comparison with SAS only
	int isTieGroup;
	int maxStepHalving;
	int fixedCoefIndex;
	int ageIndex; // index of age in timeVarIndex
	int diggNumIndex; // index of diggNum in timeVarIndex
	int isAbsoluteAge; // time or time-imprtime
	int isFixedEffect; // fixed effect (1), random effect (-1), both (0)
	int randStarts; //starting from this (theta) group, it will be random effect
	int isDiagFrailty;
	int isForceCheck;
	int isDebug; // only save data when debug enabled

	int nThetaGroups; // number of theta groups
	ivec dropFlag;
	ivec timeVarIndex;
	uvec keepIndex;
	vec offset;
	vec delta; // not used if beta_last is given
	vec beta_last; // use this beta to estimate delta
	uvec validFrailty; //For each fixed effect, one degree of freedom needs to be disabled
	double fixedCoef;
	double eigScalor;
	double eigPct;
	double minstop, maxstop; // an alternative way for censoring

	//Results
	vec beta;
	vec theta;
	double theta_se;
	vec se;
	vec se_robust;
	mat R;
	vec U;
	mat I;
	mat V_conv; // inv(I))
	mat V_robust; // Robust variance
	mat D;
	int converged;
	double L;
	double E; // entropy of weights distribution, for EM algorithm
	double scoreTest, LRTest, waldTest;
	List res;
	List meta;

public:
	int deltaType; // work together with beta_last to estimate delta, -1 min, 1 max, 0 proportional

	CoxReg() {
	}

	CoxReg(std::string path) {
		loadInput(path);
	}

	CoxReg(NumericMatrix head, NumericMatrix Xi, List par) {
		init(num2mat(head), num2mat(Xi), par);
	}

	void saveInput() {
		Rcout << "======Saving Input======" << endl;
		//system("exec rm -r save/*"); // removed previous files in the folder
		start.save("save/start.mat");
		stop.save("save/stop.mat");
		event.save("save/event.mat");
		weight.save("save/weight.mat");
		X.save("save/X.mat");
		strataBounds.save("save/strataBounds.mat");
		if (!interMap.is_empty()) interMap.save("save/interMap.mat");
		if (!fgroups.is_empty()) conv_to<imat>::from(fgroups).save("save/fgroups_imat.mat");
		if (!thetagroups.is_empty()) thetagroups.save("save/thetagroups.mat");
		if (!egroups.is_empty()) conv_to<ivec>::from(egroups).save("save/egroups_ivec.mat");

		ivec par = zeros<ivec>(23);
		par[0] = n;
		par[1] = p;
		par[2] = nf;
		par[3] = maxIter;
		par[4] = maxOutIter;
		par[5] = thetaIterMax;
		par[6] = method;
		par[7] = verbose;
		par[8] = isREML;
		par[9] = isCutBeta;
		par[10] = isCutStep;
		par[11] = isUpdateDiggNum;
		par[12] = isTieGroup;
		par[13] = maxStepHalving;
		par[14] = fixedCoefIndex;
		par[15] = ageIndex;
		par[16] = diggNumIndex;
		par[17] = isAbsoluteAge;
		par[18] = isFixedEffect;
		par[19] = randStarts;
		par[20] = nThetaGroups;
		par[21] = isDiagFrailty;
		par[22] = isForceCheck;
		par.save("save/par.mat");

		if (!dropFlag.is_empty()) dropFlag.save("save/dropFlag.mat");
		if (!timeVarIndex.is_empty()) timeVarIndex.save("save/timeVarIndex.mat");
		if (!keepIndex.is_empty()) conv_to<ivec>::from(keepIndex).save("save/keepIndex_ivec.mat");
		if (!offset.is_empty()) offset.save("save/offset.mat");
		if (!delta.is_empty()) delta.save("save/delta.mat");
		if (!beta_last.is_empty()) beta_last.save("save/beta_last.mat");
		if (!validFrailty.is_empty()) conv_to<ivec>::from(validFrailty).save("save/validFrailty_ivec.mat");
		beta.save("save/beta.mat");
		if (!theta.is_empty()) theta.save("save/theta.mat");

		vec parf = zeros(6);
		parf[0] = eps;
		parf[1] = tol;
		parf[2] = tol_par;
		parf[3] = fixedCoef;
		parf[4] = eigScalor;
		parf[5] = eigPct;
		parf[6] = minstop;
		parf[7] = maxstop;
		parf.save("save/parf.mat");
	}

	void loadInput(std::string path) {
		Rcout << "======Loading Input from " << path << " ======" << endl;
		start.load(path + "start.mat");
		stop.load(path + "stop.mat");
		event.load(path + "event.mat");
		weight.load(path + "weight.mat");
		X.load(path + "X.mat");
		strataBounds.load(path + "strataBounds.mat");
		if (exist(path + "interMap.mat")) interMap.load(path + "interMap.mat");
		if (exist(path + "fgroups_imat.mat")) {
			imat fgroups_imat;
			fgroups_imat.load(path + "fgroups_imat.mat");
			fgroups = conv_to<umat>::from(fgroups_imat);
		}
		if (exist(path + "thetagroups.mat")) thetagroups.load(path + "thetagroups.mat");
		if (exist(path + "egroups_ivec.mat")) {
			ivec egroups_ivec;
			egroups_ivec.load(path + "egroups_ivec.mat");
			egroups = conv_to<uvec>::from(egroups_ivec);
		}

		ivec par;
		par.load(path + "par.mat");
		n = par[0];
		p = par[1];
		nf = par[2];
		maxIter = par[3];
		maxOutIter = par[4];
		thetaIterMax = par[5];
		method = par[6];
		verbose = par[7];
		isREML = par[8];
		isCutBeta = par[9];
		isCutStep = par[10];
		isUpdateDiggNum = par[11];
		isTieGroup = par[12];
		maxStepHalving = par[13];
		fixedCoefIndex = par[14];
		ageIndex = par[15];
		diggNumIndex = par[16];
		isAbsoluteAge = par[17];
		isFixedEffect = par[18];
		randStarts = par[19];
		nThetaGroups = par[20];
		isDiagFrailty = par[21];
		isForceCheck = par[22];

		if (exist(path + "dropFlag.mat")) dropFlag.load(path + "dropFlag.mat");
		if (exist(path + "timeVarIndex.mat")) timeVarIndex.load(path + "timeVarIndex.mat");
		if (exist(path + "keepIndex_ivec.mat")) {
			ivec keepIndex_ivec;
			keepIndex_ivec.load(path + "keepIndex_ivec.mat");
			keepIndex = conv_to<uvec>::from(keepIndex_ivec);
		}

		if (exist(path + "offset.mat")) offset.load(path + "offset.mat");
		if (exist(path + "delta.mat")) delta.load(path + "delta.mat");
		if (exist(path + "beta_last.mat")) beta_last.load(path + "beta_last.mat");
		if (exist(path + "validFrailty_ivec.mat")) {
			ivec validFrailty_ivec;
			validFrailty_ivec.load(path + "validFrailty_ivec.mat");
			validFrailty = conv_to<uvec>::from(validFrailty_ivec);
		}

		beta.load(path + "beta.mat");
		if (exist(path + "theta.mat")) theta.load(path + "theta.mat");

		vec parf;
		parf.load(path + "parf.mat");
		eps = parf[0];
		tol = parf[1];
		tol_par = parf[2];
		fixedCoef = parf[3];
		eigScalor = parf[4];
		eigPct = parf[5];
		minstop = parf[6];
		maxstop = parf[7];
	}

	inline bool exist(const std::string& name) {
		std::ifstream file(name.c_str());
		if (!file) //if the file was not found, then file is 0, i.e. !file=1 or true
			return false;
		else
			return true;
	}

	void init(mat m, mat Xi, List par) {
		Rcout << "=====Initializing=====\n";
		start = m.col(0);
		stop = m.col(1);
		event = m.col(2);
		weight = m.col(3);
		X = Xi;
		this->par = par;
		n = start.n_elem;
		isDebug = (par.containsElementNamed("isDebug")) ? as<int> (par["isDebug"]) : 0;
		if (par.containsElementNamed("delta")) delta = num2vec(as<NumericVector> (par["delta"]));
		if (par.containsElementNamed("beta_last")) beta_last = num2vec(as<NumericVector> (par["beta_last"]));
		nf = (par.containsElementNamed("nf")) ? as<int> (par["nf"]) : 0;
		if (par.containsElementNamed("egroups")) egroups = int2uvec(as<IntegerVector> (par["egroups"]));
		if (nf > 0) {
			fgroups = conv_to<umat>::from(int2imat(as<IntegerMatrix> (par["fgroups"])));
			theta = (par.containsElementNamed("theta")) ? num2vec(as<NumericVector> (par["theta"])) : ones(nf);
			if (par.containsElementNamed("thetaScale")) theta *= as<double> (par["thetaScale"]);
			thetagroups = int2ivec(as<IntegerVector> (par["thetagroups"]));
			isFixedEffect = as<int> (par["isFixedEffect"]);
			randStarts = as<int> (par["randStarts"]);
			nThetaGroups = max(thetagroups);
			maxOutIter = (par.containsElementNamed("maxOutIter")) ? as<int> (par["maxOutIter"]) : 50;
			thetaIterMax = (par.containsElementNamed("thetaIterMax")) ? as<int> (par["thetaIterMax"]) : 100;
			isREML = (par.containsElementNamed("isREML")) ? as<int> (par["isREML"]) : 0;
			isDiagFrailty = (par.containsElementNamed("isDiagFrailty")) ? as<int> (par["isDiagFrailty"]) : 1;
		}
		if (par.containsElementNamed("strataBounds")) {
			strataBounds = int2imat(as<IntegerMatrix> (par["strataBounds"]));
		} else {
			strataBounds = Mat<int> (1, 2);
			strataBounds(0, 0) = 0;
			strataBounds(0, 1) = n - 1;
		}
		maxIter = (par.containsElementNamed("maxIter")) ? as<int> (par["maxIter"]) : 50;
		method = (par.containsElementNamed("method")) ? as<int> (par["method"]) : 1; // default efron
		if (method > 3 || method < 0) {
			Rcout << "Error: invalid method type " << method << "!!!\n";
			throw std::invalid_argument( "Invalid method type!" );
		}
		verbose = (par.containsElementNamed("verbose")) ? as<int> (par["verbose"]) : 0;
		deltaType = (par.containsElementNamed("deltaType")) ? as<int> (par["deltaType"]) : 0;
		isCutBeta = (par.containsElementNamed("isCutBeta")) ? as<int> (par["isCutBeta"]) : 0;
		isCutStep = (par.containsElementNamed("isCutStep")) ? as<int> (par["isCutStep"]) : 0;
		isAbsoluteAge = (par.containsElementNamed("isAbsoluteAge")) ? as<int> (par["isAbsoluteAge"]) : 1;
		isUpdateDiggNum = (par.containsElementNamed("isUpdateDiggNum")) ? as<int> (par["isUpdateDiggNum"]) : 1;
		ageIndex = (par.containsElementNamed("ageIndex")) ? as<int> (par["ageIndex"]) : -1;
		diggNumIndex = (par.containsElementNamed("diggNumIndex")) ? as<int> (par["diggNumIndex"]) : -1;
		isTieGroup = (par.containsElementNamed("isTieGroup")) ? as<int> (par["isTieGroup"]) : 0;
		maxStepHalving = (par.containsElementNamed("maxStepHalving")) ? as<int> (par["maxStepHalving"]) : 10;
		fixedCoefIndex = (par.containsElementNamed("fixedCoefIndex")) ? as<int> (par["fixedCoefIndex"]) : -1;
		fixedCoef = (par.containsElementNamed("fixedCoef")) ? as<double> (par["fixedCoef"]) : 0;
		isForceCheck = (par.containsElementNamed("isForceCheck")) ? as<int> (par["isForceCheck"]) : 0;
		if (par.containsElementNamed("dropFlag")) dropFlag = int2ivec(as<IntegerVector> (par["dropFlag"]));
		if (par.containsElementNamed("timeVarIndex")) timeVarIndex = int2ivec(as<IntegerVector> (par["timeVarIndex"]));
		if (par.containsElementNamed("interMap")) interMap = int2imat(as<IntegerMatrix> (par["interMap"]));
		if (par.containsElementNamed("offset")) offset = num2vec(as<NumericVector> (par["offset"]));
		eps = (par.containsElementNamed("eps")) ? as<double> (par["eps"]) : 1e-9; // The maximal distinguishable difference is around 1e-15
		tol = (par.containsElementNamed("tol")) ? as<double> (par["tol"]) : 1e-6;
		tol_par = (par.containsElementNamed("tol_par")) ? as<double> (par["tol_par"]) : 1e-4;
		eigScalor = (par.containsElementNamed("eigScalor")) ? as<double> (par["eigScalor"]) : 2;
		eigPct = (par.containsElementNamed("eigPct")) ? as<double> (par["eigPct"]) : 1e-2;
		minstop = (par.containsElementNamed("minstop")) ? as<double> (par["minstop"]) : -DBL_MAX;
		maxstop = (par.containsElementNamed("maxstop")) ? as<double> (par["maxstop"]) : DBL_MAX;
		p = X.n_rows;
		if (!dropFlag.is_empty()) {
			keepIndex = find(dropFlag == 0);
			p = keepIndex.n_elem;
		}
		if (par.containsElementNamed("beta")) beta = num2vec(as<NumericVector> (par["beta"]));
		else beta = zeros(p + nf);
		if (nf > 0) {
			fgroups = fgroups + p; // stores index of each frailty in beta
			if (isFixedEffect >= 0) { // 0 - both, 1 - fixed
				ivec vf = ones<ivec> (p + nf);
				vf(span(p, p + nf - 1)) = int2ivec(as<IntegerVector> (par["validFrailty"]));
				validFrailty = find(vf);
				find(vf - 1).t().print(Rcout, "Fixed effects reference levels = ");
			}
			Rcout << "n=" << n << ", p=" << p << ", nf=" << nf << ", isFixedEffect=" << isFixedEffect
					<< ", randStarts=" << randStarts << ", unique theta groups=" << nThetaGroups + 1 << endl;
		}
		if(isDebug==1) saveInput();
	}

	//vec(aux_mem*, number_of_elements, copy_aux_mem = true, strict = true)

	vec num2vec(NumericVector v) {
		return vec(v.begin(), v.size(), false);
	}

	uvec int2uvec(IntegerVector v) {
		uvec res = zeros<uvec>(v.size());
		for (int i = 0; i < v.size(); i++) res[i] = (uword) (v[i]);
		return res;
	}

	ivec int2ivec(IntegerVector v) {
		ivec res = zeros<ivec>(v.size());
		for (int i = 0; i < v.size(); i++) res[i] = v[i];
		return res;
	}

	//mat(aux_mem*, n_rows, n_cols, copy_aux_mem = true, strict = true)

	mat num2mat(NumericMatrix m) {
		return mat(m.begin(), m.nrow(), m.ncol(), false);
	}

	Mat<int> int2imat(IntegerMatrix m) {
		return Mat<int> (m.begin(), m.nrow(), m.ncol(), false);
	}

	inline void centerX(mat &m) {
		m.each_col() -= mean(m, 1);
	}

	void compileResult() {
		#if defined(DEBUG_MODE)
			beta(span(0, fmin(9, beta.n_elem - 1))).t().print(Rcout, "Final beta=");
			//se(span(0, fmin(9, beta.n_elem - 1))).t().print(Rcout, "Final se=");
		#endif
		if (nf == 0)
			res = List::create(_["p"] = p, _["nf"] = nf, _["beta"] = beta, _["se"] = se, _["rse"] = se_robust, _["wald"] = waldTest, _["LR"] = LRTest, _["score"] = scoreTest, _["converged"] = converged,
					_["R"] = R, _["U"] = U, _["I"] = I, _["D"] = D, _["V"] = V_conv, _["Vr"] = V_robust, _["meta"] = meta, _["likelihood"] = L, _["entropy"] = E);
		else
			res = List::create(_["p"] = p, _["nf"] = nf, _["beta"] = beta, _["se"] = se, _["wald"] = waldTest, _["LR"] = LRTest, _["score"] = scoreTest, _["converged"] = converged,
					_["R"] = R, _["U"] = U, _["I"] = I, _["D"] = D, _["V"] = V_conv, _["theta"] = theta, _["meta"] = meta, _["likelihood"] = L, _["entropy"] = E, _["theta_se"] = theta_se);
		Rcout << "Results compiled!\n";
	}

	inline void setMax(vec &v, double t) {
		uvec ix = find(v > t);
		if (verbose >= 3) v(ix).print(Rcout, "overflow values:");
		v(ix).fill(t);
	}

	inline void setMin(vec &v, double t) {
		uvec ix = find(v < t);
		if (verbose >= 3) v(ix).print(Rcout, "underflow values:");
		v(ix).fill(t);
	}

	inline bool isFinite(double x) {
		return (x <= DBL_MAX && x >= -DBL_MAX);
	}

	inline double getEntropy(vec v) {
		double e = 0;
		for (int i = 0; i < (int) v.n_elem; i++) e += v[i] > 0 ? -v[i] * log(v[i]) : 0;
		return e;
	}

	void estimateNormal() {
		converged = 0;
		double Ln = -DBL_MAX;
		double L0 = 0;
		vec beta0 = beta; // "=" results in clone by default
		vec betaOld = beta;
		vec Un;
		mat In;
		uword* riskset = new uword[n];
		int nrisk = 0;
		int ntied = 0;
		double time, time_lb, time_ub;
		vec S1 = zeros(p);
		mat S2 = zeros(p, p);
		int a, b;
		double tmp;
		int ntrials = 0;
		vec fixedCoefOffset;
		// Strata by strata, data ordered by -strata, event, stop, and start in descending order
		for (int iter = 1; iter <= maxIter; iter++) {
			L = 0;
			E = 0;
			U = zeros(p);
			I = zeros(p, p);
			for (int s = 0; s < (int) (strataBounds.n_rows); s++) {
				// Traverse events within a strata
				for (int i = strataBounds(s, 0); i <= strataBounds(s, 1); i++) {
					if (event[i] == 0) break;
					time = stop[i];
					if(time<minstop - eps) continue;
					if(time>maxstop + eps) break;
					time_lb = time * (1 - eps); // numerical underflow may cause problem; that happens, change the inequality below
					time_ub = time * (1 + eps);
					nrisk = 0;
					ntied = 0;
					for (int j = i; j <= strataBounds(s, 1); j++) {
						if (stop[j] > time_lb && start[j] < time) {
							riskset[nrisk++] = j;
							if (event[j] == 1 && stop[j] < time_ub) ntied++;
						} else if (event[j] == 0 && stop[j] < time) break;
					}
					for (int j = strataBounds(s, 0); j < i; j++) {
						if (stop[j] >= time && start[j] < time) riskset[nrisk++] = j;
					}
					if (ntied <= 0) {
						Rcout << "Error!!! ntied=" << ntied << ", start=" << start[i] << ", stop=" << stop[i] << "\n";
						throw std::invalid_argument( "Logical error: ntied<=0!" );
					}
					i += ntied - 1;
					// Update X
					uvec v = uvec(riskset, nrisk, false, true);
					vec wr = weight(v);
					vec wt = wr(span(0, ntied - 1));
					mat Xr = X.cols(v);
					// Deal with time-varying variables: Digg number and time, if subtracting mean before interaction, the results may be different with SAS
					if (!timeVarIndex.is_empty()) {
						if (ageIndex >= 0) {
							if (isAbsoluteAge == 1) Xr.row(timeVarIndex[ageIndex]).fill(time / 86400.0);
							else Xr.row(timeVarIndex[ageIndex]) = (time - Xr.row(timeVarIndex[ageIndex])) / 86400.0;
						}
						if (diggNumIndex >= 0 && isUpdateDiggNum == 1) {
							Xr.row(timeVarIndex[diggNumIndex]).fill(X(timeVarIndex[diggNumIndex], i - ntied + 1));
						}
						for (a = 0; (uword) a < interMap.n_rows; a++) Xr.row(interMap(a, 0)) = Xr.row(interMap(a, 1)) % Xr.row(interMap(a, 2));
					}
					if (fixedCoefIndex >= 0) fixedCoefOffset = fixedCoef * Xr.row(fixedCoefIndex);
					if (!dropFlag.is_empty()) Xr = Xr.rows(keepIndex);
					centerX(Xr);
					mat tXr = trans(Xr);
					mat Xt = Xr.cols(0, ntied - 1);
					// compute XB, eXB, Q, Q1,Q2, S,S1,S2,Eq,Es
					vec xbeta = vec(tXr * beta);
					if (!offset.is_empty()) xbeta += offset(v);
					if (fixedCoefIndex >= 0) xbeta += fixedCoefOffset;
					xbeta -= mean(xbeta);
					if (isCutBeta == 1) {
						setMax(xbeta, 22);
						setMin(xbeta, -200);
					}
					vec exp_xbeta = exp(xbeta);
					vec w_exp_xbeta = wr % exp_xbeta;
					vec tied_exp_xbeta = w_exp_xbeta(span(0, ntied - 1));
					if (method == 2) tied_exp_xbeta = exp_xbeta(span(0, ntied - 1));
					double S = sum(w_exp_xbeta);
					double Q = sum(tied_exp_xbeta);
					for (a = 0; (uword) a < tXr.n_cols; a++) {
						vec cur_col = tXr.col(a) % w_exp_xbeta;
						S1[a] = sum(cur_col);
						for (b = 0; b <= a; b++) {
							tmp = dot(cur_col, tXr.col(b));
							S2.at(a, b) = tmp;
							S2.at(b, a) = tmp;
						}
					}
					vec Q1 = vec(Xt * tied_exp_xbeta);
					mat Q2 = Xt * diagmat(tied_exp_xbeta) * trans(Xt);
					vec Es = S1 / S;
					vec Eq = Q1 / Q;
					// Update L, U, I, V
					if (method == 2) { // collective cause
						// assume weights for tied obs are the same, use mean in case not equal by mistake
						L += mean(wt) * (log(Q) - log(S));
					U += mean(wt)*(Eq - Es);
					I += mean(wt) *(S2 / S - Q2 / Q + Eq * trans(Eq) - Es * trans(Es));
					} else if (method == 0) { // Breslow
						L += dot(wt, xbeta(span(0, ntied - 1))) - sum(wt) * log(S);
						U += vec(Xt * wt) - Es * sum(wt);
						I += sum(wt) *(S2 / S - Es * trans(Es));
					} else if (method == 1) { // Efron
						L += dot(wt, xbeta(span(0, ntied - 1)));
						U += vec(Xt * wt);
						for (int d = 0; d < ntied; d++) {
							double ratio = d / (double) ntied;
							double newS = S - Q * ratio;
							vec newS1 = S1 - Q1 * ratio;
							mat newS2 = S2 - Q2 * ratio;
							vec newEs = newS1 / newS;
							L -= mean(wt) * log(newS);
							U -= mean(wt) * newEs;
							I += mean(wt) *(newS2 / newS - newEs * trans(newEs));
						}
					} else if (method == 3) { // Single cause
						vec dt = ones(ntied) / (double) ntied;
						if (!beta_last.is_empty()) {
							vec delta_est = zeros(ntied);
							vec xbeta_est = vec(trans(Xt) * beta_last);
							if (!offset.is_empty()) xbeta_est += offset(v(span(0, ntied - 1)));
							if (fixedCoefIndex >= 0) xbeta_est += fixedCoefOffset(span(0, ntied - 1));
							uword index;
							if (deltaType == 1) {
								xbeta_est.max(index);
								delta_est[index] = 1;
							} else if (deltaType == -1) {
								xbeta_est.min(index);
								delta_est[index] = 1;
							} else if (deltaType == 0) {
								vec exp_xbeta_est = exp(xbeta_est - max(xbeta_est));
								delta_est = exp_xbeta_est / sum(exp_xbeta_est);
							} else Rcout << "@@@delta type not specified, beta_last info will not be used\n";
							dt = delta_est;
						} else if (!delta.is_empty()) dt = delta(v(span(0, ntied - 1)));
						wt %= dt;
						E += getEntropy(dt);
						L += getEntropy(dt);
						L += dot(wt, xbeta(span(0, ntied - 1))) - sum(wt) * log(S);
						U += vec(Xt * wt) - Es * sum(wt);
						I += sum(wt) *(S2 / S - Es * trans(Es));
					}
				} // end of strata
			} // end of all observations
			// Update beta
			vec dbeta = getDbeta(U, I);
			if (verbose >= 3 && (!isFinite(L) || !dbeta.is_finite())) {
				I.print(Rcout, "I=");
				U.print(Rcout, "U=");
				beta.print(Rcout, "beta=");
				dbeta.print(Rcout, "dbeta=");
			}
			if (fabs((L - Ln) / L) < tol && isFinite(L)) {
				if (verbose >= 1)
					Rcout << "Converged after " << iter << " iterations, relative change in L =" << fabs((L - Ln) / L) << "\n";
				if (L < Ln) { // go back if L<Ln
					L = Ln;
					beta = betaOld;
					U = Un;
					I = In;
				}
				beta(span(0, p - 1)).t().print(Rcout);
				converged = 2;
				break;
			}
			if (L < Ln || !isFinite(L) || !dbeta.is_finite()) {
				if (verbose >= 2) Rcout << "--Over reach in the " << iter << " iterations after " << ntrials << " trials, L="
						<< L << ", dbeta finite: " << dbeta.is_finite() << "\n";
				if (all(abs(beta - betaOld) < abs(beta) * tol_par) || ntrials > maxStepHalving) {
					Rcout << "---Recovered to last best because no better results can be found in the suggested direction\n";
					L = Ln;
					beta = betaOld;
					U = Un;
					I = In;
					converged = 1;
					break;
				}
				beta = (beta + betaOld) / 2;
				iter--; // not converging
				ntrials++;
			} else {
				ntrials = 1;
				if (isCutStep == 1) {
					setMax(dbeta, 1000);
					setMin(dbeta, -1000);
				}
				// Score Test and Likelihood ratio test, be aware of the order of the adjacent lines
				if (iter == 1) {
					scoreTest = dot(U, dbeta); // =trans(U) * dbeta;
					L0 = L;
				}
				betaOld = beta;
				beta += dbeta;
				if (verbose >= 2) {
					Rcout << "-- The " << iter << "/" << maxIter << " iteration finished, L=" << L << ", Ln=" << Ln << ", next beta =\n";
					beta(span(0, p - 1)).t().print(Rcout);
				}
				Ln = L;
				Un = U;
				In = I;
			}
			if (iter == maxIter) Rcout << "Warning: the model is not converged after the maximal " << maxIter << " iterations\n";
		} // end of loop
		delete [] riskset;
		//Likelihood ratio test
		LRTest = 2 * (L - L0);
		//standard errors
		V_conv = getInverse(I);
		se = sqrt(diagvec(V_conv));
		//Robust variance
		getRobustV();
		waldTest = as_scalar(trans(beta - beta0) * I * (beta - beta0)); // if replace I with getInverse(V_robust) for robust Wald, then sensitive to singularity problems
		beta_last = beta; // clone beta for potential two-stage or multi-stage
	}

	List survivalNormal() {
		beta = num2vec(as<NumericVector> (par["beta"]));
		ivec receiver = int2ivec(as<IntegerVector> (par["receiver"]));
		ivec urr = unique(receiver);
		List base; // cumulative baseline hazard by strata
		mat ind = zeros(urr.n_elem, strataBounds.n_rows); // cumulative individual receiver hazard by strata
		mat indvar = zeros(urr.n_elem, strataBounds.n_rows);

		uword* riskset = new uword[n];
		int nrisk = 0;
		int ntied = 0;
		double time, time_lb, time_ub;
		vec fixedCoefOffset;
		// Strata by strata, data ordered by -strata, event, stop, and start in descending order
		for (int s = 0; s < (int) (strataBounds.n_rows); s++) {
			// count number of events
			ivec subReceiver = receiver(span(strataBounds(s, 0),strataBounds(s, 1)));
			vec subEvent = event(span(strataBounds(s, 0),strataBounds(s, 1)));
			ivec rset = unique(subReceiver(find(subEvent)));
			mat db = zeros(rset.n_elem + 1, 3);  // cumulative baseline hazard
			int k = 0; // number of unique event time points
			// Traverse events within a strata
			for (int i = strataBounds(s, 0); i <= strataBounds(s, 1); i++) {
				if (event[i] == 0) break;
				time = stop[i];
				time_lb = time * (1 - eps); // numerical underflow may cause problem; that happens, change the inequality below
				time_ub = time * (1 + eps);
				nrisk = 0;
				ntied = 0;
				for (int j = i; j <= strataBounds(s, 1); j++) {
					if (stop[j] > time_lb && start[j] < time) {
						riskset[nrisk++] = j;
						if (event[j] == 1 && stop[j] < time_ub) ntied++;
					}
				}
				if (ntied <= 0) {
					Rcout << "Error!!! ntied=" << ntied << ", start=" << start[i] << ", stop=" << stop[i] << "\n";
					throw std::invalid_argument( "Logical error: ntied<=0!" );
				}
				i += ntied - 1;
				// Update X
				uvec v = uvec(riskset, nrisk, false, true);
				vec wr = weight(v);
				vec wt = wr(span(0, ntied - 1));
				mat Xr = X.cols(v);
				// Deal with time-varying variables: Digg number and time, if subtracting mean before interaction, the results may be different with SAS
				if (!timeVarIndex.is_empty()) {
					if (ageIndex >= 0) {
						if (isAbsoluteAge == 1) Xr.row(timeVarIndex[ageIndex]).fill(time / 86400.0);
						else Xr.row(timeVarIndex[ageIndex]) = (time - Xr.row(timeVarIndex[ageIndex])) / 86400.0;
					}
					if (diggNumIndex >= 0 && isUpdateDiggNum == 1) {
						Xr.row(timeVarIndex[diggNumIndex]).fill(X(timeVarIndex[diggNumIndex], i - ntied + 1));
					}
					for (int a = 0; (uword) a < interMap.n_rows; a++) Xr.row(interMap(a, 0)) = Xr.row(interMap(a, 1)) % Xr.row(interMap(a, 2));
				}
				if (fixedCoefIndex >= 0) fixedCoefOffset = fixedCoef * Xr.row(fixedCoefIndex);
				if (!dropFlag.is_empty()) Xr = Xr.rows(keepIndex);
				centerX(Xr);
				mat tXr = trans(Xr);
				mat Xt = Xr.cols(0, ntied - 1);
				// compute XB, eXB, Q, S
				vec xbeta = vec(tXr * beta);
				if (!offset.is_empty()) xbeta += offset(v);
				if (fixedCoefIndex >= 0) xbeta += fixedCoefOffset;
				xbeta -= mean(xbeta);
				if (isCutBeta == 1) {
					setMax(xbeta, 22);
					setMin(xbeta, -200);
				}
				vec exp_xbeta = exp(xbeta);
				double S = sum(exp_xbeta);
				// Number of events
				double nEvent = 1;
				if(method==0 || method==1) nEvent = sum(wt);
				// Compute cumulative baseline hazard, based on event time
				k++;
				db(k,0) = time;
				db(k,1) += nEvent / S;
				db(k,2) += nEvent / (S*S);
				// Compute individual cumulative hazard, based on event time
				/* S doesn't change in between events. Two assumptions:
				 * 1) Ties have negligible effect on the value of S (often true for large datasets)
				 * 2) w_exp_xbeta doesn't change in between events (true if split based on adoptions)
				 * 3) Doesn't support general case weight so far, but support real ties with unit case weight!!!
				 */
				for (int j = 0; j < nrisk; j++) {
					int r = riskset[j];
					ind(receiver[r], s) += nEvent * exp_xbeta[j] / S;
					indvar(receiver[r], s) += nEvent * exp_xbeta[j] * exp_xbeta[j] / (S*S);
				}
			} // end of events
			base.push_back(db);
		} // end of all strata
		delete [] riskset;
		indvar = sqrt(indvar);
		base.push_back(ind);
		base.push_back(indvar);
		return base;
	}

	void getRobustV() {
		R = zeros(p, n);
		uword* riskset = new uword[n];
		int nrisk = 0;
		int ntied = 0;
		double time, time_lb, time_ub;
		int k, d, j;
		double tmp;
		int a;
		uvec tiegroup = zeros<uvec>(n);
		int tieID = 0;
		vec fixedCoefOffset;
		// Strata by strata, data ordered by -strata, event, stop, and start in descending order
		for (int s = 0; s < (int) (strataBounds.n_rows); s++) {
			// Traverse events within a strata
			for (int i = strataBounds(s, 0); i <= strataBounds(s, 1); i++) {
				if (event[i] == 0) break;
				time = stop[i];
				time_lb = time * (1 - eps);
				time_ub = time * (1 + eps);
				nrisk = 0;
				ntied = 0;
				for (j = i; j <= strataBounds(s, 1); j++) {
					//                    Rcout << j << "," << stop[j] << "," << start[j] << "," << time << "," << time_ub << "," << time_lb << "," << event[j] << "\n";
					if (stop[j] > time_lb && start[j] < time) {
						riskset[nrisk++] = j;
						if (event[j] == 1 && stop[j] < time_ub) ntied++;
					} else if (event[j] == 0 && stop[j] < time) break;
				}
				for (j = strataBounds(s, 0); j < i; j++) {
					if (stop[j] >= time && start[j] < time) riskset[nrisk++] = j;
				}
				for (k = i; k <= i + ntied - 1; k++) tiegroup[k] = tieID;
				//                Rcout << i << "," << ntied << "," << tieID << "\n";
				tieID++;
				i += ntied - 1;
				// Update X
				uvec v = uvec(riskset, nrisk, false, true);
				vec wr = weight(v);
				vec wt = wr(span(0, ntied - 1));
				mat Xr = X.cols(v);
				// Deal with time-varying variables: Digg number and time, if subtracting mean before interaction, the results may be different with SAS
				if (!timeVarIndex.is_empty()) {
					if (ageIndex >= 0) {
						if (isAbsoluteAge == 1) Xr.row(timeVarIndex[ageIndex]).fill(time / 86400.0);
						else Xr.row(timeVarIndex[ageIndex]) = (time - Xr.row(timeVarIndex[ageIndex])) / 86400.0;
					}
					if (diggNumIndex >= 0 && isUpdateDiggNum == 1) {
						Xr.row(timeVarIndex[diggNumIndex]).fill(X(timeVarIndex[diggNumIndex], i - ntied + 1));
					}
					for (a = 0; (uword) a < interMap.n_rows; a++) Xr.row(interMap(a, 0)) = Xr.row(interMap(a, 1)) % Xr.row(interMap(a, 2));
				}
				if (fixedCoefIndex >= 0) fixedCoefOffset = fixedCoef * Xr.row(fixedCoefIndex);
				if (!dropFlag.is_empty()) Xr = Xr.rows(keepIndex);
				centerX(Xr);
				mat Xt = Xr.cols(0, ntied - 1);
				// compute XB, eXB, Q, Q1,Q2, S,S1,S2,Eq,Es
				vec xbeta = vec(trans(trans(beta) * Xr));
				if (!offset.is_empty()) xbeta += offset(v);
				if (fixedCoefIndex >= 0) xbeta += fixedCoefOffset;
				xbeta -= mean(xbeta);
				if (isCutBeta == 1) {
					setMax(xbeta, 22);
					setMin(xbeta, -200);
				}
				vec exp_xbeta = exp(xbeta);
				vec w_exp_xbeta = wr % exp_xbeta;
				vec tied_exp_xbeta = w_exp_xbeta(span(0, ntied - 1));
				if (method == 2) tied_exp_xbeta = exp_xbeta(span(0, ntied - 1));
				double S = sum(w_exp_xbeta);
				double Q = sum(tied_exp_xbeta);
				vec S1 = vec(Xr * w_exp_xbeta);
				vec Q1 = vec(Xt * tied_exp_xbeta);
				vec Es = S1 / S;
				vec Eq = Q1 / Q;
				// Update R
				// 0-Breslow, 1-Efron, 2-Collective Cause, 3-Single Cause
				if (method == 0) {
					tmp = sum(wt) / S;
					for (k = 0; k < ntied; k++) R.col(v[k]) += (1 - tmp * exp_xbeta[k]) * (Xt.col(k) - Es);
					for (k = ntied; (uword) k < v.n_elem; k++) R.col(v[k]) -= tmp * exp_xbeta[k] * (Xr.col(k) - Es);
				} else if (method == 1) {
					for (k = 0; k < ntied; k++) R.col(v[k]) += Xr.col(k);
					for (d = 0; d < ntied; d++) {
						double ratio = d / (double) ntied;
						double newS = S - Q * ratio;
						vec newS1 = S1 - Q1 * ratio;
						vec newEs = newS1 / newS;
						tmp = sum(wt) * (1 - ratio) / newS;
						for (k = 0; k < ntied; k++) R.col(v[k]) -= (newEs + exp_xbeta[k] * tmp * (Xt.col(k) - newEs)) / ntied;
						tmp = mean(wt) / newS;
						for (j = ntied; (uword) j < v.n_elem; j++) R.col(v[j]) -= tmp * exp_xbeta[j] * (Xr.col(j) - newEs);
					}
				} else if (method == 2) {
					vec part1 = (1 - mean(wt) * Q / S) / ntied * (Eq - Es);
					for (k = 0; k < ntied; k++) R.col(v[k]) += part1;
					tmp = mean(wt) / S;
					for (k = ntied; (uword) k < v.n_elem; k++) R.col(v[k]) -= tmp * exp_xbeta[k] * (Xr.col(k) - Es);
				} else if (method == 3) {
					vec dt = ones(ntied) / (double) ntied;
					if (!beta_last.is_empty()) {
						vec delta_est = zeros(ntied);
						vec xbeta_est = vec(trans(Xt) * beta_last);
						if (!offset.is_empty()) xbeta_est += offset(v(span(0, ntied - 1)));
						if (fixedCoefIndex >= 0) xbeta_est += fixedCoefOffset(span(0, ntied - 1));
						uword index;
						if (deltaType == 1) {
							xbeta_est.max(index);
							delta_est[index] = 1;
						} else if (deltaType == -1) {
							xbeta_est.min(index);
							delta_est[index] = 1;
						} else if (deltaType == 0) {
							vec exp_xbeta_est = exp(xbeta_est - max(xbeta_est));
							delta_est = exp_xbeta_est / sum(exp_xbeta_est);
						} else Rcout << "@@@delta type not specified, beta_last info will not be used\n";
						dt = delta_est;
					} else if (!delta.is_empty()) dt = delta(v(span(0, ntied - 1)));
					wt %= dt;
					tmp = sum(wt) / S;
					for (k = 0; k < ntied; k++) R.col(v[k]) += (dt[k] - tmp * exp_xbeta[k]) * (Xt.col(k) - Es);
					for (k = ntied; (uword) k < v.n_elem; k++) R.col(v[k]) -= tmp * exp_xbeta[k] * (Xr.col(k) - Es);
				}
			} // end of strata
		} // end of all observations
		delete [] riskset;
		for (int i = 0; i < n; i++) R.col(i) *= weight[i];
		D = trans(R) * V_conv;
		if (isTieGroup == 1 && method == 2) {
			if (!egroups.is_empty()) Rcout << "Warning: error group has been overrided by tie groups\n";
			egroups = tiegroup;
		}
		if (!egroups.is_empty()) {
			uword ng = max(egroups) + 1;
			mat Dt = zeros(ng, p);
			for (uword i = 0; i < ng; i++) Dt.row(i) = sum(D.rows(find(egroups - i == 0)), 0);
			V_robust = trans(Dt) * Dt;
		} else V_robust = trans(D) * D;
		se_robust = sqrt(diagvec(V_robust));
	}

	inline vec getS1(mat m, vec w_exp_xbeta, umat groups) {
		vec res = zeros(p + nf);
		res(span(0, p - 1)) = vec(trans(m) * w_exp_xbeta);
		for (uword i = 0; i < groups.n_rows; i++) res(groups.row(i)) += w_exp_xbeta[i];
		return res;
	}

	inline mat getS2(mat m, vec w_exp_xbeta, umat groups, double Q) {
		mat res = zeros(p + nf, p + nf);
		res(span(0, p - 1), span(0, p - 1)) = trans(m) * diagmat(w_exp_xbeta / Q) * m;
		for (uword i = 0; i < groups.n_rows; i++) {
			for (uword j = 0; j < groups.n_cols; j++) {
				res(groups.at(i, j), span(0, p - 1)) += w_exp_xbeta[i] * m.row(i) / Q;
				//res(span(0, p - 1), groups.at(i, j)) = res(groups.at(i, j), span(0, p - 1)).t();
				for (uword k = j; k < groups.n_cols; k++) {
					res.at(groups.at(i, k), groups.at(i, j)) += w_exp_xbeta[i] / Q;
					//res.at(groups.at(i, j), groups.at(i, k)) = res.at(groups.at(i, k), groups.at(i, j));
				}
			}
		}
		return res;
	}

	inline vec getFrailtyDiag(vec theta) {
		vec eta = zeros(p + nf);
		if (isFixedEffect != 1) eta(span(p + randStarts, p + nf - 1)) = 1 / theta(span(randStarts, nf - 1));
		return eta;
	}

	inline mat getVecOuter(vec v) {
		int m = (int) v.n_elem;
		mat res = zeros(m, m);
		for (int i = 0; i < m; i++) {
			if (fabs(v[i]) > 1e-14) {
				res.col(i) = v[i] * v;
				res.row(i) = res.col(i).t();
			}
		}
		return res;
	}

	double getSparsity(mat m, std::string head) {
		uvec v = find(abs(m) < 1e-14);
		Rcout << "Matrix " << head << " sparsity level = " << (double) v.n_elem / m.n_elem << endl;
		return (double) v.n_elem / m.n_elem;
	}

	double getSparsity(vec m, std::string head) {
		uvec v = find(abs(m) < 1e-14);
		Rcout << "Vector " << head << " sparsity level = " << (double) v.n_elem / m.n_elem << endl;
		return (double) v.n_elem / m.n_elem;
	}

	void estimateFrailty() {
		converged = 0;
		Rcout << "n=" << n << ", p=" << p << ", nf=" << nf << endl;
		double Ln = -DBL_MAX;
		double L0 = 0;
		vec beta0 = beta; // "=" results in clone by default
		vec betaOld = beta;
		vec thetaOld = theta;
		vec Un;
		mat I = zeros(p + nf, (isDiagFrailty == 1) ? p : (p + nf));
		mat In;
		vec D, Dn; // diag of frailty
		vec diagH22;
		mat H22;
		vec S1 = zeros(p + nf);
		uword* riskset = new uword[n];
		int nrisk = 0;
		int ntied = 0;
		uword a, b, c;
		double time, time_lb, time_ub;
		vec fixedCoefOffset;
		// Strata by strata, data ordered by -strata, event, stop, and start in descending order
		for (int iter_out = 1; iter_out <= maxOutIter; iter_out++) { // update theta in outer loop
			if (isFixedEffect != 1) Rcout << "======== The " << iter_out << " outer loop started, theta=" << mean(theta(span(randStarts, nf - 1))) << ", L=" << L << "========\n";
			int ntrials = 1;
			vec diagA = getFrailtyDiag(theta);
			for (int iter = 1; iter <= maxIter; iter++) { // update beta in inner loop
				E = 0;
				L = -dot(square(beta), diagA) / 2;
				U = -diagA % beta;
				I.zeros();
				if (isDiagFrailty == 1) D = diagA(span(p, diagA.n_elem - 1));
				else I.diag() = diagA;

				for (int s = 0; s < (int) (strataBounds.n_rows); s++) {
					// Traverse events within a strata
					for (int i = strataBounds(s, 0); i <= strataBounds(s, 1); i++) {
						if (event[i] == 0) break;
						time = stop[i];
						if(time<minstop - eps) continue;
						if(time>maxstop + eps) break;
						time_lb = time * (1 - eps);
						time_ub = time * (1 + eps);
						nrisk = 0;
						ntied = 0;
						for (int j = i; j <= strataBounds(s, 1); j++) {
							if (stop[j] > time_lb && start[j] < time) {
								riskset[nrisk++] = j;
								if (event[j] == 1 && stop[j] < time_ub) ntied++;
							} else if (event[j] == 0 && stop[j] < time) break;
						}
						for (int j = strataBounds(s, 0); j < i; j++) {
							if (stop[j] >= time && start[j] < time) riskset[nrisk++] = j;
						}
						if (ntied <= 0) {
							Rcout << "Error!!! ntied=" << ntied << ", start=" << start[i] << ", stop=" << stop[i] << "\n";
							throw std::invalid_argument( "Logical error: ntied<=0!" );
						}
						i += ntied - 1;
						// Update X
						uvec v = uvec(riskset, nrisk, false, true);
						vec wr = weight(v);
						vec wt = wr(span(0, ntied - 1));
						umat groups = fgroups.rows(v);
						umat tied_groups = groups.rows(0, ntied - 1);
						mat Xr = X.cols(v);
						// Deal with time-varying variables: Digg number and time, if subtracting mean before interaction, the results may be different with SAS
						if (!timeVarIndex.is_empty()) {
							if (ageIndex >= 0) {
								if (isAbsoluteAge == 1) Xr.row(timeVarIndex[ageIndex]).fill(time / 86400.0);
								else Xr.row(timeVarIndex[ageIndex]) = (time - Xr.row(timeVarIndex[ageIndex])) / 86400.0;
							}
							if (diggNumIndex >= 0 && isUpdateDiggNum == 1) {
								Xr.row(timeVarIndex[diggNumIndex]).fill(X(timeVarIndex[diggNumIndex], i - ntied + 1));
							}
							for (a = 0; a < interMap.n_rows; a++) Xr.row(interMap(a, 0)) = Xr.row(interMap(a, 1)) % Xr.row(interMap(a, 2));
						}
						if (fixedCoefIndex >= 0) fixedCoefOffset = fixedCoef * Xr.row(fixedCoefIndex);
						if (!dropFlag.is_empty()) Xr = Xr.rows(keepIndex);
						centerX(Xr);
						mat tXr = trans(Xr);
						mat Xt = Xr.cols(0, ntied - 1);
						mat tXt = trans(Xt);
						// compute XB, eXB, Q, Q1,Q2, S,S1,S2,Eq,Es
						vec xbeta = vec(tXr * beta(span(0, p - 1)));
						for (a = 0; a < groups.n_cols; a++) xbeta += beta(groups.col(a));
						if (!offset.is_empty()) xbeta += offset(v);
						if (fixedCoefIndex >= 0) xbeta += fixedCoefOffset;
						xbeta -= mean(xbeta);
						if (isCutBeta == 1) {
							setMax(xbeta, 22);
							setMin(xbeta, -200);
						}
						vec exp_xbeta = exp(xbeta);
						vec w_exp_xbeta = wr % exp_xbeta;
						vec tied_exp_xbeta = w_exp_xbeta(span(0, ntied - 1));
						if (method == 2) tied_exp_xbeta = exp_xbeta(span(0, ntied - 1));
						double S = sum(w_exp_xbeta);
						double Q = sum(tied_exp_xbeta);
						//S1 = getS1(tXr, w_exp_xbeta, groups);
						//S2 = getS2(tXr, w_exp_xbeta, groups);
						// The code below were used to compute S1, S2 for better speed
						S1.zeros();
						vec Q1 = getS1(tXt, tied_exp_xbeta, tied_groups);

						// I += S2
						mat w_tXr = tXr;
						w_tXr.each_col() %= w_exp_xbeta;
						mat w_tXr_dS = w_tXr / S;
						vec w_exp_xbeta_dS = w_exp_xbeta / S;
						for (a = 0; (int) a < p; a++) {
							S1[a] = sum(w_tXr.col(a));
							for (b = 0; b <= a; b++) {
								I.at(a, b) += dot(w_tXr_dS.col(a), tXr.col(b));
							}
						}
						// below are identifical if only one frailty term
						if (isDiagFrailty == 1) {
							for (a = 0; a < groups.n_rows; a++) {
								S1(groups.row(a)) += w_exp_xbeta[a];
								for (b = 0; b < groups.n_cols; b++) {
									I(groups.at(a, b), span(0, p - 1)) += w_tXr_dS.row(a);
									D[groups.at(a, b) - p] += w_exp_xbeta_dS[a]; // ignore off diagonal elements
								}
							}
						} else {
							for (a = 0; a < groups.n_rows; a++) {
								S1(groups.row(a)) += w_exp_xbeta[a];
								for (b = 0; b < groups.n_cols; b++) {
									I(groups.at(a, b), span(0, p - 1)) += w_tXr_dS.row(a);
									for (c = b; c < groups.n_cols; c++) {
										I.at(groups.at(a, c), groups.at(a, b)) += w_exp_xbeta_dS[a];
									}
								}
							}
						}

						// I += Q2
						//mat Q2 = getS2(tXt, tied_exp_xbeta, tied_groups, Q);
						mat w_tXt = tXt;
						w_tXt.each_col() %= tied_exp_xbeta;
						mat w_tXt_dS = w_tXt / Q;
						vec tied_exp_xbeta_dS = tied_exp_xbeta / Q;
						for (a = 0; (int) a < p; a++) {
							for (b = 0; b <= a; b++) {
								I.at(a, b) -= dot(w_tXt_dS.col(a), tXt.col(b));
							}
						}
						if (isDiagFrailty == 1) {
							for (a = 0; a < tied_groups.n_rows; a++) {
								for (b = 0; b < tied_groups.n_cols; b++) {
									I(tied_groups.at(a, b), span(0, p - 1)) -= w_tXt_dS.row(a);
									D[tied_groups.at(a, b) - p] -= tied_exp_xbeta_dS[a];
								}
							}
						} else {
							for (a = 0; a < tied_groups.n_rows; a++) {
								for (b = 0; b < tied_groups.n_cols; b++) {
									I(tied_groups.at(a, b), span(0, p - 1)) -= w_tXt_dS.row(a);
									for (c = b; c < tied_groups.n_cols; c++) {
										I.at(tied_groups.at(a, c), tied_groups.at(a, b)) -= tied_exp_xbeta_dS[a];
									}
								}
							}
						}


						vec P = getS1(tXt, wt, tied_groups); // vec(Xt * wt) with frailty;
						vec Es = S1 / S;
						vec Eq = Q1 / Q;
						// Update L, U, I, V
						// 0-Breslow, 1-Efron, 2-Collective Cause, 3-Single Cause
						if (method == 2) { // assume weights for tied obs are the same, use mean in case not equal by mistake
							L += mean(wt) * (log(Q) - log(S));
						U += mean(wt)*(Eq - Es);
						if (isDiagFrailty == 1) {
							for (a = 0; a < (uword) p; a++)
								I(span(a, I.n_rows - 1), a) += Eq(span(a, I.n_rows - 1)) * Eq[a] - Es(span(a, I.n_rows - 1)) * Es[a];
							D += square(Eq(span(p, Eq.n_elem - 1))) - square(Es(span(p, Es.n_elem - 1)));
						} else {//the code below results in off diagonal elements in I for frailty part
							I += getVecOuter(Eq);
							I -= Es * trans(Es);
						}
						if (fabs(mean(wt) - 1) > 1e-14) I *= mean(wt);
						} else {
							Rcout << "This tuned version only support method 2 (collective cause model) for now!" << endl;
							throw std::invalid_argument( "This tuned version only support method 2 (collective cause model) for now!" );
						}
					} // end of strata
				} // end of all observations
				vec dbeta;
				if (isDiagFrailty == 1) {
					I(span(0, p - 1), span(0, p - 1)) = symmatl(I(span(0, p - 1), span(0, p - 1)));
					dbeta = getDbeta(U, I, D);
				} else {
					I = symmatl(I);
					dbeta = getDbeta(U, I);
				}
				if (fabs(L - Ln) < fabs(L) * tol && isFinite(L)) {
					if (verbose >= 1)
						Rcout << "- Converged after " << iter << " iterations, relative change in L =" << fabs((L - Ln) / L) << "\n";
					if (L < Ln) { // go back if L<Ln
						L = Ln;
						beta = betaOld;
						U = Un;
						I = In;
						D = Dn;
					} else Ln = L; // record for later use in theta updating
					beta(span(0, p - 1)).t().print(Rcout);
					beta(span(p, p + fmin(9, nf - 1))).t().print(Rcout);
					break;
				}
				// If decreases by a lot, over reach problem
				if (L < Ln || !isFinite(L) || !dbeta.is_finite()) {
					if (verbose >= 2) Rcout << "--Over reach in the " << iter << " iterations after " << ntrials << " trials, L="
							<< L << ", dbeta finite: " << dbeta.is_finite() << "\n";
					if (all(abs(beta - betaOld) < abs(beta) * tol_par) || ntrials > maxStepHalving) {
						Rcout << "---Recovered to last best because no better results can be found in the suggested direction\n";
						L = Ln;
						beta = betaOld;
						U = Un;
						I = In;
						D = Dn;
						break;
					}
					if(nf>0 && eigPct<0.02 && ntrials>3){
						eigPct *= 3;
						Rcout << "Recalculating dbeta with eigPct increased to " << eigPct << endl;
						beta = betaOld + getDbeta(U, I, isDiagFrailty == 1 ? D : vec());
					} else  beta = (beta + betaOld) / 2;
					iter--; // not converging
					ntrials++;
				} else {
					ntrials = 1;
					if (isCutStep == 1) {
						setMax(dbeta, 1000);
						setMin(dbeta, -1000);
					}
					// Score Test and Likelihood ratio test, be aware of the order of the adjacent lines
					if (iter == 1) {
						scoreTest = dot(U, dbeta); // =trans(U) * dbeta;
						L0 = L;
					}
					betaOld = beta;
					beta += dbeta;
					if (verbose >= 2) {
						Rcout << "-- The " << iter << "/" << maxIter << " iteration finished, L=" << L << ", Ln=" << Ln << ", next beta =\n";
						beta(span(0, p - 1)).t().print(Rcout);
						beta(span(p, p + fmin(9, nf - 1))).t().print(Rcout);
					}
					Ln = L;
					Un = U;
					In = I;
					Dn = D;
					if(nf>0 && mean( abs( dbeta(span(0,p-1)) ) ) < 0.02 && eigPct > 1e-3 + 1e-6) {
						eigPct /= sqrt(3.0);
						Rcout << "Last successful step is too small. Reducing eigPct to " << eigPct << endl;
					}
				}
				if (iter == maxIter) Rcout << "Warning: the model is not converged after the maximal " << maxIter << " inner iterations\n";
			} // end of inner loop
			if (isFixedEffect == 1) break;
			int x1 = p + randStarts;
			int x2 = p + nf - 1;
			// Update theta
			thetaOld = theta;
			// Note that L may decrease in the following iteration
			for (int iter = 1; iter <= thetaIterMax; iter++) {
				// always assume isREML = 0 if isDiagFrailty = 1
				if (isDiagFrailty == 1) {
					diagH22 = 1 / D(span(randStarts, nf - 1));
					// D(span(randStarts, randStarts + 19)).t().print(Rcout, "D=");
				}
				else {
					if (isREML == 0) {
						H22 = getInverse(I(span(x1, x2), span(x1, x2)));
						// I(span(x1, x1+19), span(x1, x1+19)).diag().t().print(Rcout, "diagI=");
						diagH22 = diagvec(H22);
					} else {
						H22 = getInverse(I, true)(span(x1, x2), span(x1, x2));
						diagH22 = diagvec(H22);
					}
				}
				vec thetaLast = theta;
				theta(span(randStarts, nf - 1)) = square(beta(span(x1, x2))) + diagH22;
				for (int g = 0; g <= nThetaGroups; g++) {
					uvec ix = find(thetagroups - g == 0); // thetagroups==g can change the value of thetagroups
					double mean_value = mean(theta(ix));
					theta(ix).fill(mean_value);
				}
				// Update L,U,I, may need halving strategy for theta as well
				vec eta = zeros(p + nf);
				eta(span(x1, x2)) = 1 / thetaLast(span(randStarts, nf - 1)) - 1 / theta(span(randStarts, nf - 1));
				//double dL = as_scalar(trans(beta) * diffA * beta) / 2;
				double dL = dot(square(beta), eta) / 2;
				L += dL;
				U += eta % beta;
				if (isDiagFrailty == 1) D -= eta(span(p, eta.n_elem - 1));
				else I.diag() -= eta;
				if (all(abs(theta - thetaLast) < abs(theta) * tol_par)) break;
			}
			if (all(abs(theta - thetaOld) < abs(theta) * tol_par)) {
				if (verbose >= 0)
					Rcout << "- Outer loop converged after " << iter_out << " iterations, theta=" << mean(theta) << ", L=" << L << "\n";
				converged = 1;
				break;
			}
			betaOld = beta;
			Ln = L;
			Un = U;
			In = I;
			Dn = D;
			if (isDiagFrailty == 1) beta += getDbeta(U, I, D);
			else beta += getDbeta(U, I);
			if (iter_out == maxOutIter) Rcout << "Warning: the model is not converged after the maximal " << maxOutIter << " outer iterations\n";
		}// end of outer loop
		delete [] riskset;
		//Likelihood ratio test
		LRTest = 2 * (L - L0);
		/* The codes below needs to be updated */
		if (isDiagFrailty == 1) {
			vec dbeta = beta - beta0;
			waldTest = as_scalar(trans(dbeta(span(0, p - 1))) * I(span(0, p - 1), span(0, p - 1)) * dbeta(span(0, p - 1)));
			waldTest += dot(square(dbeta(span(p, p + nf - 1))), D);
			V_conv = getInverse(I(span(0, p - 1), span(0, p - 1)));
			se = zeros(p + nf);
			se(span(0, p - 1)) = sqrt(diagvec(V_conv));
			se(span(p, p + nf - 1)) = sqrt(1.0 / D);
			vec et = theta(span(randStarts, nf - 1));
			double meanThetha = mean(et);
			theta_se = 2 * meanThetha / (et.n_elem + sum(square(diagH22)) / (meanThetha * meanThetha) - 2 * sum(diagH22) / meanThetha );
		} else {
			waldTest = as_scalar(trans(beta - beta0) * I * (beta - beta0));
			V_conv = getInverse(I, true, true);
			se = zeros(p + nf);
			se(span(0, V_conv.n_rows - 1)) = sqrt(diagvec(V_conv));
			// variance of theta, may not be correct if theta are allowed to be different for individuals
			vec et = theta(span(randStarts, nf - 1));
			double meanThetha = mean(et);
			theta_se = 2 * meanThetha / (et.n_elem + trace(H22 * H22) / (meanThetha * meanThetha) - 2 * trace(H22) / meanThetha);
		}
		beta_last = beta; // clone beta for potential two-stage or multi-stage
		theta_se = sqrt(theta_se);
		Rcout << "Estimation finished!" << endl;
	}

	static inline sp_mat genSPI(mat I, vec D) {
		uword d1 = I.n_rows; //p + nf - degree of supressed freedom (fixed effects)
		uword d2 = I.n_cols; //p
		if (d1 <= d2 || d1 - d2 != D.n_elem) {
			Rcout << "Dimensions doesn't match, d1=" << d1 << ", d2=" << d2 << ", D.n_elem=" << D.n_elem << endl;
			throw std::invalid_argument( "Incorrect parameters for genSPI!" );
		}
		uword len = 2 * I.n_elem - d2 * d2 + D.n_elem;
		umat loc = zeros<umat>(2, len);
		vec values = zeros(len);
		for (uword i = 0; i < d2; i++) {
			loc(0, span(i * d1, i * d1 + d1 - 1)) = linspace<urowvec>(0, d1 - 1, d1);
			loc(1, span(i * d1, i * d1 + d1 - 1)).fill(i);
			values(span(i * d1, i * d1 + d1 - 1)) = I.col(i);
		}
		for (uword i = 0; i < d1 - d2; i++) {
			loc(0, span(I.n_elem + i * (d2 + 1), I.n_elem + i * (d2 + 1) + d2 - 1)) = linspace<urowvec>(0, d2 - 1, d2);
			loc(1, span(I.n_elem + i * (d2 + 1), I.n_elem + i * (d2 + 1) + d2 - 1)).fill(i + d2);
			values(span(I.n_elem + i * (d2 + 1), I.n_elem + i * (d2 + 1) + d2 - 1)) = I.row(i + d2).t();
			loc(0, I.n_elem + i * (d2 + 1) + d2) = i + d2;
			loc(1, I.n_elem + i * (d2 + 1) + d2) = i + d2;
			values[I.n_elem + i * (d2 + 1) + d2] = D[i];
		}
		return sp_mat(loc, values, d1, d1, false, false);
	}

	/* It is yet not clear why sometimes some frailty's beta stays zero forever */
	vec getDbeta(vec U, mat I, vec D = vec()) {
		uword dim = U.n_elem;
		if (!validFrailty.is_empty()) {
			U = U(validFrailty);
			if (!D.is_empty()) {
				uvec tmp = validFrailty(find(validFrailty >= (uword) p)) - (uword) p;
				D = D(tmp);
				I = I.rows(validFrailty);
			} else I = I(validFrailty, validFrailty);
		}
		vec dbeta = zeros(U.n_elem);
		try {
			// check singularity first, otherwise often over reach
			double minDiag = min(I.diag());
			double maxDiag = max(I.diag());
			if (!D.is_empty()) {
				minDiag = fmin(minDiag, min(D));
				maxDiag = fmax(maxDiag, max(D));
			}else if (nf == 0) { // use eig decomposition for normal estimate
				vec eigval = eig_sym(I);
				minDiag = min(eigval);
				maxDiag = max(eigval);
			}
			if (minDiag / maxDiag < eigPct) {
				if (verbose >= 2)
					Rcout << "@@@ Negative entries in information matrix diagonal! minDiag="
					<< minDiag << ", max=" << maxDiag << ", min/max=" << minDiag / maxDiag << endl;
				I.diag() += fmax(maxDiag * eigPct, -eigScalor * minDiag) * ones(I.n_cols);
				if (!D.is_empty()) D += fmax(maxDiag * eigPct, -eigScalor * minDiag) * ones(D.n_elem);
			}
			superlu_opts settings;
			settings.symmetric = true;
			if (!D.is_empty()) {
				dbeta = spsolve(genSPI(I, D), U, "superlu", settings);
			} else if (nf > 0) {
				dbeta = spsolve(sp_mat(I), U, "superlu", settings);
			} else dbeta = solve(I, U);
		} catch (...) {
			Rcout << "@@@ Caution: unalbe to fixed the exception by adding scalor to diagnoal of I, stop updating here!!!" << endl;
		}
		if(nf>0 && any(abs(dbeta(span(0, p-1)))>1) && eigPct<0.02) {
			eigPct *= 3;
			Rcout << "Increasing eigPct to " << eigPct << endl;
			return getDbeta(U, I, D);
		}
		if (!validFrailty.is_empty() && dim > dbeta.n_elem) {
			vec newdbeta = zeros(dim);
			newdbeta(validFrailty) = dbeta;
			return newdbeta;
		} else return dbeta;
	}

	mat getInverse(mat m, bool isCheckFixedEffect = false, bool isFirstP = false) {
		mat inv = zeros(m.n_rows, m.n_cols);
		mat tmp;
		if (isCheckFixedEffect && !validFrailty.is_empty()) m = m(validFrailty, validFrailty);
		try {
			tmp = m.i();
		} catch (...) {
			if (isFirstP) { // focus on the first p when I, for example, is too large to be non-singular
				tmp = zeros(m.n_rows, m.n_cols);
				try {
					Rcout << ">>>Trying to ignore frailty terms while estimating standard errors\n";
					tmp(span(0, p - 1), span(0, p - 1)) = m(span(0, p - 1), span(0, p - 1)).i();
				} catch (...) {
					Rcout << "@@@ Caution: unable to make I non-singular by focusing on the first p elements!!!" << endl;
				}
			} else {
				Rcout << "Exception while calculating inverse!" << endl;
				diagvec(m(span(0, fmin(9, m.n_rows)), span(0, fmin(9, m.n_rows)))).t().print(Rcout,"diag of matrix=");
			}
		}
		if (isCheckFixedEffect && !validFrailty.is_empty()) inv(validFrailty, validFrailty) = tmp;
		else inv = tmp;
		return inv;
	}

	void printDeltaType() {
		Rcout << "The current delta type is " << deltaType << "\n";
	}

	List estimate() {
		try {
			if (nf > 0) estimateFrailty();
			else estimateNormal();
		} catch (...) {
			Rcout << "@@@ Caution: Unknown exception captured!!! Returning current state of parameter estimates..." << endl;
		}
		compileResult();
		return res;
	}

	~CoxReg() {
		;
	}
};

RCPP_MODULE(cox_module) {
	class_<CoxReg> ("CoxReg")
            		.constructor<NumericMatrix, NumericMatrix, List >()
            		.field("deltaType", &CoxReg::deltaType)
            		.method("printDeltaType", &CoxReg::printDeltaType)
            		.method("survivalNormal", &CoxReg::survivalNormal)
            		.method("estimate", &CoxReg::estimate);
}
