// [[Rcpp::depends(RcppArmadillo)]] 
#include <RcppArmadillo.h>
#include <string>
#include <cmath>
using namespace arma;
using namespace std;
using namespace Rcpp;
#include "DJPTtools.h"
#include "optim.h"
#include "stats.h"
#include "SSpace.h"
#include "ARMAmodel.h"
#include "BSMmodel.h"

// [[Rcpp::export]]
SEXP UCompC(SEXP commands, SEXP ys, SEXP us, SEXP models, SEXP periodss, SEXP rhoss,
            SEXP hs, SEXP tTests, SEXP criterions, SEXP ps, SEXP rubbish2s, SEXP rubbishs,
            SEXP verboses, SEXP stepwises, SEXP estimOks,
            SEXP p0s, SEXP vs, SEXP yFitVs, SEXP nonStationaryTermss,
            SEXP rubbish3s, SEXP harmonicss, SEXP criterias, SEXP cycleLimitss, 
            SEXP betass, SEXP typeOutlierss){
    
    // setbuf(stdout, NULL);
    // Converting R inputs to C++
    string command = CHAR(STRING_ELT(commands, 0));
    NumericVector yr(ys);
    NumericMatrix ur(us);
    string model = CHAR(STRING_ELT(models, 0));
    NumericVector periodsr(periodss);
    NumericVector rhosr(rhoss);
    int h = as<int>(hs);
    bool tTest = as<bool>(tTests);
    string criterion = CHAR(STRING_ELT(criterions, 0));
    NumericVector pr(ps);
    NumericVector rubbishr(rubbishs);
    bool verbose = as<bool>(verboses);
    bool stepwise = as<bool>(stepwises);
    string estimOk = CHAR(STRING_ELT(estimOks, 0));
    NumericVector p0r(p0s);
    NumericVector vr(vs);
    NumericVector yFitVr(yFitVs);
    int nonStationaryTerms = as<int>(nonStationaryTermss);
    NumericVector harmonicsr(harmonicss);
    NumericVector criteriar(criterias);
    NumericMatrix rubbish2r(rubbish2s);
    NumericMatrix rubbish3r(rubbish3s);
    NumericMatrix betar(betass);
    NumericMatrix typeOutliersr(typeOutlierss);

    vec y(yr.begin(), yr.size(), false);
    mat u(ur.begin(), ur.nrow(), ur.ncol(), false);
    vec periods(periodsr.begin(), periodsr.size(), false);
    vec rhos(rhosr.begin(), rhosr.size(), false);
    vec p(pr.begin(), pr.size(), false);
    vec p0(p0r.begin(), p0r.size(), false);
    vec v(vr.begin(), vr.size(), false);
    vec yFitV(yFitVr.begin(), yFitVr.size(), false);
    vec harmonics(harmonicsr.begin(), harmonicsr.size(), false);
    vec criteria(criteriar.begin(), criteriar.size(), false);
    vec rubbish(rubbishr.begin(), rubbishr.size(), false);
    mat rubbish2(rubbish2r.begin(), rubbish2r.nrow(), rubbish2r.ncol(), false);
    mat rubbish3(rubbish3r.begin(), rubbish3r.nrow(), rubbish3r.ncol(), false);
    mat betas(betar.begin(), betar.nrow(), betar.ncol(), false);
    mat typeOutliers(typeOutliersr.begin(), typeOutliersr.nrow(), typeOutliersr.ncol(), false);
    // Correcting dimensions of u (k x n)
    size_t k = u.n_rows;
    size_t n = u.n_cols;
    if (k > n){
        u = u.t();
    }
    if (k == 1 && n == 2){
        u.resize(0);
    }
    if (typeOutliers(0, 0) == -1){
        typeOutliers.reset();
    }
    // Setting inputs
    SSinputs inputsSS;
    BSMinputs inputsBSM;
    inputsSS.y = y;
    inputsSS.u = u;
    inputsBSM.model = model;
    inputsBSM.periods = periods;
    inputsBSM.rhos = rhos;
    inputsSS.h = h;
    inputsBSM.tTest = tTest;
    inputsBSM.criterion = criterion;
    inputsSS.grad = rubbish2.col(0);
    inputsSS.p = p;
    inputsSS.p0 = p0;
    inputsSS.v = v;
    inputsSS.F = yFitV;
    inputsSS.d_t = rubbish(0);
    inputsSS.innVariance = rubbish(1);
    inputsSS.objFunValue = rubbish(2);
    inputsSS.cLlik = rubbish(3);
    inputsSS.outlier = rubbish(4);
    vec aux(1); aux(0) = inputsSS.outlier;
    if (aux.has_nan()){
        inputsSS.outlier = 0;
    }
    inputsSS.verbose = verbose;
    inputsSS.estimOk = estimOk;
    inputsSS.nonStationaryTerms = nonStationaryTerms;
    inputsSS.criteria = criteria;
    inputsSS.betaAug = betas.col(0);
    inputsSS.betaAugVar = betas.col(1);

    inputsBSM.stepwise = stepwise;
    inputsBSM.ns = rubbish3.col(0);
    inputsBSM.nPar = rubbish3.col(1);
    if (harmonics.has_nan()){
        inputsBSM.harmonics.resize(1);
        inputsBSM.harmonics(0) = 0;
    } else {
        inputsBSM.harmonics = conv_to<uvec>::from(harmonics);
    }
    inputsBSM.constPar = rubbish2.col(1);
    inputsBSM.typePar = rubbish2.col(2);
    inputsBSM.typeOutliers = typeOutliers;
    inputsBSM.arma = rubbish(5);
    // Building model
    BSMmodel sysBSM = BSMmodel(inputsSS, inputsBSM);
    // Commands
    SSinputs inputs;
    BSMinputs inputs2;
    if (command == "estimate"){
        // Estimating and Forecasting
        sysBSM.estim();
        sysBSM.forecast();
        // Values to return
        inputs = sysBSM.SSmodel::getInputs();
        inputs2 = sysBSM.getInputs();
        vec harmonicsVec = conv_to<vec>::from(inputs2.harmonics);
        vec rubbish(4);
        mat rubbish2(inputs.p.n_elem, 3),
            rubbish3(5, 2),
            betas(inputs.betaAug.n_rows, 2);
        rubbish(0) = inputs.d_t;
        rubbish(1) = inputs.innVariance;
        rubbish(2) = inputs.objFunValue;
        rubbish2.col(0) = inputs.grad;
        rubbish2.col(1) = inputs2.constPar;
        rubbish2.col(2) = inputs2.typePar;
        // if (inputs2.ns.n_elem == 5){
            rubbish3.col(0) = inputs2.ns;
            rubbish3.col(1) = inputs2.nPar;
        // }
        inputsBSM.harmonics = conv_to<uvec>::from(harmonics);
        betas.col(0) = inputs.betaAug;
        betas.col(1) = inputs.betaAugVar;
        
        // inputs.yFor.print("yFor");
        // Converting back to R
        return List::create(Named("p") = inputs.p,
                            Named("p0") = inputs.p0,
                            Named("model") = inputs2.model,
                            Named("yFor") = inputs.yFor,
                            Named("periods") = inputs2.periods,
                            Named("rhos") = inputs2.rhos,
                            Named("yForV") = inputs.FFor,
                            Named("estimOk") = inputs.estimOk,
                            Named("rubbish") = rubbish,
                            Named("harmonics") = harmonicsVec,
                            Named("rubbish2") = rubbish2,
                            Named("rubbish3") = rubbish3,
                            Named("cycleLimits") = inputs2.cycleLimits,
                            Named("nonStationaryTerms") = inputs.nonStationaryTerms,
                            Named("betas") = betas,
                            Named("u") = inputs.u,
                            Named("typeOutliers") = inputs2.typeOutliers,
                            Named("criteria") = inputs.criteria);
    } else if (command == "validate"){
        sysBSM.validate();
        // Values to return
        inputs = sysBSM.SSmodel::getInputs();
        // Converting back to R
        return List::create(Named("table") = inputs.table,
                            Named("v") = inputs.v);
    } else if (command == "filter" || command == "smooth" || command == "disturb"){
        sysBSM.setSystemMatrices();
        if (command == "filter"){
            sysBSM.filter();
        } else if (command == "smooth") {
            sysBSM.smooth(false);
        } else {
            sysBSM.disturb();
        }
        // Values to return
        inputs = sysBSM.SSmodel::getInputs();
        inputs2 = sysBSM.getInputs();
        return List::create(Named("a") = inputs.a,
                            Named("P") = inputs.P,
                            Named("v") = inputs.v,
                            Named("yFitV") = inputs.F,
                            Named("yFit") = inputs.yFit,
                            Named("eps") = inputs2.eps,
                            Named("eta") = inputs.eta);
    } else if (command == "components"){
        sysBSM.setSystemMatrices();
        sysBSM.components();
        // Values to return
        inputs2 = sysBSM.getInputs();
        // Converting back to R
        return List::create(Named("comp") = inputs2.comp,
                            Named("compV") = inputs2.compV,
                            Named("m") = inputs2.comp.n_rows);
    }
    return List::create(Named("void") = datum::nan);
}

