#include "sc.h"
#include "core.h"
#include <set>

using namespace Rcpp;

std::vector<int> matrixToVector(IntegerMatrix& yy){
    int n = yy.nrow();
    int m = yy.ncol();
    
    // create y as single vector
    std::vector<int> y(n,0);
    int maxDom = 1;
    for(int j = 0; j < m; j++){
        std::set<int> setj;
        for(int i = 0; i < n; i++){
            int curr = yy(i,j);
            setj.insert(curr);
            if(curr > maxDom){
                maxDom = curr;
            }
        }
    }
    maxDom++;
    int currentF = 1;
    for(int j = 0; j < m; j++){
        for(int i = 0; i < n; i++){
            int curr = (yy(i,j) + 1) * currentF;
            y[i] += curr;
        }
        currentF *= (maxDom + 1);
    }
    std::vector<int> yN = getNiceCategories(y);
    return yN;
}

std::vector<int> joinVectors(std::vector<int> yy, std::vector<int>& zz){
    int n = yy.size();
    
    // create y as single vector
    int maxDom = 1;
    std::set<int> setj;
    for(int i = 0; i < n; i++){
        int curr = yy[i];
        setj.insert(curr);
        if(curr > maxDom){
            maxDom = curr;
        }
    }
    for(int i = 0; i < n; i++){
        int curr = zz[i];
        setj.insert(curr);
        if(curr > maxDom){
            maxDom = curr;
        }
    }
    maxDom++;
    int currentF = maxDom + 1;
    for(int i = 0; i < n; i++){
        int curr = (zz[i] + 1) * currentF;
        yy[i] += curr;
    }
    std::vector<int> yN = getNiceCategories(yy);
    return yN;
}

// Define function as extern with RcppExport
RcppExport SEXP conditionalFNML(SEXP xEXP, SEXP yEXP) {
    return(wrap(conditionalNML(xEXP,yEXP,true)));
}

// Define function as extern with RcppExport
RcppExport SEXP conditionalQNML(SEXP xEXP, SEXP yEXP) {
    return(wrap(conditionalNML(xEXP,yEXP,false)));
}

// Define function as extern with RcppExport
double conditionalNML(SEXP& xEXP, SEXP& yEXP, bool useFNML) {
    // convert input
    IntegerMatrix xx(xEXP);
    IntegerMatrix yy(yEXP);
    int n = yy.nrow();
    
    std::vector<int> xN = matrixToVector(xx);
    std::vector<int> yN = matrixToVector(yy);
    int dX = xN[xN.size()-1];
    int dY = yN[yN.size()-1];
    xN.pop_back();
    yN.pop_back();
    
    double modelCosts = 0.0;
    if(!useFNML){
        modelCosts += regret(n,dX * dY) - regret(n,dY);
    }
        
    // compute conditional sc
    double result = 0.0;
    if(useFNML){
        result = conditionalSC(xN,yN);
    }else{
        result = (double)n * conditionalEntropy(xN,yN);
    }
    result += modelCosts;
    
    // return result
    return result;
}

RcppExport SEXP indepfNML(SEXP xEXP, SEXP yEXP, SEXP xyEXP, SEXP zEXP){
    return(wrap(indepNML(xEXP,yEXP,xyEXP,zEXP,true)));
}

RcppExport SEXP indepqNML(SEXP xEXP, SEXP yEXP, SEXP xyEXP, SEXP zEXP){
    return(wrap(indepNML(xEXP,yEXP,xyEXP,zEXP,false)));
}

RcppExport SEXP indepAsymfNML(SEXP xEXP, SEXP yEXP, SEXP zEXP){
    return(wrap(indepAsymNML(xEXP,yEXP,zEXP,true)));
}

RcppExport SEXP indepAsymqNML(SEXP xEXP, SEXP yEXP, SEXP zEXP){
    return(wrap(indepAsymNML(xEXP,yEXP,zEXP,false)));
}

// Define function as extern with RcppExport
double indepNML(SEXP& xEXP, SEXP& yEXP, SEXP& xyEXP, SEXP& zEXP, bool useFNML) {
    // convert input
    IntegerMatrix xx(xEXP);
    IntegerMatrix yy(yEXP);
    IntegerMatrix xxyy(xyEXP);
    IntegerMatrix zz(zEXP);
    int n = yy.nrow();
    
    std::vector<int> xN = matrixToVector(xx);
    std::vector<int> yN = matrixToVector(yy);
    std::vector<int> xyN = matrixToVector(xxyy);
    std::vector<int> zN = matrixToVector(zz);
    int dX = xN[xN.size()-1];
    int dY = yN[yN.size()-1];
    int dZ = zN[zN.size()-1];
    xN.pop_back();
    yN.pop_back();
    xyN.pop_back();
    zN.pop_back();
    
    // additional model costs
    double result = 0.0;
    double modelCosts = 0.0;
    // calculate costs
    if(useFNML){
        result = conditionalSC(xN,zN) + conditionalSC(yN,zN) - conditionalSC(xyN,zN);
    }else{
        modelCosts += regret(n, dX * dZ) + regret(n, dY * dZ) - regret(n, dZ) - regret(n, dX * dZ * dY);
        // costs
        result = (double)n * (conditionalEntropy(xN,zN) + conditionalEntropy(yN,zN) - conditionalEntropy(xyN,zN));
    }
    result += modelCosts;
    
    // return result
    return result;
}

// Define function as extern with RcppExport
double indepAsymNML(SEXP& xEXP, SEXP& yEXP, SEXP& zEXP, bool useFNML) {
    // convert input
    IntegerMatrix xx(xEXP);
    IntegerMatrix yy(yEXP);
    IntegerMatrix zz(zEXP);
    int n = yy.nrow();
    
    std::vector<int> xN = matrixToVector(xx);
    std::vector<int> yN = matrixToVector(yy);
    std::vector<int> zN = matrixToVector(zz);
    int dX = xN[xN.size()-1];
    int dY = yN[xN.size()-1];
    int dZ = zN[zN.size()-1];
    xN.pop_back();
    yN.pop_back();
    zN.pop_back();
    
    // join z*, y
    std::vector<int> yzN = joinVectors(yN, zN);
    yzN.pop_back();
    
    // additional model costs
    double result = 0.0;
    double modelCosts = 0.0;
    // calculate costs
    if(useFNML){
        result = conditionalSC(xN,zN) - conditionalSC(xN,yzN);
    }else{
        modelCosts += regret(n, dX * dZ) - regret(n, dZ) - regret(n, dX * dZ * dY) + regret(n, dZ * dY);
        // costs
        result = (double)n * (conditionalEntropy(xN,zN) - conditionalEntropy(xN,yzN));
    }
    result += modelCosts;
    
    // return result
    return result;
}

// Define function as extern with RcppExport
RcppExport SEXP shannonEntropy(SEXP xEXP) {
    // convert input
    IntegerVector x(xEXP);
    
    std::map<int,int> counts;
    int sum = x.size();
    for(int i = 0; i < sum; i++){
        counts[x[i]]++;
    }
    double score = entropy(counts, sum);
    
    // return result
    return(wrap(score));
}

// Define function as extern with RcppExport
RcppExport SEXP conditionalShannonEntropy(SEXP xEXP, SEXP yEXP) {
    // convert input
    IntegerVector x(xEXP);
    IntegerMatrix yy(yEXP);
    int n = yy.nrow();
    int m = yy.ncol();
    
    // create y as single vector
    std::vector<int> y(n,0);
    int maxDom = 1;
    int sizeCart = 1;
    for(int j = 0; j < m; j++){
        std::set<int> setj;
        for(int i = 0; i < n; i++){
            int curr = yy(i,j);
            setj.insert(curr);
            if(curr > maxDom){
                maxDom = curr;
            }
        }
        sizeCart *= setj.size();
    }
    maxDom++;
    int currentF = 1;
    for(int j = 0; j < m; j++){
        for(int i = 0; i < n; i++){
            int curr = (yy(i,j) + 1) * currentF;
            y[i] += curr;
        }
        currentF *= (maxDom + 1);
    }
    std::vector<int> newx = as<std::vector<int> >(x);
    double score = conditionalEntropy(newx,y);
    
    // return result
    return(wrap(score));
}

// Define function as extern with RcppExport
RcppExport SEXP regret(SEXP mEXP, SEXP kEXP) {
    // convert input
    int m = as<int>(mEXP), k = as<int>(kEXP);
    double reg = regret(m,k);
    
    // return result
    return(wrap(reg));
}
