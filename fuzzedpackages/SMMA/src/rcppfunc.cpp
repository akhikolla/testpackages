//// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "auxfunc.h"

using namespace std;
using namespace arma;

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// algorithm ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
//[[Rcpp::export]]
Rcpp::List pga(arma::mat Phi1, arma::mat Phi2, arma::mat Phi3,
               Rcpp::NumericVector resp,
               std::string penalty,
               double zeta,
               double c,
               arma::vec lambda,
               int nlambda,
               int makelamb,
               double lambdaminratio,
               arma::mat penaltyfactor,
               double reltol,
               int maxiter,
               int steps,
               int btmax,
               int mem,
               double tau,
               double nu,
               int alg,
               int ll,
               double Lmin){
Rcpp::List output;
Rcpp::NumericVector vecY(resp);
Rcpp::IntegerVector YDim = vecY.attr("dim");
const arma::cube Z(vecY.begin(), YDim[0], YDim[1], YDim[2], false);

int ascent, ascentmax,
    bt, btenter = 0, btiter = 0,
    endmodelno = nlambda,
    n1 = Phi1.n_rows, n2 = Phi2.n_rows, n3 = Phi3.n_rows, Nog = Z.n_slices, ng = n1 * n2 * n3,
    p1 = Phi1.n_cols, p2 = Phi2.n_cols, p3 = Phi3.n_cols, p = p1 * p2 * p3,
    Stopconv = 0, Stopmaxiter = 0, Stopbt = 0;

double  alphamax, ascad = 3.7,
       delta, deltamax,
       L, lossBeta, lossProp, lossX,
       penBeta, penProp,
       relobj,
       val;

arma::vec df(nlambda),
          eig1, eig2, eig3,
          Iter(nlambda),  Pen(maxiter),
          obj(maxiter + 1),
          Stops(3),
          eevX;

arma::mat absBeta(p1, p2 * p3),
          Beta(p1, p2 * p3), Betaprev(p1, p2 * p3), Betas(p, nlambda), BT(nlambda, maxiter),
          Delta(maxiter, nlambda),  dpen(p1, p2 * p3),
          Gamma(p1, p2 * p3), GradlossX(p1, p2 * p3), GradlossXprev(p1, p2 * p3), GradlossX2(p1, p2 * p3),
          Obj(maxiter, nlambda),
          Phi1tPhi1, Phi2tPhi2, Phi3tPhi3, PhitPhiBeta, PhitPhiX, pospart(p1, p2 * p3),
          Prop(p1, p2 * p3), PhiBeta(n1, n2 * n3), PhiProp(n1, n2 * n3), PhiX(n1, n2 * n3),
          wGamma(p1, p2 * p3),
          R,
          S, Sumsqdiff(Nog, Nog),
          X(p1, p2 * p3), Xprev,
          Zi;

arma::cube PhitZ(p1, p2 * p3, Nog);

////fill variables
ascentmax = 4;
double scale = 0.9;

obj.fill(NA_REAL);
Betas.fill(42);
Iter.fill(0);
GradlossXprev.fill(0);
BT.fill(-1);

Delta.fill(NA_REAL);
Obj.fill(NA_REAL);
Pen.fill(0);

////precompute
Phi1tPhi1 = Phi1.t() * Phi1;
Phi2tPhi2 = Phi2.t() * Phi2;
Phi3tPhi3 = Phi3.t() * Phi3;
eig1 = arma::eig_sym(Phi1tPhi1);
eig2 = arma::eig_sym(Phi2tPhi2);
eig3 = arma::eig_sym(Phi3tPhi3);
alphamax = as_scalar(max(kron(eig1, kron(eig2 , eig3))));

////precompute
for(int j = 0; j < Nog; j++){

PhitZ.slice(j) = RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), Z.slice(j), n2, n3), n3, p1), p1, p2);

}

////proximal step size
// Sumsqdiff.fill(0);
// for(int i = 0; i < Nog; i++){
// Zi = Z.slice(i);
// for(int j = i + 1; j < Nog; j++){Sumsqdiff(i,j) = sum_square(Zi - Z.slice(j));}
// }

mat A(Nog, Nog);
mat PhitZi;
A.fill(0);
for(int i = 0; i < Nog; i++){
PhitZi = PhitZ.slice(i);
for(int j = i + 1; j < Nog; j++){
A(i,j) = sum_square(PhitZi- PhitZ.slice(j));
}
}

L =  4  / pow(ng, 2) * (max(max(A)) + alphamax*ng / 2); //upper bound on Lipschitz constant
//L =  4 * alphamax / pow(ng, 2) * (max(max(Sumsqdiff)) + ng / 2); //upper bound on Lipschitz constant
delta = nu * 1.9 / L; //stepsize scaled up by  nu
deltamax = 1.99 / L; //maximum theoretically allowed stepsize
if(Lmin == 0){Lmin = (1 / nu) * L;}

////initialize
Betaprev.fill(0);
Beta = Betaprev;
PhiBeta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta, p2, p3), p3, n1), n1, n2);
PhitPhiBeta = RHmat(Phi3tPhi3, RHmat(Phi2tPhi2, RHmat(Phi1tPhi1, Beta, p2, p3), p3, p1), p1, p2);

lossBeta = softmaxloss(-eev(PhiBeta, Z, ng), zeta, ll);

////make lambda sequence
if(makelamb == 1){

arma::mat Ze = zeros<mat>(n1, n2 * n3);
arma::mat absgradzeroall(p1, p2 * p2);

absgradzeroall = abs(gradloss(PhitZ, PhitPhiBeta, -eev(PhiBeta, Z, ng), ng, zeta, ll));

arma::mat absgradzeropencoef = absgradzeroall % (penaltyfactor > 0);
arma::mat penaltyfactorpencoef = (penaltyfactor == 0) * 1 + penaltyfactor;
double lambdamax = as_scalar(max(max(absgradzeropencoef / penaltyfactorpencoef)));
double m = log(lambdaminratio);
double M = 0;
double difflamb = abs(M - m) / (nlambda - 1);
double l = 0;

for(int i = 0; i < nlambda ; i++){

lambda(i) = lambdamax * exp(l);
l = l - difflamb;

}

}else{std::sort(lambda.begin(), lambda.end(), std::greater<int>());}

///////////start lambda loop
for (int j = 0; j < nlambda; j++){

Gamma = penaltyfactor * lambda(j);

ascent = 0;

//start MSA loop
for (int s = 0; s < steps; s++){

if(s == 0){

if(penalty != "lasso"){wGamma = Gamma / lambda(j);}else{wGamma = Gamma;}

}else{

if(penalty == "scad"){

absBeta = abs(Beta);
pospart = ((ascad * Gamma - absBeta) + (ascad * Gamma - absBeta)) / 2;
dpen = sign(Beta) % Gamma % ((absBeta <= Gamma) + pospart / (ascad - 1) % (absBeta > Gamma));
wGamma = abs(dpen) % Gamma / lambda(j) % (Beta != 0) + lambda(j) * (Beta == 0);

}

}

/////start proximal loop
if(alg == 1){/////////////////NPG algorithm from chen2016

for (int k = 0; k < maxiter; k++){

if(k == 0){

Betaprev = Beta;
Xprev = Betaprev;
obj(k) = lossBeta + l1penalty(wGamma, Beta);
Obj(k, j) = obj(k);
Delta(k, j) = delta;

}else{//if not the first iteration

//if(acc == 1){
//X = Beta + (k - 2) / (k + 1) * (Beta - Betaprev);//  nesterov update.....
//}else{
X = Beta;
Xprev = Betaprev;

//}

PhiX = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, X, p2, p3), p3, n1), n1, n2);
eevX = -eev(PhiX, Z, ng);
PhitPhiX = RHmat(Phi3tPhi3, RHmat(Phi2tPhi2, RHmat(Phi1tPhi1, X, p2, p3), p3, p1), p1, p2);
GradlossX = gradloss(PhitZ, PhitPhiX, eevX, ng, zeta, ll);
lossX = softmaxloss(eevX, zeta, ll);

if(k > 1){
S = X - Xprev;
R = GradlossX - GradlossXprev;
double tmp = as_scalar(accu(S % R )/ sum_square(S)); //is this corrert???
tmp = std::max(tmp, Lmin);
delta = 1 / std::min(tmp, pow(10,8));
}else{
delta = 1;
}
////proximal backtracking from chen2016
BT(j, k) = 0;
while(BT(j, k) < btmax){

Prop = prox_l1(X - delta * GradlossX, delta * wGamma);
PhiProp = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Prop, p2, p3), p3, n1), n1, n2);
lossProp = softmaxloss(-eev(PhiProp, Z, ng), zeta, ll);

val = as_scalar(max(obj(span(std::max(0, k - mem), k - 1))) - c / 2 * sum_square(Prop - Betaprev));
penProp = l1penalty(wGamma, Prop);

if (lossProp + penProp <= val + 0.0000001){

break;

}else{

delta = delta / tau; //tau>1, scaling delta down instead of L up...
BT(j, k) = BT(j, k) + 1;

}

}//end line search

////check if maximum number of proximal backtraking step is reached
if(BT(j, k) == btmax){Stopbt = 1;}

Betaprev = Beta;
GradlossXprev = GradlossX;
Beta = Prop;
lossBeta = lossProp;
penBeta = penProp;
obj(k) = lossBeta + penBeta;
Pen(k) = penBeta;
Obj(k, j) = obj(k);
Iter(j) = k;
Delta(k, j) = delta;

////proximal convergence check
relobj = abs(obj(k) - obj(k - 1)) / (reltol + abs(obj(k - 1)));

if(k < maxiter - 1 && relobj < reltol){

df(j) = p - accu((Beta == 0));
Betas.col(j) = vectorise(Beta);
obj.fill(NA_REAL);
Stopconv = 1;
break;

}else if(k == maxiter - 1){

df(j) = p - accu((Beta == 0));
Betas.col(j) = vectorise(Beta);
obj.fill(NA_REAL);
Stopmaxiter = 1;
break;

}

}

////break proximal loop if maximum number of proximal backtraking step is reached
if(Stopbt == 1 || Stopmaxiter == 1){

Betas.col(j) = vectorise(Beta);
break;

}

}//end proximal loop

}else{////FISTA

for (int k = 0; k < maxiter; k++){

if(k == 0){

Betaprev = Beta;
X = Beta;
obj(k) = lossBeta + l1penalty(wGamma, Beta);
Obj(k, j) = obj(k);
BT(j, k) = 1; //force initial backtracking
Delta(k, j) = delta;

}else{

X = Beta + (k - 2) / (k + 1) * (Beta - Betaprev);

PhiX = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, X, p2, p3), p3, n1), n1, n2);
eevX = -eev(PhiX, Z, ng);
PhitPhiX = RHmat(Phi3tPhi3, RHmat(Phi2tPhi2, RHmat(Phi1tPhi1, X, p2, p3), p3, p1), p1, p2);
GradlossX = gradloss(PhitZ, PhitPhiX, eevX, ng, zeta, ll);

////check if proximal backtracking occurred last iteration
if(BT(j, k - 1) > 0){bt = 1;}else{bt = 0;}

////check for divergence
if(ascent > ascentmax){bt  = 1;}

if((bt == 1 && deltamax < delta) || nu > 1){//backtrack

lossX = softmaxloss(eevX, zeta, ll);
BT(j, k) = 0;

while(BT(j, k) < btmax){//start backtracking

Prop = prox_l1(X - delta * GradlossX, delta * wGamma);
PhiProp = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Prop, p2, p3), p3, n1), n1, n2);
lossProp = softmaxloss(-eev(PhiProp, Z, ng), zeta, ll);

val = as_scalar(lossX + accu(GradlossX % (Prop - X))
                      + 1 / (2 * delta) * sum_square(Prop - X));

if(lossProp <= val + 0.0000001){ //need to add a little due to numerical issues

break;

}else{

delta = scale * delta;
BT(j, k) = BT(j, k) + 1;

//if(delta < deltamax){delta = deltamax;}

}

}//end backtracking
 ////check if maximum number of proximal backtraking step is reached
if(BT(j, k) == btmax){Stopbt = 1;}

}else{//no backtracking

Prop = prox_l1(X - delta * GradlossX, delta * wGamma);
PhiProp = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Prop, p2, p3), p3, n1), n1, n2);
lossProp = softmaxloss(-eev(PhiProp, Z, ng), zeta, ll);

}



    Betaprev = Beta;
    Beta = Prop;
    lossBeta = lossProp;
    obj(k) = lossBeta + l1penalty(wGamma, Beta);
    Iter(j) = k;
    Delta(k, j) = delta;
    Obj(k, j) = obj(k);


    ////proximal divergence check
    if(obj(k) > obj(k - 1)){ascent = ascent + 1;}else{ascent = 0;}

  relobj = abs(obj(k) - obj(k - 1)) / (reltol + abs(obj(k - 1)));

  if(k < maxiter - 1 && relobj < reltol){

    df(j) = p - accu((Beta == 0));
    Betas.col(j) = vectorise(Beta);
    obj.fill(NA_REAL);
    Stopconv = 1;
    break;

  }else if(k == maxiter - 1){

    df(j) = p - accu((Beta == 0));
    Betas.col(j) = vectorise(Beta);
    obj.fill(NA_REAL);
    Stopmaxiter = 1;
    break;

  }

}

////break proximal loop if maximum number of proximal backtraking step is reached
if(Stopbt == 1){

  Betas.col(j) = vectorise(Beta);
  break;

}

}//end proximal loop

}

//Stop msa loop if maximum number of backtracking steps or maxiter is reached
if(Stopbt == 1 || Stopmaxiter == 1){

endmodelno = j;
break;

}

}//end MSA loop

//Stop lambda loop if maximum number of backtracking steps or maxiter is reached
if(Stopbt == 1 || Stopmaxiter == 1){

  endmodelno = j;
  break;

}

}//end lambda loop

Stops(0) = Stopconv;
Stops(1) = Stopmaxiter;
Stops(2) = Stopbt;
btenter = accu((BT > -1));
btiter = accu((BT > 0) % BT);

output = Rcpp::List::create(Rcpp::Named("Beta") = Betas,
                            Rcpp::Named("df") = df,
                            Rcpp::Named("btenter") = btenter,
                            Rcpp::Named("btiter") = btiter,
                            Rcpp::Named("Obj") = Obj,
                            Rcpp::Named("Iter") = Iter,
                            Rcpp::Named("endmodelno") = endmodelno,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("BT") = BT,
                            Rcpp::Named("L") = L,
                            Rcpp::Named("Delta") = Delta,
                            Rcpp::Named("Sumsqdiff") = Sumsqdiff,
                            Rcpp::Named("Stops") = Stops);

return output;

}
