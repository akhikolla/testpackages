/*  
Rcpp functions used to perform penalized estimation in generalized linear array models (GLAM).
The first function, gdpg, contains a gradient descent and proximal gradient based algorithm that solves
the penalized (LASSO and SCAD) problem in the GLAM framework.
The second function, getobj, computes the objective values for the corresponding problem.

Intended for use with R.
Copyright (C) 2017 Adam Lund

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License  for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

//// [[Rcpp::depends(RcppArmadillo)]]
//#define TIMING
#include <RcppArmadillo.h>
#include "auxfunc.h"
//#include "/Users/adamlund/Documents/KU/Phd/Project/Computer/Vincent/timer/simple_timer.h"

using namespace std;
using namespace arma;

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// gd pg algorithm ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
//[[Rcpp::export]]
Rcpp::List gdpg(arma::mat Phi1, arma::mat Phi2, arma::mat Phi3, 
                arma::mat Psi, arma::mat Psirot, //rr
                arma::mat Y, arma::mat Yrot, //rr
                arma::mat V,//S
                arma::mat Weights,  arma::mat Weightsrot,
                arma::mat Betainit, arma::mat Beta12init, arma::mat Beta3init,  //rr
                arma::vec alphainit, 
                std::string family,
                std::string penalty,
                int nonten, 
                std::string iwls, 
                double nu,
                arma::vec lambda, 
                int makelamb,  
                int nlambda, 
                double lambdaminratio,
                arma::mat penaltyfactor,  arma::mat penaltyfactor12, arma::mat penaltyfactor3, //rr
                arma::mat penaltyfactoralpha,
                double reltolprox,
                double reltolnewt,
                double reltolalt, //rr
                int maxiter,
                int steps,
                int maxiterprox,
                int maxiternewt,
                int maxalt, //rr
                int btproxmax,
                int btnewtmax,
                int weightedgaussian,
                int S, 
                int RR, //rr
                int n1, int n2, int n3){
  
Rcpp::List output;
int blocks;
arma::mat penaltyfactortr, Weightstr, Psitr, Phi1tr, Phi2tr, Phi3tr, Ytr;
  
if(RR == 1){

blocks = 3;
Weightstr = Weights;
Ytr = Y;
Psitr = Psi;
Phi1tr = Phi1;
Phi2tr = Phi2;
Phi3tr = Phi3;
penaltyfactortr = penaltyfactor;
  
}else{

blocks = 2; 
maxalt = 1;
   
}

//declare some global variables
int p1 = Phi1.n_cols;
int p2 = Phi2.n_cols;
int p3 = Phi3.n_cols;
int p = p1 * p2 * p3;
int n = n1 * n2 * n3;

int btenterprox = 0, btiternewt = 0, btiterprox = 0,
endmodelno = nlambda - 1,
STOPmaxiter = 0,  STOPnewt = 0, STOPprox = 0,
q = Psi.n_cols;

double ascad = 3.7;

//make lambda sequence
if(makelamb == 1){
arma::mat absgradzeroall;
arma::mat Ze = zeros<mat>(n1, n2 * n3);

if(S == 1){

arma::mat tmpH = Weights % dtheta(Ze, family) % (Ze - Y);
arma::mat sumVtmpH(n1, n2);
sumVtmpH.fill(0);
for(int i = 0; i < n3; i++){

sumVtmpH = sumVtmpH + V(span::all, span(i * n2, (i + 1) * n2 - 1)) % tmpH(span::all, span(i * n2, (i + 1) * n2 - 1));

}
absgradzeroall = abs(Phi1.t() * sumVtmpH * Phi2 / n);

}else{absgradzeroall = abs(gradloglike(Y, Weights, Phi1, Phi2, Phi3, mu(Ze, family), Ze, n2, n3, p1, p2, n, family));}

arma::mat absgradzeropencoef = absgradzeroall % (penaltyfactor > 0);
arma::mat penaltyfactorpencoef = (penaltyfactor == 0) * 1 + penaltyfactor;

double lambdamaxalpha;
if(nonten == 1){ //compute the gradient in zero for non tensor component

arma::mat tmpalpha = Weights % dtheta(Ze, family) % (Ze - Y);
arma::mat absgradzeroallalpha = abs(Psi.t() * vectorise(tmpalpha)) / n;
arma::mat absgradzeropencoefalpha = absgradzeroallalpha % (penaltyfactoralpha > 0);
arma::mat penaltyfactorpencoefalpha = (penaltyfactoralpha == 0) * 1 + penaltyfactoralpha;
lambdamaxalpha = as_scalar(max(max(absgradzeropencoefalpha / penaltyfactorpencoefalpha)));

}else{lambdamaxalpha = 0;}
lambdamaxalpha = 0;
double lambdamax = as_scalar(max(max(max(absgradzeropencoef / penaltyfactorpencoef)), lambdamaxalpha));
double difflamb = abs(0 - log(lambdaminratio)) / (nlambda - 1);
double itr = 0;

for(int i = 0; i < nlambda ; i++){

lambda(i) = lambdamax * exp(itr);
itr = itr - difflamb;

}

}else{std::sort(lambda.begin(), lambda.end(), std::greater<int>());}
//}else{lambda = arma::sort(lambda, "descend");} better??

if(family == "gaussian"){//gaussian#################################################

if(weightedgaussian == 0){//no prior weights
  
////declare variables
double delta, 
       L,  
       relobjprox,
       sqlossBeta, sqlossBeta12, sqlossBeta3, sqlossProp, sumw12sq, sumw3sq;

arma::vec absalpha(q), alpha(q), alphaprev(q), 
          df(nlambda), dpenalpha(q),
          eig1, eig2, eig3,
          Gammaalpha(q), GradsqlossXalpha(q),
          Iter(nlambda),
          obj(maxiterprox + 1), objalt(maxalt),
          Psity, pospartalpha(q), Propalpha(q),
          wGammaalpha(q),
          Xalpha(q);

arma::mat absBeta(p1, p2 * p3),  alphas(q, nlambda),
          Beta(p1, p2 * p3), Betaa, Betaprev(p1, p2 * p3), Betas(p, nlambda), Beta12, Beta12prev, Beta3, Beta3prev, Betas12(p1 * p2, nlambda), Betas3(p3, nlambda),
          dpen(p1, p2 * p3),
          Eta,
          fix,
          Gamma(p1, p2 * p3), GradsqlossX(p1, p2 * p3),
          Phi1tPhi1, Phi2tPhi2, Phi3tPhi3,  PhitY, PhitY12(p1, p2), PhitY3(p3, 1), pospart(p1, p2 * p3), Prop(p1, p2 * p3), PsitPsi, PsirottPsirot, Psirottyrot, PsiXalpha, PsirotXalpha,
          sumVtmp(n1, n2),sumVPsiXalpha(n1, n2), sumVY(n1, n2), sumVsq(n1, n2), sumw12PsirotXalpha(n3, 1), sumw3PsiXalpha(n1, n2), 
          wGamma(p1, p2 * p3), wfix,
          X(p1, p2 * p3);

cube Lambdas(nlambda, maxalt, blocks - 1);

////fill variables
Lambdas.fill(42);
Betas.fill(42);
obj.fill(0);
Iter.fill(0);

////precompute
mat Phi1tYPhi2(p1 * p2, n3);
mat Phi3tYrot(p3, n1 * n2);

if(RR == 1){ 

for(int i  = 0; i < n3; i++){Phi1tYPhi2.col(i) = vectorise(Phi1.t() * Ytr.cols(i * n2, (i + 1) * n2 - 1) * Phi2);}
for(int i = 0; i < n1; i++){for(int j = 0; j < n2; j++){Phi3tYrot.col(i * n2 + j) = Phi3.t() * Yrot.col(i * n2 + j);}}

Phi1tPhi1 = Phi1.t() * Phi1;
Phi2tPhi2 = Phi2.t() * Phi2;
Phi3tPhi3 = Phi3.t() * Phi3;
PhitY = RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), Ytr, n2, n3), n3, p1), p1, p2);
eig1 = arma::eig_sym(Phi1tPhi1);
eig2 = arma::eig_sym(Phi2tPhi2);
eig3 = arma::eig_sym(Phi3tPhi3);

}

if(S == 1){

//weights for gradient
double sumij;
for(int i = 0; i < n1; i++) {
for(int j = 0; j < n2; j++) {
sumij = 0;
for(int itr = 0; itr < n3; itr++){sumij = sumij + pow(V(i, j + itr * n2), 2);}

sumVsq(i, j) = sumij;

}
}

sumVY.fill(0);
for(int i = 0; i < n3; i++){sumVY = sumVY + V(span::all, span(i * n2, (i + 1) * n2 - 1)) % Y(span::all, span(i * n2, (i + 1) * n2 - 1));}
PhitY = Phi1.t() * sumVY * Phi2;

}else{PhitY = RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), Y, n2, n3), n3, p1), p1, p2);}

Phi1tPhi1 = Phi1.t() * Phi1;
Phi2tPhi2 = Phi2.t() * Phi2;
Phi3tPhi3 = Phi3.t() * Phi3;

if(nonten == 1){
  
Psity = Psi.t() * vectorise(Y);
  
if(RR == 1){

Psirottyrot = Psirot.t() * vectorise(Yrot);
PsirottPsirot = Psirot.t() * Psirot;

}

}else{Psity.fill(0);}

PsitPsi = Psi.t() * Psi;
  
eig1 = arma::eig_sym(Phi1tPhi1);
eig2 = arma::eig_sym(Phi2tPhi2);
eig3 = arma::eig_sym(Phi3tPhi3);

L = as_scalar(max(eig1) * max(eig2) * max(eig3) + nonten * max(arma::eig_sym(PsitPsi))) / n;

if(S == 1){L = as_scalar(max(eig1) * max(eig2) * max(max(sumVsq)) + nonten * max(arma::eig_sym(PsitPsi))) / n;}

delta = 1 / L;

////initialize
GradsqlossX.fill(42);
PsirotXalpha.fill(42);

Beta = Betainit;
alpha = alphainit;
sqlossBeta = sqloss(Phi1, Phi2, Phi3, Psi, Y, V, Beta, alpha, n, p2, p3, n1, n2, nonten, S);

if(RR == 1){
  
Beta12 = Beta12init;
Beta3 = Beta3init;
Betaa = outermat(Beta12, Beta3);
sqlossBeta12 = sqlossrr(Phi1, Phi2, Phi3, Psi, Ytr, Beta12, Phi3 * Beta3, alpha, n, 1, nonten);
sqlossBeta3 = sqlossrr(Phi1, Phi2, Phi3, Psirot, Yrot, Beta3, Phi1 * Beta12 * Phi2.t(), alpha, n, 2, nonten);

}

//start lambda loop
for (int l = 0; l < nlambda; l++){

Gamma = penaltyfactor * lambda(l);
if(nonten == 1){Gammaalpha = penaltyfactoralpha * lambda(l);}

for(int a = 0; a < maxalt; a++){//#alternation loop over blocks

for(int b = 1; b < blocks ; b++){//#block loop

if(RR == 1){

if(b == 1){//######## block 1 //overwrite!!

Y = Ytr;
Psi = Psitr;
Beta = Beta12;
fix = Beta3;
wfix = Phi3 * fix; // #n3x1
penaltyfactor = penaltyfactor12;
delta = n * 1.99 / (max(eig1) * max(eig2) * sum_square(wfix) + nonten * max(arma::eig_sym(PsitPsi))); ////sum_square(wfix)) not max(..) because its scalars 
double lmax = lambmaxrr(Y, Phi1, Phi2, Phi3,  Psi, ones(n1, n2 * n3), wfix,  n,   b,  family);

if(lmax < lambda(l)){Lambdas(l, a, b - 1) = lmax * 0.9;}else{Lambdas(l, a, b - 1) = lambda(l);}
Gamma =  penaltyfactor12 * Lambdas(l, a, b - 1);
sqlossBeta = sqlossBeta12;

sumw3sq =  accu(pow(wfix, 2));
mat sumw3Phi1tYPhi2(p1 * p2, 1); //re-def?!?!?!?!
sumw3Phi1tYPhi2.fill(0);
for(int i = 0; i < n3; i++) {sumw3Phi1tYPhi2 = sumw3Phi1tYPhi2 + wfix(i) * Phi1tYPhi2.col(i);}
sumw3Phi1tYPhi2.reshape(p1, p2);
PhitY12 = sumw3Phi1tYPhi2;

}else if(b == 2){//######## block 2 //overwrite!!

Y = Yrot;
Psi = Psirot;
Beta = Beta3;
fix = Beta12;
wfix = Phi1 * fix * Phi2.t(); //#n1xn2
penaltyfactor = penaltyfactor3;
delta = n * 1.99 / (max(eig3) * sum_square(wfix) + nonten * max(arma::eig_sym(PsitPsi))); 
double lmax = lambmaxrr(Y, Phi1, Phi2, Phi3,  Psi, ones(n3, n1 * n2), wfix,  n,   b,  family);

if(lmax < lambda(l)){Lambdas(l, a, b - 1) = lmax * 0.9;}else{Lambdas(l, a, b - 1) = lambda(l);}
Gamma = penaltyfactor3 * Lambdas(l, a, b - 1);
sqlossBeta = sqlossBeta3; 

sumw12sq = accu(pow(wfix, 2));
PhitY3.fill(0);

for(int j = 0; j < n2; j++){for(int i = 0; i < n1; i++){PhitY3 = PhitY3 + wfix(i, j) * Phi3tYrot.col(j * n1 + i);}}

}
}

/////start MSA loop
for (int s = 0; s < steps; s++){

if(s == 0){//first MSA step

if(penalty != "lasso"){

wGamma =  Gamma / lambda(l);
wGammaalpha =  Gammaalpha / lambda(l);
if(RR == 1){wGamma = Gamma / min(Lambdas(l, a, b - 1), lambda(l));}//rr!!!!!!!

}else{

wGamma = Gamma;
wGammaalpha = Gammaalpha;

}

}else{

if(penalty == "scad"){

absBeta =  abs(Beta);
pospart = ((ascad * Gamma - absBeta) + (ascad * Gamma - absBeta)) / 2;
dpen = sign(Beta) % (Gamma % (absBeta <= Gamma) + pospart / (ascad - 1) % (absBeta > Gamma));
wGamma = abs(dpen) % penaltyfactor % (Beta != 0) + Gamma % (Beta == 0);

absalpha =  abs(alpha);
pospartalpha = ((ascad * Gammaalpha - absalpha) + (ascad * Gammaalpha - absalpha)) / 2;
dpenalpha = sign(alpha) % (Gammaalpha % (absalpha <= Gammaalpha) + pospartalpha / (ascad - 1) % (absalpha > Gammaalpha));
wGammaalpha = abs(dpenalpha) % penaltyfactoralpha % (alpha != 0) + Gammaalpha % (alpha == 0);

}

}

////start proximal loop
for (int k = 0; k < maxiterprox; k++){

if(k == 0){

obj(0) = sqlossBeta + l1penalty(wGamma, Beta) + nonten * l1penalty(wGammaalpha, alpha);
Betaprev = Beta;
alphaprev = alpha;

}else{

//next val
X = Beta + (k - 2) / (k + 1) * (Beta - Betaprev);
Xalpha = alpha + (k - 2) / (k + 1) * (alpha - alphaprev);

//gradient and proposed val
if(nonten == 1){

if(S == 1){

PsiXalpha = Psi * Xalpha; 
PsiXalpha.reshape(n1, n2 * n3);
sumVPsiXalpha.fill(0);

for(int i = 0; i < n3; i++){sumVPsiXalpha = sumVPsiXalpha + V(span::all, span(i * n2, (i + 1) * n2 - 1)) % PsiXalpha(span::all, span(i * n2, (i + 1) * n2 - 1));}

GradsqlossX = (Phi1.t() * (sumVsq % (Phi1 * X * Phi2.t())) * Phi2 + Phi1.t() * sumVPsiXalpha * Phi2 - PhitY) / n;
GradsqlossXalpha = (Psi.t() * vectorise(etafunc(Phi1, Phi2, Phi3, V, X, n, S)) + PsitPsi * Xalpha - Psity) / n;

}else if(RR == 1){

if(b == 1){

PsiXalpha = Psi * Xalpha; 
PsiXalpha.reshape(n1, n2 * n3); 
sumw3PsiXalpha.fill(0);
for(int i = 0; i < n3; i++){sumw3PsiXalpha = sumw3PsiXalpha + wfix(i) * PsiXalpha(span::all, span(i * n2, (i + 1) * n2 - 1));}

GradsqlossX = (Phi1tPhi1 * X * Phi2tPhi2 * sumw3sq + Phi1.t() * sumw3PsiXalpha * Phi2 - PhitY12) / n;
GradsqlossXalpha = (Psi.t() * vectorise(vectorise(Phi1 * X * Phi2.t()) * wfix.t()) + PsitPsi * Xalpha - Psity) / n;

}else if(b == 2){//b=2,rotated data/model!!!

PsirotXalpha = Psirot * Xalpha;
PsirotXalpha.reshape(n3, n1 * n2);
sumw12PsirotXalpha.fill(0);
 
for(int j = 0; j < n2; j++){for(int i = 0; i < n1; i++){sumw12PsirotXalpha = sumw12PsirotXalpha + wfix(i, j) * PsirotXalpha.col(j * n1 + i);}}

GradsqlossX = (Phi3tPhi3 * X * sumw12sq + Phi3.t() * sumw12PsirotXalpha - PhitY3) / n;
GradsqlossXalpha = (Psirot.t() * vectorise(Phi3 * X * vectorise(wfix).t()) + PsirottPsirot * Xalpha - Psirottyrot) / n;

}

}else{

PsiXalpha = Psi * Xalpha; 
PsiXalpha.reshape(n1, n2 * n3);
GradsqlossX = (RHmat(Phi3tPhi3, RHmat(Phi2tPhi2, RHmat(Phi1tPhi1, X, p2, p3), p3, p1), p1, p2) + RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), PsiXalpha, n2, n3), n3, p1), p1, p2) - PhitY) / n;
GradsqlossXalpha = (Psi.t() * vectorise(RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, X, p2, p3), p3, n1), n1, n2)) + PsitPsi * Xalpha - Psity) / n;

}

Prop = prox_l1(X - delta * GradsqlossX, delta * wGamma);
Propalpha = prox_l1(Xalpha - delta * GradsqlossXalpha, delta * wGammaalpha);

}else{//no nonten

if(S == 1){

GradsqlossX = (Phi1.t() * (sumVsq % (Phi1 * X * Phi2.t())) * Phi2  - PhitY) / n;

}else if(RR == 1){

if(b == 1){

GradsqlossX = (Phi1tPhi1 * X * Phi2tPhi2 * sumw3sq - PhitY12) / n;

}else if (b == 2){GradsqlossX = (Phi3tPhi3 * X * sumw12sq - PhitY3) / n;}

}else{GradsqlossX = (RHmat(Phi3tPhi3, RHmat(Phi2tPhi2, RHmat(Phi1tPhi1, X, p2, p3), p3, p1), p1, p2) - PhitY) / n;}

Prop = prox_l1(X - delta * GradsqlossX, delta * wGamma);

}

if(RR == 1){

sqlossProp =  sqlossrr(Phi1, Phi2, Phi3, Psi, Y, Prop, wfix, Propalpha, n, b, nonten);

}else{

sqlossProp = sqloss(Phi1, Phi2, Phi3, Psi, V, Y, Prop, Propalpha, n, p2, p3, n1, n2, nonten, S);

}

Betaprev = Beta;
Beta = Prop;
alphaprev = alpha;
alpha = Propalpha;
sqlossBeta = sqlossProp;

}

Iter(l) = k + 1;
obj(k + 1) = sqlossBeta + l1penalty(wGamma, Beta) + nonten * l1penalty(wGammaalpha, alpha);

////proximal convergence check //fista not descent!
relobjprox = abs(obj(k + 1) - obj(k)) / abs(obj(k));

if(k > 0 && k < maxiterprox &&  relobjprox < reltolprox){//go to next ...

obj.fill(0);
break;

}else if(k == maxiterprox){//go to next  ...

obj.fill(0);
break;

}

}//end proximal loop

//check if maximum number of iterations for current lambda is reached
if(accu(Iter) > maxiter){STOPmaxiter = 1;}

//stop program if maxiter is reached
if(STOPmaxiter == 1){

endmodelno = l;
break;

}

}//end MSA loop

if(RR == 1){if(b == 1){Beta12 = Beta;}else if(b == 2){Beta3 = Beta;}}

}//#blocks

if(RR == 1){

Betaa = outermat(Beta12, Beta3);
Eta =  RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Betaa, p2, p3), p3, n1), n1, n2);

if(nonten == 1){

mat tmp =   Psitr * alpha;
tmp.reshape(n1, n2 * n3);
Eta = Eta + tmp;
   
}
 
objalt(a) = loglike(Ytr, Weights, Eta, n, family) + penfunc(penaltyfactortr * lambda(l), Betaa, penalty) + penfunc(penaltyfactoralpha * lambda(l), alpha, penalty);

if(a > 0 && abs(objalt(a) - objalt(a - 1)) /  abs(objalt(a)) < reltolalt){break;}

}

}//#altenations

if(RR == 1){

df(l) = p1 * p2 + p3 + nonten *  q - accu((Beta12 == 0)) - accu((Beta3 == 0)) - nonten * accu((alpha == 0));
Betas12.col(l) = vectorise(Beta12);
Betas3.col(l) = vectorise(Beta3);
alphas.col(l) = alpha;

}else{

df(l) = p + nonten * q - accu((Beta == 0)) - nonten *  accu((alpha == 0));
Betas.col(l) = vectorise(Beta);
alphas.col(l) = alpha;

}
 
}//end lambda loop

output = Rcpp::List::create(Rcpp::Named("Beta") = Betas, 
                            Rcpp::Named("Beta12") = Betas12,
                            Rcpp::Named("Beta3") = Betas3,
                            Rcpp::Named("alpha") = alphas,
                            Rcpp::Named("df") = df,
                            Rcpp::Named("endmodelno") = endmodelno,
                            Rcpp::Named("Iter") = Iter,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("STOPmaxiter") = STOPmaxiter,
                            Rcpp::Named("STOPnewt") = STOPnewt,
                            Rcpp::Named("STOPprox") = STOPprox);

}else if (weightedgaussian == 1){//weigthed gaussian, if prior weights are used, solve (one) weighted ls problem

////declare variables
int ascentprox, ascentproxmax,
btprox;

double maxeig,
       delta, deltamin,
       L,  
       relobjprox, 
       sprox, 
       valprox, valalpha,
       wsqlossBeta, wsqlossProp, wsqlossX;  

arma::vec absalpha(q), alpha(q), alphaprevprox(q), 
          df(nlambda), dpenalpha(q),
          eig1, eig2, eig3,
          Gammaalpha(q), GradwsqlossXalpha(q),
          Iter(nlambda),
          objprox(maxiterprox + 1), objalt(maxalt),
          Psitwz, pospartalpha(q), Propalpha(q),
          wGammaalpha(q),
          Xalpha(q);
 
arma::mat absBeta,  alphas(q, nlambda),
          Beta, Beta12, Beta3, Betaa,  Beta12prevprox, Beta3prevprox, Betaprevprox, BTprox(nlambda, maxiterprox + 1), Betas(p, nlambda), Betas12(p1 * p2, nlambda), Betas3(p3, nlambda),
          dpen,
          fix,
          Eta(n1, n2 * n3), 
          GradwsqlossX, Gamma, 
          PhitWZ, Phi1tPhi1, Phi2tPhi2, Phi3tPhi3, pospart(p1, p2 * p3), Prop, PsitPsi, PsiXalpha, Psirottyrot, Psitwzrot, Psitwztr,
          SqrtW, SqrtWtr, SqrtWrot, SqrtWV, SqrtWZ, SqrtWZtr, SqrtWZrot, sumVtmp, sumWVsq(n1, n2), sumWVZ(n1, n2), sumw12PsiXalpha,
          tmp,
          W(n1, n2 * n3),  wfix, w12(n1, n2), w3(n3, 1), wGamma,  Wtr, Wrot(n3, n1 * n2), WZrot, WV, WZ, WZtr,
          X,
          Z, Zrot(n3, n1 * n2), Ztr(n3, n1 * n2);

cube Lambdas(nlambda, maxalt, blocks - 1);//rr!!!!!!!!!!!!!!!!fix
Lambdas.fill(42);//rr
Betas.fill(42);//rr

////fill variables
ascentproxmax = 4;
sprox = 0.9; 
objprox.fill(NA_REAL);
Betas.fill(NA_REAL);
Iter.fill(0);
BTprox.fill(-1);

Phi1tPhi1 = Phi1.t() * Phi1;
Phi2tPhi2 = Phi2.t() * Phi2;
Phi3tPhi3 = Phi3.t() * Phi3;
PsitPsi = Psi.t() * Psi;

////prior weight matrix
W = Weights;
Z = Y;
SqrtW = sqrt(W);
SqrtWZ = SqrtW % Z;
SqrtWV = SqrtW % V;
WV = W % V;
WZ = W % Z;

if(RR == 1){

Ztr = Ytr;
Zrot = Yrot;
Wtr = Weights;
Wrot = Weightsrot;
SqrtWtr = SqrtW;
SqrtWZtr = SqrtWZ;
WZtr = WZ;
SqrtWrot = sqrt(Wrot);
SqrtWZrot = SqrtWrot % Zrot;
WZrot = Wrot % Zrot;

}

if(S == 1){ //weighted weights for gradient

double sumij;

for(int i = 0; i < n1; i++) {
for(int j = 0; j < n2; j++) {
sumij = 0;
for(int itr = 0; itr < n3; itr++) {sumij = sumij + pow(SqrtWV(i, j +  itr * n2), 2);}

sumWVsq(i, j) = sumij;

}
}

sumWVZ.fill(0);
for(int i = 0; i < n3; i++){sumWVZ = sumWVZ + WV(span::all, span(i * n2, (i + 1) * n2 - 1)) % Z(span::all, span(i * n2, (i + 1) * n2 - 1));}
PhitWZ = Phi1.t() * sumWVZ * Phi2;

}else{PhitWZ = RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), WZ, n2, n3), n3, p1), p1, p2);}

if(nonten == 1){

Psitwz = Psi.t() * vectorise(WZ);

if(RR == 1){

Psitwztr = Psitwz;
Psitwzrot = Psirot.t() * vectorise(WZrot);

}

}else{Psitwz.fill(0);}

////proximal step size
eig1 = arma::eig_sym(Phi1tPhi1);
eig2 = arma::eig_sym(Phi2tPhi2);
eig3 = arma::eig_sym(Phi3tPhi3);
maxeig = as_scalar(max(eig1) * max(eig2) * max(eig3));

if(S == 1){

L = as_scalar(max(eig1) * max(eig2) * max(max(sumWVsq)) + nonten * max(arma::eig_sym(PsitPsi)) * max(max(W))) / n;

}else{

L = (maxeig + nonten * max(arma::eig_sym(PsitPsi))) * max(max(W)) / n;

}

deltamin = 1.99 / L; //minimum stepsize

//initial step size
if(nu > 0){delta = 1.9 / (nu * L);}else{delta = 1;}

////initialize
Beta = Betainit;
alpha = alphainit;
Betaprevprox = Beta;
alphaprevprox = alpha;

if(RR == 1){

Beta12prevprox = Beta12init;
Beta12 = Beta12init;
Beta3prevprox = Beta3init;
Beta3 = Beta3init;
Betaa = outermat(Beta12, Beta3);

}

//start lambda loop
for (int l = 0; l < nlambda; l++){

Gamma = penaltyfactor * lambda(l);
  if(nonten == 1){Gammaalpha = penaltyfactoralpha * lambda(l);}

ascentprox = 0;

for(int a = 0; a < maxalt; a++){//#alternation loop over blocks

for(int b = 1; b < blocks ; b++){//#block loop

if(RR == 1){ 

if(b == 1){//######## block 1 //overwrite!!
Z = Ztr;
W = Wtr;
SqrtW = SqrtWtr;
SqrtWZ = SqrtWZtr;
WZ = WZtr;
Psi = Psitr;
Psitwz = Psitwztr;

penaltyfactor = penaltyfactor12;
 Beta = Beta12;
Betaprevprox = Beta12prevprox;
fix = Beta3;
 wfix = Phi3 * Beta3; // #n3x1
 V.reshape(n1, n2 * n3);

for(int i = 0; i < n3; i++){V.cols(i * n2, (i + 1) * n2 - 1) = w12.fill(wfix(i));}

SqrtWV = SqrtW % V;
WV = W % V; 
sumWVsq.reshape(n1, n2);

double sumij;
for(int i = 0; i < n1; i++) {
for(int j = 0; j < n2; j++) {
sumij = 0;
for(int itr = 0; itr < n3; itr++){sumij = sumij + pow(SqrtWV(i, j + itr * n2), 2);}

sumWVsq(i, j) = sumij;

}
}

sumWVZ.reshape(n1, n2);
sumWVZ.fill(0);
for(int i = 0; i < n3; i++){sumWVZ = sumWVZ + WV.cols(i * n2, (i + 1) * n2 - 1) % Z.cols(i * n2, (i + 1) * n2 - 1);}
PhitWZ = Phi1.t() * sumWVZ * Phi2;

L = as_scalar(max(eig1) * max(eig2) * max(max(sumWVsq)) + nonten * max(arma::eig_sym(PsitPsi)) * max(max(W))) / n;
deltamin = 1.99 / L; //maximum theoretically allowed stepsize
if(nu > 0){delta = 1.9 / (nu * L);}else{delta = 1;}

double lmax = lambmaxrr(Ytr, Phi1, Phi2, Phi3,  Psi, Weights, wfix,  n,   b,  family);
if(lmax < lambda(l)){Lambdas(l, a, b - 1) = lmax * 0.9;}else{Lambdas(l, a, b - 1) = lambda(l);}
Gamma =  penaltyfactor12 * Lambdas(l, a, b - 1);

}else if(b == 2){//######## block 2 //overwrite!!

Z = Zrot;
W = Wrot;
SqrtW = SqrtWrot;
SqrtWZ = SqrtWZrot;
WZ = WZrot;
Psi = Psirot;
Psitwz = Psitwzrot;

penaltyfactor = penaltyfactor3;
 Beta = Beta3;
Betaprevprox = Beta3prevprox;
fix = Beta12;
wfix = Phi1 * fix * Phi2.t(); //#n1xn2
V.reshape(n3, n1 * n2); //V is rotated!!
for (int i = 0; i < n1; i++){for (int j = 0; j < n2; j++){V.col(i + j * n1) = w3.fill(wfix(i, j));}}

SqrtWV = SqrtW % V; //n3xn1n2
WV = Wrot % V; 

sumWVsq.reshape(n3, 1);
double sumi;
for(int i = 0; i < n3; i++) {

sumi = 0;
for(int j = 0; j < n1; j++) {for(int itr = 0; itr < n2; itr++){sumi = sumi + pow(SqrtWV(i, j + itr * n1), 2);}}

sumWVsq(i) = sumi;

//}?????????????
}

sumWVZ.reshape(n3, 1);

sumWVZ.fill(0);
for(int i = 0; i < n1; i++){for(int j = 0; j < n2; j++){sumWVZ = sumWVZ + WV.col(i + j * n1) % Z.col(i + j * n1);}}
PhitWZ = Phi3.t() * sumWVZ ;

L = as_scalar(max(eig3) * max(max(sumWVsq)) + nonten * max(arma::eig_sym(PsitPsi)) * max(max(W))) / n;
deltamin = 1.99 / L; //maximum theoretically allowed stepsize
if(nu > 0){delta = 1.9 / (nu * L);}else{delta = 1;}

double lmax = lambmaxrr(Yrot, Phi1, Phi2, Phi3,  Psi, Weightsrot, wfix,  n,   b,  family);

if(lmax < lambda(l)){Lambdas(l, a, b - 1) = lmax * 0.9;}else{Lambdas(l, a, b - 1) = lambda(l);}
Gamma = penaltyfactor3 * Lambdas(l, a, b - 1);

}

}

//start MSA loop
for (int s = 0; s < steps; s++){

if(s == 0){//first MSA step

if(penalty != "lasso"){

wGamma =  Gamma / lambda(l);
wGammaalpha =  Gammaalpha / lambda(l);

}else{

wGamma = Gamma;
wGammaalpha = Gammaalpha;

}

}else{

  if(penalty == "scad"){

absBeta =  abs(Beta);
pospart = ((ascad * Gamma - absBeta) + (ascad * Gamma - absBeta)) / 2;
dpen = sign(Beta) % (Gamma % (absBeta <= Gamma) + pospart / (ascad - 1) % (absBeta > Gamma));
wGamma = abs(dpen) % penaltyfactor % (Beta != 0) + Gamma % (Beta == 0);

absalpha =  abs(alpha);
pospartalpha = ((ascad * Gammaalpha - absalpha) + (ascad * Gammaalpha - absalpha)) / 2;
dpenalpha = sign(alpha) % (Gammaalpha % (absalpha <= Gammaalpha) + pospartalpha / (ascad - 1) % (absalpha > Gammaalpha));
wGammaalpha = abs(dpenalpha) % penaltyfactoralpha % (alpha != 0) + Gammaalpha % (alpha == 0);

}

}

/////start proximal loop
for (int k = 0; k < maxiterprox; k++){

if(k == 0){

Betaprevprox = Beta;
alphaprevprox = alpha;

  if(RR == 1){

objprox(0) =  wsqlossrr(Phi1, Phi2, Phi3, Psi, SqrtWZ, Beta, SqrtWV, alpha, n, b, nonten) + l1penalty(wGamma, Beta) + nonten * l1penalty(wGammaalpha, alpha);

  }else{

objprox(0) = wsqloss(SqrtW, SqrtWV, Phi1, Phi2, Phi3, Psi, SqrtWZ, Beta, alpha, n, p2, p3, n1, n2, nonten, S)
+ l1penalty(wGamma, Beta) +  nonten * l1penalty(wGammaalpha, alpha);

  }

  BTprox(l, k) = 1; //force initial backtracking (if deltamin < delta)

}else{

X = Beta + (k - 2) / (k + 1) * (Beta - Betaprevprox);
Xalpha = alpha + (k - 2) / (k + 1) * (alpha - alphaprevprox);

//gradient and proposed val
if(nonten == 1){

PsiXalpha = Psi * Xalpha;

if(S == 1){

  PsiXalpha.reshape(n1, n2 * n3);

tmp =  W % (etafunc(Phi1, Phi2, Phi3, V, X, n, S) + PsiXalpha);//remove W% by precomputing W%Psi and using V=sqrtWV
sumVtmp.reshape(n1, n2);
sumVtmp.fill(0);
for(int i = 0; i < n3; i++){sumVtmp = sumVtmp + V(span::all, span(i * n2, (i + 1) * n2 - 1)) % tmp(span::all, span(i * n2, (i + 1) * n2 - 1));}//span to cols!!!!!!!

GradwsqlossX =  (Phi1.t() * sumVtmp * Phi2 - PhitWZ) / n;
GradwsqlossXalpha = (Psi.t() * vectorise(tmp) - Psitwz) / n;

}else if(RR == 1){

  if(b == 1){

PsiXalpha.reshape(n1, n2 * n3);
sumVtmp.reshape(n1, n2);
tmp =  W % (etafunc(Phi1, Phi2, Phi3, V, X, n, b) + PsiXalpha);//remove W% by precomputing W%Psi and using V=sqrtWV
sumVtmp.fill(0);///???????dims!!!
for(int i = 0; i < n3; i++){sumVtmp = sumVtmp + V.cols(i * n2, (i + 1) * n2 - 1) % tmp.cols(i * n2, (i + 1) * n2 - 1);}
GradwsqlossX =  (Phi1.t() * sumVtmp * Phi2 - PhitWZ) / n;
GradwsqlossXalpha = (Psi.t() * vectorise(tmp) - Psitwz) / n;

}else if(b == 2){//b=2,rotated data/model

PsiXalpha.reshape(n3, n1 * n2);
sumVtmp.reshape(n3, 1);

tmp =  W % (etafunc(Phi1, Phi2, Phi3, V, X, n, b) + PsiXalpha);//remove W% by precomputing W%Psi and using V=sqrtWV
sumVtmp.fill(0);///???????dims!!!
for(int j = 0; j < n2; j++){for(int i = 0; i < n1; i++){sumVtmp = sumVtmp + V.col(i + j * n1) % tmp.col(i + j * n1);}}
GradwsqlossX =  (Phi3.t() * sumVtmp - PhitWZ) / n;
GradwsqlossXalpha = (Psi.t() * vectorise(tmp) - Psitwz) / n;

}

}else{

PsiXalpha.reshape(n1, n2 * n3);
tmp =  W % (RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, X, p2, p3), p3, n1), n1, n2) + PsiXalpha);//remove W% by precomputing W%Psi and using V=sqrtWV
GradwsqlossX = (RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), tmp, n2, n3), n3, p1), p1, p2) - PhitWZ) / n;
GradwsqlossXalpha = (Psi.t() * vectorise(tmp) - Psitwz) / n;

}

}else{//no nonten

if(S == 1){

  GradwsqlossX = (Phi1.t() * (sumWVsq % (Phi1 * X * Phi2.t())) * Phi2 - PhitWZ) / n;

}else if(RR == 1){//rr

if(b == 1){

GradwsqlossX = (Phi1.t() * (sumWVsq % (Phi1 * X * Phi2.t())) * Phi2 - PhitWZ) / n;

}else if(b == 2){//rotated data/model!!!b=2

GradwsqlossX = (Phi3.t() * (sumWVsq % (Phi3 * X)) - PhitWZ) / n;

}

  }else{

GradwsqlossX = (winprod(W, Phi1, Phi2, Phi3, X, n1, n2, n3, p1, p2, p3) - PhitWZ) / n;

}

}

////check if proximal backtracking occurred last iteration
if(BTprox(l, k - 1) > 0){btprox = 1;}else{btprox = 0;}

////check for divergence
if(ascentprox > ascentproxmax){btprox = 1;}

if((btprox == 1 && deltamin < delta) || nu == 0){//backtrack

  if(RR == 1){

   wsqlossX = wsqlossrr(Phi1, Phi2, Phi3, Psi, SqrtWZ, X, SqrtWV, alpha, n, b, nonten);//  fix !!!!!!!ok??

  }else{

wsqlossX = wsqloss(SqrtW, SqrtWV, Phi1, Phi2, Phi3, Psi, SqrtWZ, X, Xalpha, n, p2, p3, n1, n2, nonten, S);

  }

////proximal backtracking line search
BTprox(l, k) = 0;

while(BTprox(l, k) < btproxmax){//start backtracking

Prop = prox_l1(X - delta * GradwsqlossX, delta * wGamma);

if(nonten == 1){

Propalpha = prox_l1(Xalpha - delta * GradwsqlossXalpha, delta * wGammaalpha);
valalpha = accu(GradwsqlossXalpha % (Propalpha - Xalpha)) + 1 / (2 * delta) * sum_square(Propalpha - Xalpha);

}else{valalpha = 0;}

if(RR == 1){

  wsqlossProp = wsqlossrr(Phi1, Phi2, Phi3, Psi, SqrtWZ, Prop, SqrtWV, alpha, n, b, nonten);//  fix !!!!!!!ok??

  }else{

wsqlossProp = wsqloss(SqrtW, SqrtWV, Phi1, Phi2, Phi3, Psi, SqrtWZ, Prop, Propalpha, n, p2, p3, n1, n2, nonten, S);

  }

valprox = as_scalar(wsqlossX + accu(GradwsqlossX % (Prop - X)) + 1 / (2 * delta) * sum_square(Prop - X) + valalpha);

if (wsqlossProp <= valprox + 0.0000001){ //need to add a little due to numerical issues

break;

}else{

delta = sprox * delta;
BTprox(l, k) = BTprox(l, k) + 1;

if(delta < deltamin){delta = deltamin;}

}

}//end backtracking

////check if maximum number of proximal backtraking step is reached
if(BTprox(l, k) == btproxmax){STOPprox = 1;}

}else{//no backtracking

Prop = prox_l1(X - delta * GradwsqlossX, delta * wGamma);
if(nonten == 1){Propalpha = prox_l1(Xalpha - delta * GradwsqlossXalpha, delta * wGammaalpha);}else{Propalpha = alpha;}

if( RR == 1){//fix!!!

 wsqlossProp =  wsqlossrr(Phi1, Phi2, Phi3, Psi, SqrtWZ, Prop, SqrtWV, alpha, n, b, nonten);//  fix !!!!!!!ok??

}else{

  wsqlossProp = wsqloss(SqrtW, SqrtWV, Phi1, Phi2, Phi3, Psi, SqrtWZ, Prop, Propalpha, n, p2, p3, n1, n2, nonten, S);

  }

}

Betaprevprox = Beta;
alphaprevprox = alpha;
Beta = Prop;
alpha = Propalpha;
wsqlossBeta = wsqlossProp;
objprox(k) = wsqlossBeta + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);
Iter(l) = k;

////proximal divergence check
if(objprox(k) > objprox(k - 1)){ascentprox = ascentprox + 1;}else{ascentprox = 0;}

////proximal convergence check
relobjprox = abs(objprox(k) - objprox(k - 1)) / (reltolprox + abs(objprox(k - 1)));

if(k < maxiterprox && relobjprox < reltolprox){

objprox.fill(NA_REAL);
break;

}else if(k == maxiterprox){

objprox.fill(NA_REAL);
break;

}

}

////break proximal loop if maximum number of proximal backtraking step is reached
if(STOPprox == 1){break;}

}//end proximal loop

if(RR == 1){

  if(b == 1){

Beta12 = Beta;
Beta12prevprox = Betaprevprox;

  }

  if(b == 2){

Beta3 = Beta;
Beta3prevprox = Betaprevprox;

  }

}

//stop program if maximum number of backtracking steps or maxiter is reached
if(STOPprox == 1){

endmodelno = l;
break;

}

}//end MSA loop

}//#blocks

if(RR == 1){

Betaa = outermat(Beta12, Beta3);
Eta =  RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Betaa, p2, p3), p3, n1), n1, n2);

if(nonten == 1){

mat  Psialpha =   Psitr * alpha;
Psialpha.reshape(n1, n2 * n3);
Eta = Eta + Psialpha;

}

objalt(a) = loglike(Ytr, Weights, Eta, n, family) + penfunc(penaltyfactortr * lambda(l), Betaa, penalty) + penfunc(penaltyfactoralpha * lambda(l), alpha, penalty);

if(a > 0 && abs(objalt(a) - objalt(a - 1)) /  abs(objalt(a)) < reltolalt){break;}

}

}//#altenations

if(RR == 1){
  
df(l) = p1 * p2 + p3 + nonten *  q - accu((Beta12 == 0)) - accu((Beta3 == 0)) - nonten * accu((alpha == 0));
Betas12.col(l) = vectorise(Beta12);
Betas3.col(l) = vectorise(Beta3);
alphas.col(l) = alpha;

  }else{

df(l) = p + nonten * q - accu((Beta == 0)) - nonten *  accu((alpha == 0));
Betas.col(l) = vectorise(Beta);
alphas.col(l) = alpha;

  }


}//end lambda loop

btenterprox = accu((BTprox > -1));
btiterprox = accu((BTprox > 0) % BTprox);

output = Rcpp::List::create(Rcpp::Named("Beta") = Betas,
Rcpp::Named("Beta12") = Betas12, 
Rcpp::Named("Beta3") = Betas3,
Rcpp::Named("alpha") = alphas,
Rcpp::Named("btenterprox") = btenterprox,
Rcpp::Named("btiterprox") = btiterprox,
Rcpp::Named("df") = df,
Rcpp::Named("endmodelno") = endmodelno,
Rcpp::Named("Iter") = Iter,
Rcpp::Named("lambda") = lambda,
Rcpp::Named("STOPmaxiter") = STOPmaxiter,
Rcpp::Named("STOPnewt") = STOPnewt,
Rcpp::Named("STOPprox") = STOPprox);  
  
}


}else{//general#################################################

if (iwls != "exact"){ //use a kronecker approximation and solve subproblems as pure ls

////declare variables

double constnewt = 0.2, 
 delta, 
 L, loglikeBeta, logliketmp,  lmax,
 relobjnewt, relobjprox, 
 sqlossBeta, sqlossProp, snewt = 0.5, sumw12sq, sumw3sq, 
 tmp, tnewt, 
 valalphanewt, valnewt;
 
arma::uvec idx(n3), tmpidx(n3);   
 
arma::vec absalpha(q), alpha(q), alphaprevprox(q), alphaprevnewt(q), 
Deltaalpha(q), df(nlambda), dpenalpha(q),
eig1, eig2, eig3, 
Gammaalpha(q), GradsqlossXalpha(q),  gradloglikealpha(q),
mwtrue(n3), 
objalt(maxalt), objnewt(maxiternewt + 1), objprox(maxiterprox + 1),
pospartalpha(q), Propalpha(q), Psitwz(q),
vhat1(n1), vhat2(n2), vhat3(n3),
w(n), wGammaalpha(q),
Xalpha(q);

arma::mat absBeta(p1, p2 * p3), alphas(q, nlambda),
Beta ,  Betaa, Beta12, Beta3, Betaprevnewt, Betaprevnewt12 , Betaprevnewt3, Betaprevprox , Betaprevprox12, Betaprevprox3 ,  Betas(p, nlambda), Betas12(p1 * p2, nlambda), Betas3(p3, nlambda), BTnewt(nlambda, maxiternewt + 1), 
DeltaBeta,  dpen(p1, p2 * p3), 
Eta, Etatmp, Etaalpha(n, 1),
fix,
Gamma, GradloglikeBeta, GradsqlossX, 
Iter(nlambda, maxiternewt), 
MuEta, MuEtatmp, 
Phi1tPhi1, Phi2tPhi2, Phi3tPhi3, Phi1tW1Phi1, Phi2tW2Phi2, Phi3tW3Phi3, Phi3tYrot(p3, n1 * n2), Phi1tYPhi2(p1 * p2, n3), PhitY12, PhitY3(p3, 1), PhitWZ, pospart(p1, p2 * p3), Prop(p1, p2 * p3), PsitPsi(q, q), PsitWPsi(q, q), PsiXalpha, PsirotXalpha, Psity,
SqrtW, SqrtW1(n1, n1), SqrtW2(n2, n2), SqrtW3(n3, n3), SqrtW1Phi1, SqrtW2Phi2, SqrtW3Phi3, SqrtWV, SqrtWZ, SqrtWPsi(n, q), Submat(n1, n2), sumVtmp(n1, n2), sumWVZ(n1, n2), sumw3Phi1tYPhi2(p1 * p2, 1), sumw3PsiXalpha(n1, n2), sumw12PsiXalpha(n3, 1), 
Tmp,
U,
W(n1, n2 * n3), W1(n1, n1), W2(n2, n2), W3(n3, n3), wfix, wapp, wtmp, wtrue, wGamma(p1, p2 * p3), WV, sumWVsq(n1, n2), Wtrue(n1, n2 * n3), Wtruetr,WZ,
X(p1, p2 * p3), 
Z(n1, n2 * n3);

cube  Lambdas(nlambda, maxalt, blocks - 1);//rr

////fill variables
W1.eye();
SqrtW1 = W1;
W2.eye();
SqrtW2 = W2;
W3.eye();
SqrtW3 = W3;
W.fill(1);

Iter.fill(0);
BTnewt.fill(-1);
Betas.fill(NA_REAL);
objnewt.fill(NA_REAL);
objprox.fill(NA_REAL);

//idx = 0 * n2, 1 * n2,  2 * n2, ...,  (n3 - 1) * n2
for (int i = 0; i < n3; i++){idx(i) = i * n2;}

////initialize
Beta = Betainit; 
Betaprevprox = Beta; 
Betaprevnewt = Beta;
alpha  = alphainit; 
alphaprevprox = alpha; 
alphaprevnewt = alpha;

if(RR == 1){//rr
  
  Betaprevprox12 = Beta12init;
  Betaprevnewt12 = Beta12init;
  Beta12 = Beta12init;
  Betaprevprox3 = Beta3init;
  Betaprevnewt3 = Beta3init;
  Beta3 = Beta3init;
  Betaa = outermat(Beta12, Beta3);
  
}

  if(S == 1){

Eta = etafunc(Phi1, Phi2, Phi3, V, Beta, n, S);

  }else if (RR == 1){

Eta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Betaa, p2, p3), p3, n1), n1, n2);

  } else{

Eta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta, p2, p3), p3, n1), n1, n2);

  }
  
if(nonten == 1){

Etaalpha = Psi * alpha;
Etaalpha.reshape(n1, n2 * n3);
Eta = Eta + Etaalpha;

}

loglikeBeta = loglike(Y, Weightstr, Eta, n, family);
MuEta = mu(Eta, family);

/////lambda loop
for (int l = 0; l < nlambda; l++){

Gamma = penaltyfactor * lambda(l);
if(nonten == 1){Gammaalpha = penaltyfactoralpha * lambda(l);}

for(int a = 0; a < maxalt; a++){//#alternation loop over blocks
  
  for(int b = 2; b < blocks ; b++){//#block loop

    if(RR == 1){

if(b == 1){//######## block 1 //overwrite!!

  Weights = Weightstr;
  Y = Ytr;
  Psi = Psitr;
   Beta = Beta12;
   Betaprevnewt = Beta;
  fix = Beta3;
  wfix = Phi3 * fix; // #n3x1
  sumw3sq =  accu(pow(wfix, 2));
   penaltyfactor = penaltyfactor12;
  lmax = lambmaxrr(Y, Phi1, Phi2, Phi3, Psi,  Weights,  wfix,  n,   b, family);
 if(lmax < lambda(l)){Lambdas(l, a, b - 1) = lmax * 0.9;}else{Lambdas(l, a, b - 1) = lambda(l);}
  Gamma =  penaltyfactor12 * Lambdas(l, a, b - 1);
   Eta = vectorise(Phi1 * Beta * Phi2.t()) * wfix.t();
   Eta.reshape(n1, n2 * n3);//ok????
   if(nonten == 1){

 Etaalpha = Psi * alpha;
 Etaalpha.reshape(n1, n2 * n3);
 Eta = Eta + Etaalpha;

   }
}

if(b == 2){//######## block 2 //overwrite!!

  Weights = Weightsrot;
  Y = Yrot;
   Psi = Psirot;
  Beta = Beta3;
   Betaprevnewt = Beta;
  fix = Beta12;
  wfix = Phi1 * fix * Phi2.t(); //#n1xn2
  sumw12sq = accu(pow(wfix, 2));
  penaltyfactor = penaltyfactor3;
 lmax = lambmaxrr(Y, Phi1, Phi2, Phi3, Psi,  Weightsrot,  wfix,  n,   b, family);
  if(lmax < lambda(l)){Lambdas(l, a, b - 1) = lmax * 0.9;}else{Lambdas(l, a, b - 1) = lambda(l);}
  Gamma = penaltyfactor3 * Lambdas(l, a, b - 1);

   Eta =  Phi3 * (Beta * vectorise(wfix).t());
   
  if(nonten == 1){

Etaalpha = Psi * alpha;
Etaalpha.reshape(n3, n1 * n2);
Eta = Eta + Etaalpha;

  }
 }

 loglikeBeta = loglike(Y, Weights, Eta, n, family);
 MuEta = mu(Eta, family);

}

///start MSA loop
for (int s = 0; s < steps; s++){

  if(s == 0){//first MSA step

if(penalty != "lasso"){

wGamma =  Gamma / lambda(l);
wGammaalpha =  Gammaalpha / lambda(l);

}else{

wGamma = Gamma;
wGammaalpha = Gammaalpha;

}

  }else{

if(penalty == "scad"){

absBeta =  abs(Beta);
pospart = ((ascad * Gamma - absBeta) + (ascad * Gamma - absBeta)) / 2;
dpen = sign(Beta) % (Gamma % (absBeta <= Gamma) + pospart / (ascad - 1) % (absBeta > Gamma));
wGamma = abs(dpen) % penaltyfactor % (Beta != 0) + Gamma % (Beta == 0);

absalpha =  abs(alpha);
pospartalpha = ((ascad * Gammaalpha - absalpha) + (ascad * Gammaalpha - absalpha)) / 2;
dpenalpha = sign(alpha) % (Gammaalpha % (absalpha <= Gammaalpha) + pospartalpha / (ascad - 1) % (absalpha > Gammaalpha));
wGammaalpha = abs(dpenalpha) % penaltyfactoralpha % (alpha != 0) + Gammaalpha % (alpha == 0);

}

  }

objnewt(0) = loglikeBeta + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);

/////outer newton loop
for (int m = 0; m < maxiternewt; m++){

////true iwls weights W
Wtrue = Weights % dmu(Eta, family) % dtheta(Eta, family); //a * (mu’)^2 * theta’/mu’ = a * theta'mu'
 Wtruetr = Wtrue;
  if(RR == 1){
  if(b == 2){
wtrue = vectorise(Wtrue);
wtmp.reshape(n, 1);
   int  ind = 0;
int n1n2 = (n1 * n2);

for (int i = 0; i < n3; i++){
for (int j = 0; j < n1n2; j++){
wtmp(ind) = wtrue(i + j * n3);
ind = ind + 1;
}
  }

wtmp.reshape(n1, n2 * n3);
Wtrue = wtmp;
  }
  W.reshape(n1, n2 * n3);

}
////approximate newton weights
if(iwls == "kron1"){

tmp = exp(1.0 / (n1 * n2 * n3) * sum(sum(log(Wtrue))));
//pow(prod(prod(Wtrue)), 1 / (n1 * n2 * n3));

for (int i1 = 0; i1 < n1; i1++){

vhat1(i1) = exp(1.0 / (n2 * n3) * sum(log(Wtrue.row(i1))));
//pow(prod(Wtrue.row(i1)), 1 / (n2 * n3));

}

for (int i2 = 0; i2 < n2; i2++){

tmpidx = i2 + idx;
vhat2(i2) = exp(1.0 / (n1 * n3) * sum(sum(log(Wtrue.cols(tmpidx)))));
//pow(prod(prod(Wtrue.cols(tmpidx))), 1 / (n1 * n3));

}

for (int i3 = 0; i3 < n3; i3++){

vhat3(i3) = exp(1.0 / (n1 * n2) * sum(sum(log(Wtrue.cols(i3 * n2, (i3 + 1) * n2 - 1)))));
//pow(prod(prod(Wtrue.cols(i3 * n2, (i3 + 1) * n2 - 1))), 1 / (n1 * n2));

}

W1.diag() = vhat1 / pow(tmp, 2);
W2.diag() = vhat2;
W3.diag() = vhat3;

for (int i1 = 0; i1 < n1; i1++){
for (int i2 = 0; i2 < n2; i2++){
for (int i3 = 0; i3 < n3; i3++){

W(i1, i2 + i3 * n2) = vhat1(i1) * vhat2(i2) * vhat3(i3);

}
}
}

W = W / pow(tmp, 2);

}else if (iwls == "kron2"){

//kron approx to W ie W approx W3 kron W2 kron W1. How to pick W3,W2,W1?
//let W2 = I, W1 = I and W3 such that diag(W3) = mwtrue and offdiag(W3) = 0 then
//each of the first n1n2 vals in diag(W) are approx by their average i.e. mwtrue(1),
//each of the next n1n2 vals in diag(W) are approx by their average i.e. mwtrue(2), etc

for (int i3 = 0; i3 < n3; i3++){

mwtrue(i3) = mean(mean(Wtrue.cols(i3 * n2, (i3 + 1) * n2 - 1)));
Submat.fill(mwtrue(i3));
W.cols(i3 * n2, (i3 + 1) * n2 - 1) = Submat;

}

W3.diag() = mwtrue;

}

if(RR == 1){
  if(b == 2){
wapp =  vectorise(W);
wtmp.reshape(n, 1);
int ind = 0;
int n1n2 = (n1 * n2);
for (int i = 0; i < n1n2; i++){
for (int j = 0; j < n3; j++){
  wtmp(ind) = wapp(i + j * n1n2);
  ind = ind + 1;
}
}

wtmp.reshape(n3, n1 * n2);
W = wtmp; //overwrite!!!
  }

}

////working response with approx weights
U = dtheta(Eta, family) % (Y - mu(Eta, family));  //theta'(Eta) * (Y - mu(Eta)) / psi
Z = pow(W, -1) % U + Eta;

//precompute using that  now W has tensor structure
SqrtW1 = sqrt(W1);
SqrtW2 = sqrt(W2);
SqrtW3 = sqrt(W3);

SqrtW1Phi1 = SqrtW1 * Phi1;
SqrtW2Phi2 = SqrtW2 * Phi2;
if(S == 1){SqrtW3Phi3 = Phi3;}else{SqrtW3Phi3 = SqrtW3 * Phi3;}
if(nonten == 1){

w = vectorise(W);
for (int i = 0; i < n; i++){SqrtWPsi.row(i) = sqrt(w(i)) * Psi.row(i);}

}

Phi1tW1Phi1 = Phi1.t() * W1 * Phi1;
Phi2tW2Phi2 = Phi2.t() * W2 * Phi2;
if(S == 1){Phi3tW3Phi3 = Phi3.t() * Phi3;}else{Phi3tW3Phi3 = Phi3.t() * W3 * Phi3;}

if(nonten == 1){PsitWPsi = SqrtWPsi.t() * SqrtWPsi;}else{PsitWPsi.fill(0);}

SqrtW = sqrt(W);
SqrtWZ = SqrtW % Z;
WZ = W % Z;
if(S == 1){
SqrtWV = SqrtW % V;
WV = W % V;
}


if(S == 1){//weighted weights

double sumij;
for(int i = 0; i < n1; i++){
for(int j = 0; j < n2; j++){
sumij = 0;
for(int itr = 0; itr < n3; itr++){sumij = sumij + pow(SqrtWV(i, j + itr * n2), 2);}

sumWVsq(i, j) = sumij;

}
}

sumWVZ.fill(0);
for(int i = 0; i < n3; i++){sumWVZ = sumWVZ + WV.cols(i * n2, (i + 1) * n2 - 1) % Z.cols(i * n2, (i + 1) * n2 - 1);}
PhitWZ = Phi1.t() * sumWVZ * Phi2;

}else if(RR == 0){//PhitY is  used when RR==1

PhitWZ = RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), WZ, n2, n3), n3, p1), p1, p2);

}

if(nonten == 1){Psitwz = Psi.t() * vectorise(WZ);}else{Psitwz.fill(0);}

//proximal step size
eig1 = arma::eig_sym(Phi1tW1Phi1);
eig2 = arma::eig_sym(Phi2tW2Phi2);
eig3= arma::eig_sym(Phi3tW3Phi3);
if(RR != 1){
if(S == 1){

L = as_scalar(max(eig1) * max(eig2) * max(max(sumWVsq)) + nonten * max(arma::eig_sym(PsitWPsi))) / n;

}else {L = (max(eig1) * max(eig2) * max(eig3) * max(max(W)) + nonten * max(arma::eig_sym(PsitWPsi)) ) / n;}

delta = 1 / L; //can go up to 2 / L!

}

if(RR == 1){////overwrite to get sqloss problem, no tr or rot here !!!

Y = SqrtWZ;
Phi1 = SqrtW1Phi1;//???????????  overwrite like that!!!!
Phi2 = SqrtW2Phi2;//???????????  overwrite like that!!!!
Phi3 = SqrtW3Phi3;//???????????  overwrite like that!!!!
Psi = SqrtWPsi;//???????????  overwrite like that!!!!
Phi1tPhi1 = Phi1tW1Phi1;
Phi2tPhi2 = Phi2tW2Phi2;
Phi3tPhi3 = Phi3tW3Phi3;
PsitPsi = PsitWPsi;
Psity = Psi.t() * vectorise(Y);

if(b == 1){//######## block 1

for(int i  = 0; i < n3; i++){Phi1tYPhi2.col(i) = vectorise(Phi1.t() * Y.cols(i * n2, (i + 1) * n2 - 1) * Phi2);}
delta = n * 1.99 / (max(eig1) * max(eig2) * sum_square(wfix) + nonten * max(arma::eig_sym(PsitWPsi))  ); ////sum_square(wfix)) not max(..) because its scalars

sumw3Phi1tYPhi2.zeros(p1 * p2, 1);
for(int i = 0; i < n3; i++) {sumw3Phi1tYPhi2 = sumw3Phi1tYPhi2 + wfix(i) * Phi1tYPhi2.col(i);}
sumw3Phi1tYPhi2.reshape(p1, p2);
PhitY12 = sumw3Phi1tYPhi2;

}else if(b == 2){//######## block 2 //overwrite!!

for(int i = 0; i < n1; i++){for(int j = 0; j < n2; j++){Phi3tYrot.col(i * n2 + j) = Phi3.t() * Yrot.col(i * n2 + j);}}
delta = n * 1.99 / (max(eig3) * sum_square(wfix) + nonten * max(arma::eig_sym(PsitWPsi)) ); //fix nonten!!!!!!!!!! ok?
PhitY3.fill(0);
for(int j = 0; j < n2; j++){for(int i = 0; i < n1; i++){PhitY3 = PhitY3 + wfix(i, j) * Phi3tYrot.col(j * n1 + i);}}

}

}

// ///proximal loop
for (int k = 0; k < maxiterprox; k++){

if(k == 0){

Betaprevprox = Beta;
alphaprevprox = alpha;
if(RR == 1){

objprox(0) = sqlossrr(Phi1, Phi2, Phi3, Psi, Y, Beta, wfix, alpha, n, b, nonten) + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);

}else{

objprox(0) = sqloss(SqrtW1Phi1, SqrtW2Phi2, SqrtW3Phi3, SqrtWPsi, V, SqrtWZ, Beta, alpha, n, p2, p3, n1, n2, nonten, S) + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);

}

}else{

X = Beta + (k - 2) / (k + 1) * (Beta - Betaprevprox);

if(nonten == 1){

Xalpha = alpha + (k - 2) / (k + 1) * (alpha - alphaprevprox);

if(S == 1){
  
mat tmp = SqrtWPsi * Xalpha;
tmp.reshape(n1, n2 * n3);
sumVtmp.fill(0);
for(int i = 0; i < n3; i++){sumVtmp = sumVtmp + V.cols(i * n2, (i + 1) * n2 - 1) % tmp.cols(i * n2, (i + 1) * n2 - 1);}

GradsqlossX = (SqrtW1Phi1.t() * (sumWVsq % (SqrtW1Phi1 * X * SqrtW2Phi2.t())) * SqrtW2Phi2 + SqrtW1Phi1.t() * sumVtmp * SqrtW2Phi2 - PhitWZ) / n;
GradsqlossXalpha = (SqrtWPsi.t() * vectorise(etafunc(SqrtW1Phi1, SqrtW2Phi2, Phi3, V, X, n, S)) + PsitWPsi * Xalpha - Psitwz) / n;

}else if(RR == 1){

PsiXalpha = Psi * Xalpha;

if(b == 1){

PsiXalpha.reshape(n1, n2 * n3);
sumw3PsiXalpha.fill(0);
for(int i = 0; i < n3; i++){sumw3PsiXalpha = sumw3PsiXalpha + wfix(i) * PsiXalpha(span::all, span(i * n2, (i + 1) * n2 - 1));}

GradsqlossX = (Phi1tPhi1 * X * Phi2tPhi2 * sumw3sq + Phi1.t() * sumw3PsiXalpha * Phi2 - PhitY12) / n;
GradsqlossXalpha = (Psi.t() * vectorise(vectorise(Phi1 * X * Phi2.t()) * wfix.t()) + PsitPsi * Xalpha - Psity) / n;

}else if (b == 2){//b=2,rotated data/model!!!

PsiXalpha.reshape(n3, n1 * n2);
sumw12PsiXalpha.fill(0);

for(int j = 0; j < n2; j++){
for(int i = 0; i < n1; i++){

sumw12PsiXalpha = sumw12PsiXalpha + wfix(i, j) * PsiXalpha.col(j * n1 + i);

}
}

GradsqlossX = (Phi3tPhi3 * X * sumw12sq + Phi3.t() * sumw12PsiXalpha - PhitY3) / n;
GradsqlossXalpha = (Psi.t() * vectorise(Phi3 * X * vectorise(wfix).t()) + PsitPsi * Xalpha - Psity) / n;

}

}else{

mat tmp = SqrtWPsi * Xalpha;
tmp.reshape(n1, n2 * n3);

GradsqlossX = (RHmat(Phi3tW3Phi3, RHmat(Phi2tW2Phi2, RHmat(Phi1tW1Phi1, X, p2, p3), p3, p1), p1, p2) + RHmat(SqrtW3Phi3.t(), RHmat(SqrtW2Phi2.t(), RHmat(SqrtW1Phi1.t(), tmp, n2, n3), n3, p1), p1, p2) - PhitWZ) / n;
GradsqlossXalpha = (SqrtWPsi.t() * vectorise(RHmat(SqrtW3Phi3, RHmat(SqrtW2Phi2, RHmat(SqrtW1Phi1, X, p2, p3), p3, n1), n1, n2)) + PsitWPsi * Xalpha - Psitwz) / n;

}

Prop = prox_l1(X - delta * GradsqlossX, delta * wGamma);
Propalpha = prox_l1(Xalpha - delta * GradsqlossXalpha, delta * wGammaalpha);

}else{//if no noten

if(S == 1){

GradsqlossX = (SqrtW1Phi1.t() * (sumWVsq % (SqrtW1Phi1 * X * SqrtW2Phi2.t())) * SqrtW2Phi2 - PhitWZ) / n;

}else if(RR == 1){

if(b == 1){

GradsqlossX = (Phi1tPhi1 * X * Phi2tPhi2 * sumw3sq -  PhitY12) / n;

}else if(b == 2){GradsqlossX = (Phi3tPhi3 * X * sumw12sq - PhitY3) / n;}

}else{GradsqlossX = (RHmat(Phi3tW3Phi3, RHmat(Phi2tW2Phi2, RHmat(Phi1tW1Phi1, X, p2, p3), p3, p1), p1, p2) - PhitWZ) / n;}

Prop = prox_l1(X - delta * GradsqlossX, delta * wGamma);

}

if(RR == 1){

sqlossProp = sqlossrr(Phi1, Phi2, Phi3, Psi, Y, Prop, wfix, Propalpha, n, b, nonten);

}else{

sqlossProp = sqloss(SqrtW1Phi1, SqrtW2Phi2, SqrtW3Phi3, SqrtWPsi, V, SqrtWZ, Prop, Propalpha, n, p2, p3, n1, n2, nonten, S);

}

Betaprevprox = Beta;
Beta = Prop;
alphaprevprox = alpha;
alpha = Propalpha;

sqlossBeta = sqlossProp;

Iter(l, m) = k;
objprox(k) = sqlossBeta + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);
relobjprox = abs(objprox(k) - objprox(k - 1)) / abs(objprox(k - 1));

if(k < maxiterprox && relobjprox < reltolprox){

objprox.fill(NA_REAL);
break;

}else if(k == maxiterprox){

objprox.fill(NA_REAL);
break;

}

}

}//end proximal loop

if(RR == 1){

Phi1 = Phi1tr;//???????????  overwrite like that!!!!
Phi2 = Phi2tr;//???????????  overwrite like that!!!!
Phi3 = Phi3tr;//???????????  overwrite like that!!!!

if(b == 1){

Psi = Psitr;
Y = Ytr;

}else if(b == 2){

Psi = Psirot;
Y = Yrot;

}

}

/////newton line search (backtracking line search in tseng yun 2009 with gamma = 0)

if(S == 1){

Eta = etafunc(Phi1, Phi2, Phi3, V, Beta, n, S);

}else if(RR == 1){

if(b == 1){

Eta = vectorise(Phi1 * Beta * Phi2.t()) * wfix.t();
Eta.reshape(n1, n2 * n3);

}else if (b == 2){Eta =  (Phi3 * Beta) * vectorise(wfix).t();}

}else{Eta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta, p2, p3), p3, n1), n1, n2);}

if(nonten == 1){

Etaalpha = Psi * alpha;

if(b == 2){Etaalpha.reshape(n3, n1 * n2);}else{Etaalpha.reshape(n1, n2 * n3);}

Eta = Eta + Etaalpha;

}

loglikeBeta = loglike(Y, Weights, Eta, n, family);
MuEta = mu(Eta, family);

if(S == 1){

Tmp = Weights % dtheta(Eta, family) % (MuEta - Y);
sumVtmp.fill(0);
for(int i = 0; i < n3; i++){sumVtmp = sumVtmp + V.cols(i * n2, (i + 1) * n2 - 1) % Tmp.cols(i * n2, (i + 1) * n2 - 1);}
GradloglikeBeta = Phi1.t() * sumVtmp * Phi2 / n;

}else if(RR == 1){

Tmp = Weights % dtheta(Eta, family) % (MuEta - Y);

if(b ==1){
  
sumVtmp.reshape(n1, n2);//necessary??
sumVtmp.fill(0);
for(int i = 0; i < n3; i++){sumVtmp = sumVtmp + wfix(i) * Tmp.cols(i * n2, (i + 1) * n2 - 1);}

GradloglikeBeta = Phi1.t() * sumVtmp * Phi2 / n;

}else if(b == 2){
  
sumVtmp.reshape(n3, 1); //necessary??
sumVtmp.fill(0);

for(int j = 0; j < n2; j++){for(int i = 0; i < n1; i++){sumVtmp = sumVtmp + wfix(i, j) * Tmp.col(j * n1 + i);}}

GradloglikeBeta = Phi3.t() * sumVtmp / n;

}

}else{GradloglikeBeta = gradloglike(Y, Weights, Phi1, Phi2, Phi3, MuEta, Eta, n2, n3, p1, p2, n, family);}

if(nonten == 1){ //compute the gradient in alpha for non tensor component

gradloglikealpha = Psi.t() * vectorise(Weights % dtheta(Eta, family) % (MuEta - Y)) / n;
Deltaalpha = alpha - alphaprevnewt;
valalphanewt = accu(gradloglikealpha % Deltaalpha) + l1penalty(wGammaalpha, alpha) - l1penalty(wGammaalpha, alphaprevnewt);

}else{valalphanewt = 0;}

DeltaBeta = Beta - Betaprevnewt;
valnewt = accu(GradloglikeBeta % DeltaBeta) + l1penalty(wGamma, Beta) - l1penalty(wGamma, Betaprevnewt) + valalphanewt;

tnewt = 1;
BTnewt(l, m) = 0;

while (BTnewt(l, m) < btnewtmax) {

DeltaBeta = tnewt * DeltaBeta;
  if(S == 1){

Etatmp = etafunc(Phi1, Phi2, Phi3, V, Beta + DeltaBeta, n, S);

  }else if(RR == 1){

if(b == 1){

Etatmp = vectorise(Phi1 * (Beta + DeltaBeta) * Phi2.t()) * wfix.t();
Etatmp.reshape(n1, n2 * n3);

}else{

Etatmp =  (Phi3 * (Beta + DeltaBeta)) * vectorise(wfix).t();

}

}else{

Etatmp =  RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta + DeltaBeta, p2, p3), p3, n1), n1, n2);

  }

if(nonten == 1){

Deltaalpha = tnewt * Deltaalpha;
mat tmp = Psi * (alpha + Deltaalpha);
if(b == 2){tmp.reshape(n3, n1 *n2);}else{tmp.reshape(n1, n2 *n3);}
Etatmp = Etatmp + tmp;

}

MuEtatmp = mu(Etatmp, family);
logliketmp = loglike(Y, Weights, Etatmp, n, family);

if(logliketmp <= loglikeBeta + constnewt * tnewt * valnewt){

Beta = (1 - tnewt) * Betaprevnewt + tnewt * Beta;
if(nonten == 1){alpha = (1 - tnewt) * alphaprevnewt + tnewt * alpha;}
break;

}else{ 

tnewt = snewt * tnewt;
BTnewt(l, m) = BTnewt(l, m) + 1;

}

}

if(tnewt < 1){//Beta has changed

if(S == 1){

Eta = etafunc(Phi1, Phi2, Phi3, V, Beta, n, S);

}else if(RR == 1){

if(b == 1){

Eta = vectorise(Phi1 * Beta * Phi2.t()) * wfix.t();
Eta.reshape(n1, n2 * n3);

}else{

Eta = (Phi3 * Beta) * vectorise(wfix).t(); //n3xn1n2

}

}else{

Eta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta, p2, p3), p3, n1), n1, n2);

}

if(nonten == 1){

mat tmp = Psi * alpha;
if(b == 2){tmp.reshape(n3, n1 *n2);}else{tmp.reshape(n1, n2 *n3);}
Eta = Eta + tmp;

}

MuEta = mu(Eta, family);
loglikeBeta = loglike(Y, Weights, Eta, n, family);
objnewt(m + 1) = loglikeBeta + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);

}else{objnewt(m + 1) = loglikeBeta + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);}

relobjnewt = abs(objnewt(m + 1) - objnewt(m)) / (reltolnewt + abs(objnewt(m)));
Betaprevnewt = Beta;
alphaprevnewt = alpha;

/////newton convergence check
if(relobjnewt < reltolnewt){//go to next ...

objnewt.fill(NA_REAL);
break;

}else if(m + 1 == maxiternewt){//go to next ...

objnewt.fill(NA_REAL);
break;

}

////check if maximum number of newton backtraking step is reached
if(BTnewt(l, m) >= btnewtmax){STOPnewt = 1;}

//check if maximum number of iterations for current lambda is reached
if(accu(Iter.row(l)) > maxiter){STOPmaxiter = 1;}

//break newton loop if maxiter or btnewtmax is reached
if(STOPmaxiter == 1 || STOPnewt == 1){break;}

} //end newton loop

//stop program if maximum number of backtracking steps or maxiter is reached
if(STOPmaxiter == 1 || STOPnewt == 1){

endmodelno = l;
break;

}

}//end MSA loop

if(RR == 1){if(b == 1){Beta12 = Beta;}else if(b == 2){Beta3 = Beta;}}

}//#blocks
  
if(RR == 1){

Betaa = outermat(Beta12, Beta3);
Eta =  RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Betaa, p2, p3), p3, n1), n1, n2);

if(nonten == 1){

mat tmp =   Psitr * alpha;
tmp.reshape(n1, n2 * n3);
Eta = Eta + tmp;

}

objalt(a) = loglike(Ytr, Weightstr, Eta, n, family) + penfunc(penaltyfactortr * lambda(l), Betaa, penalty) + penfunc(penaltyfactoralpha * lambda(l), alpha, penalty);

if(a > 0 && abs(objalt(a) - objalt(a - 1)) /  abs(objalt(a)) < reltolalt){break;}

}

  
}//#altenations

if(RR == 1){

df(l) = p1 * p2 + p3 + nonten *  q - accu((Beta12 == 0)) - accu((Beta3 == 0)) - nonten * accu((alpha == 0));
Betas12.col(l) = vectorise(Beta12);
Betas3.col(l) = vectorise(Beta3);

}else{
 
df(l) = p + nonten * q - accu((Beta == 0)) - nonten * accu((alpha == 0));
Betas.col(l) = vectorise(Beta);

}

if(nonten == 1){alphas.col(l) = alpha;}


} //end lambda loop

btiternewt = accu((BTnewt > 0) % BTnewt);

output = Rcpp::List::create(Rcpp::Named("Beta") = Betas, 
                            Rcpp::Named("Beta12") = Betas12, 
                            Rcpp::Named("Beta3") = Betas3,
Rcpp::Named("alpha") = alphas,
Rcpp::Named("btiternewt") = btiternewt,
Rcpp::Named("btiterprox") = btiterprox,
Rcpp::Named("endmodelno") = endmodelno,
Rcpp::Named("Iter") = Iter,
Rcpp::Named("lambda") = lambda,
Rcpp::Named("STOPmaxiter") = STOPmaxiter,
Rcpp::Named("STOPnewt") = STOPnewt,
Rcpp::Named("STOPprox") = STOPprox);

}else if(iwls == "exact"){//solve wls suproblems/////////////////////////////////

////declare variables
int ascentprox, ascentproxmax,
     btprox;

double constnewt,
 delta, deltamin,
 L, lmax, loglikeBeta, logliketmp, 
 relobjnewt, relobjprox, 
 snewt, sprox, 
 tnewt,  
 valalpha, valalphanewt, valnewt, valprox, 
 wsqlossBeta, wsqlossProp, wsqlossX;  
 
arma::vec absalpha(q), alpha(q), alphaprevprox(q), alphaprevnewt(q), 
Deltaalpha(q), df(nlambda), dpenalpha(q),
eig1, eig2, eig3, 
Gammaalpha(q), GradwsqlossXalpha(q),  gradloglikealpha(q),
mwtrue(n3), 
objnewt(maxiternewt + 1), objprox(maxiterprox + 1), objalt(maxalt),
pospartalpha(q), Propalpha(q), Psitwz(q),
vhat1(n1), vhat2(n2), vhat3(n3),
w(n), wGammaalpha(q),
Xalpha(q);
  
arma::mat absBeta(p1, p2 * p3), alphas(q, nlambda),
Beta(p1, p2 * p3),  Betaa, Beta12, Beta3, Betaprevnewt(p1, p2 * p3), Betaprevprox(p1, p2 * p3), Betaprevprox12, Betaprevprox3 ,  Betas(p, nlambda), Betas12(p1 * p2, nlambda), Betas3(p3, nlambda), BTnewt(nlambda, maxiternewt + 1), 
DeltaBeta(p1, p2 * p3),  dpen(p1, p2 * p3), 
Eta(n1, n2 * n3), Etatmp(n1, n2 * n3), Etaalpha(n, 1),
fix,
Gamma(p1, p2 * p3), GradloglikeBeta(p1, p2 * p3), GradwsqlossX(p1, p2 * p3), 
Iter(nlambda, maxiternewt), 
MuEta(n1, n2 * n3), MuEtatmp(n1, n2 * n3), 
Phi1tPhi1, Phi2tPhi2, Phi3tPhi3,  PhitWZ, pospart(p1, p2 * p3), Prop(p1, p2 * p3), PsitPsi(q, q), PsitWPsi(q, q), PsiXalpha,
SqrtW, SqrtWV, SqrtWZ, SqrtWPsi(n, q), sumVTmp(n1, n2), sumVtmp, sumWVZ(n1, n2),
tmp, Tmp,
U, what,
W(n1, n2 * n3), wfix,  w12(n1, n2), w3(n3, 1), wGamma(p1, p2 * p3), WV, sumWVsq(n1, n2), WZ,
X(p1, p2 * p3), 
Z(n1, n2 * n3);
  
arma::cube BTprox(nlambda, maxiternewt, maxiterprox + 1),
           Lambdas(nlambda, maxalt, blocks - 1);//rr


////fill variables
ascentproxmax = 4;

constnewt = 0.2;
snewt = 0.5; 
sprox= 0.9; 
 
objnewt.fill(NA_REAL);
objprox.fill(NA_REAL);

Betas.fill(NA_REAL);
BTnewt.fill(-1); 
Iter.fill(0);

BTprox.fill(-1); //negative vals?

////precompute
Phi1tPhi1 = Phi1.t() * Phi1;
Phi2tPhi2 = Phi2.t() * Phi2;
Phi3tPhi3 = Phi3.t() * Phi3;

//proximal step size
eig1 = arma::eig_sym(Phi1tPhi1);
eig2 = arma::eig_sym(Phi2tPhi2);
eig3 = arma::eig_sym(Phi3tPhi3);

////initialize
Beta = Betainit; 
Betaprevprox = Beta; 
Betaprevnewt = Beta;
alpha = alphainit; 
alphaprevprox = alpha; 
alphaprevnewt = alpha;

if(RR == 1){//rr
  
Beta12 = Beta12init;
Beta3 = Beta3init;
Betaa = outermat(Beta12, Beta3);
  
}


if(S == 1){

Eta = etafunc(Phi1, Phi2, Phi3, V, Beta, n, S);

}else if(RR != 1){

Eta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta, p2, p3), p3, n1), n1, n2);

}

if(nonten == 1){

Etaalpha = Psi * alpha;
Etaalpha.reshape(n1, n2 * n3);
Eta = Eta + Etaalpha;

}

loglikeBeta = loglike(Y, Weights, Eta, n, family);
MuEta = mu(Eta, family);

////lambda loop
for (int l = 0; l < nlambda; l++){
  
Gamma = penaltyfactor * lambda(l);
if(nonten == 1){Gammaalpha = penaltyfactoralpha * lambda(l);}
  
for(int a = 0; a < maxalt; a++){//#alternation loop over blocks

for(int b = 1; b < 3 ; b++){//#block loop
 
if(RR == 1){////compute Eta
  
if(b == 1){//######## block 1 //overwrite!!

Weights = Weightstr;
Y = Ytr;
Psi = Psitr;
Beta = Beta12;
Betaprevnewt = Beta;
fix = Beta3;
wfix = Phi3 * fix; // #n3x1
penaltyfactor = penaltyfactor12;
lmax = lambmaxrr(Y, Phi1, Phi2, Phi3, Psi,  Weights,  wfix,  n,   b, family);
if(lmax < lambda(l)){Lambdas(l, a, b - 1) = lmax * 0.9;}else{Lambdas(l, a, b - 1) = lambda(l);}
Gamma =  penaltyfactor12 * Lambdas(l, a, b - 1);

V.reshape(n1, n2 * n3);
for (int i = 0; i < n3; i++){V.cols(i * n2, (i + 1) * n2 - 1) = w12.fill(wfix(i));}
Eta = vectorise(Phi1 * Beta * Phi2.t()) * wfix.t();
Eta.reshape(n1, n2 * n3);
   
if(nonten== 1 ){  
 
Etaalpha = Psi * alpha;
Etaalpha.reshape(n1, n2 * n3);
Eta = Eta + Etaalpha;
 
}

}else if(b == 2){//######## block 2 //overwrite!!

Weights = Weightsrot;
Y = Yrot;
Psi = Psirot;
Beta = Beta3;
Betaprevnewt = Beta;
fix = Beta12;
wfix = Phi1 * fix * Phi2.t(); //#n1xn2
penaltyfactor = penaltyfactor3;
lmax = lambmaxrr(Y, Phi1, Phi2, Phi3, Psi,  Weightsrot,  wfix,  n,   b, family);
if(lmax < lambda(l)){Lambdas(l, a, b - 1) = lmax * 0.9;}else{Lambdas(l, a, b - 1) = lambda(l);}
Gamma = penaltyfactor3 * Lambdas(l, a, b - 1);

V.reshape(n3, n1 * n2); //V is rotated!!
for (int i = 0; i < n1; i++){for (int j = 0; j < n2; j++){V.col(i + j * n1) = w3.fill(wfix(i, j));}} ///?????index!!!!!
Eta =  (Phi3 * Beta) * vectorise(wfix).t();

if(nonten== 1 ){  

Etaalpha = Psi * alpha;
Etaalpha.reshape(n3, n1 * n2);
Eta = Eta + Etaalpha;

}

}
  
}

if(nonten== 1 ){   PsitPsi = Psi.t() * Psi;}else{PsitPsi.fill(0);}
  
loglikeBeta = loglike(Y, Weights, Eta, n, family);
MuEta = mu(Eta, family);
  
/////start MSA loop
for (int s = 0; s < steps; s++){

if(s == 0){//first MSA step

if(penalty != "lasso"){
  
wGamma =  Gamma / lambda(l);
wGammaalpha =  Gammaalpha / lambda(l);
  
}else{
  
wGamma = Gamma;
wGammaalpha = Gammaalpha;
  
}

}else{ // not first MSA step

if(penalty == "scad"){
  
absBeta =  abs(Beta);
pospart = ((ascad * Gamma - absBeta) + (ascad * Gamma - absBeta)) / 2;
dpen = sign(Beta) % (Gamma % (absBeta <= Gamma) + pospart / (ascad - 1) % (absBeta > Gamma));
wGamma = abs(dpen) % penaltyfactor % (Beta != 0) + Gamma % (Beta == 0);
  
absalpha =  abs(alpha);
pospartalpha = ((ascad * Gammaalpha - absalpha) + (ascad * Gammaalpha - absalpha)) / 2;
dpenalpha = sign(alpha) % (Gammaalpha % (absalpha <= Gammaalpha) + pospartalpha / (ascad - 1) % (absalpha > Gammaalpha));
wGammaalpha = abs(dpenalpha) % penaltyfactoralpha % (alpha != 0) + Gammaalpha % (alpha == 0);
  
}

}

objnewt(0) = loglikeBeta + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);

/////start newton loop
for (int m = 0; m < maxiternewt; m++){

////iwls weights
W = Weights % dmu(Eta, family) % dtheta(Eta, family); //a * (mu’)^2 * theta’/mu’
////working responses
Z = (Y - MuEta) % dg(MuEta, family) + Eta;

////precompute
SqrtW = sqrt(W);
SqrtWZ = SqrtW % Z;
WZ = W % Z;

if(S == 1){//weighted weights for gradient
  
SqrtWV = SqrtW % V;
WV = W % V;
  
double sumij;
for(int i = 0; i < n1; i++) {
for(int j = 0; j < n2; j++) {
sumij = 0;
for(int itr = 0; itr < n3; itr++) {sumij = sumij + pow(SqrtWV(i, j + itr * n2), 2);}

sumWVsq(i, j) = sumij;

}
} 
  
sumWVZ.fill(0);
for(int i = 0; i < n3; i++){sumWVZ = sumWVZ + WV(span::all, span(i * n2, (i + 1) * n2 - 1)) % Z(span::all, span(i * n2, (i + 1) * n2 - 1));}
PhitWZ = Phi1.t() * sumWVZ * Phi2;
  
L = as_scalar(max(eig1) * max(eig2) * max(max(sumWVsq)) + nonten * max(arma::eig_sym(PsitPsi)) * max(max(W))) / n;
  
}else if(RR == 1){
  
if(b == 1){//######## block 1 //overwrite!!

SqrtWV = SqrtW % V;
WV = W % V; 
sumWVsq.reshape(n1, n2);

double sumij;
for(int i = 0; i < n1; i++) {
for(int j = 0; j < n2; j++) {
sumij = 0;

for(int itr = 0; itr < n3; itr++) {sumij = sumij + pow(SqrtWV(i, j + itr * n2), 2);}

sumWVsq(i, j) = sumij;

}
}

sumWVZ.reshape(n1, n2);
sumWVZ.fill(0);
for(int i = 0; i < n3; i++){sumWVZ = sumWVZ + WV.cols(i * n2, (i + 1) * n2 - 1) % Z.cols(i * n2, (i + 1) * n2 - 1);}
PhitWZ = Phi1.t() * sumWVZ * Phi2;

L = as_scalar(max(eig1) * max(eig2) * max(max(sumWVsq)) + nonten * max(arma::eig_sym(PsitPsi)) * max(max(W))) / n;

}else if(b == 2){//######## block 2 //overwrite!!

SqrtWV = SqrtW % V; //n3xn1n2
WV = W % V; 

sumWVsq.reshape(n3, 1);
double sumi;
for(int i = 0; i < n3; i++){

sumi = 0;
for(int j = 0; j < n1; j++){
  
  for(int itr = 0; itr < n2; itr++) {sumi = sumi + pow(SqrtWV(i, j + itr * n1), 2);}}
  
sumWVsq(i) = sumi;
  
//}???????

}

sumWVZ.reshape(n3, 1);
sumWVZ.fill(0);
for(int i = 0; i < n1; i++){for(int j = 0; j < n2; j++){sumWVZ = sumWVZ + WV.col(i + j * n1) % Z.col(i + j * n1);}}

PhitWZ = Phi3.t() * sumWVZ ;

L = as_scalar(max(eig3) * max(max(sumWVsq)) + nonten * max(arma::eig_sym(PsitPsi)) * max(max(W))) / n;
 
}
  
}else{
  
PhitWZ = RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), WZ, n2, n3), n3, p1), p1, p2);
L = (max(eig1) * max(eig2) * max(eig3) + nonten * max(arma::eig_sym(PsitPsi))) * max(max(W)) / n;
  
}

if(nonten == 1){Psitwz = Psi.t() * vectorise(WZ);}

////initial step size
deltamin = 1.99 / L; //minimum stepsize
if(nu > 0){delta = 1.9 / (nu * L);}else{delta = 1;}

ascentprox = 0;
  
/////start proximal loop
for (int k = 0; k < maxiterprox; k++){
  
if(k == 0){

Betaprevprox = Beta;
alphaprevprox = alpha;
  
if(RR == 1){

objprox(0) =  wsqlossrr(Phi1, Phi2, Phi3, Psi, SqrtWZ, Beta, SqrtWV, alpha, n, b, nonten) + l1penalty(wGamma, Beta) + nonten * l1penalty(wGammaalpha, alpha);

}else{
  
objprox(0) = wsqloss(SqrtW, SqrtWV, Phi1, Phi2, Phi3, Psi, SqrtWZ, Beta, alpha, n, p2, p3, n1, n2, nonten, S) + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);
  
}

if(nu > 0 && nu < 1){BTprox(l, m, k) = 1;} //force initial backtracking for deltamin < delta

}else{

X  = Beta + (k - 2) / (k + 1) * (Beta - Betaprevprox);
Xalpha = alpha + (k - 2) / (k + 1) * (alpha - alphaprevprox);
  
////gradient
if(nonten == 1){

PsiXalpha = Psi * Xalpha;
PsiXalpha.reshape(n1, n2 * n3);

if(S == 1){

Tmp =  W % (etafunc(Phi1, Phi2, Phi3, V, X, n, S) + PsiXalpha);//remove W% by precomputing
sumVTmp.fill(0);
for(int i = 0; i < n3; i++){sumVTmp = sumVTmp + V(span::all, span(i * n2, (i + 1) * n2 - 1)) % Tmp(span::all, span(i * n2, (i + 1) * n2 - 1));}

GradwsqlossX =  (Phi1.t() * sumVTmp * Phi2 - PhitWZ) / n;
GradwsqlossXalpha = (Psi.t() * vectorise(Tmp) - Psitwz) / n;

}else if(RR == 1){

if(b == 1){
  
PsiXalpha.reshape(n1, n2 * n3);
sumVtmp.reshape(n1, n2);
tmp =  W % (etafunc(Phi1, Phi2, Phi3, V, X, n, b) + PsiXalpha);//remove W% by precomputing W%Psi and using V=sqrtWV
sumVtmp.fill(0);///???????dims!!!
for(int i = 0; i < n3; i++){sumVtmp = sumVtmp + V.cols(i * n2, (i + 1) * n2 - 1) % tmp.cols(i * n2, (i + 1) * n2 - 1);}
GradwsqlossX =  (Phi1.t() * sumVtmp * Phi2 - PhitWZ) / n;
GradwsqlossXalpha = (Psi.t() * vectorise(tmp) - Psitwz) / n;
  
}else if(b == 2){//b=2,rotated data/model
  
PsiXalpha.reshape(n3, n1 * n2);
sumVtmp.reshape(n3, 1);
tmp =  W % (etafunc(Phi1, Phi2, Phi3, V, X, n, b) + PsiXalpha);//remove W% by precomputing W%Psi and using V=sqrtWV
sumVtmp.fill(0);///???????dims!!!
for(int j = 0; j < n2; j++){for(int i = 0; i < n1; i++){sumVtmp = sumVtmp + V.col(i + j * n1) % tmp.col(i + j * n1);}}
GradwsqlossX =  (Phi3.t() * sumVtmp - PhitWZ) / n;
GradwsqlossXalpha = (Psi.t() * vectorise(tmp) - Psitwz) / n;
  
}

}else{//glamlasso

Tmp =  W % (RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, X, p2, p3), p3, n1), n1, n2) + PsiXalpha);
GradwsqlossX = (RHmat(Phi3.t(), RHmat(Phi2.t(), RHmat(Phi1.t(), Tmp, n2, n3), n3, p1), p1, p2) - PhitWZ) / n;
GradwsqlossXalpha = (Psi.t() * vectorise(Tmp) - Psitwz) / n;

}

}else{//no nonten

if(S == 1){

GradwsqlossX = (Phi1.t() * (sumWVsq % (Phi1 * X * Phi2.t())) * Phi2  - PhitWZ) / n;

}else if(RR == 1){//rr

if(b == 1){
  
GradwsqlossX = (Phi1.t() * (sumWVsq % (Phi1 * X * Phi2.t())) * Phi2 - PhitWZ) / n;
  
}else{//rotated data/model!!!b=2
  
GradwsqlossX = (Phi3.t() * (sumWVsq % (Phi3 * X)) - PhitWZ) / n;
  
}

}else{

GradwsqlossX = (winprod(W, Phi1, Phi2, Phi3, X, n1, n2, n3, p1, p2, p3) - PhitWZ) / n;

}

}
  
////check if proximal backtracking occurred last iteration
if(BTprox(l, m, k - 1) > 0){btprox = 1;}else{btprox = 0;}

////check for divergence
if(ascentprox > ascentproxmax){btprox = 1;}

if((btprox == 1 && deltamin < delta) || nu == 0){//backtrack

  
if(RR == 1){

wsqlossX = wsqlossrr(Phi1, Phi2, Phi3, Psi, SqrtWZ, X, SqrtWV, alpha, n, b, nonten);

}else{

wsqlossX = wsqloss(SqrtW, SqrtWV, Phi1, Phi2, Phi3, Psi, SqrtWZ, X, Xalpha, n, p2, p3, n1, n2, nonten, S);

}
  
BTprox(l, m, k) = 0;

while(BTprox(l, m, k) < btproxmax){////proximal line search

Prop = prox_l1(X - delta * GradwsqlossX, delta * wGamma); 

if(nonten == 1){

Propalpha = prox_l1(Xalpha - delta * GradwsqlossXalpha, delta * wGammaalpha);
valalpha =  accu(GradwsqlossXalpha % (Propalpha - Xalpha)) + 1 / (2 * delta) * sum_square(Propalpha - Xalpha);
  
}else{valalpha = 0;}
  
if(RR == 1){

wsqlossProp = wsqlossrr(Phi1, Phi2, Phi3, Psi, SqrtWZ, Prop, SqrtWV, alpha, n, b, nonten);

}else{
    
wsqlossProp = wsqloss(SqrtW, SqrtWV, Phi1, Phi2, Phi3, Psi, SqrtWZ, Prop, Propalpha, n, p2, p3, n1, n2, nonten, S);
  
}
  
valprox = as_scalar(wsqlossX + accu(GradwsqlossX % (Prop - X)) + 1 / (2 * delta) * sum_square(Prop - X)) + valalpha;

if (wsqlossProp <= valprox + 0.0000001){ //need to add a little due to numerical issues

break;

}else{

delta = sprox * delta;
BTprox(l, m, k) = BTprox(l, m, k) + 1;

if(delta < deltamin){delta = deltamin;}

}

}//end backtracking

////check if maximum number of proximal backtraking step is reached
if(BTprox(l, m, k) == btproxmax){STOPprox = 1;}

}else{//no backtracking

Prop = prox_l1(X - delta * GradwsqlossX, delta * wGamma);
if(nonten == 1){Propalpha = prox_l1(Xalpha - delta * GradwsqlossXalpha, delta * wGammaalpha);}

if(RR == 1){
  
wsqlossProp =  wsqlossrr(Phi1, Phi2, Phi3, Psi, SqrtWZ, Prop, SqrtWV, alpha, n, b, nonten);
  
}else{
  
wsqlossProp = wsqloss(SqrtW, SqrtWV, Phi1, Phi2, Phi3, Psi, SqrtWZ, Prop, Propalpha, n, p2, p3, n1, n2, nonten, S);

}

}

Betaprevprox = Beta;
alphaprevprox = alpha;
Beta = Prop;
alpha = Propalpha;
wsqlossBeta = wsqlossProp;
objprox(k) = wsqlossBeta + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha); 

////check if objective has increased
if(objprox(k) > objprox(k - 1)){ascentprox = ascentprox + 1;}else{ascentprox = 0;}

Iter(l, m) = k;

////proximal convergence check
relobjprox = abs(objprox(k) - objprox(k - 1)) / (reltolprox + abs(objprox(k - 1))); 

if(k < maxiterprox && relobjprox < reltolprox){

objprox.fill(NA_REAL);
break;

}else if(k == maxiterprox){

objprox.fill(NA_REAL);
break;

}

}

////break proximal loop if maximum number of proximal backtraking step is reached
if(STOPprox == 1){break;}

////break proximal loop if maximum number of iterations is reached
if(accu(Iter.row(l)) > maxiter){
  
STOPmaxiter = 1;
break;
  
}

} //end proximal loop

if(S == 1){

Eta = etafunc(Phi1, Phi2, Phi3, V, Beta, n, S);

}else if(RR == 1){

if(b == 1){

Eta = vectorise(Phi1 * Beta * Phi2.t()) * wfix.t();
Eta.reshape(n1, n2 * n3);
  
}else if(b == 2){Eta =  (Phi3 * Beta) * vectorise(wfix).t();}

}else{Eta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta, p2, p3), p3, n1), n1, n2);}

if(nonten == 1){

Etaalpha = Psi * alpha;
if(b == 2){Etaalpha.reshape(n3, n1 * n2);}else{Etaalpha.reshape(n1, n2 * n3);}
Eta = Eta + Etaalpha;

}

loglikeBeta = loglike(Y, Weights, Eta, n, family);
MuEta = mu(Eta, family);
if(S == 1){

Tmp = Weights % dtheta(Eta, family) % (MuEta - Y);
sumVTmp.fill(0);
for(int i = 0; i < n3; i++){sumVTmp = sumVTmp + V.cols(i * n2, (i + 1) * n2 - 1) % Tmp.cols(i * n2, (i + 1) * n2 - 1);}

GradloglikeBeta = Phi1.t() * sumVTmp * Phi2 / n;

}else if(RR == 1){
Tmp = Weights % dtheta(Eta, family) % (MuEta - Y);

if(b ==1){

sumVtmp.reshape(n1, n2);//necessary??
sumVtmp.fill(0);

for(int i = 0; i < n3; i++){sumVtmp = sumVtmp + wfix(i) * Tmp.cols(i * n2, (i + 1) * n2 - 1);}

GradloglikeBeta = Phi1.t() * sumVtmp * Phi2 / n;

}else if(b == 2){
    
sumVtmp.reshape(n3, 1); //necessary??

sumVtmp.fill(0);

for(int j = 0; j < n2; j++){for(int i = 0; i < n1; i++){sumVtmp = sumVtmp + wfix(i, j) * Tmp.col(j * n1 + i); }}

GradloglikeBeta = Phi3.t() * sumVtmp / n;

}

}else{GradloglikeBeta = gradloglike(Y, Weights, Phi1, Phi2, Phi3, MuEta, Eta, n2, n3, p1, p2, n, family);}

if(nonten == 1){ //compute the gradient in alpha for non tensor component

gradloglikealpha = Psi.t() * vectorise(Weights % dtheta(Eta, family) % (MuEta - Y)) / n;
Deltaalpha = alpha - alphaprevnewt;
valalphanewt = accu(gradloglikealpha % Deltaalpha) + l1penalty(wGammaalpha, alpha) - l1penalty(wGammaalpha, alphaprevnewt);

}else{valalphanewt = 0;}

DeltaBeta = Beta - Betaprevnewt;
valnewt = accu(GradloglikeBeta % DeltaBeta) + l1penalty(wGamma, Beta) - l1penalty(wGamma, Betaprevnewt) + valalphanewt;

tnewt = 1;
BTnewt(l, m) = 0;

while(BTnewt(l, m) < btnewtmax){/////newton line search (backtracking line search in tseng yun 2009 with gamma=0)

DeltaBeta = tnewt * DeltaBeta;
if(S == 1){

Etatmp = etafunc(Phi1, Phi2, Phi3, V, Beta + DeltaBeta, n, S);

}else if(RR == 1){

if(b == 1){

Etatmp = vectorise(Phi1 * (Beta + DeltaBeta) * Phi2.t()) * wfix.t();
Etatmp.reshape(n1, n2 * n3);

}else{

Etatmp =  (Phi3 * (Beta + DeltaBeta)) * vectorise(wfix).t();

}

}else{

Etatmp =  RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta + DeltaBeta, p2, p3), p3, n1), n1, n2);

}

if(nonten == 1){

Deltaalpha = tnewt * Deltaalpha;
mat tmp = Psi * (alpha + Deltaalpha);
if(b == 2){tmp.reshape(n3, n1 *n2);}else{tmp.reshape(n1, n2 *n3);}
Etatmp = Etatmp + tmp;

}

MuEtatmp = mu(Etatmp, family);
logliketmp = loglike(Y, Weights, Etatmp, n, family);

if(logliketmp <= loglikeBeta + constnewt * tnewt * valnewt){

Beta = (1 - tnewt) * Betaprevnewt + tnewt * Beta;
if(nonten == 1){alpha = (1 - tnewt) * alphaprevnewt + tnewt * alpha;}

break;

}else{

tnewt = snewt * tnewt;
BTnewt(l, m) = BTnewt(l, m) + 1;

}

}

if(tnewt < 1){//Beta has changed from backtracking

if(S == 1){

Eta = etafunc(Phi1, Phi2, Phi3, V, Beta, n, S);

}else if(RR == 1){
 
if(b == 1){

Eta = vectorise(Phi1 * Beta * Phi2.t()) * wfix.t();
Eta.reshape(n1, n2 * n3);

}else if(b == 2){Eta = (Phi3 * Beta) * vectorise(wfix).t();}

}else{Eta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta, p2, p3), p3, n1), n1, n2);}

if(nonten == 1){

mat tmp = Psi * alpha;
if(b == 2){tmp.reshape(n3, n1 *n2);}else{tmp.reshape(n1, n2 *n3);}

Eta = Eta + tmp;

}

MuEta = mu(Eta, family);
loglikeBeta = loglike(Y, Weights, Eta, n, family);
objnewt(m + 1) = loglikeBeta + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);

}else{objnewt(m + 1) = loglikeBeta + l1penalty(wGamma, Beta) + l1penalty(wGammaalpha, alpha);}

relobjnewt = abs(objnewt(m + 1) - objnewt(m)) / (reltolnewt + abs(objnewt(m)));
Betaprevnewt = Beta;
alphaprevnewt = alpha;

/////newton convergence check
if(relobjnewt < reltolnewt){//go to next lambda

objnewt.fill(NA_REAL);
break;

}else if(m + 1 == maxiternewt){//go to next lambda

objnewt.fill(NA_REAL);
break;

}

////check if maximum number of newton backtraking step is reached
if(BTnewt(l, m) >= btnewtmax){STOPnewt = 1;}

////break newton loop if maximum number of backtracking steps or maxiter is reached
if(STOPprox == 1 || STOPmaxiter == 1 || STOPnewt == 1){break;}

} //end newton loop

//stop program if maximum number of backtracking steps or maxiter is reached
if(STOPnewt == 1 || STOPprox == 1 || STOPmaxiter == 1){
  
endmodelno = l;
break;
  
}
  
//convergence check?? or just run steps times.......?????????????
  
} //end MSA loop
  
if(RR == 1){ 

if(b == 1){Beta12 = Beta;}else if(b == 2){Beta3 = Beta;}

}

}//#blocks

if(RR == 1){

Betaa = outermat(Beta12, Beta3);
Eta =  RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Betaa, p2, p3), p3, n1), n1, n2);

if(nonten == 1){
  
mat tmp = Psitr * alpha;
tmp.reshape(n1, n2 * n3);
Eta = Eta + tmp;
  
}

objalt(a) = loglike(Ytr, Weightstr, Eta, n, family) + penfunc(penaltyfactortr * lambda(l), Betaa, penalty) + penfunc(penaltyfactoralpha * lambda(l), alpha, penalty);
  
if(a > 0 && abs(objalt(a) - objalt(a - 1)) /  abs(objalt(a)) < reltolalt){break;}
  
}

}//#altenations
  
if(RR == 1){

df(l) = p1 * p2 + p3 + nonten *  q - accu((Beta12 == 0)) - accu((Beta3 == 0)) - nonten * accu((alpha == 0));
Betas12.col(l) = vectorise(Beta12);
Betas3.col(l) = vectorise(Beta3);

}else{

df(l) = p + nonten * q - accu((Beta == 0)) - nonten * accu((alpha == 0));
Betas.col(l) = vectorise(Beta);

}
  
if(nonten == 1){alphas.col(l) = alpha;}
  
} //end lambda loop

btiternewt = accu((BTnewt > 0) % BTnewt);
btenterprox = accu((BTprox > -1));
btiterprox = accu((BTprox > 0) % BTprox);

output = Rcpp::List::create(Rcpp::Named("Beta") = Betas,
                            Rcpp::Named("Beta12") = Betas12,
                            Rcpp::Named("Beta3") = Betas3,
                            Rcpp::Named("alpha") = alphas,
                            Rcpp::Named("btenterprox") = btenterprox,
                            Rcpp::Named("btiternewt") = btiternewt,
                            Rcpp::Named("btiterprox") = btiterprox,
                            Rcpp::Named("df") = df,
                            Rcpp::Named("endmodelno") = endmodelno,
                            Rcpp::Named("Iter") = Iter,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("STOPmaxiter") = STOPmaxiter,
                            Rcpp::Named("STOPnewt") = STOPnewt,
                            Rcpp::Named("STOPprox") = STOPprox);  

}

}

return output;

} //end function

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// objective values  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
//[[Rcpp::export]]
Rcpp::List getobj(arma::mat Y, arma::mat Weights,
  arma::mat Phi1, arma::mat Phi2, arma::mat Phi3,
  Rcpp::NumericVector beta,
  arma::vec lambda,
  arma::mat penaltyfactor,
  std::string family,
  std::string penalty){

Rcpp::NumericVector vecbeta(beta);
Rcpp::IntegerVector BetaDim = vecbeta.attr("dim");
arma::cube Beta(vecbeta.begin(), BetaDim[0], BetaDim[1], BetaDim[2], false);

int p2 = Phi2.n_cols;
int p3 = Phi3.n_cols;
int n1 = Phi1.n_rows;
int n2 = Phi2.n_rows;
int n3 = Phi3.n_rows;
int n = n1 * n2 * n3;
int nlambda = lambda.n_elem;

arma::mat Eta, MuEta;

arma::vec Obj(nlambda), Loss(nlambda), Pen(nlambda);

for (int j = 0; j < nlambda; j++){

if(penalty == "lasso"){

Pen(j) = l1penalty(penaltyfactor * lambda(j), Beta.slice(j));

}

if(penalty == "scad"){

double ascad = 3.7;

Pen(j) = scadpenalty(penaltyfactor * lambda(j), ascad, Beta.slice(j));

}

Eta = RHmat(Phi3, RHmat(Phi2, RHmat(Phi1, Beta.slice(j), p2, p3), p3, n1), n1, n2);
MuEta = mu(Eta, family);
Loss(j) = loglike(Y, Weights, Eta, n, family);
Obj(j) = Loss(j) + Pen(j);

}

Rcpp::List output = Rcpp::List::create(Rcpp::Named("Obj") = Obj,
                                       Rcpp::Named("Loss") = Loss,
                                       Rcpp::Named("Pen") = Pen);
return output;

}
