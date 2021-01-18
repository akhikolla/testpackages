#include<RcppArmadillo.h>

using namespace Rcpp;

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP nselprcpp(SEXP y, SEXP x, SEXP sam, SEXP frho, SEXP yfrac) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix samr(sam);
Rcpp::NumericVector frhor(frho);
Rcpp::NumericVector yfracr(yfrac);

int B = samr.nrow(), p = xr.ncol(), n = xr.nrow();

arma::mat Y(yr.begin(), B, n, false);
arma::mat X(xr.begin(), n, p, false);
arma::mat SAM(samr.begin(), B, p, false);
arma::vec a(frhor.begin(), B, false);
arma::vec b(yfracr.begin(), B, false);

arma::mat SAMt = SAM.t();

arma::mat suff = X.t()*Y.t();
arma::mat out = arma::zeros(p,B);
double temp;
double LL;

for (int i=0; i<B; i++) {

LL = 0;
for (int j=0; j<B; j++) {
temp = dot(suff.unsafe_col(i),SAMt.unsafe_col(j));
temp -= b.at(i);
temp -= a.at(j);
temp = exp(temp);
LL += temp;
out.unsafe_col(i) += temp*SAMt.unsafe_col(j);}
out.unsafe_col(i) /= LL;}

return as<NumericMatrix>(wrap(out.t()));
}



// NEW FUNCTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP PRLAP2cpp(SEXP y, SEXP sam, SEXP X, SEXP pm, SEXP pv) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix samr(sam);
Rcpp::NumericMatrix Xr(X);
Rcpp::NumericVector pmr(pm);
Rcpp::NumericVector pvr(pv);

int n = Xr.nrow();
int k = Xr.ncol();
int B = yr.nrow();

arma::mat Y(yr.begin(), B, n, false);
arma::mat BETA(samr.begin(), B, k, false);
arma::mat x(Xr.begin(), n, k, false);
arma::vec PM(pmr.begin(), k, false);
arma::vec PV(pvr.begin(), k, false);

arma::mat BETAT = BETA.t();

double eps;
eps = 0.0001;

double scale;
scale = 0.25;

arma::mat xt = x.t();
arma::mat Yt = Y.t();

arma::vec iPV = 1/PV;
arma::vec miPV = -1/PV;
arma::mat invV = arma::zeros(k,k);
for (int j=0; j<k; j++) {
invV(j,j) = iPV(j);}

arma::vec beta(k);
arma::vec eta(n);
arma::vec mu(n);
arma::mat XWX(k,k);

arma::vec temp1(n);
arma::vec temp2(k);
arma::mat Sigma(k,k);
double crit;
int counter;
double dett;

double cons;
cons =  sum(log(PV));
cons *= 0.5;
arma::vec out(B);
out.fill(cons);

for (int i=0; i<B; i++) {

beta = PM;
crit = eps+1;
counter = 0;
while ((crit>eps)&(counter<50)) {

eta = x*beta;
mu = exp(eta);
XWX.zeros();
for (int ii=0; ii<k; ii++) {
for (int jj=ii; jj<k; jj++) {
for (int kk=0; kk<n; kk++) {
XWX.at(ii,jj) += mu.at(kk)*x.at(kk,ii)*x.at(kk,jj);}
XWX.at(jj,ii) = XWX.at(ii,jj);}}

for (int ii=0; ii<k; ii++) {
XWX.at(ii,ii) += iPV.at(ii);}

XWX = inv_sympd(XWX);

temp2 = beta;
temp2 -= PM;
temp2 %= miPV;

temp1 = Yt.unsafe_col(i);
temp1 -= mu;
temp2 += xt*temp1;

temp2 = XWX*temp2;

beta += scale*temp2;
counter += 1;
crit = dot(temp2,temp2);}

eta = x*BETAT.col(i);
out(i) += dot(Y.row(i),eta);
eta = exp(eta);
out(i) -= sum(eta);

eta = x*beta;
out(i) -= dot(Y.row(i),eta);
eta = exp(eta);
out(i) += sum(eta);

beta -= PM;
temp2 = beta;
beta %= temp2;
beta %= iPV;
out(i) += 0.5*sum(beta);
log_det(dett,crit,XWX);
out(i) -= 0.5*dett;

}

return as<NumericVector>(wrap(out));

}

// NEW FUNCTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP PRNSELLAP2cpp(SEXP y, SEXP sam, SEXP X, SEXP pm, SEXP pv) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix samr(sam);
Rcpp::NumericMatrix Xr(X);
Rcpp::NumericVector pmr(pm);
Rcpp::NumericVector pvr(pv);

int n = Xr.nrow();
int k = Xr.ncol();
int B = yr.nrow();

arma::mat Y(yr.begin(), B, n, false);
arma::mat BETA(samr.begin(), B, k, false);
arma::mat x(Xr.begin(), n, k, false);
arma::vec PM(pmr.begin(), k, false);
arma::vec PV(pvr.begin(), k, false);

arma::mat BETAT = BETA.t();

double eps;
eps = 0.0001;

double scale;
scale = 0.25;

arma::mat xt = x.t();
arma::mat Yt = Y.t();

arma::vec iPV = 1/PV;
arma::vec miPV = -1/PV;
arma::mat invV = arma::zeros(k,k);
for (int j=0; j<k; j++) {
invV(j,j) = iPV(j);}

arma::vec beta(k);
arma::vec eta(n);
arma::vec mu(n);
arma::mat XWX(k,k);

arma::vec temp1(n);
arma::vec temp2(k);
double crit;
int counter;

arma::vec out = arma::zeros(B);

for (int i=0; i<B; i++) {

beta = PM;
crit = eps+1;
counter = 0;
while ((crit>eps)&(counter<50)) {

eta = x*beta;
mu = exp(eta);
XWX.zeros();
for (int ii=0; ii<k; ii++) {
for (int jj=ii; jj<k; jj++) {
for (int kk=0; kk<n; kk++) {
XWX.at(ii,jj) += mu.at(kk)*x.at(kk,ii)*x.at(kk,jj);}
XWX.at(jj,ii) = XWX.at(ii,jj);}}

for (int ii=0; ii<k; ii++) {
XWX.at(ii,ii) += iPV.at(ii);}

XWX = inv_sympd(XWX);

temp2 = beta;
temp2 -= PM;
temp2 %= miPV;

temp1 = Yt.unsafe_col(i);
temp1 -= mu;
temp2 += xt*temp1;

temp2 = XWX*temp2;

beta += scale*temp2;
counter += 1;
crit = dot(temp2,temp2);}

beta -= BETAT.col(i);
temp2 = beta;
beta %= temp2;

for (int ii=0; ii<k; ii++) {
out.at(i) -= beta.at(ii);}

}

return as<NumericVector>(wrap(out));

}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP sigprcpp(SEXP y, SEXP x, SEXP sam, SEXP frho, SEXP yfrac) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix samr(sam);
Rcpp::NumericVector frhor(frho);
Rcpp::NumericVector yfracr(yfrac);

int B = samr.nrow(), p = xr.ncol(), n = xr.nrow();

arma::mat Y(yr.begin(), B, n, false);
arma::mat X(xr.begin(), n, p, false);
arma::mat SAM(samr.begin(), B, p, false);
arma::vec a(frhor.begin(), B, false);
arma::vec b(yfracr.begin(), B, false);

arma::mat SAMt = SAM.t();

arma::mat suff = X.t()*Y.t();
arma::vec out = arma::zeros(B);
double temp;

for (int i=0; i<B; i++) {

for (int j=0; j<B; j++) {
temp = dot(suff.unsafe_col(i),SAMt.unsafe_col(j));
temp -= b.at(i);
temp -= a.at(j);
out.at(i) += exp(temp);}}

out /= B;

out = log(out);

return as<NumericVector>(wrap(out));
}


// NEW FUNCTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP LRNSELLAP2cpp(SEXP y, SEXP sam, SEXP X, SEXP pm, SEXP pv) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix samr(sam);
Rcpp::NumericMatrix Xr(X);
Rcpp::NumericVector pmr(pm);
Rcpp::NumericVector pvr(pv);

int n = Xr.nrow();
int k = Xr.ncol();
int B = yr.nrow();

arma::mat Y(yr.begin(), B, n, false);
arma::mat BETA(samr.begin(), B, k, false);
arma::mat x(Xr.begin(), n, k, false);
arma::vec PM(pmr.begin(), k, false);
arma::vec PV(pvr.begin(), k, false);

arma::mat BETAT = BETA.t();

double eps;
eps = 0.0001;

double scale;
scale = 0.25;

arma::mat xt = x.t();
arma::mat Yt = Y.t();

arma::vec iPV = 1/PV;
arma::vec miPV = -1/PV;
arma::mat invV = arma::zeros(k,k);
for (int j=0; j<k; j++) {
invV(j,j) = iPV(j);}

arma::vec beta(k);
arma::vec eta(n);
arma::vec prob(n);
arma::vec W(n);
arma::mat XWX(k,k);

arma::vec temp1(n);
arma::vec temp2(k);
//arma::mat Sigma(k,k);
double crit;
int counter;
//double dett;

arma::vec out = arma::zeros(B);

for (int i=0; i<B; i++) {

beta = PM;
crit = eps+1;
counter = 0;
while ((crit>eps)&(counter<50)) {

eta = x*beta;
prob = exp(-eta);
prob += 1;
prob = 1/prob;
W = 1 - prob;
W %= prob;
XWX.zeros();
for (int ii=0; ii<k; ii++) {
for (int jj=ii; jj<k; jj++) {
for (int kk=0; kk<n; kk++) {
XWX.at(ii,jj) += W.at(kk)*x.at(kk,ii)*x.at(kk,jj);}
XWX.at(jj,ii) = XWX.at(ii,jj);}}

for (int ii=0; ii<k; ii++) {
XWX.at(ii,ii) += iPV.at(ii);}

XWX = inv_sympd(XWX);

temp2 = beta;
temp2 -= PM;
temp2 %= miPV;

temp1 = Yt.unsafe_col(i);
temp1 -= prob;
temp2 += xt*temp1;

temp2 = XWX*temp2;

beta += scale*temp2;
counter += 1;
crit = dot(temp2,temp2);}

beta -= BETAT.col(i);
temp2 = beta;
beta %= temp2;

for (int ii=0; ii<k; ii++) {
out.at(i) -= beta.at(ii);}

}

return as<NumericVector>(wrap(out));

}

// NEW FUNCTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP LRLAP2cpp(SEXP y, SEXP sam, SEXP X, SEXP pm, SEXP pv) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix samr(sam);
Rcpp::NumericMatrix Xr(X);
Rcpp::NumericVector pmr(pm);
Rcpp::NumericVector pvr(pv);

int n = Xr.nrow();
int k = Xr.ncol();
int B = yr.nrow();

arma::mat Y(yr.begin(), B, n, false);
arma::mat BETA(samr.begin(), B, k, false);
arma::mat x(Xr.begin(), n, k, false);
arma::vec PM(pmr.begin(), k, false);
arma::vec PV(pvr.begin(), k, false);

arma::mat BETAT = BETA.t();

double eps;
eps = 0.0001;

double scale;
scale = 0.25;

arma::mat xt = x.t();
arma::mat Yt = Y.t();

arma::vec iPV = 1/PV;
arma::vec miPV = -1/PV;
arma::mat invV = arma::zeros(k,k);
for (int j=0; j<k; j++) {
invV(j,j) = iPV(j);}

arma::vec beta(k);
arma::vec eta(n);
arma::vec prob(n);
arma::vec W(n);
arma::mat XWX(k,k);

arma::vec temp1(n);
arma::vec temp2(k);
arma::mat Sigma(k,k);
double crit;
int counter;
double dett;

double cons;
cons =  sum(log(PV));
cons *= 0.5;
arma::vec out(B);
out.fill(cons);

for (int i=0; i<B; i++) {

beta = PM;
crit = eps+1;
counter = 0;
while ((crit>eps)&(counter<50)) {

eta = x*beta;
prob = exp(-eta);
prob += 1;
prob = 1/prob;
W = 1 - prob;
W %= prob;
XWX.zeros();
for (int ii=0; ii<k; ii++) {
for (int jj=ii; jj<k; jj++) {
for (int kk=0; kk<n; kk++) {
XWX.at(ii,jj) += W.at(kk)*x.at(kk,ii)*x.at(kk,jj);}
XWX.at(jj,ii) = XWX.at(ii,jj);}}

for (int ii=0; ii<k; ii++) {
XWX.at(ii,ii) += iPV.at(ii);}

XWX = inv_sympd(XWX);

temp2 = beta;
temp2 -= PM;
temp2 %= miPV;

temp1 = Yt.unsafe_col(i);
temp1 -= prob;
temp2 += xt*temp1;

temp2 = XWX*temp2;

beta += scale*temp2;
counter += 1;
crit = dot(temp2,temp2);}

eta = x*BETAT.col(i);
out(i) += dot(Y.row(i),eta);
eta = exp(eta);
eta += 1;
eta = log(eta);
out(i) -= sum(eta);

eta = x*beta;
out(i) -= dot(Y.row(i),eta);
eta = exp(eta);
eta += 1;
eta = log(eta);
out(i) += sum(eta);

beta -= PM;
temp2 = beta;
beta %= temp2;
beta %= iPV;
out(i) += 0.5*sum(beta);
log_det(dett,crit,XWX);
out(i) -= 0.5*dett;

}

return as<NumericVector>(wrap(out));

}



// NEW FUNCTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP Enlmcpp(SEXP jac, SEXP dims) {

Rcpp::NumericMatrix jacr(jac);
Rcpp::NumericVector dimsr(dims);

int p = jacr.ncol();
int r = jacr.nrow();

arma::mat JAC(jacr.begin(), r, p, false);
arma::vec dimsrr(dimsr.begin(), 2, false);

arma::uvec DIMS = arma::conv_to<arma::uvec>::from(dimsrr);
int n = DIMS(0);
int B = DIMS(1);

int nm = n - 1;
int F;
int L;

arma::vec out = arma::zeros(B);
arma::mat J = arma::zeros(n,p);
arma::mat FIM = arma::zeros(p,p);
arma::vec eigs = arma::zeros(p);

for (int i=0; i<B; i++) {

F = i*n;
L = F + nm;
J = JAC.rows(F,L);
FIM = J.t()*J;

eig_sym(eigs,FIM);

out.at(i) = eigs.min();

}

return as<NumericVector>(wrap(out));

}

// NEW FUNCTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP Anlmcpp(SEXP jac, SEXP dims) {

Rcpp::NumericMatrix jacr(jac);
Rcpp::NumericVector dimsr(dims);

int p = jacr.ncol();
int r = jacr.nrow();

arma::mat JAC(jacr.begin(), r, p, false);
arma::vec dimsrr(dimsr.begin(), 2, false);

arma::uvec DIMS = arma::conv_to<arma::uvec>::from(dimsrr);
int n = DIMS(0);
int B = DIMS(1);

int nm = n - 1;
int F;
int L;

arma::vec out = arma::zeros(B);
arma::mat J = arma::zeros(n,p);
arma::mat FIM = arma::zeros(p,p);

for (int i=0; i<B; i++) {

F = i*n;
L = F + nm;
J = JAC.rows(F,L);
FIM = J.t()*J;
FIM = inv_sympd(FIM);

for (int j=0; j<p; j++) {
out.at(i) -= FIM.at(j,j);}

}

return as<NumericVector>(wrap(out));

}

// NEW FUNCTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP Dnlmcpp(SEXP jac, SEXP dims) {

Rcpp::NumericMatrix jacr(jac);
Rcpp::NumericVector dimsr(dims);

int p = jacr.ncol();
int r = jacr.nrow();

arma::mat JAC(jacr.begin(), r, p, false);
arma::vec dimsrr(dimsr.begin(), 2, false);

arma::uvec DIMS = arma::conv_to<arma::uvec>::from(dimsrr);
int n = DIMS(0);
int B = DIMS(1);

int nm = n - 1;
int F;
int L;

arma::vec out = arma::zeros(B);
arma::mat J = arma::zeros(n,p);
arma::mat FIM = arma::zeros(p,p);
double temp;

for (int i=0; i<B; i++) {

F = i*n;
L = F + nm;
J = JAC.rows(F,L);
FIM.zeros();
for (int j1=0; j1<p; j1++) {
for (int j2=j1; j2<p; j2++) {
FIM.at(j1,j2) = dot(J.col(j1),J.col(j2));}}
FIM = symmatu(FIM);

log_det(out(i),temp,FIM);}

return as<NumericVector>(wrap(out));

}


// NEW FUNCTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP SIGnlmcpp(SEXP y, SEXP mu1, SEXP mu2, SEXP sig1, SEXP sig2) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix mu1r(mu1);
Rcpp::NumericMatrix mu2r(mu2);
Rcpp::NumericVector sig1r(sig1);
Rcpp::NumericVector sig2r(sig2);

int B = mu2r.nrow();
int n = mu2r.ncol();

arma::mat Y(yr.begin(), B, n, false);
arma::mat MU1(mu1r.begin(), B, n, false);
arma::mat MU2(mu2r.begin(), B, n, false);
arma::vec SIG1(sig1r.begin(), B, false);
arma::vec SIG2(sig2r.begin(), B, false);

arma::vec iSIG1 = -0.5/SIG1;
arma::vec iSIG2 = -0.5/SIG2;
arma::vec LSIG1 = -0.5*n*log(SIG1);
arma::vec LSIG2 = -0.5*n*log(SIG2);

arma::mat Yt = Y.t();
arma::mat MU1t = MU1.t();
arma::mat MU2t = MU2.t();

arma::vec diff(n);
double temp;
double LL;
arma::vec yi(n);

arma::vec out = arma::zeros(B);
for (int i=0; i<B; i++) {

yi = Yt.unsafe_col(i);
diff = yi;
diff -= MU1t.unsafe_col(i);
out(i) += dot(diff,diff);
out(i) *= iSIG1.at(i);
out(i) += LSIG1.at(i);

LL = 0;
for (int j=0; j<B; j++) {
diff = yi;
diff -= MU2t.unsafe_col(j);
temp = dot(diff,diff);
temp *= iSIG2.at(j);
temp += LSIG2.at(j);
temp = exp(temp);
LL += temp;}
LL /= B;

out.at(i) -= log(LL);}

return as<NumericVector>(wrap(out));

}

// NEW FUNCTION /////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP NSELnlmcpp(SEXP y, SEXP mu2, SEXP sig2, SEXP theta1, SEXP theta2) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix mu2r(mu2);
Rcpp::NumericVector sig2r(sig2);
Rcpp::NumericMatrix theta1r(theta1);
Rcpp::NumericMatrix theta2r(theta2);

int B = mu2r.nrow();
int n = mu2r.ncol();
int p = theta1r.ncol();

arma::mat Y(yr.begin(), B, n, false);
arma::mat MU2(mu2r.begin(), B, n, false);
arma::vec SIG2(sig2r.begin(), B, false);
arma::mat THETA1(theta1r.begin(), B, p, false);
arma::mat THETA2(theta2r.begin(), B, p, false);

arma::vec iSIG2 = -0.5/SIG2;
arma::vec LSIG2 = -0.5*n*log(SIG2);

arma::mat Yt = Y.t();
arma::mat MU2t = MU2.t();
arma::mat THETA1t = THETA1.t();
arma::mat THETA2t = THETA2.t();

arma::vec diff(n);
double temp;
double LL;
arma::vec yi(n);
arma::vec dcr = arma::zeros(p);

arma::vec out = arma::zeros(B);
for (int i=0; i<B; i++) {

yi = Yt.unsafe_col(i);

LL = 0;
dcr.zeros();
for (int j=0; j<B; j++) {
diff = yi;
diff -= MU2t.unsafe_col(j);
temp = dot(diff,diff);
temp *= iSIG2.at(j);
temp += LSIG2.at(j);
temp = exp(temp);
dcr += temp*THETA2t.unsafe_col(j);
LL += temp;}

dcr /= LL;

dcr -= THETA1t.unsafe_col(i);

out.at(i) -= dot(dcr,dcr);}

return as<NumericVector>(wrap(out));

}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP LRDcpp(SEXP x, SEXP beta) {

Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix betar(beta);

int n = xr.nrow();
int p = xr.ncol();
int B = betar.nrow();

arma::mat xs(xr.begin(), n, p, false);
arma::mat betas(betar.begin(), B, p, false);

arma::mat eta(B,n);
arma::mat W(B,n);
eta = betas*xs.t();
W = arma::trunc_exp(eta);
W += 1;
W = arma::trunc_log(W);
W *= -2;
W += eta;
W = exp(W);

arma::mat xwx(p,p);
arma::vec DETS = arma::zeros(B);

for(int j=0; j<B; j++){

for(int i1=0; i1<p; i1++){
for(int i2=i1; i2<p; i2++){
xwx(i1,i2) = 0;
for(int i3=0; i3<n; i3++){
xwx(i1,i2) += W(j,i3)*xs(i3,i1)*xs(i3,i2);
}
xwx(i2,i1) = xwx(i1,i2);
}}

xwx = arma::chol(xwx);

for (int i4=0; i4<p; i4++){
DETS(j) += log(xwx(i4,i4));}
DETS(j) += DETS(j);

}

return as<NumericVector>(wrap(DETS));

}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP HLRDcpp(SEXP x, SEXP z, SEXP beta, SEXP gam, SEXP S) {

Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix zr(z);
Rcpp::NumericMatrix betar(beta);
Rcpp::NumericMatrix gamr(gam);
Rcpp::NumericMatrix Sr(S);

int n = xr.nrow();
int p = xr.ncol();
int B = betar.nrow();
int G = gamr.ncol();

arma::mat X(xr.begin(), n, p, false);
arma::mat Z(zr.begin(), n, zr.ncol(), false);
arma::mat BETA(betar.begin(), B, p, false);
arma::mat GAMMA(gamr.begin(), B, G, false);
arma::mat s(Sr.begin(), B, p, false);

G /= p;
int m = n/G;
int mm = m - 1;

arma::mat eta(B,n);
arma::mat eta2(B,n);
arma::mat iW(B,n);
eta = BETA*X.t();
eta2 = GAMMA*Z.t();
eta += eta2;
iW = arma::trunc_exp(eta);
iW += 1;
iW = arma::trunc_log(iW);
iW *= 2;
iW -= eta;
iW = exp(iW);

arma::mat xwx(p,p);
arma::mat temp(m,m);
arma::mat littlex(m,p);
arma::mat littletx(p,m);
int F;
int L;
arma::vec DETS = arma::zeros(B);
double temp2;

for(int i=0; i<B; i++){
xwx.zeros();

for(int j=0; j<G; j++){
F = j*m;
L = F + mm;
littlex = X(arma::span(F,L),arma::span::all);
littletx = littlex.t();

for(int i1=0; i1<m; i1++){
for(int i2=i1; i2<m; i2++){
temp(i1,i2) = 0;
for(int i3=0; i3<p; i3++){
temp(i1,i2) += s(i,i3)*littlex(i1,i3)*littlex(i2,i3);
}
temp(i2,i1) = temp(i1,i2);
}}
for(int i1=0; i1<m; i1++){
temp(i1,i1) += iW(i,F+i1);}

temp = inv_sympd(temp);

littletx = littletx*temp;
xwx += littletx*littlex;}

log_det(DETS(i), temp2, xwx);

}

return as<NumericVector>(wrap(DETS));

}


// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP LRAcpp(SEXP x, SEXP beta) {

Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix betar(beta);

int n = xr.nrow();
int p = xr.ncol();
int B = betar.nrow();

arma::mat xs(xr.begin(), n, p, false);
arma::mat betas(betar.begin(), B, p, false);

arma::mat eta(B,n);
arma::mat W(B,n);
eta = betas*xs.t();
W = arma::trunc_exp(eta);
W += 1;
W = arma::trunc_log(W);
W *= -2;
W += eta;
W = exp(W);

arma::mat xwx(p,p);
arma::vec DETS = arma::zeros(B);
arma::vec eigs(p);
double temp;

for(int j=0; j<B; j++){

for(int i1=0; i1<p; i1++){
for(int i2=i1; i2<p; i2++){
xwx(i1,i2) = 0;
for(int i3=0; i3<n; i3++){
xwx(i1,i2) += W(j,i3)*xs(i3,i1)*xs(i3,i2);
}
xwx(i2,i1) = xwx(i1,i2);
}}

eig_sym(eigs,xwx);

for (int i4=0; i4<p; i4++){
temp = 1/eigs(i4);
DETS(j) -= temp;}

}

return as<NumericVector>(wrap(DETS));

}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP HLRAcpp(SEXP x, SEXP z, SEXP beta, SEXP gam, SEXP S) {

Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix zr(z);
Rcpp::NumericMatrix betar(beta);
Rcpp::NumericMatrix gamr(gam);
Rcpp::NumericMatrix Sr(S);

int n = xr.nrow();
int p = xr.ncol();
int B = betar.nrow();
int G = gamr.ncol();

arma::mat X(xr.begin(), n, p, false);
arma::mat Z(zr.begin(), n, zr.ncol(), false);
arma::mat BETA(betar.begin(), B, p, false);
arma::mat GAMMA(gamr.begin(), B, G, false);
arma::mat s(Sr.begin(), B, p, false);

G /= p;
int m = n/G;
int mm = m - 1;

arma::mat eta(B,n);
arma::mat eta2(B,n);
arma::mat iW(B,n);
eta = BETA*X.t();
eta2 = GAMMA*Z.t();
eta += eta2;
iW = arma::trunc_exp(eta);
iW += 1;
iW = arma::trunc_log(iW);
iW *= 2;
iW -= eta;
iW = exp(iW);

arma::mat xwx(p,p);
arma::mat temp(m,m);
arma::mat littlex(m,p);
arma::mat littletx(p,m);
int F;
int L;
arma::vec DETS = arma::zeros(B);
arma::vec eigs(p);
double temp2;

for(int i=0; i<B; i++){
xwx.zeros();

for(int j=0; j<G; j++){
F = j*m;
L = F + mm;
littlex = X(arma::span(F,L),arma::span::all);
littletx = littlex.t();

for(int i1=0; i1<m; i1++){
for(int i2=i1; i2<m; i2++){
temp(i1,i2) = 0;
for(int i3=0; i3<p; i3++){
temp(i1,i2) += s(i,i3)*littlex(i1,i3)*littlex(i2,i3);
}
temp(i2,i1) = temp(i1,i2);
}}
for(int i1=0; i1<m; i1++){
temp(i1,i1) += iW(i,F+i1);}

temp = inv_sympd(temp);

littletx = littletx*temp;
xwx += littletx*littlex;}

eig_sym(eigs,xwx);

for (int i4=0; i4<p; i4++){
temp2 = 1/eigs(i4);
DETS(i) -= temp2;}

}

return as<NumericVector>(wrap(DETS));

}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP Dcpp(SEXP x, SEXP W) {

Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix Wr(W);

int n = xr.nrow();
int p = xr.ncol();
int B = Wr.nrow();

arma::mat xs(xr.begin(), n, p, false);
arma::mat Ws(Wr.begin(), B, n, false);

arma::mat xwx(p,p);
arma::vec DETS = arma::zeros(B);

for(int j=0; j<B; j++){

for(int i1=0; i1<p; i1++){
for(int i2=i1; i2<p; i2++){
xwx(i1,i2) = 0;
for(int i3=0; i3<n; i3++){
xwx(i1,i2) += Ws(j,i3)*xs(i3,i1)*xs(i3,i2);
}
xwx(i2,i1) = xwx(i1,i2);
}}

xwx = arma::chol(xwx);

for (int i4=0; i4<p; i4++){
DETS(j) += log(xwx(i4,i4));}
DETS(j) += DETS(j);

}

return as<NumericVector>(wrap(DETS));

}


// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP Acpp(SEXP x, SEXP W) {

Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix Wr(W);

int n = xr.nrow();
int p = xr.ncol();
int B = Wr.nrow();

arma::mat xs(xr.begin(), n, p, false);
arma::mat Ws(Wr.begin(), B, n, false);

arma::mat xwx(p,p);
arma::vec DETS = arma::zeros(B);
arma::vec eigs(p);
double temp;

for(int j=0; j<B; j++){

for(int i1=0; i1<p; i1++){
for(int i2=i1; i2<p; i2++){
xwx(i1,i2) = 0;
for(int i3=0; i3<n; i3++){
xwx(i1,i2) += Ws(j,i3)*xs(i3,i1)*xs(i3,i2);
}
xwx(i2,i1) = xwx(i1,i2);
}}

eig_sym(eigs,xwx);

for (int i4=0; i4<p; i4++){
temp = 1/eigs(i4);
DETS(j) -= temp;}

}

return as<NumericVector>(wrap(DETS));

}


// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP Ecpp(SEXP x, SEXP W) {

Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix Wr(W);

int n = xr.nrow();
int p = xr.ncol();
int B = Wr.nrow();

arma::mat xs(xr.begin(), n, p, false);
arma::mat Ws(Wr.begin(), B, n, false);

arma::mat xwx(p,p);
arma::vec DETS = arma::zeros(B);
arma::vec eigs(p);

for(int j=0; j<B; j++){

for(int i1=0; i1<p; i1++){
for(int i2=i1; i2<p; i2++){
xwx(i1,i2) = 0;
for(int i3=0; i3<n; i3++){
xwx(i1,i2) += Ws(j,i3)*xs(i3,i1)*xs(i3,i2);
}
xwx(i2,i1) = xwx(i1,i2);
}}

eig_sym(eigs,xwx);

DETS(j) = eigs.min();

}

return as<NumericVector>(wrap(DETS));

}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP pvalcpp(SEXP oldeval, SEXP neweval) {

Rcpp::NumericVector oldevalr(oldeval);
Rcpp::NumericVector newevalr(neweval);

int B = oldevalr.size();

arma::vec oldv(oldevalr.begin(), B, false);
arma::vec newv(newevalr.begin(), B, false);

int n = 2*B;

double omean;
double nmean;
omean=0;
nmean=0;
for(int i=0; i<B; i++){
omean += oldv(i)/B;
nmean += newv(i)/B;}

double a;
double aa;
a = 0;
for(int i=0; i<B; i++){
a += oldv(i)*(oldv(i) - omean);
a += newv(i)*(newv(i) - nmean);}
aa = ::sqrt(a);

arma::vec ans(2);

ans(0) = -B*(nmean-omean)/aa;
ans(1) = n;

return as<NumericVector>(wrap(ans));

}

// NEW FUNCTION;

RcppExport SEXP distcpp(SEXP Dij) {

Rcpp::NumericVector xr(Dij);

arma::vec xs(xr.begin(), xr.size(), false);

int Q = xs.n_elem;

arma::mat Aarr(Q,Q);

for (int i=0; i<Q; i++){
for (int j=i; j<Q; j++){
Aarr(i,j) = xs(i)-xs(j);
Aarr(i,j) *= Aarr(i,j);
Aarr(j,i) = Aarr(i,j);}}

return as<NumericMatrix>(wrap(Aarr));

}

// NEW FUNCTION;

RcppExport SEXP GPpredcpp(SEXP paras, SEXP dist, SEXP z, SEXP newDij, SEXP Dij) {

Rcpp::NumericVector thetar(paras);
Rcpp::NumericMatrix Aar(dist);
Rcpp::NumericVector zzzr(z);
Rcpp::NumericVector xxxr(newDij);
Rcpp::NumericVector xr(Dij);

int Q = zzzr.size();

arma::mat ARR(Aar.begin(), Q, Q, false);
arma::vec THETA(thetar.begin(), thetar.size(), false);
arma::mat Z(zzzr.begin(), Q , 1 , false);
arma::vec xp(xxxr.begin(), xxxr.size(), false);
arma::vec x(xr.begin(), xr.size(), false);

int m = xp.n_elem;

double RR;
double ET;
RR = exp(THETA(1));
ET = exp(THETA(0));

arma::mat B = exp(-RR*ARR);
for (int i=0; i<Q; i++){
B(i,i) += ET;}
B = inv_sympd(B);
arma::mat ibz = B*Z;

arma::mat ART(m,Q);
for(int i=0; i<m; i++){
for(int j=0; j<Q; j++){
ART(i,j) = xp(i)-x(j);
ART(i,j) *= ART(i,j);
ART(i,j) *= -RR;
ART(i,j) = exp(ART(i,j));}}

arma::mat ans = ART*ibz;

return as<NumericMatrix>(wrap(ans));

}

// NEW FUNCTION;

RcppExport SEXP solve2Dcpp(SEXP x) {

Rcpp::NumericMatrix xr(x);

arma::mat xs(xr.begin(), 2 , 2 , false);

arma::mat iB(2,2);

double deti;

deti = xs(0,0);
deti *= xs(1,1);
deti -= xs(1,0)*xs(1,0);

iB(0,0) = xs(1,1);
iB(0,1) = -xs(0,1);
iB(1,0) = iB(0,1);
iB(1,1) = xs(0,0);

iB /= deti;

return as<NumericMatrix>(wrap(iB));

}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;


RcppExport SEXP utilcomp18badcpp(SEXP d, SEXP sam) {

Rcpp::NumericMatrix dr(d);
Rcpp::NumericMatrix samr(sam);

int n = dr.nrow();
int BB = samr.nrow();

arma::mat ds(dr.begin(), n, 1, false);
arma::mat sams(samr.begin(), BB, 2, false);

arma::vec times = vectorise(ds);
times += 1;
times *= 12;
arma::vec times2 = times;
times2 %= times;

double ltheta34 = log(21.8);
ltheta34 *= 4;

arma::vec theta = arma::zeros(2);
arma::vec e1(n);
arma::vec e2(n);
arma::vec e11(n);
arma::vec e22(n);
arma::vec e12(n);
arma::vec ec(n);

arma::vec c11 = arma::zeros(BB);
arma::vec c22 = arma::zeros(BB);
arma::vec c33 = arma::zeros(BB);
arma::vec c12 = arma::zeros(BB);
arma::vec c13 = arma::zeros(BB);
arma::vec c23 = arma::zeros(BB);

double temp;
arma::vec tempv(n);

for(int j=0; j<BB; j++){

theta = vectorise(sams.row(j));

e1 = times;
e2 = times;
e1 *= -theta(0);
e2 *= -theta(1);
e11 = e1;
e22 = e2;
ec = e1;
e11 += e1;
e22 += e2;
ec += e2;
e1 = exp(e1);
e2 = exp(e2);
e11 = exp(e11);
e22 = exp(e22);
ec = exp(ec);
e12 = e1;
e12 -= e2;

temp = dot(times2,e11);
c11(j) = temp;
temp = dot(times2,e22);
c22(j) = temp;
temp = dot(times2,ec);
c12(j) = temp;
tempv = times%e1%e12;
temp = sum(tempv);
c13(j) = temp;
tempv = times%e2%e12;
temp = sum(tempv);
c23(j) = temp;
tempv = e12;
tempv %= e12;
temp = sum(tempv);
c33(j) = temp;

}

arma::vec deter = c11%c22%c33+2*c12%c23%c13-c13%c13%c22-c12%c12%c33-c23%c23%c11;
deter = log(deter);
deter += ltheta34;

return as<NumericVector>(wrap(deter));

}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;


RcppExport SEXP utilcomp15badcpp(SEXP d2, SEXP theta) {

Rcpp::NumericMatrix dr(d2);
Rcpp::NumericMatrix samr(theta);

int n = dr.nrow();
int B = samr.nrow();

arma::mat D(dr.begin(), n, 1, false);
arma::mat THETA(samr.begin(), B, 3, false);

int DDD = 400;
double tau = 0.01;
double sig = 0.1;
double tausig = tau/sig;
double diff;
double rat;

arma::vec ANS = arma::zeros(B);
arma::mat dmu(n,3);
arma::mat mu(n,1);
arma::mat t1(n,1);
arma::mat t2(n,1);
double c;
arma::mat dc(1,3);
arma::vec phi(n);
arma::vec v(n);
arma::mat dphi(n,3);
arma::mat dv(n,3);
arma::vec iv(n);
arma::mat I1(3,3);
arma::mat I2(3,3);

for (int i=0; i<B; i++) {
dmu.zeros();
mu.zeros();

t1 = D;
t2 = D;
t1 *= THETA(i,0);
t2 *= THETA(i,1);
t1 = exp(-t1);
t2 = exp(-t2);
mu = t1;
mu -= t2;
t1 %= D;
t2 %= D;
t1 *= -1;
dmu.col(0) = t1;
dmu.col(1) = t2;

diff = THETA(i,1) - THETA(i,0);
rat = DDD;
rat /= THETA(i,2);

c = rat;
c *= THETA(i,1);
c /= diff;

dc(0,2) = THETA(i,1);
dc(0,2) /= diff;
dc(0,2) *= rat;
dc(0,2) /= THETA(i,2);
dc(0,2) *= -1;
diff *= diff;
dc(0,1) = THETA(i,0);
dc(0,1) *= rat;
dc(0,1) /= diff;
dc(0,0) = THETA(i,1);
dc(0,0) *= rat;
dc(0,0) /= diff;
dc(0,1) *= -1;

dphi = dmu;
dphi *= c;
dphi += mu*dc;
dv = dphi;

mu *= c;
phi = vectorise(mu);
v = phi;
//v = v%v;
v %= phi;
v *= tausig;
v += 1;

dv.each_col() %= phi;
dv *= tausig;
dv *= 2;

iv = 1/v;
dv.each_col() %= iv;

I1 = dv.t()*dv;
I1 *= 0.5;

v = sqrt(v);
iv = 1/v;
dphi.each_col() %= iv;

I2 = dphi.t()*dphi;
I2 /= sig;
I1 = I1 + I2;

log_det(diff, rat, I1);

ANS(i) = diff;

}

return as<NumericVector>(wrap(ANS));

}


// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP rowSumscpp(SEXP Z) {

Rcpp::NumericMatrix xr(Z);

int nrows = xr.nrow(), ncols = xr.ncol();

arma::mat x(xr.begin(), nrows, ncols, false);

arma::vec sums(nrows);
for (int row=0; row<nrows; row++) {
    sums(row) = 0;
    for (int col=0; col<ncols; col++) {
        sums(row) += x(row,col);
    }
}
return as<NumericVector>(wrap(sums));
}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP siglrcpp(SEXP y, SEXP x, SEXP sam, SEXP frho) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix samr(sam);
Rcpp::NumericVector frhor(frho);

int B = samr.nrow(), p = xr.ncol(), n = xr.nrow();

arma::mat Y(yr.begin(), B, n, false);
arma::mat X(xr.begin(), n, p, false);
arma::mat SAM(samr.begin(), B, p, false);
arma::vec FRHO(frhor.begin(), B, false);

arma::vec ANS = arma::zeros(B);
double LL;
arma::vec suffi= arma::zeros(p);;

for (int i=0; i<B; i++){
for (int k=0; k<p; k++){
suffi.at(k) = dot(Y.row(i),X.col(k));}
for (int j=0; j<B; j++){
LL = dot(SAM.row(j),suffi);
LL += FRHO.at(j);
ANS.at(i) += exp(LL);}}

return as<NumericVector>(wrap(ANS));
}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP nsellrcpp(SEXP y, SEXP x, SEXP sam, SEXP frho) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix xr(x);
Rcpp::NumericMatrix samr(sam);
Rcpp::NumericVector frhor(frho);

int B = samr.nrow(), p = xr.ncol(), n = xr.nrow();

arma::mat Y(yr.begin(), B, n, false);
arma::mat X(xr.begin(), n, p, false);
arma::mat SAM(samr.begin(), B, p, false);
arma::vec FRHO(frhor.begin(), B, false);

arma::mat ANS = arma::zeros(B,p);
double LL;
double SS;
arma::vec suffi= arma::zeros(p);

for (int i=0; i<B; i++){
for (int k=0; k<p; k++){
suffi.at(k) = dot(Y.row(i),X.col(k));}
SS = 0;
for (int j=0; j<B; j++){
LL = dot(SAM.row(j),suffi);
LL += FRHO.at(j);
LL = exp(LL);
SS += LL;
ANS.row(i) += LL*SAM.row(j);}
ANS.row(i) /= SS;}

return as<NumericMatrix>(wrap(ANS));
}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP LMcpp(SEXP x) {

Rcpp::NumericMatrix xr(x);

int p = xr.ncol(), n = xr.nrow();

arma::mat X(xr.begin(), n, p, false);

arma::mat XTX = arma::zeros(p,p);
double ANS = 0;
double temp;

for (int i=0; i<p; i++){
for (int j=i; j<p; j++){
for (int k=0; k<n; k++){
temp = X(k,i);
temp *= X(k,j);
XTX(i,j) += temp;}
XTX(j,i) = XTX(i,j);}}

arma::log_det(ANS, temp, XTX);


return as<NumericVector>(wrap(ANS));
}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP utilcomp15sigcpp(SEXP y, SEXP mu, SEXP vv, SEXP frho) {

Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericMatrix vvr(vv);
Rcpp::NumericVector frhor(frho);

int B = yr.nrow(), n = yr.ncol();

arma::mat Y(yr.begin(), B, n, false);
arma::mat MU(mur.begin(), B, n, false);
arma::mat VV(vvr.begin(), B, n, false);
arma::vec FRHO(frhor.begin(), B, false);

arma::vec ANS = arma::zeros(B);
arma::rowvec temp1(n);
arma::rowvec temp2(n);
arma::rowvec temp3(n);
double LL;
double harris;
harris = arma::datum::pi;
harris *= 2;
harris = log(harris);
harris *= 0.5*n;

for (int i=0; i<B; i++){
temp1 = Y.row(i);
for (int j=0; j<B; j++){
temp2 = temp1;
temp2 -= MU.row(j);
//temp2 = temp2%temp2;
temp3 = temp2;
temp2 %= temp3; 
temp2 /= VV.row(j);
LL = sum(temp2);
LL += FRHO(j);
LL *= -0.5;
LL -= harris;
LL = exp(LL);
ANS(i) += LL;}}

return as<NumericVector>(wrap(ANS));

}

// NEW FUNCTION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////;


RcppExport SEXP beetlecpp(SEXP phi, SEXP y, SEXP mu, SEXP frho, SEXP ncr) {

Rcpp::NumericVector phir(phi);
Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix mur(mu);
Rcpp::NumericVector frhor(frho);
Rcpp::NumericVector ncrr(ncr);

int B = yr.nrow(), n = yr.ncol();

arma::vec PHI(phir.begin(), B, false);
arma::mat Y(yr.begin(), B, n, false);
arma::mat MU(mur.begin(), B, n, false);
arma::vec FRHO(frhor.begin(), B, false);
arma::vec NCR(ncrr.begin(), B, false);

arma::rowvec yy(n);
arma::rowvec mm(n);

double LL;
double SS;
double PP;
double temp;
double temp2;

arma::vec ANS = arma::zeros(B);

for(int i=0; i<B; i++){

yy = Y.row(i);
SS = 0;
PP = 0;
temp = NCR(i);

for(int j=0; j<B; j++){
mm = MU.row(j);
LL = temp;
temp2 = dot(yy,mm);
LL += temp2;
LL += FRHO(j);
LL = exp(LL);
SS += LL;
LL *= PHI(j);
PP += LL;}

ANS(i) = PP;
ANS(i) /= SS;

}

return as<NumericVector>(wrap(ANS));

}

// NEW FUNCTION

RcppExport SEXP nselhlrcpp(SEXP Fd, SEXP y, SEXP fam, SEXP frho) {

Rcpp::NumericMatrix Fdr(Fd);
Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix famr(fam);
Rcpp::NumericVector frhor(frho);

int n = Fdr.nrow(), p = Fdr.ncol(), B = famr.nrow();

arma::mat F(Fdr.begin(), n, p, false);
arma::mat Y(yr.begin(), B, n, false);
arma::mat FAM(famr.begin(), B, p, false);
arma::vec FRHO(frhor.begin(), B, false);

int k = 5;

arma::mat ANS = arma::zeros(B,k);
double LL;
double SS;
arma::vec suffi= arma::zeros(p);

for (int i=0; i<B; i++){
for (int k=0; k<p; k++){
suffi.at(k) = dot(Y.row(i),F.col(k));}
SS = 0;
for (int j=0; j<B; j++){
LL = dot(FAM.row(j),suffi);
LL += FRHO.at(j);
LL = exp(LL);
SS += LL;
ANS.row(i) += LL*FAM(j,arma::span(0,4));}
ANS.row(i) /= SS;}

return as<NumericMatrix>(wrap(ANS));

}

// NEW FUNCTION

RcppExport SEXP sighlrcpp(SEXP Fd, SEXP y, SEXP fam, SEXP frho, SEXP sam, SEXP etax, SEXP etaz) {

Rcpp::NumericMatrix Fdr(Fd);
Rcpp::NumericMatrix yr(y);
Rcpp::NumericMatrix famr(fam);
Rcpp::NumericVector frhor(frho);
Rcpp::NumericMatrix samr(sam);
Rcpp::NumericMatrix etaxr(etax);
Rcpp::NumericMatrix etazr(etaz);

int n = Fdr.nrow(), p = samr.ncol(), p2 = Fdr.ncol(), B = samr.nrow();

int B3 = famr.nrow();
int B2 = B3/2;

arma::mat F(Fdr.begin(), n, p2, false);
arma::mat Y(yr.begin(), B, n, false);
arma::mat FAM(famr.begin(), B3, p2, false);
arma::vec FRHO(frhor.begin(), B2, false);
arma::mat SAM(samr.begin(), B, p, false);
arma::mat ETAX(etaxr.begin(), B, n, false);
arma::mat ETAZ(etazr.begin(), B2, n, false);

int p3 = p2 - p;
int temp = p2-1;

arma::mat ANS = arma::zeros(B,2);
double AA;
double LL;
double LL2;
arma::vec suffi= arma::zeros(p2);;
arma::vec suffiA= arma::zeros(p);;
arma::vec suffiB= arma::zeros(p3);;
arma::mat tempmat = arma::zeros(B2,n);
arma::mat tempvec = arma::zeros(B2);

for (int i=0; i<B; i++){

tempmat = ETAZ;
tempmat.each_row() %= ETAX.row(i);
tempmat += 1;
tempmat = log(tempmat);
for (int j=0; j<B2; j++){
tempvec.at(j) = 0;
for (int k=0; k<n; k++){
tempvec.at(j) += tempmat.at(j,k);}}

for (int k=0; k<p2; k++){
suffi.at(k) = dot(Y.row(i),F.col(k));}

for (int k=0; k<p; k++){
suffiA.at(k) = dot(Y.row(i),F.col(k));}

for (int k=0; k<p3; k++){
suffiB.at(k) = dot(Y.row(i),F.col(p+k));}

AA = dot(SAM.row(i),suffiA);

for (int j=0; j<B2; j++){
LL = dot(FAM.row(j),suffi);
LL += FRHO.at(j);
ANS.at(i,0) += exp(LL);
LL2 = dot(FAM(B2+j,arma::span(5,temp)),suffiB);
LL2 += AA;
LL2 -= tempvec.at(j);
ANS.at(i,1) += exp(LL2);}}

return as<NumericMatrix>(wrap(ANS));

}
