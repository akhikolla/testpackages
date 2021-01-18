/*
    A number of functions utilized by rcppfunc.cpp.

    Intended for use with R.
    Copyright (C) 2015 Adam Lund

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
using namespace std;
//using namespace arma;

////////////////////////////////// Auxiliary functions
//////////////////// Direct RH-transform of a flat 3d array (matrix) M by a matrix X
arma::mat RHmat(arma::mat const& X, arma::mat const& M,int col, int sli){

int rowx = X.n_rows;

////matrix multiply
arma::mat XM = X * M;

////make matrix into rotated (!) cube (on matrix form)
arma::mat Mnew(col, sli * rowx);
for (int s = 0; s < sli; s++) {

for (int c = 0; c < col; c++) {

for (int r = 0; r < rowx; r++) {

Mnew(c, s + r * sli) = XM(r, c + s * col);

}

}

}

return Mnew;

}

////////////////// empirical explained variance function
arma::vec eev(arma::mat XBeta, arma::cube Z, int ng){

arma::vec eevar(Z.n_slices);
double sumXBeta2 = accu(pow(XBeta, 2));

for(int j = 0; j < Z.n_slices; j++) {eevar(j) = (2 * accu(XBeta  %  Z.slice(j))  - sumXBeta2) / ng;}

return eevar;

}

//////////////////// softmax loss
double softmaxloss(arma::vec h, double c, int ll){

if(ll == 1){

double k =  max(h);
return log(accu(exp(c * (h - k)))) / c +  k;

}else{
return accu(exp(c * h));

}

}
//////////////////// gradloss
arma::mat gradloss(arma::cube const& PhitZ, arma::mat const& XtXb, arma::vec const& h, int ng, double c, int ll){

arma::mat gradout(XtXb.n_rows, XtXb.n_cols);
gradout.fill(0);

if(ll == 1){
double  k = max(h);
double tmp =  accu(exp(c * (h - k)));
 for(int j = 0; j < PhitZ.n_slices; j++){
gradout = exp(c * (h(j) - k)) * (XtXb - PhitZ.slice(j)) + gradout;
}
 return 2 * gradout / (tmp * ng);

}else{

for(int j = 0; j < PhitZ.n_slices; j++){gradout = exp(c * h(j)) * (XtXb - PhitZ.slice(j)) + gradout;}

return 2 * c * gradout / ng;

}

}

//////////////////// Sum of squares function
double sum_square(arma::mat const& x){return accu(x % x);}

//////////////////// Proximal operator for the l1-penalty (soft threshold)
arma::mat prox_l1(arma::mat const& zv, arma::mat const& gam){

return (zv >= gam) % (zv - gam) + (zv <= -gam) % (zv + gam);

}

//////////////////// The weighted (gam = penaltyfactor * lambda) l1-penalty function
double l1penalty(arma::mat const& gam, arma::mat const& zv){return accu(gam % abs(zv));}

//////////////////// The weighted (gam = penaltyfactor * lambda) scad-penalty function
double scadpenalty(arma::mat const& gam, double a, arma::mat const& zv){

arma::mat absbeta = abs(zv);

return accu(gam % absbeta % (absbeta <= gam) - (pow(zv, 2) - 2 * a * gam % absbeta + pow(gam, 2)) / (2 * (a - 1)) % (gam < absbeta && absbeta <= a * gam) + (a + 1) * pow(gam, 2) / 2 % (absbeta > a * gam));

}
