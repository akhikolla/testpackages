//    dpparmadillo.cpp: C++ functions implementing determinantal point process 
//         based sampling of (approximately) optimal designs for GP regression.
//    Copyright (C) 2018  Matthew T. Pratola <mpratola@stat.osu.edu>
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as published
//    by the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.
//
//    You should have received a copy of the GNU Affero General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
sp_mat sqrt_(sp_mat X) 
{ 
	return sqrt(X);
}

// [[Rcpp::export]]
mat subspace_(mat V)
{
	double scal1,scal2;

	if(V.is_empty()) { Rprintf("Empty matrix in subspace_!\n"); return(V); }

/* Modified Gram-Schmidt procedure */
	for(uword i=0;i<V.n_cols;i++) {
		scal1=((V.col(i)).t()*V.col(i)).eval()(0,0);
		scal1=sqrt(scal1);
		V.col(i)=V.col(i)/scal1;
		for(uword j=i+1;j<V.n_cols;j++) {
			scal2=((V.col(i)).t()*V.col(j)).eval()(0,0);
			V.col(j)=V.col(j)-scal2*V.col(i);
		}
	}

/* Gram-Schmidt procedure
	// Orthogonal basis for Vnew orthogonal to v.orthto
	for(uword i=0;i<V.n_cols;i++) {
		for(uword j=0;j<i;j++) {
			scal1=((V.col(i)).t()*V.col(j)).eval()(0,0);
			scal2=((V.col(j)).t()*V.col(j)).eval()(0,0);
			V.col(i)=V.col(i)-scal1/scal2*V.col(j);
		}
	}

	// Normalize
	for(uword i=0;i<V.n_cols;i++){
		scal1=((V.col(i)).t()*V.col(i)).eval()(0,0);
		scal1=sqrt(scal1);
		V.col(i)=V.col(i)/scal1;
	} 
*/
	return(V);
 }

// [[Rcpp::export]]
vec simDppModal_(sp_mat R,uword n)
{
	vec vals,p;
	mat V;
//	vec pts(0);

	if(n<=0) { Rprintf("Invalid n!\n"); return(NULL); }
	if(R.is_empty()) { Rprintf("Invalid R!\n"); return(NULL); }

	// sparse eigendecomposition for first n evecs/evals
	eigs_sym(vals, V, R, n);

	// Is this annoying? Why yes, yes it is.
	// (eigs_sym returns things in reverse order...)
	V=fliplr(V);
	vals=flipud(vals);
	
	// calculate p
	p=vals/(1.0+vals);
	// initialize N and pts
	uword N=V.n_rows;	
	vec pts(n);
	vec pvec(N);
	vec v_orthto(N);
	uword p_dx;

	// main loop
	for(uword i=0;i<(n-1);i++)
	{
		pvec=sum(square(V),1)/(n-i+1.0);
		p_dx=pvec.index_max();
		pts(i)=p_dx;

		v_orthto=V.col(0);
		V.shed_col(0);
		V=V-v_orthto*V.row(p_dx)/v_orthto(p_dx);
		V=subspace_(V);
	}
	pvec=sum(square(V),1)/2.0; //(n-i+1.0);
	p_dx=pvec.index_max();
	pts(n-1)=p_dx;

	return(pts);
}


// // [[Rcpp::export]]
// vec simDppModalEig_(vec vals, mat& V,uword n)
// {
// 	vec p;

// 	if(n<=0) Rprintf("Invalid n!\n");
	
// 	// calculate p
// 	p=vals/(1.0+vals);
// 	// initialize N and pts
// 	uword N=V.n_rows;
// 	vec pts(n);
// 	vec pvec(N);
// 	vec v_orthto(N);
// 	uword p_dx;

// 	// main loop
// 	for(uword i=0;i<n;i++)
// 	{
// 		pvec=sum(square(V),1)/(n-i+1.0);
// 		p_dx=pvec.index_max();
// 		pts(i)=p_dx;

// 		v_orthto=V.col(0);
// 		V.shed_col(0);
// 		V=V-v_orthto*V.row(p_dx)/v_orthto(p_dx);

// 		if(i<(n-1))
// 			V=subspace_(V);
// 	}

// 	return(pts);
// }
