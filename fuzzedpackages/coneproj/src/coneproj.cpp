// includes from the plugin
//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP coneACpp( SEXP y, SEXP amat, SEXP face) ;
}

// definition

SEXP coneACpp( SEXP y, SEXP amat, SEXP face){
BEGIN_RCPP

    Rcpp::NumericVector y_(y);
    Rcpp::NumericMatrix amat_(amat);
    int n = y_.size(), m = amat_.nrow();
    arma::mat namat(amat_.begin(), m, n, false);
    arma::colvec ny(y_.begin(), n, false);
    float sm = 1e-8;
    arma::colvec h(m); h.fill(0);
    arma::colvec obs = arma::linspace(0, m-1, m);
    int check = 0;
    arma::mat amat_in = namat;
//new:
    arma::colvec face_ = Rcpp::as<arma::colvec>(face);
    int nf = face_.n_rows;    

    for(int im = 0; im < m; im ++){
        arma::mat nnamat = namat.row(im) * namat.row(im).t();
        amat_in.row(im) = namat.row(im) / sqrt(nnamat(0,0));
    } 

    arma::mat delta = -amat_in;
    arma::colvec b2 = delta * ny;
    arma::colvec theta(n); theta.fill(0);

    if(max(b2) > 2 * sm){
        int i = min(obs.elem(find(b2 == max(b2))));
        h(i) = 1;
//new:
        if(!face_.is_empty()) {
            for (int i = 0; i < nf; i ++) {
                int posi = face_(i) - 1;
                h(posi) = 1;
            }
        }
    }
    else{check = 1;}

    int nrep = 0;
//new:
    arma::mat xmat_use;
    while((check==0) & (nrep<(n*n))){
        nrep ++ ;
        arma::colvec indice = arma::linspace(0, delta.n_rows - 1, delta.n_rows);
        indice = indice.elem(find(h == 1));
        arma::mat xmat(indice.n_elem, delta.n_cols); xmat.fill(0);
  
        for(int k = 0; k < indice.n_elem; k ++){
            xmat.row(k) = delta.row(indice(k));
        }

        arma::colvec a = solve(xmat * xmat.t(), xmat * ny);
        arma::colvec avec(m); avec.fill(0);

        if(min(a) < (-sm)){
            avec.elem(find(h == 1)) = a;
            int i = min(obs.elem(find(avec == min(avec))));
            h(i) = 0;
            check = 0;
        }
      
        else{
            check = 1;
            theta = xmat.t() * a;
            b2 = delta * (ny - theta)/n;

            if(max(b2) > 2 * sm){
                int  i = min(obs.elem(find(b2 == max(b2))));
                h(i) = 1;
                check = 0;
            }
        }
//new:
       xmat_use = xmat;
    }

    //if(nrep > (n * n - 1)){Rcpp::Rcout << "Fail to converge in coneproj!Too many steps! Number of steps:" << nrep << std::endl;}

    return wrap(Rcpp::List::create(Rcpp::Named("thetahat") = ny - theta, Named("xmat") = xmat_use, Named("dim") = n - sum(h), Named("nrep") = nrep, Named("h") = h));

END_RCPP
}


// declarations
extern "C" {
SEXP coneBCpp( SEXP y, SEXP delta, SEXP vmat, SEXP face) ;
}

// definition

SEXP coneBCpp( SEXP y, SEXP delta, SEXP vmat, SEXP face){
BEGIN_RCPP

    Rcpp::NumericVector y_(y);
    Rcpp::NumericMatrix delta_(delta);
    arma::mat nvmat = Rcpp::as<arma::mat>(vmat);
    int n = y_.size(), m = delta_.nrow(), p = nvmat.n_cols;
    arma::colvec ny(y_.begin(), n, false);
    arma::mat ndelta(delta_.begin(), m, n, false);
    arma::mat a;
    arma::mat sigma;
    arma::colvec h;
//new:
    arma::colvec face_ = Rcpp::as<arma::colvec>(face);
    int nf = face_.n_rows;
    arma::colvec obs;
    arma::mat theta(n, 1);

    float sm = 1e-8;
//new: test!
    //float sm = 1e-5;
    int check = 0;

    arma::colvec scalar(m);
    arma::mat delta_in = ndelta;

    for(int im = 0; im < m; im ++){
        arma::mat nndelta = ndelta.row(im) * ndelta.row(im).t();
        scalar(im) = sqrt(nndelta(0,0));
        delta_in.row(im) = ndelta.row(im) / scalar(im);
    } 

    if(nvmat.is_empty()){
       p = p - 1;
       sigma.set_size(m, n);
       sigma = delta_in;
       h.set_size(m); 
       h.fill(0);
       obs.set_size(m);
       obs = arma::linspace(0, m - 1, m);
       theta.fill(0);
    }

    if(!nvmat.is_empty()){ 
        sigma.set_size(m + p, n);
        sigma.rows(0, p - 1) = nvmat.t(); sigma.rows(p, m + p - 1) = delta_in;
        h.set_size(m + p); 
        h.fill(0);
        for(int i = 0; i < p; i ++){
          h(i) = 1;
        }
        obs.set_size(m + p);
        obs = arma::linspace(0, m + p - 1, m + p);
        theta = nvmat * solve(nvmat.t() * nvmat, nvmat.t() * ny);
    } 

    arma::colvec b2 = sigma * (ny - theta) / n;

    if(max(b2) > 2 * sm){
        int i = min(obs.elem(find(b2 == max(b2))));
        h(i) = 1;
//new:
        if(!face_.is_empty()) {
            for (int i = 0; i < nf; i ++) {
                int posi = face_(i) - 1;
                h(posi) = 1;
            }
        }
    }
    
    int nrep = 0;  

    if(max(b2) <= 2 * sm){
        check = 1;
        theta.fill(0);

        if(nvmat.is_empty()){
           a.set_size(m, 1); a.fill(0);
        }

        if(!nvmat.is_empty()){
           a.set_size(p, 1); 
           a = solve(nvmat.t() * nvmat, nvmat.t() * ny);
           theta = nvmat * solve(nvmat.t() * nvmat, nvmat.t() * ny);
        }
        arma::colvec avec(m + p); avec.fill(0);
        if(!nvmat.is_empty()){
          avec.elem(find(h == 1)) = a;
        }
        return wrap(Rcpp::List::create(Named("yhat") = theta, Named("coefs") = avec, Named("nrep") = nrep, Named("dim") = sum(h)));
    }
//new: 
    //int upper = n*n - 1;
    //int upper = 1000;
   // while(check == 0 & nrep < (n * n)){
//double sc = 0;
   while((check==0) & (nrep<1e+6)){
        nrep ++;
        //if(nrep > (n * n)){
          // throw (Rcpp::exception("Fail to converge in coneproj! nrep > n^2 !"));
        //}
        arma::colvec indice = arma::linspace(0, sigma.n_rows-1, sigma.n_rows); 
        indice = indice.elem(find(h == 1));
        arma::mat xmat(indice.n_elem, sigma.n_cols); xmat.fill(0);

        for(int k = 0; k < indice.n_elem; k ++){
            xmat.row(k) = sigma.row(indice(k));
        }
 	//double sc = arma::norm(xmat * xmat.t(), 2);
        a = solve(xmat * xmat.t(), xmat * ny);
//new: 
       if (a.n_elem > p) {
            arma::colvec a_sub(a.n_elem - p);

            for(int i = p; i <= a.n_elem - 1; i ++){
                a_sub(i-p) = a(i);
            }

            if(min(a_sub) < (- sm)){
                arma::colvec avec(m + p); avec.fill(0);
                avec.elem(find(h == 1)) = a;
                arma::colvec avec_sub(m);

                for(int i = p; i <= p + m - 1; i ++){
                    avec_sub(i-p) = avec(i);
                }

                int i = max(obs.elem(find(avec == min(avec_sub))));
                h(i) = 0;
                check = 0;
            }
 
            if(min(a_sub) > (-sm)){
                check = 1;
                theta = xmat.t() * a;
                b2 = sigma * (ny - theta) / n;
//arma::mat sc0 = sqrt(b2.t() * b2);
//sc = sqrt(sc0(0,0));
//sc = arma::as_scalar(b2.t() * b2);
                //if(max(b2) > 2 * sc * sm){
		//if((max(b2) * sc) > 2 * sm){
                if(max(b2) > 2 * sm){
                    int i = min(obs.elem(find(b2 == max(b2))));
                    check = 0;
                    h(i) = 1;
                }
            }
        } else {
            check = 1;
       }
//new: avoid the mismatch problem
       if (nrep == 1e+6) {
            arma::colvec indiceEnd = arma::linspace(0, sigma.n_rows-1, sigma.n_rows); 
       	    indiceEnd = indiceEnd.elem(find(h == 1));
            arma::mat xmat(indiceEnd.n_elem, sigma.n_cols); xmat.fill(0);
            for(int k = 0; k < indiceEnd.n_elem; k ++){
                xmat.row(k) = sigma.row(indiceEnd(k));
            }
            //sc = norm(xmat * xmat.t(), 2);
            a = solve(xmat * xmat.t(), xmat * ny);
            theta = xmat.t() * a;
       }
   }

   arma::colvec avec(m + p); avec.fill(0);
   avec.elem(find(h == 1)) = a;
   arma::colvec avec_orig(m + p); avec_orig.fill(0);
 
   for(int i = 0; i < p; i ++){
      avec_orig(i) = avec(i);
   }
	
   for(int i = p; i < (m + p); i ++){
      avec_orig(i) = avec(i) / scalar(i - p);
   }
   // if(nrep > (n * n - 1)){Rcpp::Rcout << "Fail to converge in coneproj!Too many steps! Number of steps:" << nrep << std::endl;}
   return wrap(Rcpp::List::create(Named("yhat") = theta, Named("coefs") = avec_orig, Named("nrep") = nrep, Named("dim") = sum(h)));

END_RCPP
}


// declarations
extern "C" {
SEXP qprogCpp( SEXP q, SEXP c, SEXP amat, SEXP b, SEXP face) ;
}

// definition

SEXP qprogCpp( SEXP q, SEXP c, SEXP amat, SEXP b, SEXP face){
BEGIN_RCPP

    Rcpp::NumericVector c_(c);
    Rcpp::NumericMatrix q_(q);
    Rcpp::NumericMatrix amat_(amat);
    Rcpp::NumericVector nb(b);
    int n = c_.size(), m = amat_.nrow();
    arma::colvec nc(c_.begin(), n, false);
    arma::mat namat(amat_.begin(), m, n, false);
    arma::mat nq(q_.begin(), n, n, false);
    bool constr = is_true(any( nb != 0 ));
    arma::colvec theta0(n);
    arma::colvec nnc(n);
//new:
    arma::colvec face_ = Rcpp::as<arma::colvec>(face);
    int nf = face_.n_rows;
    
    if(constr){
        arma::colvec b_(nb.begin(), m, false);
        theta0 = solve(namat, b_);
        nnc = nc - nq * theta0;
    } 

    else{nnc = nc;}

    arma::mat preu = chol(nq);
    arma::mat u = trimatu(preu); 
    arma::colvec z = inv(u).t() * nnc;
    arma::mat atil = namat * inv(u);

    float sm = 1e-8;
    arma::colvec h(m); h.fill(0);
    arma::colvec obs = arma::linspace(0, m-1, m);
    int check = 0;

    for(int im = 0; im < m; im ++){
        arma::mat atilnorm = atil.row(im) * atil.row(im).t();
        atil.row(im) = atil.row(im) / sqrt(atilnorm(0,0));
    } 

    arma::mat delta = -atil;
    arma::colvec b2 = delta * z;
    arma::colvec phi(n); phi.fill(0);

    if(max(b2) > 2 * sm){
        int i = min(obs.elem(find(b2 == max(b2))));
        h(i) = 1;
//new:
        if(!face_.is_empty()) {
            for (int i = 0; i < nf; i ++) {
                int posi = face_(i) - 1;
                h(posi) = 1;
            }
        }
    }

    else{check = 1;}

    int nrep = 0;
//new:
    arma::mat xmat_use;
    while((check==0) & (nrep<(n*n))){
        nrep ++ ;
       // if(nrep > (n * n)){
         //  throw (Rcpp::exception("Fail to converge in coneproj! nrep > n^2 !"));}
        arma::colvec indice = arma::linspace(0, delta.n_rows - 1, delta.n_rows);
        indice = indice.elem(find(h == 1));
        arma::mat xmat(indice.n_elem, delta.n_cols); xmat.fill(0);
  
        for(int k = 0; k < indice.n_elem; k ++){
        xmat.row(k) = delta.row(indice(k));
        }

        arma:: colvec a = solve(xmat * xmat.t(), xmat * z);
        arma:: colvec avec(m); avec.fill(0);

        if(min(a) < (-sm)){
            avec.elem(find(h == 1)) = a;
            int i = min(obs.elem(find(avec == min(avec))));
            h(i) = 0;
            check = 0;
        }
      
        else{
            check = 1;
            phi = xmat.t() * a;
            b2 = delta * (z - phi)/n;
    
            if(max(b2) > 2 * sm){
                int  i = min(obs.elem(find(b2 == max(b2))));
                h(i) = 1;
                check = 0;
            }
        }
	//new:
    	xmat_use = xmat;
    }

    arma::colvec thetahat = solve(u, z - phi);

    if(constr){
        thetahat = thetahat + theta0;
    }

    // if(nrep > (n * n - 1)){Rcpp::Rcout << "Fail to converge in qprog!Too many steps! Number of steps:" << nrep << std::endl;}

    return wrap(Rcpp::List::create(Rcpp::Named("thetahat") = thetahat, Named("xmat") = xmat_use, Named("dim") = n - sum(h), Named("nrep") = nrep, Named("h") = h));

END_RCPP
}


