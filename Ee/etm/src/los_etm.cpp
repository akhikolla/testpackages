// diff start at .begin() + 1

#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

using namespace arma;

vec my_diff(const vec bouh);

/*
  Length of stay without competing outcomes 
*/
RcppExport SEXP los_nocp(SEXP __times,
			 SEXP __tr_mat,
			 SEXP __tau) {

    Rcpp::NumericVector _times(__times), _tr_mat(__tr_mat);
    Rcpp::IntegerVector Dims = _tr_mat.attr("dim");
    const double tau = Rcpp::as<double>(__tau);
    const int lt = _times.size();
    const int nstate = Dims[0];
    cube tr_mat(_tr_mat.begin(), nstate, nstate, lt);
    vec times(_times.begin(), _times.size(), false);
    
    vec T(times);
    T.resize(lt+1); T[lt] = tau;

    // new stuffs we'll need
    mat::fixed<3, 3> I;
    I.eye();
    
    cube aj(nstate, nstate, lt);

    for (int i = 0; i < lt; ++i) {
	tr_mat.slice(i) = tr_mat.slice(i) + I;
	aj.slice(i).eye();
    }

    vec los0(lt), los1(lt), a00(lt), a11(lt), a01(lt);
    a00.zeros();
    a01.zeros();
    a11.zeros();
    
    los0.fill(tau);
    los1.fill(tau);

    for (int t = (times.size() - 2); t >= 0; --t) {

	vec dd = my_diff(T(span(t+1, lt)));
	
	for (int j = t; j < lt; ++j) {
	    aj.slice(j) = aj.slice(j) * tr_mat.slice(t+1);
	    a00(j) = aj(0, 0, j);
	    a01(j) = aj(0, 1, j);
	    a11(j) = aj(1, 1, j);
	}

	colvec a = a00(span(t, lt - 2)) + a01(span(t, lt - 2));
	colvec b = a11(span(t, lt - 2));
	
	los0[t] = T[t+1] + as_scalar(dd.t() * a);
	los1[t] = T[t+1] + as_scalar(dd.t() * b);
    }

    return Rcpp::List::create(Rcpp::Named("los0") = los0,
			      Rcpp::Named("los1") = los1);

}



/* 
   Length of stay with competing outcomes
*/
RcppExport SEXP los_cp(SEXP __times,
		       SEXP __tr_mat,
		       SEXP __tau) {

    Rcpp::NumericVector _times(__times), _tr_mat(__tr_mat);
    Rcpp::IntegerVector Dims = _tr_mat.attr("dim");
    const double tau = Rcpp::as<double>(__tau);
    const int lt = _times.size();
    const int nstate = Dims[0];
    cube tr_mat(_tr_mat.begin(), nstate, nstate, lt);
    vec times(_times.begin(), _times.size(), false);
    
    vec T(times);
    T.resize(lt+1); T[lt] = tau;

    // new stuffs we'll need
    mat::fixed<4, 4> I;
    I.eye();
    
    cube aj(nstate, nstate, lt);

    for (int i = 0; i < lt; ++i) {
	tr_mat.slice(i) = tr_mat.slice(i) + I;
	aj.slice(i).eye();
    }

    vec los0(lt), los1(lt), a00(lt), a11(lt), a01(lt), a12(lt), a13(lt),
	phi2case(lt), phi2control(lt),
	phi3case(lt), phi3control(lt);
    
    a00.zeros();
    a01.zeros();
    a11.zeros();
    phi2case.zeros();
    phi2control.zeros();
    phi3case.zeros();
    phi3control.zeros();
    
    los0.fill(tau);
    los1.fill(tau);
    
    for (int t = (times.size() - 2); t >= 0; --t) {
	
	vec dd = my_diff(T(span(t+1, lt)));
	
	for (int j = t; j < lt; ++j) {
	    aj.slice(j) = tr_mat.slice(t+1) * aj.slice(j);
	    a00(j) = aj(0, 0, j);
	    a01(j) = aj(0, 1, j);
	    a11(j) = aj(1, 1, j);
	    a12(j) = aj(1, 2, j);
	    a13(j) = aj(1, 3, j);
	}
	
	colvec a = a00(span(t, lt - 2)) + a01(span(t, lt - 2));
	colvec b = a11(span(t, lt - 2));
	colvec c = a12(span(t, lt - 2));
	colvec d = a13(span(t, lt - 2));
	
	los0[t] = T[t+1] + as_scalar(dd.t() * a);
	los1[t] = T[t+1] + as_scalar(dd.t() * b);

	phi2case[t] = T[lt - 1]  * aj(1, 2, lt - 1) - as_scalar(dd.t() * c);
	phi3case[t] = T[lt - 1]  * aj(1, 3, lt - 1) - as_scalar(dd.t() * d);

	phi2control[t] = aj(1, 2, lt - 1) * los0[t];	
	phi3control[t] = aj(1, 3, lt - 1) * los0[t];
    }

    return Rcpp::List::create(Rcpp::Named("los0") = los0,
			      Rcpp::Named("los1") = los1,
			      Rcpp::Named("phi2control") = phi2control,
			      Rcpp::Named("phi3control") = phi3control,
			      Rcpp::Named("phi2case") = phi2case,
			      Rcpp::Named("phi3case") = phi3case);
    
}


    
