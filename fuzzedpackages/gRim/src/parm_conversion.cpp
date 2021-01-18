#include <RcppArmadillo.h>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
List parm_ghk2pms_ (List parms){

  int ii;
  
  List out = clone( parms );
  CharacterVector out_names = out.names();
  out_names[0] = "p";
  out_names[1] = "mu";
  out_names[2] = "Sigma";
  out.names() = out_names;
  //Rf_PrintValue( out_names );
  
  std::string gentype = as<std::string>(parms["gentype"]); 
  
  if (strcmp(gentype.c_str(),"discrete")==0){
    //Rcpp::Rcout << "generator: discrete" << std::endl;
    
    NumericVector gg3 = parms[0]; 
    NumericVector gg2 = clone(gg3);
    
    gg3 = exp( gg3 - mean( gg3 ) );
    gg3 = gg3 / sum( gg3 );
    
    for (ii=0; ii<gg2.size(); ii++){ 
      gg2(ii) = gg3(ii); 
    }
    
    out[0] = wrap(gg2);
    out[1] = R_NilValue;
    out[2] = R_NilValue;
    
  } else {
    arma::mat QQ = parms["K"];
    arma::mat LL = parms["h"];
    arma::mat QQ_out = inv_sympd( QQ );
    arma::mat LL_out = QQ_out * LL;
    int ndd = LL_out.n_cols, Q =QQ_out.n_rows;
    rowvec DD_out  = rowvec(ndd);
    out[1] = wrap(LL_out);
    out[2] = wrap(QQ_out);
    
    if (strcmp(gentype.c_str(),"mixed")==0){
      //Rcpp::Rcout << "generator: mixed" << std::endl;      
      vec dd = vec(Q);
      arma::vec DD = as<arma::vec>(parms["g"]);
      arma::mat tmp = ( QQ_out * LL ); 
      arma::mat qmat = LL % tmp; 
      rowvec quad2 = sum(qmat, 0); // same as colSums
      
      for (ii=0; ii<ndd; ii++){
	quad2(ii) = DD(ii) + quad2(ii)/2;
      }
      double mm2 = mean(quad2);
      
      DD_out = exp(quad2-mm2);
      DD_out = DD_out / sum(DD_out);
      
      NumericVector gg1 = parms["g"];
      NumericVector gg2 = clone(gg1);
      for (ii=0; ii<ndd; ii++){ 
	gg2(ii) = DD_out(ii); 
      }
      
      out[0] = wrap(gg2);
    } else {
      out[0] = 1;
    } 
  }
  return out;
}


//[[Rcpp::export]]
List parm_pms2ghk_ ( List parms ){
  using namespace arma;
  using namespace Rcpp;
  // Rcpp::List parms(parms_);

  List out = clone( parms );
  CharacterVector out_names = out.names();
  out_names[0] = "g";
  out_names[1] = "h";
  out_names[2] = "K";
  out.names() = out_names;
  //Rf_PrintValue( out_names );
  
  std::string gentype = as<std::string>(parms["gentype"]); 
  if (strcmp(gentype.c_str(),"discrete")==0){
    //Rcpp::Rcout << "generator: discrete" << std::endl;
    vec DD = as<vec>(parms["p"]);
    int ii, ndd = DD.n_elem;
    NumericVector gg1 = parms["p"];
    NumericVector gg2 = clone(gg1);
    for (ii=0; ii<ndd; ii++){ 
      gg2(ii) = log(DD(ii)); 
    }
    
    out[0] = wrap(gg2);
    out[1] = R_NilValue;
    out[2] = R_NilValue;
    
  } else {
    arma::mat QQ = as<arma::mat>(parms["Sigma"]);
    arma::mat LL = as<arma::mat>(parms["mu"]);
    //arma::mat QQ_out  = inv( sympd(QQ) );
    arma::mat QQ_out  = inv_sympd( QQ );
    arma::mat LL_out  = QQ_out * LL;
    int ii, ndd = LL_out.n_cols, Q =QQ.n_rows;
    double cst  = log(det(QQ_out)) - Q*log(2*datum::pi);
    
    if (strcmp(gentype.c_str(),"mixed")==0){
      //Rcpp::Rcout << "generator: mixed" << std::endl;      
      vec dd = vec(Q);
      vec DD = as<arma::vec>(parms["p"]);
      vec gg  = vec(ndd);
      
      for (ii=0; ii<ndd; ii++){
	dd = LL.col(ii);
	gg(ii)   = log(DD(ii)) + (cst - as_scalar(dd.t() * QQ_out * dd)) / 2;
      }
      NumericVector gg1 = parms["p"];
      NumericVector gg2 = clone(gg1);
      
      // FIXME: Not necessary :
      for (ii=0; ii<ndd; ii++){ gg2(ii) = gg(ii); }
      out[0] = wrap(gg2);
      
      out[1] = wrap(LL_out);
      out[2] = wrap(QQ_out);			
    } else { 
      //Rcpp::Rcout << "generator: continuous" << std::endl;      	
      double gg2;
      gg2 = (cst - as_scalar(LL.t() * QQ_out * LL))/2;
      
      out[0] = wrap(gg2);
      out[1] = wrap(LL_out);
      out[2] = wrap(QQ_out);
    }
  }
  return out;
}

//[[Rcpp::export]]
arma::mat updateA ( arma::mat A, arma::mat E, arma::uvec row, arma::uvec col){
  
  Rf_PrintValue(wrap(A));	Rf_PrintValue(wrap(E));
  Rf_PrintValue(wrap(row)); 	Rf_PrintValue(wrap(col));
  
  //A[row-1,col-1] = 	A[row-1,col-1] + E;
  A(row-1,col-1) = A(row-1,col-1) + E;
  return A;
}


//[[Rcpp::export]]
List parm_update_ghk_ (List Cparms, IntegerVector dgen_idx, IntegerVector cgen_idx, 
		       List ghk_obs, List pms_obs, List ghk_fit, List pms_fit, double scale, double details=0) {
  
  Environment stats("package:stats");
  Function rnorm = stats["rnorm"];
  
  Environment gRbase("package:gRbase");
  Function R_tabAdd=gRbase["tabAdd"];
  List out = clone( Cparms );
  
  std::string gentype = as<std::string>(ghk_fit["gentype"]); 
  if (strcmp(gentype.c_str(),"discrete")==0){
    Rcout << "gentype=discrete" << std::endl;
    NumericVector upd_g  = scale * ((NumericVector) ghk_obs[ 0 ] - (NumericVector) ghk_fit[ 0 ]);
    NumericVector g_new  = R_tabAdd(Cparms[ 0 ], upd_g);
    //NumericVector g.new  = tableOp2(Cparms[[ g.idx ]], upd.g, `+`, restore=TRUE);
    //List res = List::create(g=g_new, h=Cparms[1], K=Cparms[2], gentype="discrete");
    // List res = List::create(Named("g", g_new), Named("h", Cparms[1]), 
    // 												Named("K", Cparms[2]), Named("gentype", "discrete"));
    
    // Not necessary; just stick into the clone
    List res = List::create(_["g"] = g_new, _["h"] = Cparms[1], _["K"]=Cparms[2], _["gentype"]="discrete");
    
  } else {
    if (strcmp(gentype.c_str(),"continuous")==0){
      Rcout << "gentype=continuous" << std::endl;
      
      arma::mat h_new    = Cparms[ 1 ];
      arma::mat ghk_obs_ = ghk_obs[1];
      arma::mat ghk_fit_ = ghk_fit[1];
      arma::uvec cgen_idx_ = as<arma::uvec>(cgen_idx);
      
      mat upd_h    = scale * (ghk_obs_ - ghk_fit_);
      Rf_PrintValue(wrap(upd_h));
      Rf_PrintValue(wrap(cgen_idx));
      
      // Rf_PrintValue(wrap(h_new));
      // for (size_t j=0; j<h_new.n_cols; ++j){
      // 	h_new(cgen_idx_-1, j) =  h_new(cgen_idx_-1, j) + upd_h;
      // }
      // Rf_PrintValue(wrap(h_new));
      
      //NumericMatrix upd_h = scale * ((NumericMatrix) ghk_obs[1] - (NumericMatrix) ghk_fit[1]);
      // for (j in 1:ncol(h.new))
      //     h.new[cgen.idx, j] <- Cparms[[ h.idx ]][cgen.idx, j, drop=FALSE] + upd.h
      
      // upd.k   <- scale * (ghk.obs[[ K.idx ]] - ghk.fit[[ K.idx ]])
      // K.new   <- Cparms[[ K.idx ]]
      // K.new[cgen.idx, cgen.idx] <- K.new[cgen.idx, cgen.idx] + upd.k
      // ## cat("cont: upd.h:\n"); print(cbind(ghk.obs[[ h.idx ]], ghk.fit[[ h.idx ]], upd.h))
      // ## cat("cont: upd.k:\n"); print(cbind(ghk.obs[[ K.idx ]], ghk.fit[[ K.idx ]], upd.k))
      // res <- list(g=Cparms[[ g.idx ]], h=h.new, K=K.new, gentype="continuous")
      
      
    } else {
      Rcout << "gentype=mixed" << std::endl;
    }
  }
  return out;
}





//[[Rcpp::export]]
RcppExport SEXP C_pms2ghk ( SEXP parms_ ){
  using namespace arma;
  using namespace Rcpp;
  Rcpp::List parms(parms_);
  //Rcout << parms_ << std::endl;
  std::string gentype = as<std::string>(parms["gentype"]); 
  if (strcmp(gentype.c_str(), "discrete")==0){
    //Rcpp::Rcout << "generator: discrete" << std::endl;
    vec DD = as<vec>(parms["p"]);
    int ii, ndd = DD.n_elem;
    NumericVector gg1 = parms["p"];
    NumericVector gg2 = clone(gg1);
    for (ii=0; ii<ndd; ii++){ gg2(ii) = log(DD(ii)); }
    
    return List::create(Named("g", wrap(gg2)),  Named("h", R_NilValue),
			Named("K", R_NilValue), Named("gentype", "discrete"));
  } else {
    arma::mat QQ = as<arma::mat>(parms["Sigma"]);
    arma::mat LL = as<arma::mat>(parms["mu"]);
    //arma::mat QQ_out  = inv( sympd(QQ) );
    arma::mat QQ_out = inv_sympd( QQ );
    arma::mat LL_out  = QQ_out * LL;
    int ii, ndd = LL_out.n_cols, Q =QQ.n_rows;
    double cst  = log(det(QQ_out)) - Q*log(2*datum::pi);
    
    if (strcmp(gentype.c_str(), "mixed")==0){
      //Rcpp::Rcout << "generator: mixed" << std::endl;      
      vec dd = vec(Q);
      vec DD = as<arma::vec>(parms["p"]);
      vec gg  = vec(ndd);
      
      for (ii=0; ii<ndd; ii++){
	dd = LL.col(ii);
	gg(ii)   = log(DD(ii)) + (cst - as_scalar(dd.t() * QQ_out * dd)) / 2;
      }
      NumericVector gg1 = parms["p"];
      NumericVector gg2 = clone(gg1);
      for (ii=0; ii<ndd; ii++){ gg2(ii) = gg(ii); }
      //      Rcout << "gg1 :" << gg1 << std::endl;
      //Rcout << "gg2 :" << gg2 << std::endl;
      return List::create(Named("g", wrap(gg2)),    Named("h", wrap(LL_out)),
			  Named("K", wrap(QQ_out)), Named("gentype", gentype.c_str()));	  	
    } else { 
      //Rcpp::Rcout << "generator: continuous" << std::endl;      	
      double gg2;
      gg2 = (cst - as_scalar(LL.t() * QQ_out * LL))/2;
      return List::create(Named("g", wrap(gg2)),    Named("h", wrap(LL_out)),
			  Named("K", wrap(QQ_out)), Named("gentype", gentype.c_str()));	  
    }
  }
}


//[[Rcpp::export]]
RcppExport SEXP C_ghk2pms ( SEXP parms_ ){
  using namespace arma;
  using namespace Rcpp;
  Rcpp::List parms(parms_);
  std::string gentype = as<std::string>(parms["gentype"]); 
  if (strcmp(gentype.c_str(),"discrete")==0){
    //Rcpp::Rcout << "generator: discrete" << std::endl;
    vec DD = as<arma::vec>(parms["g"]);
    int ii, ndd = DD.n_elem;
    rowvec DD_out = rowvec(ndd);
    double mm = mean(DD);
    for (ii=0; ii<ndd; ii++) {DD_out(ii) = exp(DD(ii)-mm);}
    double ss = sum(DD_out);
    for (ii=0; ii<ndd; ii++) {DD_out(ii) = DD_out(ii)/ss;}
    
    NumericVector gg1 = parms["g"];
    NumericVector gg2 = clone(gg1);
    for (ii=0; ii<ndd; ii++){ gg2(ii) = DD_out(ii); }

    return List::create(Named("p", wrap(gg2)),      Named("mu", R_NilValue),
			Named("Sigma", R_NilValue), Named("gentype", "discrete"));
  } else {
    arma::mat QQ = as<arma::mat>(parms["K"]);
    arma::mat LL = as<arma::mat>(parms["h"]);
    //arma::mat QQ_out = inv( sympd(QQ) );
		arma::mat QQ_out = inv_sympd( QQ );
    arma::mat LL_out = QQ_out * LL;
    int ii, ndd = LL_out.n_cols, Q =QQ_out.n_rows;
    rowvec DD_out  = rowvec(ndd);

    if (strcmp(gentype.c_str(),"mixed")==0){
      //Rcpp::Rcout << "generator: mixed" << std::endl;      
      vec dd = vec(Q);
      arma::vec DD = as<arma::vec>(parms["g"]);

      vec quad = vec(ndd);
      
      for (ii=0; ii<ndd; ii++){
	dd   = LL.col(ii);
	quad(ii) = DD(ii) + as_scalar(dd.t() * QQ_out * dd)/2;
      }
      double mm = mean(quad);
      for (ii=0; ii<ndd; ii++){ DD_out(ii) = exp(quad(ii)-mm); }
      double ss = sum(DD_out);
      for (ii=0; ii<ndd; ii++) {DD_out(ii) = DD_out(ii)/ss;}

      NumericVector gg1 = parms["g"];
      NumericVector gg2 = clone(gg1);
      for (ii=0; ii<ndd; ii++){ gg2(ii) = DD_out(ii); }
      return Rcpp::List::create(Rcpp::Named("p", wrap(gg2)),
				Rcpp::Named("mu", wrap(LL_out)),
				Rcpp::Named("Sigma", wrap(QQ_out)),
				Rcpp::Named("gentype", gentype.c_str())	);	
    } else {
      //Rcpp::Rcout << "generator: continuous" << std::endl;      	
      return Rcpp::List::create(Rcpp::Named("p", 1),
				Rcpp::Named("mu", wrap(LL_out)),
				Rcpp::Named("Sigma", wrap(QQ_out)),
				Rcpp::Named("gentype", gentype.c_str())	);	
    } 
  }
}























//[[Rcpp::export]]
List parm_normalize_ghk_(List parms){

	NumericVector gg= parms[0]; // g

	mat h2 = parms[1];
	mat K2 = parms[2];
	int Q  = K2.n_rows;
	double logdetK, sign;
	log_det(logdetK, sign, K2);

	mat mu = inv_sympd( K2 ) * h2;
	mat qmat = h2 % mu ;
	rowvec quad = sum(qmat, 0); // same as colSums

	NumericVector quad2 = wrap(quad);
	NumericVector zzz   = gg + quad2/2;
	NumericVector ppp   = exp(zzz - mean(zzz));
	NumericVector pppn  = ppp / sum( ppp );
	NumericVector gnew = log(pppn) + (logdetK - Q*log(2*3.141593) - quad2)/2;
	
	gnew.attr("dim")      = gg.attr("dim");
	gnew.attr("dimnames") = gg.attr("dimnames");

	List out=clone(parms);
	out[0] = gnew;
	return out;
} 






