// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
 
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do 
// 
// [[Rcpp::depends(RcppArmadillo)]] 

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
#include "losses.h"
#include "transformations.h"

#define TYPE_PROJ_NO  0
#define TYPE_PROJ_YES  1
#define TYPE_PROJ_AUTO  2
#define TYPE_PREDICT_DEFAULT 0
#define TYPE_PREDICT_RESPONSE 1
#define TYPE_PREDICT_PROB 2
#define TYPE_PREDICT_CLASS 3

#define FISTA_MAX_ITER_INNER 50000
#define FISTA_TOL 1e-7
#define FISTA_T0 2
#define FISTA_STEP 0.1

// Soft-threshold operator
arma::vec soft_thresh(const arma::vec & z, double l){
  arma::vec S = arma::vec(z.n_elem);
  for(unsigned int i=0; i<z.n_elem; i++){
    if(std::abs(z[i]) <= l){
      S(i) = 0;
    }else{
      if(z[i]<=0){
        S[i] = z[i] + l;
      }else{
        S[i] = z[i] - l;
      }
    }
  }
  return S;
}

class s2net
{
private:
    arma::mat xL;
    arma::mat xU;
    arma::mat yL;
    arma::mat T;
    double lambda1;
    double lambda2;
    double gamma1;
    double gamma2;
    double gamma3;
    int p;
    int nL;
    int nU;
    Loss R;
    Loss Rkern;
    arma::vec beta;
    double intercept;
    double mean_yL;
    int loss;
    int proj;
    int frame;
    arma::vec s_scale;
    arma::vec s_center;
    arma::vec rm_cols;

    bool use_warmstart;
    int Arg_FISTA_MAX_ITER_INNER;
    double Arg_FISTA_TOL;
    double Arg_FISTA_T0;
    double Arg_FISTA_STEP;
    bool use_default;

public:
    s2net(const Rcpp::List s2Data, int loss);
    void fit_fast();
    void fit(const arma::vec & params, int frame, int proj);
    void setupFista(const Rcpp::List & s2Fista);
    void optimizeFista();
    void optimizeFista_user();
    arma::vec predict_response(const arma::mat & newX);
    arma::vec predict_probability(const arma::mat & newX);
    arma::vec predict_class(const arma::mat & newX);
    arma::vec predict(const arma::mat & newX, int type);

    double L(arma::vec beta){return R.R(beta, xL, yL) + gamma1*Rkern.R2(beta, T, mean_yL);}
    arma::vec gradL(arma::vec beta){return R.gradR(beta, xL, yL) + gamma1*Rkern.gradR2(beta, T, mean_yL);}
    arma::vec Update(arma::vec beta, arma::vec gradL_beta, double t);

    // get/set 
    arma::vec get_beta(){return beta;}
    void set_beta(arma::vec beta){this->beta = beta;}
    double get_intercept(){return intercept;}
    void set_intercept(double intercept){this->intercept = intercept;}
};


s2net::s2net(const Rcpp::List s2Data, int loss)
{
    this->xL = Rcpp::as<arma::mat>(s2Data["xL"]);
    this->yL = Rcpp::as<arma::mat>(s2Data["yL"]);
    if(Rf_isNull(s2Data["xU"])){
        this->xU = xL.rows(1,2);
    }else{
        this->xU = Rcpp::as<arma::mat>(s2Data["xU"]);
    }
    this->s_scale = Rcpp::as<arma::vec>(s2Data.attr("pr:scale"));
    this->s_center = Rcpp::as<arma::vec>(s2Data.attr("pr:center"));
    this->rm_cols = Rcpp::as<arma::vec>(s2Data.attr("pr:rm_cols"));
    // this->xL = s2Data["xL"];
    // this->yL = s2Data["yL"];
    // this->xU = s2Data["xU"];
    // this->s_scale = s2Data.attr("pr:scale");
    // this->s_center = s2Data.attr("pr:center");
    // this->rm_cols = s2Data.attr("pr:rm_cols");

    this->loss = loss;


    p = xL.n_cols;
    nL = xL.n_rows;
    nU = xU.n_rows;

    this->R = Loss(loss);
    this->Rkern = Loss(loss);

    // Default values
    lambda1 = 0;
    lambda2 = 0;
    gamma1 = 0;
    gamma2 = 0;
    gamma3 = 0;

    use_warmstart = FALSE;
    use_default = TRUE;

    beta = arma::zeros(p);

    frame = TYPE_TRANSFORM_ExtJT;
    proj = TYPE_PROJ_AUTO;
    
    // Compute intercept... or not
    switch (loss)
    {
    case LOSS_LINEAR:
        intercept = arma::mean(arma::mean(yL));
        yL = yL - intercept;
        mean_yL = 0;
        break;
    case LOSS_LOGIT:
        intercept = 0;
        mean_yL = arma::mean(arma::mean(yL));
        break;
    default:
        break;
    }


}

arma::vec s2net::Update(arma::vec beta, arma::vec gradL_beta, double t){
    // arma::vec S1 = soft_thresh(beta - t*gradL_beta, t*lambda1);
    // return S1/(1 + 2*t*lambda2);
    return soft_thresh(beta - t*gradL_beta, t*lambda1)/(1 + 2*t*lambda2);
}

void s2net::setupFista(const Rcpp::List & s2Fista){
    Arg_FISTA_MAX_ITER_INNER = s2Fista["MAX_ITER_INNER"];
    Arg_FISTA_TOL = s2Fista["TOL"];
    Arg_FISTA_T0 = s2Fista["t0"];
    Arg_FISTA_STEP = s2Fista["step"];
    use_warmstart = s2Fista["use_warmstart"];
    use_default = FALSE;
}

void s2net::fit(const arma::vec & params, int frame, int proj){
    // Update params
    lambda1 = params[0];
    lambda2 = params[1];
    gamma1 = params[2];
    gamma2 = params[3];
    gamma3 = params[4];

    this->frame = frame;
    this->proj = proj;

    // Should we transform the unlabeled data before projecting?
    if(proj != TYPE_PROJ_NO){
            arma::rowvec u = arma::mean(xU, 0);
            double angle = 1.0;
            arma::vec projection = - R.gradR(arma::zeros(p), xL, yL);
            projection = projection/arma::norm(projection, 2);

            if(proj == TYPE_PROJ_AUTO){
                angle = std::acos(arma::as_scalar(u*projection/arma::norm(u, 2)));
            }

            if(proj == TYPE_PROJ_YES || std::abs(std::cos(angle)) > std::sqrt(0.5)){
                xU = xU - arma::ones(nU)*u*projection*projection.t();
            }            
          }
    // Compute the projection
    switch (frame)
          {
          case TYPE_TRANSFORM_JT:
            T = transform_JT(xU, gamma2);
            break;
          default:
            T = transform_ExtJT(xU, gamma2, gamma3);
            break;
          }

    // is warm start used?
    if(!use_warmstart){
        beta = arma::zeros(p);
    }
    // call FISTA
    // Very important!!!
    // Should we use the default version or not?
    if (use_default)
    {
        // This one is compiled with #define stuff that makes it... faster?
        optimizeFista(); 
    }else{
        optimizeFista_user();
    }
}

void s2net::optimizeFista(){
    double t = FISTA_T0;
    double l_new = 1;
    double l_old = 1;
    //int iter = 0;
    
    arma::vec theta_new = beta;
    arma::vec theta_old;
    arma::vec g;
    
    
    double L_beta_new = L(beta);
    double L_beta = L_beta_new;
    // Null model's loss
    double L_null = L(arma::zeros(p));
    
    // Check that the optimal beta is not zero (global)
    if (arma::abs(gradL(arma::zeros(p))).max() <= lambda1)
    {
        beta = arma::zeros(p);
    }else{
        for (int iter = 0; iter < FISTA_MAX_ITER_INNER; iter++)
        {
        // This is U_{t_k-1}(\bbeta_{k-1})
        theta_old = theta_new;  
        l_old = l_new;
        
        // Compute the gradient in beta_k
        g = gradL(beta);
        
        // Compute the risk in beta_k
        L_beta = L_beta_new;
        //std::printf("Risk: %.6f\n", L_beta);
        
        //Find t such that R_updated <= R + t(g)*(beta_updated - beta) + 1/2t||beta_updated - beta||_2^2
        theta_new = Update(beta, g, t);
        while (  L(theta_new) > L_beta + arma::as_scalar(g.t()*(theta_new - beta)) + 1/(2*t)*arma::sum(arma::square(theta_new - beta))){
            t = FISTA_STEP*t;
            theta_new = Update(beta, g, t);
        }
        
        // Compute the aceleration term
        l_new = (1 + sqrt(1 + 4*pow(l_old, 2)))/2;
        
        // Update beta
        beta = theta_new + (l_old - 1)*(theta_new - theta_old)/l_new;
        iter ++;
        
        // Compute the new risk
        L_beta_new = L(beta);

        if(std::abs(L_beta_new - L_beta) < FISTA_TOL * L_null){
            break;
        }
        }
    }
     
}

void s2net::optimizeFista_user(){
    double t = Arg_FISTA_T0;
    double l_new = 1;
    double l_old = 1;
    //int iter = 0;
    
    arma::vec theta_new = beta;
    arma::vec theta_old;
    arma::vec g;
    
    
    double L_beta_new = L(beta);
    double L_beta = L_beta_new;
    // Null model's loss
    double L_null = L(arma::zeros(p));
    
    // Check that the optimal beta is not zero (global)
    if (arma::abs(gradL(arma::zeros(p))).max() <= lambda1)
    {
        beta = arma::zeros(p);
    }else{
    for (int iter = 0; iter < Arg_FISTA_MAX_ITER_INNER; iter++)
    {
        // This is U_{t_k-1}(\bbeta_{k-1})
        theta_old = theta_new;  
        l_old = l_new;
        
        // Compute the gradient in beta_k
        g = gradL(beta);
        
        // Compute the risk in beta_k
        L_beta = L_beta_new;
        //std::printf("Risk: %.6f\n", L_beta);
        
        //Find t such that R_updated <= R + t(g)*(beta_updated - beta) + 1/2t||beta_updated - beta||_2^2
        theta_new = Update(beta, g, t);
        while (arma::as_scalar(L(theta_new)) > 
                          arma::as_scalar(L_beta + g.t()*(theta_new - beta) + 1/(2*t)*arma::sum(arma::square(theta_new - beta)))){
            t = Arg_FISTA_STEP*t;
            theta_new = Update(beta, g, t);
        }
        
        // Compute the aceleration term
        l_new = (1 + sqrt(1 + 4*pow(l_old, 2)))/2;
        
        // Update beta
        beta = theta_new + (l_old - 1)*(theta_new - theta_old)/l_new;
        iter ++;
        
        // Compute the new risk
        L_beta_new = L(beta);

        if(std::abs(L_beta_new - L_beta) < Arg_FISTA_TOL * L_null){
            break;
        }
    }
    }
}

arma::vec s2net::predict_response(const arma::mat & newX){
    // Asume that newX in the same space as xL, xU
    arma::vec eta = newX*beta + intercept;
    return eta;
}

arma::vec s2net::predict_probability(const arma::mat & newX){
    arma::vec eta = predict_response(newX);
    arma::vec prob = 1/(1 + arma::exp(-eta));
    return prob;
}

arma::vec s2net::predict_class(const arma::mat & newX){
    arma::vec prob = predict_probability(newX);
    prob.for_each([](double & p_i){
        if(p_i > 0.5){ // or not?
            p_i = 1;
        }else{
            p_i = 0;
        }});
    return prob;
}

arma::vec s2net::predict(const arma::mat & newX, int type){
    switch (type)
    {
    case TYPE_PREDICT_RESPONSE:
        return predict_response(newX);
        break;
    case TYPE_PREDICT_PROB:
        return predict_probability(newX);
        break;
    case TYPE_PREDICT_CLASS:
        return predict_class(newX);
        break;
    default:    
        switch (loss)
        {
        case LOSS_LOGIT:
            return predict_probability(newX);
            break;
        default:
            return predict_response(newX);
            break;
        }
        break;
    }
}

// Expose class s2net
RCPP_MODULE(Rcpp_s2net_export){
    Rcpp::class_<s2net>("s2net")
    
    .constructor<const Rcpp::List , int>()

    .method("fit", &s2net::fit, "Computes beta using FISTA")
    
    .method("setupFista", &s2net::setupFista, "Sets the hyperparameters for the FISTA algorithm")

    .method("predict", &s2net::predict, "Predicts response vector")

    .property("beta", &s2net::get_beta, &s2net::set_beta)
    .property("intercept", &s2net::get_intercept, &s2net::set_intercept)

    ;
}
