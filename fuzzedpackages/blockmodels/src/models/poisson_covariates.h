
/*
 * This is package skeleton for poisson_covariates model
 *
 * You should delete this header when the poisson_covariates model is ready for use to avoid
 * any confusion
 */

class poisson_covariates
{
    public:

    class network
    {
        public:
         
        mat adj; // adjacency matrix
        cube covariates; // cube of covariates, two first dimention have the
                         // same size than the ajacency matrix, the 3rd is the
                         // index of covariate.

        // precomputed matrix
        mat Mones;
        mat adjZD;
        mat adjZDt;
        mat MonesZD;
        mat adjt;

        // precomuted values for the liklihood
        double accu_log_fact_XZD;
        double accu_log_fact_X;

        network(Rcpp::List & network_from_R)
        {
            adj = Rcpp::as<mat>(network_from_R["adjacency"]);

            Rcpp::List covariates_list = network_from_R["covariates"];

            covariates.set_size(adj.n_rows,adj.n_cols,covariates_list.size());
            for(int k=0; k<covariates_list.size(); k++)
                covariates.slice(k) = Rcpp::as<mat>(covariates_list[k]);

            // precomputation
            Mones = ones<mat>(adj.n_rows,adj.n_cols);
            adjZD = fill_diag(adj,0);
            adjZDt = adjZD.t();
            MonesZD = fill_diag(Mones,0);
            accu_log_fact(adj,accu_log_fact_X,accu_log_fact_XZD);
            adjt = adj.t();
        }
    };

    // parameters
    unsigned int n_parameters;
    bool symmetric;

    /* Here you must put the poisson_covariates model parameters */
    mat lambda;
    colvec beta;

    poisson_covariates(SBM & membership, poisson_covariates::network & net)
    {
        lambda = (membership.Z.t() * net.adjZD * membership.Z)
                  /
                 (membership.Z.t() * net.MonesZD * membership.Z);
        
        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = lambda.n_elem + beta.n_elem;
        symmetric = false;
    }
    
    poisson_covariates(SBM_sym & membership, poisson_covariates::network & net)
    {
        lambda = (membership.Z.t() * net.adjZD * membership.Z)
                  /
                 (membership.Z.t() * net.MonesZD * membership.Z);
        
        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = lambda.n_rows*(lambda.n_rows+1)/2 + beta.n_elem;

        symmetric = true;
    }
    
    poisson_covariates(LBM & membership, poisson_covariates::network & net)
    {
        lambda = (membership.Z1.t() * net.adj * membership.Z2)
                  /
                 (membership.Z1.t() * net.Mones * membership.Z2);
        
        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = lambda.n_rows * lambda.n_cols + net.covariates.n_slices;

        symmetric = false;
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;
        values["lambda"] = lambda;
        values["beta"] = beta;

        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    poisson_covariates(SBM &, const vec &);
    poisson_covariates(SBM_sym &, const vec &);
    poisson_covariates(LBM &, const vec &);
};





/* Function for the model.
 */

/******************************************************************************
 *
 * Notations :
 *
 * i a row node
 * j a col node
 * q a row class
 * l a col class
 *
 * (for SBM there are no distinctions between row and col nodes and classes)
 *
 * f(i,j,q,l) the prob or density of the variable
 * X_{ij} conditionnally to i belong to class q and j to class l
 *
 * logf(i,j,q,l) is a quicker way to write log(f(i,j,q,l))
 *
 ******************************************************************************/


/* models usefull function
 */



template<>
inline
void e_fixed_step(SBM & membership,
                  poisson_covariates & model,
                  poisson_covariates::network & net,
                  mat & lZ)
{
    mat eBZD = fill_diag( exp(compute_B(model.beta,net.covariates)), 0);

    lZ += -eBZD * membership.Z * model.lambda.t()
          +net.adjZD * membership.Z * log(model.lambda).t()
          -eBZD.t() * membership.Z * model.lambda
          +net.adjZDt * membership.Z * log(model.lambda);
}



template<>
inline
void e_fixed_step(SBM_sym & membership,
                  poisson_covariates & model,
                  poisson_covariates::network & net,
                  mat & lZ)
{
    mat eBZD = fill_diag( exp(compute_B(model.beta,net.covariates)), 0);

    lZ += -eBZD * membership.Z * model.lambda
          +net.adjZD * membership.Z * log(model.lambda);
}



template<>
inline
void e_fixed_step(LBM & membership,
                  poisson_covariates & model,
                  poisson_covariates::network & net,
                  mat & lZ1,
                  mat & lZ2)
{
    mat eB = exp(compute_B(model.beta,net.covariates));

    lZ1 += -eB * membership.Z2 * model.lambda.t()
           +net.adj * membership.Z2 * log(model.lambda).t();

    lZ2 += -eB.t() * membership.Z1 * model.lambda
           + net.adjt * membership.Z1 * log(model.lambda);
}


/* If you have defined m_step(SBM &, poisson_covariates &, poisson_covariates::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have direct formula to compute PL, without loops over (i,j,q,l), in
 * the case of SBM, use it.
 *
 * See poisson_covariate for example.
 */

template<>
inline
double PL(poisson_covariates & model,
          SBM & membership,
          poisson_covariates::network & net)
{
    mat B = compute_B(model.beta,net.covariates);
    mat BZD = fill_diag(B,0);
    mat eBZD = fill_diag(exp(B),0);

    return(
        accu(
                -
                model.lambda
                %
                ( membership.Z.t() * eBZD * membership.Z )
            +
                log(model.lambda)
                %
                ( membership.Z.t() * net.adjZD * membership.Z )
            )
        +
        accu( net.adj%B )
        - net.accu_log_fact_X
        );
    
}




template<>
inline
double PL(poisson_covariates & model,
          SBM_sym & membership,
          poisson_covariates::network & net)
{
    mat B = compute_B(model.beta,net.covariates);
    mat BZD = fill_diag(B,0);
    mat eBZD = fill_diag(exp(B),0);

    return((
        accu(
                -
                model.lambda
                %
                ( membership.Z.t() * eBZD * membership.Z )
            +
                log(model.lambda)
                %
                ( membership.Z.t() * net.adjZD * membership.Z )
            )
        +
        accu( net.adj%B )
        - net.accu_log_fact_X
        ) /2 );
    
}




/* If you have defined m_step(LBM &, poisson_covariates &, poisson_covariates::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have direct formula to compute PL, without loops over (i,j,q,l), in
 * the case of LBM, use it.
 *
 * See poisson_covariate for example.
 */

template<>
inline
double PL(poisson_covariates & model,
          LBM & membership,
          poisson_covariates::network & net)
{
    mat B = compute_B(model.beta,net.covariates);

    return(
        accu(
                -
                model.lambda
                %
                ( membership.Z1.t() * exp(B) * membership.Z2 )
            +
                log(model.lambda)
                %
                ( membership.Z1.t() * net.adj * membership.Z2 )
            )
        +
        accu( net.adj%B )
        - net.accu_log_fact_X
        );
}





/* If one of the following specialization
 *     - m_step(LBM &, poisson_covariates &, poisson_covariates::network)
 *     - m_step(SBM &, poisson_covariates &, poisson_covariates::network)
 *   is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * Export the model parameters as a vector.
 *
 * See naive_bernoulli for example.
 */

inline
vec poisson_covariates::to_vector()
{
    vec out(n_parameters);
    vec vlambda = (symmetric) ? vech(lambda) : reshape(lambda,lambda.n_elem,1);

    out.subvec(0,vlambda.n_elem-1) = vlambda;
    out.subvec(vlambda.n_elem,n_parameters-1) = beta;
    return(out);
}





/* If m_step(SBM &, poisson_covariates &, poisson_covariates::network) is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 * 
 * Constructor of the model, from a vector of parameter named vectorized with
 * the given membership of type SBM.
 *
 * See naive_bernoulli for example.
 */

poisson_covariates::poisson_covariates(SBM & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    lambda = reshape(vectorized.subvec(0,Q*Q-1),Q,Q);
    beta = vectorized.subvec(Q*Q,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = false;
}





poisson_covariates::poisson_covariates(SBM_sym & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    lambda = unvech(vectorized.subvec(0,Q*(Q+1)/2-1));
    beta = vectorized.subvec(Q*(Q+1)/2,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = true;
}

/* If m_step(LBM &, poisson_covariates &, poisson_covariates::network) is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 * 
 * Constructor of the model, from a vector of parameter named vectorized with
 * the given membership of type LBM.
 *
 * See naive_bernoulli for example.
 */

poisson_covariates::poisson_covariates(LBM & membership, const vec & vectorized)
{
    unsigned int Q1=membership.Z1.n_cols;
    unsigned int Q2=membership.Z2.n_cols;
    lambda = reshape(vectorized.subvec(0,Q1*Q2-1),Q1,Q2);
    beta = vectorized.subvec(Q1*Q2,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric=false;
}





/* If one of the following specialization
 *     - m_step(LBM &, poisson_covariates &, poisson_covariates::network)
 *     - m_step(SBM &, poisson_covariates &, poisson_covariates::network)
 *   is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * From the given model, in the given direction for model parameters, return the
 * maximum betwenn 1 and the maximum step the descend algorithm can do.
 *
 * See naive_bernoulli for example.
 *
 * If you don't have any constraint on you parameter vector, return 1, else, use
 * your brain.
 */

template<class membership_type>
inline
double maximum_step_in_direction(membership_type & membership, poisson_covariates & model, poisson_covariates::network & net, vec & direction)
{
    vec vlambda = (model.symmetric) ? 
                        vech(model.lambda)
                        :
                        reshape(model.lambda,model.lambda.n_elem,1);
    double amax=1;
    for(unsigned int i=0;i<vlambda.n_elem;i++)
    {
        if(direction(i)!=0 && vlambda(i)*direction(i)<0)
        {
            double a=-vlambda(i)/direction(i);
            if(a<amax)
                amax=a;
        }
    }
    // no constraints on part of vector concering beta

    return(amax);
}





/* If you have defined m_step(SBM &, poisson_covariates &, poisson_covariates::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have a direct formula to compute the gradient of PL in respect to your
 * parameters vectors without loops over (i,j,q,l), in the case of SBM, use it.
 *
 * See poisson_covariates for example.
 */

template<>
inline
vec grad(poisson_covariates & model,
         SBM & membership,
         poisson_covariates::network & net)
{
    mat B = compute_B(model.beta,net.covariates);
    mat eBZD=fill_diag(exp(B),0);

    mat dPLodlambda = -membership.Z.t() * eBZD * membership.Z
                      + (membership.Z.t() * net.adjZD * membership.Z) / model.lambda;

    mat coef = net.adjZD - ((membership.Z * model.lambda * membership.Z.t()) % eBZD);

    vec dPLodbeta(model.beta.n_elem);

    for(unsigned int k=0;k<dPLodbeta.n_elem;k++)
        dPLodbeta(k) = accu(coef % net.covariates.slice(k));

    vec out(model.n_parameters);
    out.subvec(0,dPLodlambda.n_elem-1) = reshape(dPLodlambda,dPLodlambda.n_elem,1);
    out.subvec(dPLodlambda.n_elem,model.n_parameters-1) = dPLodbeta;
    return(out);
}



template<>
inline
vec grad(poisson_covariates & model,
         SBM_sym & membership,
         poisson_covariates::network & net)
{
    mat B = compute_B(model.beta,net.covariates);
    mat eBZD=fill_diag(exp(B),0);

    mat dPLodlambda = -membership.Z.t() * eBZD * membership.Z
                      + (membership.Z.t() * net.adjZD * membership.Z) / model.lambda;

    mat coef = net.adjZD - ((membership.Z * model.lambda * membership.Z.t()) % eBZD);

    vec dPLodbeta(model.beta.n_elem);

    for(unsigned int k=0;k<dPLodbeta.n_elem;k++)
        dPLodbeta(k) = accu(coef % net.covariates.slice(k));

    vec out(model.n_parameters);
    unsigned int Q = model.lambda.n_rows;
    out.subvec(0,Q*(Q+1)/2-1) = vech(dPLodlambda)/2;
    out.subvec(Q*(Q+1)/2,model.n_parameters-1) = dPLodbeta/2;
    return(out);
}




/* If you have defined m_step(LBM &, poisson_covariates &, poisson_covariates::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have a direct formula to compute the gradient of PL in respect to your
 * parameters vectors without loops over (i,j,q,l), in the case of LBM, use it.
 *
 * See poisson_covariates for example.
 */

template<>
inline
vec grad(poisson_covariates & model,
         LBM & membership,
         poisson_covariates::network & net)
{
    mat B = compute_B(model.beta,net.covariates);
    mat eB=exp(B);

    mat dPLodlambda = -membership.Z1.t() * eB * membership.Z2
                      + (membership.Z1.t() * net.adj * membership.Z2) / model.lambda;

    mat coef = net.adj - ((membership.Z1 * model.lambda * membership.Z2.t()) % eB);

    vec dPLodbeta(model.beta.n_elem);

    for(unsigned int k=0;k<dPLodbeta.n_elem;k++)
        dPLodbeta(k) = accu(coef % net.covariates.slice(k));

    vec out(model.n_parameters);
    out.subvec(0,dPLodlambda.n_elem-1) = reshape(dPLodlambda,dPLodlambda.n_elem,1);
    out.subvec(dPLodlambda.n_elem,model.n_parameters-1) = dPLodbeta;
    return(out);
}





/* If one of the following specialization
 *     - grad(poisson_covariates &, LBM &, poisson_covariates::network)
 *     - grad(poisson_covariates &, SBM &, poisson_covariates::network)
 *   is Usefull and not defined
 * Then
 *     Usefull, Optional
 * Else
 *     Useless
 *
 * The gradient of logf(i,j,q,l) in respect to the parameter vector.
 */

// template<>
// inline
// vec grad_logf(poisson_covariates & model,
//               poisson_covariates::network & net,
//               unsigned int i,
//               unsigned int j,
//               unsigned int q,
//               unsigned int l)
// {
// }





/* If grad(poisson_covariates &,
 *         poisson_covariates::network &,
 *         unsigned int,
 *         unsigned int,
 *         unsigned int,
 *         unsigned int)
 *    is Usefull and not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 *
 * The derivative of logf(i,j,q,l) in respect to the k-th parameters of the
 * parameter vector.
 *
 * See naive_bernoulli for example.
 */

// inline
// double grad_logf(poisson_covariates & model,
//                  poisson_covariates::network & net,
//                  unsigned int i,
//                  unsigned int j,
//                  unsigned int q,
//                  unsigned int l,
//                  unsigned int k)
// {
// }





/* If one of the following specialization
 *     - PL(poisson_covariates &, SBM &, poisson_covariates::network)
 *     - PL(poisson_covariates &, LBM &, poisson_covariates::network)
 *   is Usefull and not definied
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * logf(i,j,q,l) see notations.
 */

// inline
// double logf(poisson_covariates & model,
//             poisson_covariates::network & net,
//             unsigned int i,
//             unsigned int j,
//             unsigned int q,
//             unsigned int l)
// {
// }
