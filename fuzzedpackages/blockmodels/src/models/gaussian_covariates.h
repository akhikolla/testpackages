
/*
 * This is package skeleton for gaussian_covariates model
 *
 * You should delete this header when the gaussian_covariates model is ready for use to avoid
 * any confusion
 */

class gaussian_covariates
{
    public:

    class network
    {
        public:
        /* Here you should put all the variable which define the network
         * for scalar model, a mat field is good choice.
         * See bernoulli, poisson model for example.
         *
         * Covariates must also stocked here.
         * See poisson_covariates model for example
         */

        mat adj;
        cube covariates;


        /* Here you should add all precomputed values which depends only on the
         * network usefull in various functions, to avoid computing these value
         * many time.
         */
        mat adjZD;
        mat MonesZD;
        mat Mones;
        mat Monest;


        network(Rcpp::List & network_from_R)
        {
            /* Here you must define how initialize the fields describing the
             * network. Precomputed value must be computed here.
             * See bernoulli, poisson or poisson_covariates for example.
             *
             * For scalar network, the provided list have a adjacency field
             * which contains the adjacency matrix.
             * For scalar network with covariates vectors on edges (existent or
             * not), the provided list have a covariates field which contain a
             * list of matrices. This matrices have the same size than the
             * adjacency matrix, and the i-th matrix is the matrix of the i-th
             * covariate on all edges.
             */
            adj = Rcpp::as<mat>(network_from_R["adjacency"]);

            Rcpp::List covariates_list = network_from_R["covariates"];

            covariates.set_size(adj.n_rows,adj.n_cols,covariates_list.size());
            for(int k=0; k<covariates_list.size(); k++)
                covariates.slice(k) = Rcpp::as<mat>(covariates_list[k]);

            adjZD = fill_diag(adj,0);
            Mones = ones<mat>(adj.n_rows,adj.n_cols);
            MonesZD = fill_diag(Mones,0);
            Monest = Mones.t();
        }
    };

    // parameters
    unsigned int n_parameters;
    bool symmetric;

    /* Here you must put the gaussian_covariates model parameters */
    mat mu;
    colvec beta;
    double sigma2;

    gaussian_covariates(SBM & membership, gaussian_covariates::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM case
         */
        mu = (membership.Z.t() * net.adjZD * membership.Z)
              /
             (membership.Z.t() * net.MonesZD * membership.Z);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        mat diff = fill_diag(membership.Z * mu * membership.Z.t() - net.adj,0);
        sigma2 = accu(diff % diff)/((membership.Z.n_rows-1)*membership.Z.n_rows);

        n_parameters = mu.n_elem + beta.n_elem + 1;
        symmetric = false;

    }
    
    gaussian_covariates(SBM_sym & membership, gaussian_covariates::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM_sym case
         */
        mu = (membership.Z.t() * net.adjZD * membership.Z)
              /
             (membership.Z.t() * net.MonesZD * membership.Z);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);
        
        mat diff = fill_diag(membership.Z * mu * membership.Z.t() - net.adj,0);
        sigma2 = accu(diff % diff)/((membership.Z.n_rows-1)*membership.Z.n_rows);

        n_parameters = mu.n_rows*(mu.n_rows+1)/2 + beta.n_elem + 1;
        symmetric = true;

    }
    
    gaussian_covariates(LBM & membership, gaussian_covariates::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * LBM case 
         */
        mu = (membership.Z1.t() * net.adj * membership.Z2)
              /
             (membership.Z1.t() * net.MonesZD * membership.Z2);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);
        
        mat diff = membership.Z1 * mu * membership.Z2.t() - net.adj;
        sigma2 = accu(diff % diff)/(membership.Z1.n_rows*membership.Z2.n_rows);

        n_parameters = mu.n_elem + beta.n_elem + 1;
        symmetric = false;

    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;
        values["mu"] = mu;
        values["beta"] = beta;
        values["sigma2"] = sigma2;

        /* Here you must define the export way of your model parameters in the
         * list values. This is this list which is returned to R. */

        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    gaussian_covariates(SBM &, const vec &);
    gaussian_covariates(SBM_sym &, const vec &);
    gaussian_covariates(LBM &, const vec &);
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


inline
mat gaussian_covariates_compute_B(colvec & beta, cube & covariates)
{
    // B definied as \sum \beta_k C_{:,:,k}
    mat B=zeros<mat>(covariates.n_rows,covariates.n_cols);

    for(unsigned int k=0;k<covariates.n_slices;k++)
    {
        B+=beta(k)*covariates.slice(k);
    }
    return B;
}


/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ_{iq} += \sum_{jl} Z_{jl} (logf(i,j,q,l) + logf(j,i,l,q)
 *
 * This is usefull for SBM. See bernoulli and poisson for example.
 */

template<>
inline
void e_fixed_step(SBM & membership,
                  gaussian_covariates & model,
                  gaussian_covariates::network & net,
                  mat & lZ)
{
    mat adjmBZD = fill_diag(net.adj-gaussian_covariates_compute_B(model.beta,net.covariates),0);
    lZ += 1.0/(2*model.sigma2) * (
            - net.MonesZD * membership.Z * (model.mu.t() % model.mu.t())
            + 2 * adjmBZD * membership.Z * model.mu.t()
            - net.MonesZD * membership.Z * (model.mu % model.mu)
            + 2 * adjmBZD.t() * membership.Z * model.mu
        );
}





/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ_{iq} += \sum_{jl} Z_{jl} logf(i,j,q,l)
 *
 * This is usefull for SBM_sym. See bernoulli and poisson for example.
 */

template<>
inline
void e_fixed_step(SBM_sym & membership,
                  gaussian_covariates & model,
                  gaussian_covariates::network & net,
                  mat & lZ)
{
    mat adjmBZD = fill_diag(net.adj-gaussian_covariates_compute_B(model.beta,net.covariates),0);
    lZ  += 1.0/(2*model.sigma2) * (
            - adjmBZD * membership.Z * (model.mu % model.mu)
            + 2 * net.adjZD * membership.Z * model.mu
        );

}





/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ1_{iq} += \sum_{jl} Z2_{jl} logf(i,j,q,l)
 * forall j,l ; lZ2_{jl} += \sum_{iq} Z1_{iq} logf(i,j,q,l)
 *
 * This is usefull for LBM. See bernoulli and poisson for example.
 */

template<>
inline
void e_fixed_step(LBM & membership,
                  gaussian_covariates & model,
                  gaussian_covariates::network & net,
                  mat & lZ1,
                  mat & lZ2)
{
    mat adjmB = net.adj-gaussian_covariates_compute_B(model.beta,net.covariates);
    lZ1 += 1.0/(2*model.sigma2) * (
            - net.Mones * membership.Z2 * (model.mu.t() % model.mu.t())
            + 2 * adjmB * membership.Z2 * model.mu.t()
        );
    lZ2 += 1.0/(2*model.sigma2) * (
            - net.Monest * membership.Z1 * (model.mu % model.mu)
            + 2 * adjmB.t() * membership.Z1 * model.mu
        );
}





/* Usefull, Optional
 *
 * If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of SBM.
 *
 * This function must return PL in the maximum.
 * 
 * This is usefull for SBM. See bernoulli and poisson for example.
 */

// template<>
// inline
// double m_step(SBM & membership,
//               gaussian_covariates & model,
//               gaussian_covariates::network & net)
// {
// }





/* Usefull, Optional
 *
 * If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of SBM_sym.
 *
 * This function must return PL in the maximum.
 * 
 * This is usefull for SBM_sym. See bernoulli and poisson for example.
 */

// template<>
// inline
// double m_step(SBM_sym & membership,
//               gaussian_covariates & model,
//               gaussian_covariates::network & net)
// {
// }





/* Usefull, Optional
 *
 * If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of LBM.
 *
 * This function must return PL in the maximum.
 * 
 * This is usefull for LBM. See bernoulli and poisson for example.
 */

// template<>
// inline
// double m_step(LBM & membership,
//               gaussian_covariates & model,
//               gaussian_covariates::network & net)
// {
// }





/* If you have defined m_step(SBM &, gaussian_covariates &, gaussian_covariates::network &)
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
double PL(gaussian_covariates & model,
          SBM & membership,
          gaussian_covariates::network & net)
{
    mat adjmBZD = fill_diag(net.adj-gaussian_covariates_compute_B(model.beta,net.covariates),0);
    return(
            -.5*(membership.Z.n_rows * (membership.Z.n_rows-1))*log(2*PI*model.sigma2)
            -1.0/(2*model.sigma2)*(
                accu(
                    adjmBZD % adjmBZD
                )
                +
                accu(
                    (
                        (model.mu % model.mu)
                        %
                        (membership.Z.t() * net.MonesZD * membership.Z)
                    )
                    -
                    (
                        2 * model.mu
                        %
                        (membership.Z.t() * adjmBZD * membership.Z)
                    )
                )
            )
        );

}





/* If you have defined m_step(SBM_sym &, gaussian_covariates &, gaussian_covariates::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have direct formula to compute PL, without loops over (i,j,q,l), in
 * the case of SBM_sym, use it.
 *
 * See poisson_covariate for example.
 */

template<>
inline
double PL(gaussian_covariates & model,
          SBM_sym & membership,
          gaussian_covariates::network & net)
{
    mat adjmBZD = fill_diag(net.adj-gaussian_covariates_compute_B(model.beta,net.covariates),0);
    return( .5*(
            -.5*(membership.Z.n_rows * (membership.Z.n_rows-1))*log(2*PI*model.sigma2)
            -1.0/(2*model.sigma2)*(
                accu(
                    adjmBZD % adjmBZD
                )
                +
                accu(
                    (
                        (model.mu % model.mu)
                        %
                        (membership.Z.t() * net.MonesZD * membership.Z)
                    )
                    -
                    (
                        2 * model.mu
                        %
                        (membership.Z.t() * adjmBZD * membership.Z)
                    )
                )
            )
        ));
}





/* If you have defined m_step(LBM &, gaussian_covariates &, gaussian_covariates::network &)
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
double PL(gaussian_covariates & model,
          LBM & membership,
          gaussian_covariates::network & net)
{
    mat adjmB = net.adj-gaussian_covariates_compute_B(model.beta,net.covariates);
    return(
            -.5*(membership.Z1.n_rows * membership.Z2.n_rows)*log(2*PI*model.sigma2)
            -1.0/(2*model.sigma2)*(
                accu(
                    adjmB % adjmB
                )
                +
                accu(
                    (
                        (model.mu % model.mu)
                        %
                        (membership.Z1.t() * net.Mones * membership.Z2)
                    )
                    -
                    (
                        2 * model.mu
                        %
                        (membership.Z1.t() * adjmB * membership.Z2)
                    )
                )
            )
        );
}





/* If one of the following specialization
 *     - m_step(LBM &, gaussian_covariates &, gaussian_covariates::network)
 *     - m_step(SBM &, gaussian_covariates &, gaussian_covariates::network)
 *     - m_step(SBM_sym &, gaussian_covariates &, gaussian_covariates::network)
 *   is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * Export the model parameters as a vector.
 *
 * See naive_bernoulli and poisson_covariates for example.
 */

inline
vec gaussian_covariates::to_vector()
{
    vec out(n_parameters);
    vec vmu = (symmetric) ? vech(mu) : reshape(mu,mu.n_elem,1);

    out.subvec(0,vmu.n_elem-1) = vmu;
    out.subvec(vmu.n_elem,n_parameters-2) = beta;
    out(n_parameters-1) = sigma2;
    return(out);
}





/* If m_step(SBM &, gaussian_covariates &, gaussian_covariates::network) is not defined
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

gaussian_covariates::gaussian_covariates(SBM & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    mu = reshape(vectorized.subvec(0,Q*Q-1),Q,Q);
    beta = vectorized.subvec(Q*Q,vectorized.n_elem-2);
    sigma2 = vectorized(vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = false;
}





/* If m_step(SBM_sym &, gaussian_covariates &, gaussian_covariates::network) is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 * 
 * Constructor of the model, from a vector of parameter named vectorized with
 * the given membership of type SBM_sym.
 *
 * See naive_bernoulli for example.
 */

gaussian_covariates::gaussian_covariates(SBM_sym & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    mu = unvech(vectorized.subvec(0,Q*(Q+1)/2-1));
    beta = vectorized.subvec(Q*(Q+1)/2,vectorized.n_elem-2);
    sigma2 = vectorized(vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = true;
}





/* If m_step(LBM &, gaussian_covariates &, gaussian_covariates::network) is not defined
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

gaussian_covariates::gaussian_covariates(LBM & membership, const vec & vectorized)
{
    unsigned int Q1=membership.Z1.n_cols;
    unsigned int Q2=membership.Z2.n_cols;
    mu = reshape(vectorized.subvec(0,Q1*Q2-1),Q1,Q2);
    beta = vectorized.subvec(Q1*Q2,vectorized.n_elem-2);
    sigma2 = vectorized(vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric=false;
}





/* If one of the following specialization
 *     - m_step(LBM &, gaussian_covariates &, gaussian_covariates::network)
 *     - m_step(SBM &, gaussian_covariates &, gaussian_covariates::network)
 *     - m_step(SBM_sym &, gaussian_covariates &, gaussian_covariates::network)
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
double maximum_step_in_direction(membership_type & membership, gaussian_covariates & model, gaussian_covariates::network & net, vec & direction)
{
    double ds2 = direction(direction.n_elem-1);
    if(ds2>=0)
        return 1;

    double a = -model.sigma2/ds2;

    if(a<1)
        return a;
    
    return 1;
}





/* If you have defined m_step(SBM &, gaussian_covariates &, gaussian_covariates::network &)
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
vec grad(gaussian_covariates & model,
         SBM & membership,
         gaussian_covariates::network & net)
{
    mat adjmBZD = fill_diag(net.adj-gaussian_covariates_compute_B(model.beta,net.covariates),0);
    mat dmu = 1.0/model.sigma2 *
        (
            membership.Z.t() * adjmBZD * membership.Z
            - ( (membership.Z.t() * net.MonesZD * membership.Z) % model.mu)
        );

    mat M = 1.0/model.sigma2 * fill_diag(
                -membership.Z * model.mu * membership.Z.t()
                + ((
                        membership.Z 
                        * ones<mat>(membership.Z.n_cols,membership.Z.n_cols) 
                        * membership.Z.t()
                   ) % adjmBZD)
            ,0);

    vec dbeta(model.beta.n_elem);
    for(unsigned int k=0;k<dbeta.n_elem;k++)
    {
        dbeta(k) = accu( M % net.covariates.slice(k));
    }

    double dsigma2 =
            -.5*(membership.Z.n_rows * (membership.Z.n_rows-1))/model.sigma2
            +1.0/(2*model.sigma2*model.sigma2)*(
                accu(
                    adjmBZD % adjmBZD
                )
                +
                accu(
                    (
                        (model.mu % model.mu)
                        %
                        (membership.Z.t() * net.MonesZD * membership.Z)
                    )
                    -
                    (
                        2 * model.mu
                        %
                        (membership.Z.t() * adjmBZD * membership.Z)
                    )
                )
            );
    
    vec out(model.n_parameters);
    vec vdmu = reshape(dmu,dmu.n_elem,1);

    out.subvec(0,vdmu.n_elem-1) = vdmu;
    out.subvec(vdmu.n_elem,model.n_parameters-2) = dbeta;
    out(model.n_parameters-1) = dsigma2;
    return(out);
}





/* If you have defined m_step(SBM_sym &, gaussian_covariates &, gaussian_covariates::network &)
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
vec grad(gaussian_covariates & model,
         SBM_sym & membership,
         gaussian_covariates::network & net)
{
    mat adjmBZD = fill_diag(net.adj-gaussian_covariates_compute_B(model.beta,net.covariates),0);
    mat dmu = 1.0/model.sigma2 *
        (
            membership.Z.t() * adjmBZD * membership.Z
            - ( (membership.Z.t() * net.MonesZD * membership.Z) % model.mu)
        );

    mat M = 1.0/model.sigma2 * fill_diag(
                -membership.Z * model.mu * membership.Z.t()
                + ((
                        membership.Z 
                        * ones<mat>(membership.Z.n_cols,membership.Z.n_cols) 
                        * membership.Z.t()
                   ) % adjmBZD)
            ,0);

    vec dbeta(model.beta.n_elem);
    for(unsigned int k=0;k<dbeta.n_elem;k++)
    {
        dbeta(k) = accu( M % net.covariates.slice(k));
    }

    double dsigma2 =
            -.5*(membership.Z.n_rows * (membership.Z.n_rows-1))/model.sigma2
            +1.0/(2*model.sigma2*model.sigma2)*(
                accu(
                    adjmBZD % adjmBZD
                )
                +
                accu(
                    (
                        (model.mu % model.mu)
                        %
                        (membership.Z.t() * net.MonesZD * membership.Z)
                    )
                    -
                    (
                        2 * model.mu
                        %
                        (membership.Z.t() * adjmBZD * membership.Z)
                    )
                )
            );
    
    vec out(model.n_parameters);
    vec vdmu = vech(dmu);

    out.subvec(0,vdmu.n_elem-1) = vdmu;
    out.subvec(vdmu.n_elem,model.n_parameters-2) = dbeta;
    out(model.n_parameters-1) = dsigma2;
    return(out);
}





/* If you have defined m_step(LBM &, gaussian_covariates &, gaussian_covariates::network &)
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
vec grad(gaussian_covariates & model,
         LBM & membership,
         gaussian_covariates::network & net)
{
    mat adjmB = net.adj-gaussian_covariates_compute_B(model.beta,net.covariates);
    mat dmu = 1.0/model.sigma2 *
        (
            membership.Z1.t() * adjmB * membership.Z2
            - ( (membership.Z1.t() * net.Mones * membership.Z2) % model.mu)
        );

    mat M = 1/model.sigma2 * (-membership.Z1 * model.mu * membership.Z2.t()
                + ((
                        membership.Z1
                        * ones<mat>(membership.Z1.n_cols,membership.Z2.n_cols) 
                        * membership.Z2.t()
                   ) % adjmB));

    vec dbeta(model.beta.n_elem);
    for(unsigned int k=0;k<dbeta.n_elem;k++)
    {
        dbeta(k) = accu( M % net.covariates.slice(k));
    }

    double dsigma2 =
            -.5*(membership.Z1.n_rows * membership.Z2.n_rows)/model.sigma2
            +1.0/(2*model.sigma2*model.sigma2)*(
                accu(
                    adjmB % adjmB
                )
                +
                accu(
                    (
                        (model.mu % model.mu)
                        %
                        (membership.Z1.t() * net.Mones * membership.Z2)
                    )
                    -
                    (
                        2 * model.mu
                        %
                        (membership.Z1.t() * adjmB * membership.Z2)
                    )
                )
            );
    
    vec out(model.n_parameters);
    vec vdmu = reshape(dmu,dmu.n_elem,1);

    out.subvec(0,vdmu.n_elem-1) = vdmu;
    out.subvec(vdmu.n_elem,model.n_parameters-2) = dbeta;
    out(model.n_parameters-1) = dsigma2;
    return(out);
}





/* If one of the following specialization
 *     - grad(gaussian_covariates &, LBM &, gaussian_covariates::network)
 *     - grad(gaussian_covariates &, SBM &, gaussian_covariates::network)
 *     - grad(gaussian_covariates &, SBM_sym &, gaussian_covariates::network)
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
// vec grad_logf(gaussian_covariates & model,
//               gaussian_covariates::network & net,
//               unsigned int i,
//               unsigned int j,
//               unsigned int q,
//               unsigned int l)
// {
// }





/* If grad_logf(gaussian_covariates &,
 *         gaussian_covariates::network &,
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
// double grad_logf(gaussian_covariates & model,
//                  gaussian_covariates::network & net,
//                  unsigned int i,
//                  unsigned int j,
//                  unsigned int q,
//                  unsigned int l,
//                  unsigned int k)
// {
// }





/* If one of the following specialization
 *     - PL(gaussian_covariates &, SBM &, gaussian_covariates::network)
 *     - PL(gaussian_covariates &, SBM_sym &, gaussian_covariates::network)
 *     - PL(gaussian_covariates &, LBM &, gaussian_covariates::network)
 *   is Usefull and not definied
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * logf(i,j,q,l) see notations.
 */

// inline
// double logf(gaussian_covariates & model,
//             gaussian_covariates::network & net,
//             unsigned int i,
//             unsigned int j,
//             unsigned int q,
//             unsigned int l)
// {
// }
