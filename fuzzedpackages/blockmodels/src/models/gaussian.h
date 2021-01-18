
/*
 * This is package skeleton for gaussian model
 *
 * You should delete this header when the gaussian model is ready for use to avoid
 * any confusion
 */

class gaussian
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

        //precomputed matrix for SBM
        mat adjZD;
        mat adjZDt;
        mat MonesZD;

        
        //precomputed matrix for LBM
        mat Mones;
        mat Monest;
        mat adjt;
        
        //precomputed value for likilhood
        double accu_adj_square;
        double accu_adjZD_square;
        


        /* Here you should add all precomputed values which depends only on the
         * network usefull in various functions, to avoid computing these value
         * many time.
         */

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

            adjZD = fill_diag(adj,0);
            adjt = adj.t();
            Mones = ones<mat>(adj.n_rows,adj.n_cols);
            Monest = Mones.t();
            adjZDt = adjZD.t();
            MonesZD = fill_diag(Mones,0);

            accu_adj_square = accu( adj % adj );
            accu_adjZD_square = accu( adjZD % adjZD );
        }
    };

    // parameters
    unsigned int n_parameters;
    mat mu;
    double sigma2;

    /* Here you must put the gaussian model parameters */

    gaussian(SBM & membership, gaussian::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM case
         */
        n_parameters = membership.Z.n_cols * membership.Z.n_cols+1;
        mu.set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    gaussian(SBM_sym & membership, gaussian::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM_sym case
         */
        n_parameters = membership.Z.n_cols * (membership.Z.n_cols+1)/2+1;
        mu.set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    gaussian(LBM & membership, gaussian::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * LBM case 
         */
        n_parameters = membership.Z1.n_cols * membership.Z2.n_cols+1;
        mu.set_size(membership.Z1.n_cols,membership.Z2.n_cols);
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;
        values["mu"] = mu;
        values["sigma2"] = sigma2;

        /* Here you must define the export way of your model parameters in the
         * list values. This is this list which is returned to R. */

        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    gaussian(SBM &, const vec &);
    gaussian(SBM_sym &, const vec &);
    gaussian(LBM &, const vec &);
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
                  gaussian & model,
                  gaussian::network & net,
                  mat & lZ)
{
    lZ += 1.0/(2*model.sigma2) * (
            - net.MonesZD * membership.Z * (model.mu.t() % model.mu.t())
            + 2 * net.adjZD * membership.Z * model.mu.t()
            - net.MonesZD * membership.Z * (model.mu % model.mu)
            + 2 * net.adjZDt * membership.Z * model.mu
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
                  gaussian & model,
                  gaussian::network & net,
                  mat & lZ)
{
    lZ  += 1.0/(2*model.sigma2) * (
            - net.MonesZD * membership.Z * (model.mu % model.mu)
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
                  gaussian & model,
                  gaussian::network & net,
                  mat & lZ1,
                  mat & lZ2)
{
    lZ1 += 1.0/(2*model.sigma2) * (
            - net.Mones * membership.Z2 * (model.mu.t() % model.mu.t())
            + 2 * net.adj * membership.Z2 * model.mu.t()
        );
    lZ2 += 1.0/(2*model.sigma2) * (
            - net.Monest * membership.Z1 * (model.mu % model.mu)
            + 2 * net.adjt * membership.Z1 * model.mu
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

template<>
inline
double m_step(SBM & membership,
              gaussian & model,
              gaussian::network & net)
{
    model.mu = (membership.Z.t() * net.adjZD * membership.Z)
                /
               (membership.Z.t() * net.MonesZD * membership.Z);
    
    double all_accu_except_square_adj = accu(
            (
                (model.mu % model.mu)
                %
                (membership.Z.t() * net.MonesZD * membership.Z)
            )
            -
            (
                2 * model.mu
                %
                (membership.Z.t() * net.adjZD * membership.Z)
            )
        );

    model.sigma2 = 1.0/(membership.Z.n_rows * membership.Z.n_rows) * (
            net.accu_adjZD_square + all_accu_except_square_adj
        );

    return
        (
            -.5*(membership.Z.n_rows * (membership.Z.n_rows-1))*log(2*PI*model.sigma2)
            -1.0/(2*model.sigma2)*(net.accu_adjZD_square + all_accu_except_square_adj)
        );
}





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

template<>
inline
double m_step(SBM_sym & membership,
              gaussian & model,
              gaussian::network & net)
{
    return m_step<SBM>(membership,model,net)/2;
}





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

template<>
inline
double m_step(LBM & membership,
              gaussian & model,
              gaussian::network & net)
{
    model.mu = (membership.Z1.t() * net.adj * membership.Z2)
                /
               (membership.Z1.t() * net.Mones * membership.Z2);
    
    double all_accu_except_square_adj = accu(
            (
                (model.mu % model.mu)
                %
                (membership.Z1.t() * net.Mones * membership.Z2)
            )
            -
            (
                2 * model.mu
                %
                (membership.Z1.t() * net.adj * membership.Z2)
            )
        );

    model.sigma2 = 1.0/(membership.Z1.n_rows * membership.Z2.n_rows) * (
            net.accu_adj_square + all_accu_except_square_adj
        );

    return
        (
            -.5*(membership.Z1.n_rows*membership.Z2.n_rows)*log(2*PI*model.sigma2)
            -1.0/(2*model.sigma2)*(net.accu_adj_square + all_accu_except_square_adj)
        );
}





/* If you have defined m_step(SBM &, gaussian &, gaussian::network &)
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

// template<>
// inline
// double PL(gaussian & model,
//           SBM & membership,
//           gaussian::network & net)





/* If you have defined m_step(SBM_sym &, gaussian &, gaussian::network &)
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

// template<>
// inline
// double PL(gaussian & model,
//           SBM_sym & membership,
//           gaussian::network & net)





/* If you have defined m_step(LBM &, gaussian &, gaussian::network &)
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

// template<>
// inline
// double PL(gaussian & model,
//           LBM & membership,
//           gaussian::network & net)





/* If one of the following specialization
 *     - m_step(LBM &, gaussian &, gaussian::network)
 *     - m_step(SBM &, gaussian &, gaussian::network)
 *     - m_step(SBM_sym &, gaussian &, gaussian::network)
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

// inline
// vec gaussian::to_vector()
// {
// }





/* If m_step(SBM &, gaussian &, gaussian::network) is not defined
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

// gaussian::gaussian(SBM & membership, const vec & vectorized)
// {
// }





/* If m_step(SBM_sym &, gaussian &, gaussian::network) is not defined
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

// gaussian::gaussian(SBM_sym & membership, const vec & vectorized)
// {
// }





/* If m_step(LBM &, gaussian &, gaussian::network) is not defined
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

// gaussian::gaussian(LBM & membership, const vec & vectorized)
// {
// }





/* If one of the following specialization
 *     - m_step(LBM &, gaussian &, gaussian::network)
 *     - m_step(SBM &, gaussian &, gaussian::network)
 *     - m_step(SBM_sym &, gaussian &, gaussian::network)
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





/* If you have defined m_step(SBM &, gaussian &, gaussian::network &)
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

// template<>
// inline
// vec grad(gaussian & model,
//          SBM & membership,
//          gaussian::network & net)
// {
// }





/* If you have defined m_step(SBM_sym &, gaussian &, gaussian::network &)
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

// template<>
// inline
// vec grad(gaussian & model,
//          SBM_sym & membership,
//          gaussian::network & net)
// {
// }





/* If you have defined m_step(LBM &, gaussian &, gaussian::network &)
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

// template<>
// inline
// vec grad(gaussian & model,
//          LBM & membership,
//          gaussian::network & net)
// {
// }





/* If one of the following specialization
 *     - grad(gaussian &, LBM &, gaussian::network)
 *     - grad(gaussian &, SBM &, gaussian::network)
 *     - grad(gaussian &, SBM_sym &, gaussian::network)
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
// vec grad_logf(gaussian & model,
//               gaussian::network & net,
//               unsigned int i,
//               unsigned int j,
//               unsigned int q,
//               unsigned int l)
// {
// }





/* If grad_logf(gaussian &,
 *         gaussian::network &,
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
// double grad_logf(gaussian & model,
//                  gaussian::network & net,
//                  unsigned int i,
//                  unsigned int j,
//                  unsigned int q,
//                  unsigned int l,
//                  unsigned int k)
// {
// }





/* If one of the following specialization
 *     - PL(gaussian &, SBM &, gaussian::network)
 *     - PL(gaussian &, SBM_sym &, gaussian::network)
 *     - PL(gaussian &, LBM &, gaussian::network)
 *   is Usefull and not definied
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * logf(i,j,q,l) see notations.
 */

// inline
// double logf(gaussian & model,
//             gaussian::network & net,
//             unsigned int i,
//             unsigned int j,
//             unsigned int q,
//             unsigned int l)
// {
// }
