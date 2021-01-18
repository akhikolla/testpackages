
class gaussian_multivariate
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
        cube adj;


        /* Here you should add all precomputed values which depends only on the
         * network usefull in various functions, to avoid computing these value
         * many time.
         */
        cube adjZ;
        mat Mones;
        mat MonesZ;

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
            Rcpp::List adj_list = network_from_R["adjacency"];
            mat first_mat = Rcpp::as<mat>( adj_list[0] );
            adj.set_size(first_mat.n_rows, first_mat.n_cols, adj_list.size());
            for(int k=0; k<adj_list.size(); k++)
                adj.slice(k) = Rcpp::as<mat>(adj_list[k]);

            Mones = ones<mat>(first_mat.n_rows, first_mat.n_cols);
            MonesZ = fill_diag(Mones,0);
            adjZ.set_size(adj.n_rows,adj.n_cols,adj.n_slices);
            for(unsigned int k=0; k<adj.n_slices; k++)
            {
                adjZ.slice(k) = fill_diag(adj.slice(k),0);
            }
        }
    };

    // parameters
    unsigned int n_parameters;

    /* Here you must put the gaussian_multivariate model parameters */
    cube mu;
    mat Sigma;
    mat iL; // This is not a parameter, iL^T * iL = inv(Sigma)
            // iL is chol decompo of inv(Sigma)

    gaussian_multivariate(SBM & membership, gaussian_multivariate::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM case
         */
        n_parameters = membership.Z.n_cols * membership.Z.n_cols * net.adj.n_slices + net.adj.n_slices * net.adj.n_slices;
        mu.set_size(membership.Z.n_cols, membership.Z.n_cols, net.adj.n_slices);
        Sigma.set_size(net.adj.n_slices, net.adj.n_slices);
        iL.set_size(net.adj.n_slices, net.adj.n_slices);
    }
    
    gaussian_multivariate(SBM_sym & membership, gaussian_multivariate::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM_sym case
         */
        n_parameters = membership.Z.n_cols * (membership.Z.n_cols+1)/2 * net.adj.n_slices + net.adj.n_slices * net.adj.n_slices;
        mu.set_size(membership.Z.n_cols, membership.Z.n_cols, net.adj.n_slices);
        Sigma.set_size(net.adj.n_slices, net.adj.n_slices);
        iL.set_size(net.adj.n_slices, net.adj.n_slices);
    }
    
    gaussian_multivariate(LBM & membership, gaussian_multivariate::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * LBM case 
         */
        n_parameters = membership.Z1.n_cols * membership.Z2.n_cols * net.adj.n_slices + net.adj.n_slices * net.adj.n_slices;
        mu.set_size(membership.Z1.n_cols, membership.Z2.n_cols, net.adj.n_slices);
        Sigma.set_size(net.adj.n_slices, net.adj.n_slices);
        iL.set_size(net.adj.n_slices, net.adj.n_slices);
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;
     
        values["Sigma"] = Sigma;
        Rcpp::List mu_as_list(mu.n_slices);
        for(unsigned int k=0;k<mu.n_slices;k++)
        {
            mu_as_list[k]=mu.slice(k);
        }
        values["mu"] = mu_as_list;

        /* Here you must define the export way of your model parameters in the
         * list values. This is this list which is returned to R. */

        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    gaussian_multivariate(SBM &, const vec &);
    gaussian_multivariate(SBM_sym &, const vec &);
    gaussian_multivariate(LBM &, const vec &);
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
                  gaussian_multivariate & model,
                  gaussian_multivariate::network & net,
                  mat & lZ)
{
    cube X_tilde = apply_matrix_on_tubes(model.iL, net.adjZ);
    cube mu_tilde = apply_matrix_on_tubes(model.iL, model.mu);

    for(unsigned int k=0;k<net.adj.n_slices;k++)
    {
        lZ += .5 * (
                - net.MonesZ * membership.Z * (mu_tilde.slice(k).t() % mu_tilde.slice(k).t())
                + 2 * X_tilde.slice(k) * membership.Z * mu_tilde.slice(k).t()
                - net.MonesZ.t() * membership.Z * (mu_tilde.slice(k) % mu_tilde.slice(k))
                + 2 * X_tilde.slice(k).t() * membership.Z * mu_tilde.slice(k));
    }
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
                  gaussian_multivariate & model,
                  gaussian_multivariate::network & net,
                  mat & lZ)
{
    cube X_tilde = apply_matrix_on_tubes(model.iL, net.adjZ);
    cube mu_tilde = apply_matrix_on_tubes(model.iL, model.mu);

    for(unsigned int k=0;k<net.adj.n_slices;k++)
    {
        lZ +=  .5 * (
                - net.MonesZ * membership.Z * (mu_tilde.slice(k).t() % mu_tilde.slice(k).t())
                + 2 * X_tilde.slice(k) * membership.Z * mu_tilde.slice(k).t());
    }
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
                  gaussian_multivariate & model,
                  gaussian_multivariate::network & net,
                  mat & lZ1,
                  mat & lZ2)
{
    cube X_tilde = apply_matrix_on_tubes(model.iL, net.adj);
    cube mu_tilde = apply_matrix_on_tubes(model.iL, model.mu);

    for(unsigned int k=0;k<net.adj.n_slices;k++)
    {
        lZ1 += .5 * (
                - net.Mones * membership.Z2 * (mu_tilde.slice(k).t() % mu_tilde.slice(k).t())
                + 2 * X_tilde.slice(k) * membership.Z2 * mu_tilde.slice(k).t());
        lZ2 += .5 * (
                - net.Mones.t() * membership.Z1 * (mu_tilde.slice(k) % mu_tilde.slice(k))
                + 2 * X_tilde.slice(k).t() * membership.Z1 * mu_tilde.slice(k));
    }
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
              gaussian_multivariate & model,
              gaussian_multivariate::network & net)
{
    mat provdiv = membership.Z.t() * net.MonesZ * membership.Z;
    for(unsigned int k=0; k<net.adj.n_slices;k++)
    {
        model.mu.slice(k) = (membership.Z.t() * net.adjZ.slice(k) * membership.Z)
                            /
                            provdiv;
    }

    cube residual(net.adj.n_rows,net.adj.n_cols,net.adj.n_slices);
    for(unsigned int k=0;k<net.adj.n_slices;k++)
        residual.slice(k) = fill_diag(net.adj.slice(k) - membership.Z * model.mu.slice(k) * membership.Z.t(),0);

    for(unsigned int k1=0;k1<net.adj.n_slices;k1++)
        for(unsigned int k2=0;k2<net.adj.n_slices;k2++)
            model.Sigma(k1,k2) = 1.0 / (net.adj.n_rows * (net.adj.n_cols - 1) ) *
                accu(residual.slice(k1) % residual.slice(k2));

    // just to be sure:
    model.Sigma = .5 * (model.Sigma + model.Sigma.t());
    model.iL = chol( inv_sympd( model.Sigma ));

    cube X_tilde = apply_matrix_on_tubes(model.iL, net.adjZ);
    cube mu_tilde = apply_matrix_on_tubes(model.iL, model.mu);

    double PL= -.5*(membership.Z.n_rows * (membership.Z.n_rows-1) * net.adj.n_slices)*log(2*PI);
    PL += -.5*(membership.Z.n_rows * (membership.Z.n_rows-1)) * log( det(model.Sigma));
    PL += -.5 * accu(X_tilde % X_tilde);

    for(unsigned int k=0;k<net.adj.n_slices;k++)
        PL += -.5 * accu( 
                    (mu_tilde.slice(k) % mu_tilde.slice(k))
                    %
                    (membership.Z.t() * net.MonesZ * membership.Z)
                    )
               + accu( 
                    mu_tilde.slice(k)
                    %
                    (membership.Z.t() * X_tilde.slice(k) * membership.Z)
                    );
    
    return PL;
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
              gaussian_multivariate & model,
              gaussian_multivariate::network & net)
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
              gaussian_multivariate & model,
              gaussian_multivariate::network & net)
{
    mat provdiv = membership.Z1.t() * net.Mones * membership.Z2;
    for(unsigned int k=0; k<net.adj.n_slices;k++)
    {
        model.mu.slice(k) = (membership.Z1.t() * net.adj.slice(k) * membership.Z2)
                            /
                            provdiv;
    }

    cube residual(net.adj.n_rows,net.adj.n_cols,net.adj.n_slices);
    for(unsigned int k=0;k<net.adj.n_slices;k++)
        residual.slice(k) = net.adj.slice(k) - membership.Z1 * model.mu.slice(k) * membership.Z2.t();

    for(unsigned int k1=0;k1<net.adj.n_slices;k1++)
        for(unsigned int k2=0;k2<net.adj.n_slices;k2++)
            model.Sigma(k1,k2) = 1.0 / (net.adj.n_rows * net.adj.n_cols ) *
                accu(residual.slice(k1) % residual.slice(k2));

    // just to be sure:
    model.Sigma = .5 * (model.Sigma + model.Sigma.t());
    model.iL = chol( inv_sympd( model.Sigma ));

    cube X_tilde = apply_matrix_on_tubes(model.iL, net.adj);
    cube mu_tilde = apply_matrix_on_tubes(model.iL, model.mu);

    double PL= -.5*(membership.Z1.n_rows * membership.Z2.n_rows * net.adj.n_slices)*log(2*PI);
    PL += -.5*(membership.Z1.n_rows * membership.Z2.n_rows) * log( det(model.Sigma));
    PL += -.5 * accu(X_tilde % X_tilde);

    for(unsigned int k=0;k<net.adj.n_slices;k++)
        PL += -.5 * accu( 
                    (mu_tilde.slice(k) % mu_tilde.slice(k))
                    %
                    (membership.Z1.t() * net.MonesZ * membership.Z2)
                    )
               + accu( 
                    mu_tilde.slice(k)
                    %
                    (membership.Z1.t() * X_tilde.slice(k) * membership.Z2)
                    );
    
    return PL;
}





/* If you have defined m_step(SBM &, gaussian_multivariate &, gaussian_multivariate::network &)
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
// double PL(gaussian_multivariate & model,
//           SBM & membership,
//           gaussian_multivariate::network & net)





/* If you have defined m_step(SBM_sym &, gaussian_multivariate &, gaussian_multivariate::network &)
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
// double PL(gaussian_multivariate & model,
//           SBM_sym & membership,
//           gaussian_multivariate::network & net)





/* If you have defined m_step(LBM &, gaussian_multivariate &, gaussian_multivariate::network &)
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
// double PL(gaussian_multivariate & model,
//           LBM & membership,
//           gaussian_multivariate::network & net)





/* If one of the following specialization
 *     - m_step(LBM &, gaussian_multivariate &, gaussian_multivariate::network)
 *     - m_step(SBM &, gaussian_multivariate &, gaussian_multivariate::network)
 *     - m_step(SBM_sym &, gaussian_multivariate &, gaussian_multivariate::network)
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
// vec gaussian_multivariate::to_vector()
// {
// }





/* If m_step(SBM &, gaussian_multivariate &, gaussian_multivariate::network) is not defined
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

// gaussian_multivariate::gaussian_multivariate(SBM & membership, const vec & vectorized)
// {
// }





/* If m_step(SBM_sym &, gaussian_multivariate &, gaussian_multivariate::network) is not defined
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

// gaussian_multivariate::gaussian_multivariate(SBM_sym & membership, const vec & vectorized)
// {
// }





/* If m_step(LBM &, gaussian_multivariate &, gaussian_multivariate::network) is not defined
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

// gaussian_multivariate::gaussian_multivariate(LBM & membership, const vec & vectorized)
// {
// }





/* If one of the following specialization
 *     - m_step(LBM &, gaussian_multivariate &, gaussian_multivariate::network)
 *     - m_step(SBM &, gaussian_multivariate &, gaussian_multivariate::network)
 *     - m_step(SBM_sym &, gaussian_multivariate &, gaussian_multivariate::network)
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
 *
 * You can specialize this template in three function for SBM, LBM, and SBM_sym
 * if you need
 */

// template<class membership_type>
// inline
// double maximum_step_in_direction(membership_type & membership,
//                                  gaussian_multivariate & model,
//                                  gaussian_multivariate::network & net,
//                                  vec & direction)
// {
// }





/* If you have defined m_step(SBM &, gaussian_multivariate &, gaussian_multivariate::network &)
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
// vec grad(gaussian_multivariate & model,
//          SBM & membership,
//          gaussian_multivariate::network & net)
// {
// }





/* If you have defined m_step(SBM_sym &, gaussian_multivariate &, gaussian_multivariate::network &)
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
// vec grad(gaussian_multivariate & model,
//          SBM_sym & membership,
//          gaussian_multivariate::network & net)
// {
// }





/* If you have defined m_step(LBM &, gaussian_multivariate &, gaussian_multivariate::network &)
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
// vec grad(gaussian_multivariate & model,
//          LBM & membership,
//          gaussian_multivariate::network & net)
// {
// }





/* If one of the following specialization
 *     - grad(gaussian_multivariate &, LBM &, gaussian_multivariate::network)
 *     - grad(gaussian_multivariate &, SBM &, gaussian_multivariate::network)
 *     - grad(gaussian_multivariate &, SBM_sym &, gaussian_multivariate::network)
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
// vec grad_logf(gaussian_multivariate & model,
//               gaussian_multivariate::network & net,
//               unsigned int i,
//               unsigned int j,
//               unsigned int q,
//               unsigned int l)
// {
// }





/* If grad_logf(gaussian_multivariate &,
 *         gaussian_multivariate::network &,
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
// double grad_logf(gaussian_multivariate & model,
//                  gaussian_multivariate::network & net,
//                  unsigned int i,
//                  unsigned int j,
//                  unsigned int q,
//                  unsigned int l,
//                  unsigned int k)
// {
// }





/* If one of the following specialization
 *     - PL(gaussian_multivariate &, SBM &, gaussian_multivariate::network)
 *     - PL(gaussian_multivariate &, SBM_sym &, gaussian_multivariate::network)
 *     - PL(gaussian_multivariate &, LBM &, gaussian_multivariate::network)
 *   is Usefull and not definied
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * logf(i,j,q,l) see notations.
 */

// inline
// double logf(gaussian_multivariate & model,
//             gaussian_multivariate::network & net,
//             unsigned int i,
//             unsigned int j,
//             unsigned int q,
//             unsigned int l)
// {
// }
