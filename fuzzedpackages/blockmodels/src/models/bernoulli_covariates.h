
/*
 * This is package skeleton for bernoulli_covariates model
 *
 * You should delete this header when the bernoulli_covariates model is ready for use to avoid
 * any confusion
 */

class bernoulli_covariates
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

        mat adj; // adjacency matrix
        cube covariates; // cube of covariates the 3rd dimention is the index
                         // of covariate


        /* Here you should add all precomputed values which depends only on the
         * network usefull in various functions, to avoid computing these value
         * many time.
         */

        mat adjZD;
        mat Mones;
        mat MonesZD;

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


        }
    };

    // parameters
    unsigned int n_parameters;
    bool symmetric;

    /* Here you must put the bernoulli_covariates model parameters */
    mat m;
    colvec beta;

    bernoulli_covariates(SBM & membership, bernoulli_covariates::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM case
         */

        m = (membership.Z.t() * net.adjZD * membership.Z)
            /
            (membership.Z.t() * net.MonesZD * membership.Z);

        m=logit(m);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = m.n_elem + beta.n_elem;
        symmetric=false;
    }
    
    bernoulli_covariates(SBM_sym & membership, bernoulli_covariates::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM_sym case
         */
        
        m = (membership.Z.t() * net.adjZD * membership.Z)
            /
            (membership.Z.t() * net.MonesZD * membership.Z);
        
        m=logit(m);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = m.n_rows*(m.n_rows+1)/2 + beta.n_elem;
        symmetric=true;
    }
    
    bernoulli_covariates(LBM & membership, bernoulli_covariates::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * LBM case 
         */
        m = (membership.Z1.t() * net.adj * membership.Z2)
            /
            (membership.Z1.t() * net.Mones * membership.Z2);
        
        m=logit(m);
        
        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = m.n_rows * m.n_cols + net.covariates.n_slices;
        symmetric=false;
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;

        /* Here you must define the export way of your model parameters in the
         * list values. This is this list which is returned to R. */

        values["m"] = m;
        values["beta"] = beta;

        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    bernoulli_covariates(SBM &, const vec &);
    bernoulli_covariates(SBM_sym &, const vec &);
    bernoulli_covariates(LBM &, const vec &);
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

// template<>
// inline
// void e_fixed_step(SBM & membership,
//                   bernoulli_covariates & model,
//                   bernoulli_covariates::network & net,
//                   mat & lZ)
// {
// }





/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ_{iq} += \sum_{jl} Z_{jl} logf(i,j,q,l)
 *
 * This is usefull for SBM_sym. See bernoulli and poisson for example.
 */

// template<>
// inline
// void e_fixed_step(SBM_sym & membership,
//                   bernoulli_covariates & model,
//                   bernoulli_covariates::network & net,
//                   mat & lZ)
// {
// }





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

// template<>
// inline
// void e_fixed_step(LBM & membership,
//                   bernoulli_covariates & model,
//                   bernoulli_covariates::network & net,
//                   mat & lZ1,
//                   mat & lZ2)
// {
// }





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
//               bernoulli_covariates & model,
//               bernoulli_covariates::network & net)
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
//               bernoulli_covariates & model,
//               bernoulli_covariates::network & net)
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
//               bernoulli_covariates & model,
//               bernoulli_covariates::network & net)
// {
// }





/* If you have defined m_step(SBM &, bernoulli_covariates &, bernoulli_covariates::network &)
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
// double PL(bernoulli_covariates & model,
//           SBM & membership,
//           bernoulli_covariates::network & net)





/* If you have defined m_step(SBM_sym &, bernoulli_covariates &, bernoulli_covariates::network &)
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
// double PL(bernoulli_covariates & model,
//           SBM_sym & membership,
//           bernoulli_covariates::network & net)





/* If you have defined m_step(LBM &, bernoulli_covariates &, bernoulli_covariates::network &)
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
// double PL(bernoulli_covariates & model,
//           LBM & membership,
//           bernoulli_covariates::network & net)





/* If one of the following specialization
 *     - m_step(LBM &, bernoulli_covariates &, bernoulli_covariates::network)
 *     - m_step(SBM &, bernoulli_covariates &, bernoulli_covariates::network)
 *     - m_step(SBM_sym &, bernoulli_covariates &, bernoulli_covariates::network)
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
vec bernoulli_covariates::to_vector()
{
    vec out(n_parameters);
    vec vm = (symmetric) ? vech(m) : reshape(m,m.n_elem,1);

    out.subvec(0,vm.n_elem-1) = vm;
    out.subvec(vm.n_elem,n_parameters-1) = beta;
    return(out);
}





/* If m_step(SBM &, bernoulli_covariates &, bernoulli_covariates::network) is not defined
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

bernoulli_covariates::bernoulli_covariates(SBM & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    m = reshape(vectorized.subvec(0,Q*Q-1),Q,Q);
    beta = vectorized.subvec(Q*Q,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = false;
}





/* If m_step(SBM_sym &, bernoulli_covariates &, bernoulli_covariates::network) is not defined
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

bernoulli_covariates::bernoulli_covariates(SBM_sym & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    m = unvech(vectorized.subvec(0,Q*(Q+1)/2-1));
    beta = vectorized.subvec(Q*(Q+1)/2,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = true;
}





/* If m_step(LBM &, bernoulli_covariates &, bernoulli_covariates::network) is not defined
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

bernoulli_covariates::bernoulli_covariates(LBM & membership, const vec & vectorized)
{
    unsigned int Q1=membership.Z1.n_cols;
    unsigned int Q2=membership.Z2.n_cols;
    m = reshape(vectorized.subvec(0,Q1*Q2-1),Q1,Q2);
    beta = vectorized.subvec(Q1*Q2,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric=false;
}





/* If one of the following specialization
 *     - m_step(LBM &, bernoulli_covariates &, bernoulli_covariates::network)
 *     - m_step(SBM &, bernoulli_covariates &, bernoulli_covariates::network)
 *     - m_step(SBM_sym &, bernoulli_covariates &, bernoulli_covariates::network)
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
double maximum_step_in_direction(membership_type & membership, bernoulli_covariates & model, bernoulli_covariates::network & net, vec & direction)
{
    return 1;
}





/* If you have defined m_step(SBM &, bernoulli_covariates &, bernoulli_covariates::network &)
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
// vec grad(bernoulli_covariates & model,
//          SBM & membership,
//          bernoulli_covariates::network & net)
// {
// }





/* If you have defined m_step(SBM_sym &, bernoulli_covariates &, bernoulli_covariates::network &)
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
// vec grad(bernoulli_covariates & model,
//          SBM_sym & membership,
//          bernoulli_covariates::network & net)
// {
// }





/* If you have defined m_step(LBM &, bernoulli_covariates &, bernoulli_covariates::network &)
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
// vec grad(bernoulli_covariates & model,
//          LBM & membership,
//          bernoulli_covariates::network & net)
// {
// }


/* If one of the following specialization
 *     - grad(bernoulli_covariates &, LBM &, bernoulli_covariates::network)
 *     - grad(bernoulli_covariates &, SBM &, bernoulli_covariates::network)
 *     - grad(bernoulli_covariates &, SBM_sym &, bernoulli_covariates::network)
 *   is Usefull and not defined
 * Then
 *     Usefull, Optional
 * Else
 *     Useless
 *
 * The gradient of logf(i,j,q,l) in respect to the parameter vector.
 */

template<>
inline
vec grad_logf(bernoulli_covariates & model,
              bernoulli_covariates::network & net,
              unsigned int i,
              unsigned int j,
              unsigned int q,
              unsigned int l)
{
    unsigned int m_n=(model.symmetric) ?
        model.m.n_rows*(model.m.n_rows+1)/2
        :
        model.m.n_elem;

    unsigned int m_k=model.m.n_rows*l+q;

    if(model.symmetric)
    {
        if(q>l)
        {
            unsigned int p=q;
            q=l;
            l=p;
        }
        
        m_k=(2*model.m.n_rows-1-q)*q/2+l;
    }

    colvec yij = net.covariates.tube(i,j);
    double c = model.m(q,l)+vec((model.beta.t() * yij))(0);
    double constante = net.adj(i,j)-sigmo(c);

    vec out(model.n_parameters);
    for(unsigned p=0;p<m_n;p++)
    {
        if(p==m_k)
            out(p)=constante;
        else
            out(p)=0;
    }

    out.subvec(m_n,out.n_elem-1) = constante * yij;

    return out;
}





/* If grad_logf(bernoulli_covariates &,
 *         bernoulli_covariates::network &,
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
// double grad_logf(bernoulli_covariates & model,
//                  bernoulli_covariates::network & net,
//                  unsigned int i,
//                  unsigned int j,
//                  unsigned int q,
//                  unsigned int l,
//                  unsigned int k)
// {
// }





/* If one of the following specialization
 *     - PL(bernoulli_covariates &, SBM &, bernoulli_covariates::network)
 *     - PL(bernoulli_covariates &, SBM_sym &, bernoulli_covariates::network)
 *     - PL(bernoulli_covariates &, LBM &, bernoulli_covariates::network)
 *   is Usefull and not definied
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * logf(i,j,q,l) see notations.
 */

inline
double logf(bernoulli_covariates & model,
            bernoulli_covariates::network & net,
            unsigned int i,
            unsigned int j,
            unsigned int q,
            unsigned int l)
{
    colvec yij = net.covariates.tube(i,j);
    double c = model.m(q,l)+vec((model.beta.t() * yij))(0);
    return net.adj(i,j)*c + log(1-sigmo(c));
}
